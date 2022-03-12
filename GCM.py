#!/usr/bin/env python
"""GCM.py
Implementation of the generalized spatial connectome model.
This file reproduces Fig. 3D of the manuscipt.
requires file: 
1. "journal.pcbi.1007974.s007.xlsx"
    from Fenyves et al. 2020.
    (10.1371/journal.pcbi.1007974)
outputs files:
1.  fig3_panelD.png:
    Precision rank plot of the predicted polarities
2.  predicted_pairs_sl3.csv:
    comma separated file of the predictions with rows
    corresponding to synapses indicated by neuron pairs 
    and columns corresponding to 
    1. presynaptic neuron name
    2. postsynaptic neuron name
    3. "actual" sign determined by the CM
    4. source neuron index
    5. target neuron index
    6. pseudo inverse predicted sign ($\alpha \rightarrow 0$)
    7. regularized predicted sign (optimal $\alpha$)
    8. heavily regularized sign ($\alpha \rightarrow \infty$)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.sparse import kron,coo_matrix
from scipy.sparse.linalg import lsqr,eigs,svds
from scipy.optimize import minimize_scalar
from sklearn.model_selection import train_test_split
from collections import defaultdict

plt.rcParams['svg.fonttype'] = 'none'

def read_XOY(fn='journal.pcbi.1007974.s007.xlsx',sheet_nm='5. Sign prediction'):
    df = pd.read_excel(fn,engine='openpyxl',sheet_name='1. NT expr',usecols=[0,1,2],index_col=0) ## from Fenyves et al.
    df2 = df[~df['Dominant NT'].isna()] ## remove neurons without NT expression
    X_df = pd.DataFrame(np.zeros((df2.shape[0],3)),index=df2.index,columns=['ACh','GABA','Glu']) ## Ensure consistent labeling
    nt2ind = {'ACh':0,'GABA':1,'Glu':2} ## order the NTs alphabetically
    for nm,row in df.iterrows():
        for col,elt in row[~row.isna()].iteritems():
            X_df.loc[nm,elt]=1
    ## regulation rules are from Fenyves, saved as a pkl file.
    regulation_dict = pd.read_pickle('regulation_rules2.pkl')
    recept = pd.read_csv('r_genes.tsv',sep='\t',index_col=1) 
    nz_recept = recept[recept.astype(int).iloc[:,1]>0].iloc[:,0]
    O_df = pd.DataFrame(np.zeros((3,nz_recept.shape[0])),index=X_df.columns,columns=np.sort(nz_recept.index))
    for (nt,reg),gns in regulation_dict.items():
        O_df.loc[nt,gns] = 1 if reg=='Pos' else -1
    ydf = pd.read_excel(fn,engine='openpyxl',sheet_name='3. Receptor Gene expression',skiprows=1, header=None, index_col=0)
    ydf = ydf[~ydf.index.isna()] ## remove any NAs
    Y_df = pd.DataFrame(np.zeros((O_df.shape[1],X_df.shape[0])),index=O_df.columns,columns=X_df.index)
    for neu,row in ydf.iterrows():
        for elt in row[~row.isna()].values:
            if elt=='0':
                break
            Y_df.loc[elt,neu]=1
    Y_df=Y_df.fillna(0)
    Y_df = Y_df[Y_df.sum(axis=1)>0]
    cY_df = Y_df.loc[Y_df.index.intersection(O_df.columns)] ## removes gar-2
    cO_df = O_df.loc[:,O_df.columns.intersection(cY_df.index)] ## restrict to receptors that are expressed
    cY_df = cY_df.loc[:,np.sort(cY_df.columns)] ## alphabetize the ordering of neurons
    cX_df = X_df.loc[np.sort(X_df.index.intersection(cY_df.columns))] ## alphabetize the ordering of neurons
    conn = pd.read_excel(fn,engine='openpyxl',sheet_name=sheet_nm,skiprows=1, header=0)
    conn = conn.set_index(['Neuron','Neuron.1']) ## index by neuron names
    conn = conn.loc[:,['1mary NT','2ndary NT','Edge weight','Edge type','Unnamed: 16']]
    conn.columns = ['1mary NT','2ndary NT','Edge weight','Edge type','Predicted'] ## rename columns
    conn.index.names = ['Source','Target'] #rename index columns
    conn = conn[~conn.index.duplicated()]
    return cX_df,cO_df,cY_df,conn

def min_alpha(al,*args):
    KK_trim_spv,aa_vec,evals = args
    soln= lsqr(KK_trim_spv,aa_vec.flatten(),damp=al)[0]
    r_sqr = np.linalg.norm(KK_trim_spv@soln-aa_vec)**2
    tau_sqr = np.sum(1/(evals**2/al**2+1))**2
    return r_sqr/tau_sqr

def min_alpha_95pct(inp_vals,*args):
    KK_trim_spv,aa_vec,evals,KK_spv_all,resolved_inds,actual_signs = args
    prec95max = {}
    for xx in inp_vals:
        soln, istop, itn, r1norm, r2norm, anorm, acond, arnorm, xnorm, var= lsqr(KK_trim_spv,aa_vec.flatten(),damp=xx)
        wpredicted_signs = (KK_spv_all@soln)[resolved_inds]
        ASRT_INDS = np.argsort(np.abs(wpredicted_signs))[::-1]
        TF_arr = np.sign(wpredicted_signs[ASRT_INDS])==actual_signs[ASRT_INDS]
        xvs = np.arange(1,TF_arr.shape[0]+1)
        precision = np.cumsum(TF_arr)/xvs
        gt95 = np.where(precision>.95)[0]
        if gt95.shape[0]>0:
            prec95max[xx] = np.amax(gt95)
        else:
            prec95max[xx] = 0
    return pd.Series(prec95max)

def GSCM(fn,sheet_nm,opt_reg=True,opt_damp=-1):
    ## determine which connectome we are using
    if 'Cook' in sheet_nm:
        cn = 'Cook'
    elif 'Varshney' in sheet_nm:
        cn = 'Varshney'
    else:
        cn = 'wormwiring'
    ## Read in data from Fenyves et al (.s007)
    X_df,O_df,Y_df,conn = read_XOY(fn,sheet_nm)
    ## Select only chemical synapses (exclude electrical synapses)
    conn = conn[conn['Edge type']=='chemical']
    ## Remove any postsynaptic neurons without R expressed genes
    selcY = Y_df.T[Y_df.sum(axis=0)>0].T
    ## Organize the postsynaptic neurons by labels
    NNames = np.sort(Y_df.columns)
    neuro2col = pd.Series(dict([(v,k) for k,v in enumerate(Y_df.columns)]))
    ## XOY is essentially the algorithm applied in Fenyves
    XOY = (X_df@O_df.fillna(0)@selcY)
    SS = XOY.stack() ## flatten the matrix
    SS = SS[SS.index.isin(conn.index)] ## restrict to ONLY chemical synapse edges in the connectome
    SS = SS.to_frame()
    SS.columns = ['weight']
    SS['polarity'] = conn.loc[SS.index,'Predicted'] ## add the feynves polarities
    sel_SS = SS[SS.polarity.isin(['+','-'])] ## restrict to "known" polarities
    
    sel_SS_us = sel_SS.weight.unstack()
    INDS = sel_SS_us.index
    COLS = sel_SS_us.columns
    ## construct sparse matrices
    sel_SS_s = coo_matrix(sel_SS_us.fillna(0).values)
    NROW,NCOL = sel_SS_s.shape
    ROW_inds,COL_inds = sel_SS_s.nonzero()
    aa_vec = np.sign(sel_SS_s.data)
    sel_inds = ROW_inds*NCOL+COL_inds
    sel_SS_sgn = sel_SS_s.copy()
    sel_SS_sgn.data = np.sign(sel_SS_sgn.data)
    ## calculate the Kronecker product
    KK_sl3 = kron(sel_SS_sgn,sel_SS_sgn.T)
    ## truncate to only the measured indices
    KK_trim_sl3 = KK_sl3.tocsr()[sel_inds] 
    ## identify the synapses with complex polarity, nonzero weight, and in the connectome
    ## these are the synapses whose sign is "resolved" by the model
    resolved_pairs_sl3 = SS[(SS.polarity=='complex')&(SS.weight!=0)&(SS.index.get_level_values(0).isin(INDS))&(SS.index.get_level_values(1).isin(COLS))]
    row_inds_sl3 = np.asarray(list(map(INDS.get_loc,resolved_pairs_sl3.index.get_level_values(0))))
    col_inds_sl3 = np.asarray(list(map(COLS.get_loc,resolved_pairs_sl3.index.get_level_values(1))))
    resolved_inds_sl3 = NCOL*row_inds_sl3+col_inds_sl3
    actual_signs_sl3 = np.sign(resolved_pairs_sl3.weight.values)
    presyn_neuron = resolved_pairs_sl3.index.get_level_values(0) 
    postsyn_neuron = resolved_pairs_sl3.index.get_level_values(1) 

    fig,ax = plt.subplots(1,1)
    alph_preds_d = {}
    prec95max_sl3 = {}
    for ii,dmp in enumerate([0,1,10,95,1e5]):
        ## perform regularized fit
        soln, istop, itn, r1norm, r2norm, anorm, acond, arnorm, xnorm, var= lsqr(KK_trim_sl3,aa_vec,damp=dmp)
        ## find the predicted signs
        wpredicted_signs_sl3 = (KK_sl3@soln)[resolved_inds_sl3]
        ## order the preditions by purported relevance
        ASRT_INDS = np.argsort(np.abs(wpredicted_signs_sl3))[::-1]
        ## make a precision plot
        TF_arr = np.sign(wpredicted_signs_sl3[ASRT_INDS])==actual_signs_sl3[ASRT_INDS]
        xvs = np.arange(1,TF_arr.shape[0]+1)
        precision = np.cumsum(TF_arr)/xvs
        lbl = '%.1e' % dmp**2 if dmp !=0 else str(dmp)
        ## optimality criterion for precision:
        ##    maximize the number of elements above 97 % precision
        nn95 =np.where(precision>0.97)[0]
        if nn95.shape[0]>0:
            prec95max_sl3[dmp] = np.amax(nn95)
        else:
            prec95max_sl3[dmp] = 0
        if ii==1:
            clr = 'C3'
        elif ii==3:
            clr = 'C1'
        else:
            clr = 'C%d' %ii
        ## make the plot
        ax.plot(xvs,precision,label='%s' % lbl,color=clr)
        ## store the predictions
        if dmp==0:
            alph_preds_d['pseudoinverse_predicted_sign'] = wpredicted_signs_sl3
        elif dmp==95:
            alph_preds_d['regularized_predicted_sign'] = wpredicted_signs_sl3
        elif dmp==1e5:
            alph_preds_d['heavily_regularized_sign'] = wpredicted_signs_sl3

    ax.legend(title=r'$\alpha$',prop={'size':6},ncol=2)
    ax.set_ylim(0.0,1.05)
    ax.axhline(0.95,color='C7',ls='--')
    ax.axhline(0.54,color='C7',ls=':')
    #ax.set_ylabel('Precision',size=8)
    ax.set_xlabel('Rank',size=8)
    ax.set_yticks(np.linspace(0,1,6))
    ax.set_yticklabels(['%.1f' %ee for ee in np.linspace(0,1,6)])
    ax.set_xlim(-1,273)
    ax.set_xticks(np.linspace(0,250,6))
    ax.set_xticklabels(['%d' %ee for ee in np.linspace(0,250,6)])
    
    ## organize the predictions into a DataFrame and output to a CSV
    csv_op_d = {'presynaptic_neuron':presyn_neuron,
                'postsynaptic_neuron':postsyn_neuron,
                'actual_sign':actual_signs_sl3,
                'source':neuro2col.loc[presyn_neuron].values,
                'target':neuro2col.loc[postsyn_neuron].values
               }
    csv_op_d.update(alph_preds_d)
    pred_df_sl3 = pd.DataFrame(csv_op_d)
    pred_df_sl3.to_csv('predicted_pairs_sl3.csv')

    plt.setp(ax.get_xticklabels(),size=8)
    plt.setp(ax.get_yticklabels(),size=8)
    fig.savefig('fig3_panelD.png')
    return

if __name__ == '__main__':
    GSCM('journal.pcbi.1007974.s007.xlsx','5. Sign prediction (Cook)')
    #predict_SCM('../journal.pcbi.1007974.s007.xlsx','5. Sign prediction (Varshney)')
