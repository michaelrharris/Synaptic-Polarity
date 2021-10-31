import numpy as np
import csv
import random
import pandas as pd
import matplotlib.pyplot as plt


finalprecision = []
connectionlist = []
reader=csv.reader(open("C:/Users/youve/OneDrive/Desktop/Data/ConnectomeEdgelist.csv"), delimiter=",")
connedge=np.array(list(reader)).astype("str")
reader = csv.reader(open("C:/Users/youve/OneDrive/Desktop/Data/known.csv"), delimiter = ',')
matrix = list(reader)
matrix = np.array(matrix).astype("float")

for x in range(len(matrix)):
        for y in range(len(matrix)):
            if matrix[x][y] != 0:
                connectionlist.append((x,y))
random.shuffle(connectionlist)
tenfold = [[] for j in range(10)]

for x in range(len(connectionlist)):
    if x < 175:
        tenfold[0].append((connectionlist[x][0],connectionlist[x][1]))
    if (x <= 350) and (x > 175):
        tenfold[1].append((connectionlist[x][0],connectionlist[x][1]))
    if (x <= 525) and (x > 350):
        tenfold[2].append((connectionlist[x][0],connectionlist[x][1]))
    if (x <= 700) and (x > 525):
        tenfold[3].append((connectionlist[x][0],connectionlist[x][1]))
    if (x <= 875) and (x > 700):
        tenfold[4].append((connectionlist[x][0],connectionlist[x][1]))
    if (x <= 1050) and (x > 875):
        tenfold[5].append((connectionlist[x][0],connectionlist[x][1]))
    if (x <= 1225) and (x > 1050):
        tenfold[6].append((connectionlist[x][0],connectionlist[x][1]))
    if (x <= 1400) and (x > 1225):
        tenfold[7].append((connectionlist[x][0],connectionlist[x][1]))
    if (x <= 1575) and (x > 1400):
        tenfold[8].append((connectionlist[x][0],connectionlist[x][1]))
    if (x > 1575):
        tenfold[9].append((connectionlist[x][0],connectionlist[x][1]))

for qq in range(10):  
    reader=csv.reader(open("C:/Users/youve/OneDrive/Desktop/Data/Connectome.csv"), delimiter=",")
    ta=np.array(list(reader)).astype("float")
    reader=csv.reader(open("C:/Users/youve/OneDrive/Desktop/Data/X.csv"), delimiter=",")
    tx=np.array(list(reader)).astype("float")
    reader=csv.reader(open("C:/Users/youve/OneDrive/Desktop/Data/Y.csv"), delimiter=",")
    ty=np.array(list(reader)).astype("float")
    reader=csv.reader(open("C:/Users/youve/OneDrive/Desktop/Data/known.csv"), delimiter=",")
    tb=np.array(list(reader)).astype("float")
    
    links = []
    for x in range(len(connedge)):
        if connedge[x][2] == "no pred":
            ta[int(connedge[x][0])][int(connedge[x][1])] = 0
        if connedge[x][2] == "complex":
            ta[int(connedge[x][0])][int(connedge[x][1])] = 0
                
    for x in range(len(tenfold[qq])):
        links.append(tenfold[qq][x])
        ta[tenfold[qq][x][0]][tenfold[qq][x][1]] = 0
     
    k=np.kron(tx,ty)
    taf=ta.flatten()
    kp=np.delete(k,np.nonzero(1-taf),0)
    tbf=tb.flatten()
    tbp=tbf[taf==1]
    kk=np.dot(kp.transpose(),kp) 
    alpha=31.92
    ka=np.linalg.inv(kk+alpha*np.identity(kk.shape[1]))
    ka=np.dot(ka,kp.transpose())
    ky=np.dot(kp,ka)
    tau=np.trace(np.identity(ky.shape[1])-ky)
    tof=np.dot(ka,tbp)
    r=np.linalg.norm(tbp-np.dot(kp,tof))
    obj=r*r/tau/tau
    tof=tof.reshape(tx.shape[1],ty.shape[1])   
    kpinv=np.linalg.pinv(kp)
    oMP=np.dot(kpinv,tbp)
    of=oMP.reshape(tx.shape[1],ty.shape[1])
    
    reader=csv.reader(open("C:/Users/youve/OneDrive/Desktop/Data/Connectome.csv"), delimiter=",")
    connectome = np.array(list(reader)).astype("float")
    
    B = tx@tof@ty.T
    for x in range(len(B)):
        for y in range(len(B)):
            if connectome[x][y] ==0:
                B[x][y] = 0
                 
    presults = pd.DataFrame({"value":[], 
                        "pos1":[],  
                        "pos2":[],
                        "sign":[]})
     
    for x in range(len(tb)):
        for y in range(len(tb)):
            if (x,y) in links:
                if B[x][y] < 0:
                    dfupdate1 = pd.DataFrame([[abs(B[x][y]),x,y,"neg"]],columns=['value','pos1','pos2','sign'])
                    presults = presults.append(dfupdate1, ignore_index = True)
                elif B[x][y] > 0:
                    dfupdate1 = pd.DataFrame([[abs(B[x][y]),x,y,"pos"]],columns=['value','pos1','pos2','sign'])
                    presults = presults.append(dfupdate1, ignore_index = True)
                elif B[x][y] == 0:
                    dfupdate1 = pd.DataFrame([[abs(B[x][y]),x,y,"zero"]],columns=['value','pos1','pos2','sign'])
                    presults = presults.append(dfupdate1, ignore_index = True)
            else:
                continue
            
    presults.sort_values(by=['value'], inplace=True, ascending=True)
    sortedresults = presults.sort_values(by=['value'])
    sortedresults.reset_index(inplace = True) 
    vector = []
                
    for o in range(len(sortedresults)):
        if (sortedresults.loc[o,'pos1'],sortedresults.loc[o,'pos2']) in links:
            if (matrix[int(sortedresults.loc[o,'pos1'])][int(sortedresults.loc[o,'pos2'])] < 0) and (sortedresults.loc[o,'sign'] == "neg"):
                vector.append((1,int(sortedresults.loc[o,'pos1']),int(sortedresults.loc[o,'pos2'])))
            elif (matrix[int(sortedresults.loc[o,'pos1'])][int(sortedresults.loc[o,'pos2'])] > 0) and (sortedresults.loc[o,'sign'] == "pos"):
                vector.append((1,int(sortedresults.loc[o,'pos1']),int(sortedresults.loc[o,'pos2'])))
            else:
                vector.append((0,int(sortedresults.loc[o,'pos1']),int(sortedresults.loc[o,'pos2'])))
        else:
            continue
    
    vector.reverse() 
    l3sum = 0
    precisiononly = []
    
    for oj in range(175):
        l3sum += vector[oj][0]
        l3precision = l3sum/(oj+1)
        precisiononly.append(l3precision)
    
    for x in range(175):
        if qq == 0:
            finalprecision.append(precisiononly[x])
        else:
            finalprecision[x] += precisiononly[x]

    print(qq)
    
for x in range(175):
    finalprecision[x] = finalprecision[x]/10

connectionlist = []

reader = csv.reader(open("C:/Users/youve/OneDrive/Desktop/Data/known.csv"), delimiter = ',')
refmat = list(reader)
refmat = np.array(refmat).astype("float")
finalprc1 = [0 for x in range(175)]
finalprc2 =[0 for x in range(175)]
finalprc3 =[0 for x in range(175)]

for qq in range(10):
    reader = csv.reader(open("C:/Users/youve/OneDrive/Desktop/Data/known.csv"), delimiter = ',')
    matrix = list(reader)
    matrix = np.array(matrix).astype("float")
    training = []
    testing = []
    for y in range(len(tenfold[qq])):
        matrix[int(tenfold[qq][y][0])][int(tenfold[qq][y][1])] = 0
        testing.append(tenfold[qq][y])
        
    reader = csv.reader(open("C:/Users/youve/OneDrive/Desktop/Data/Connectome.csv"), delimiter = ',')
    connectome = list(reader)
    connectome = np.array(connectome).astype("float")
    positivein = []
    positiveout = []
    negativein = []
    negativeout = []
    
    for x in range(len(matrix)):
        posin = 0
        negin = 0
        for y in range(len(matrix)):
            if matrix[x][y] > 0:
                posin += (matrix[x][y])
            elif matrix[x][y] < 0:
                negin += (matrix[x][y])
            elif matrix[x][y] == 0:
                continue
        positivein.append(posin)
        negativein.append(negin)
    for x in range(len(matrix)):
        posout = 0
        negout = 0
        for y in range(len(matrix)):
            if matrix[y][x] > 0:
                posout += (matrix[y][x])
            elif matrix[y][x] < 0:
                negout += (matrix[y][x])
            elif matrix[y][x] == 0:
                continue
        positiveout.append(posout)
        negativeout.append(negout)
        
    SPA = [[0 for x in range(len(matrix))] for y in range(len(matrix))]
    SPA = np.array(SPA).astype("float")
    presults = pd.DataFrame({"value":[], 
                        "pos1":[],  
                        "pos2":[],
                        "origin":[],
                        "sign":[]})
    presults2 = pd.DataFrame({"value":[], 
                        "pos1":[],  
                        "pos2":[],
                        "origin":[],
                        "sign":[]})
    presults3 = pd.DataFrame({"value":[], 
                        "pos1":[],  
                        "pos2":[],
                        "origin":[],
                        "sign":[]})
    
    prcL3list = []
    prcL3list2 = []
    prcL3list3 = []   
    vectorL3 = []
    vector2L3 = []
    vector3L3 = []   
    l3sum = 0
    l3sum2 = 0
    l3sum3 = 0
    
    connconnects = []
    for x in range(len(matrix)):
        for y in range(len(matrix)):
            if connectome[x][y] != 0:
                connconnects.append((x,y))
                SPA[x][y] = (positivein[x]*positiveout[y])-(negativein[x]*negativeout[y])
    SL3 = matrix@matrix.T@matrix
    SL2 = matrix.T@matrix + matrix@matrix.T
         
    for x in range(len(matrix)):
        for y in range(len(matrix)):
            if (x,y) in connconnects:
                if SL3[x][y] < 0:
                    dfupdate1 = pd.DataFrame([[abs(SL3[x][y]),x,y,0,"neg"]],columns=['value','pos1','pos2','origin','sign'])
                    presults = presults.append(dfupdate1, ignore_index = True)
                elif SL3[x][y] > 0:
                    dfupdate1 = pd.DataFrame([[abs(SL3[x][y]),x,y,0,"pos"]],columns=['value','pos1','pos2','origin','sign'])
                    presults = presults.append(dfupdate1, ignore_index = True)
                elif SL3[x][y] == 0:
                    dfupdate1 = pd.DataFrame([[abs(SL3[x][y]),x,y,0,"zero"]],columns=['value','pos1','pos2','origin','sign'])
                    presults = presults.append(dfupdate1, ignore_index = True)  
                if SPA[x][y] < 0:
                    dfupdate2 = pd.DataFrame([[abs(SPA[x][y]),x,y,0,"neg"]],columns=['value','pos1','pos2','origin','sign'])
                    presults2 = presults2.append(dfupdate2, ignore_index = True)
                elif SPA[x][y] > 0:
                    dfupdate2 = pd.DataFrame([[abs(SPA[x][y]),x,y,0,"pos"]],columns=['value','pos1','pos2','origin','sign'])
                    presults2 = presults2.append(dfupdate2, ignore_index = True)
                elif SPA[x][y] == 0:
                    dfupdate2 = pd.DataFrame([[abs(SPA[x][y]),x,y,0,"zero"]],columns=['value','pos1','pos2','origin','sign'])
                    presults2 = presults2.append(dfupdate2, ignore_index = True)            
                if SL2[x][y] < 0:
                    dfupdate3 = pd.DataFrame([[abs(SL2[x][y]),x,y,0,"neg"]],columns=['value','pos1','pos2','origin','sign'])
                    presults3 = presults3.append(dfupdate3, ignore_index = True)
                elif SL2[x][y] > 0:
                    dfupdate3 = pd.DataFrame([[abs(SL2[x][y]),x,y,0,"pos"]],columns=['value','pos1','pos2','origin','sign'])
                    presults3 = presults3.append(dfupdate3, ignore_index = True)
                elif SL2[x][y] == 0:
                    dfupdate3 = pd.DataFrame([[abs(SL2[x][y]),x,y,0,"zero"]],columns=['value','pos1','pos2','origin','sign'])
                    presults3 = presults3.append(dfupdate3, ignore_index = True)
            else:
                continue              
    presults.sort_values(by=['value'], inplace=True, ascending=True)
    sortedresults = presults.sort_values(by=['value'])
    sortedresults.reset_index(inplace = True)
    presults2.sort_values(by=['value'], inplace=True, ascending=True)
    sortedresults2 = presults2.sort_values(by=['value'])
    sortedresults2.reset_index(inplace = True)
    presults3.sort_values(by=['value'], inplace=True, ascending=True)
    sortedresults3 = presults3.sort_values(by=['value'])
    sortedresults3.reset_index(inplace = True)              
    for o in range(len(sortedresults)):
        if (sortedresults.loc[o,'pos1'],sortedresults.loc[o,'pos2']) in testing:
            if (refmat[int(sortedresults.loc[o,'pos1'])][int(sortedresults.loc[o,'pos2'])] < 0) and (sortedresults.loc[o,'sign'] == "neg"):
                vectorL3.append((1,int(sortedresults.loc[o,'pos1']),int(sortedresults.loc[o,'pos2'])))
            elif (refmat[int(sortedresults.loc[o,'pos1'])][int(sortedresults.loc[o,'pos2'])] > 0) and (sortedresults.loc[o,'sign'] == "pos"):
                vectorL3.append((1,int(sortedresults.loc[o,'pos1']),int(sortedresults.loc[o,'pos2'])))
            else:
                vectorL3.append((0,int(sortedresults.loc[o,'pos1']),int(sortedresults.loc[o,'pos2'])))
        else:
            continue
    for o in range(len(sortedresults)):
        if (sortedresults2.loc[o,'pos1'],sortedresults2.loc[o,'pos2']) in testing:
            if (refmat[int(sortedresults2.loc[o,'pos1'])][int(sortedresults2.loc[o,'pos2'])] < 0) and (sortedresults2.loc[o,'sign'] == "neg"):
                vector2L3.append((1,int(sortedresults2.loc[o,'pos1']),int(sortedresults2.loc[o,'pos2'])))
            elif (refmat[int(sortedresults2.loc[o,'pos1'])][int(sortedresults2.loc[o,'pos2'])] > 0) and (sortedresults2.loc[o,'sign'] == "pos"):
                vector2L3.append((1,int(sortedresults2.loc[o,'pos1']),int(sortedresults2.loc[o,'pos2'])))
            else:
                vector2L3.append((0,int(sortedresults2.loc[o,'pos1']),int(sortedresults2.loc[o,'pos2'])))
        else:
            continue
    for o in range(len(sortedresults)):
        if (sortedresults3.loc[o,'pos1'],sortedresults3.loc[o,'pos2']) in testing:
            if (refmat[int(sortedresults3.loc[o,'pos1'])][int(sortedresults3.loc[o,'pos2'])] < 0) and (sortedresults3.loc[o,'sign'] == "neg"):
                vector3L3.append((1,int(sortedresults3.loc[o,'pos1']),int(sortedresults3.loc[o,'pos2'])))
            elif (refmat[int(sortedresults3.loc[o,'pos1'])][int(sortedresults3.loc[o,'pos2'])] > 0) and (sortedresults3.loc[o,'sign'] == "pos"):
                vector3L3.append((1,int(sortedresults3.loc[o,'pos1']),int(sortedresults3.loc[o,'pos2'])))
            else:
                vector3L3.append((0,int(sortedresults3.loc[o,'pos1']),int(sortedresults3.loc[o,'pos2'])))
        else:
            continue
 
    vectorL3.reverse()
    vector2L3.reverse()
    vector3L3.reverse()
    
    for oj in range(175):
        l3sum += vectorL3[oj][0]
        l3precision = l3sum/(oj+1)
        prcL3list.append((l3precision,vectorL3[oj][1],vectorL3[oj][2]))   
        l3sum2 += vector2L3[oj][0]
        l3precision2 = l3sum2/(oj+1)
        prcL3list2.append((l3precision2,vector2L3[oj][1],vector2L3[oj][2]))    
        l3sum3 += vector3L3[oj][0]
        l3precision3 = l3sum3/(oj+1)
        prcL3list3.append((l3precision3,vector3L3[oj][1],vector3L3[oj][2]))
        
    for x in range(175):
        finalprc1[x] += prcL3list[x][0] 
        finalprc2[x] += prcL3list2[x][0]
        finalprc3[x] += prcL3list3[x][0]
    
    print(qq)
rank = []
for jk in range(175):
    finalprc1[jk] = finalprc1[jk]/10
    finalprc2[jk] = finalprc2[jk]/10
    finalprc3[jk] = finalprc3[jk]/10
    rank.append(jk)
plt.figure(dpi=300)
plot4 = plt.plot(rank,finalprecision, label = 'SCM', color = 'red')
plot2 = plt.plot(rank,finalprc2, label = 'SPA', color = 'teal')
plot1 = plt.plot(rank,finalprc1, label = 'SL3', color = 'orange')
plot3 = plt.plot(rank,finalprc3, label = 'SL2', color = 'green')
plt.ylabel("Precision")
plt.xlabel("Rank")
plt.xlim(0,175)
plt.ylim(0,1.1)
plt.legend(loc = 'lower right')
plt.show()
