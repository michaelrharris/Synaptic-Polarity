“Known Network.csv" - Network of known positive and negative connections from Fenyves et al., based on the Cook et al. connectome reconstruction

“Connectome.csv” - Contact adjacency matrix for the entire verified space of connections of the Cook et al. reconstruction

“ConnectomeEdgelist.csv” - The contact adjacency matrix "Connectome.csv" in edge-list format

“X.csv” - Each of the 295 neurons of interest and their accompanying neurotransmitters

“Y.csv” - Each of the 295 neurons of interest and their accompanying receptor genes

“O.csv” - Signed wiring rule network of interactions between neurotransmitters and receptor genes

"ovec.csv" - Vector format of the updated signed wiring rule network produced by the SCM

"O~.csv" - Matrix format of the updated signed wiring rule network produced by the SCM

"Kprime.csv" - Kronecker product of "X" and "Y" truncated to only include connections with known polarity

"K.csv" - Kronecker product of "X" and "Y"

"journal.pcbi.1007974.s007" - Source data from Fenyves et al. used in this work

"Neuron Tags.csv" - Labels of the 295 neurons present in "X", "Y", the known network, and connectome constructions in ascending order.

"Neurotransmitters.csv" - Labels of the 3 neurotransmitters present in "X", "O", and "O~" in ascending order

"Receptor Genes.csv" - Labels of the 42 receptor genes present in "Y", "O", and "O~" in ascending order

“SCM.py” - Implementation of the Spatial Connectome Model which produces the results "SCM_Predictions.csv", in addition to the updated signed wiring rule networks "ovec" and "O~". Running this code requires the files "X.csv", "Y.csv", "Known Network.csv", "Connectome.csv", and "ConnectomeEdgelist.csv" to be downloaded to the same directory as SCM.py.

“Network Based Prediction Cross Validation.py” - Performs a 10-fold cross validation procedure using the SCM, SL3, SPA, and SL2 methods, and produces the precision rank plot of Figure 3, panel E in this work. Running this code requires the files "X.csv", "Y.csv", "Known Network.csv", "Connectome.csv", and "ConnectomeEdgelist.csv" to be downloaded to the same directory as Network Based Prediction Cross Validation.py.

"GCM.py" - Implementation of the Generalized Connectome Model which produces the GCM predictions in the file "predicted_pairs_sl3" and the precision rank plot of Figure 3, panel D in this work. Running this code requires Python (v 3.7.6 / anaconda), scipy (v. 1.7.1), numpy (v. 1.20.1), matplotlib (v. 3.4.2), pandas (v. 1.1.3), as well as the files "journal.pcbi.1007974.s007.xlsx", "regulation_rules2.pkl", and "r_genes.tsv" to be downloaded to the same directory as GCM.py.   