“Known Network.csv" - Network of known positive and negative connections from Fenyves et al., based on the Cook et al. connectome reconstruction


“Connectome.csv” - Contact adjacency matrix for the entire verified space of connections, known and unknown, of the Cook et al. reconstruction


“ConnectomeEdgelist.csv” - The contact adjacency matrix "Connectome.csv" in edge-list format


“X.csv” - Each of the 295 neurons of interest and their accompanying neurotransmitters


“Y.csv” - Each of the 295 neurons of interest and their accompanying receptor genes


“O.csv” - Signed wiring rule network of interactions between neurotransmitters and receptor genes


"ovec.csv" - Vector format of the updated signed wiring rule network produced by the SCM


"O~.csv" - Matrix format of the updated signed wiring rule network produced by the SCM


"Kprime.csv" - Kronecker product of "X" and "Y" truncated to only include connections with known polarity 


"Source Data - Fenyves et al..xlsx" - Source data used in this work


"Neuron Tags.csv" - Labels of the 295 neurons present in "X", "Y", the known network, and connectome constructions in ascending order.


"Neurotransmitters.csv" - Labels of the 3 neurotransmitters present in "X", "O", and "O~" in ascending order


"Receptor Genes.csv" - Labels of the 42 receptor genes present in "Y", "O", and "O~" in ascending order


“SCM.py” - Implementation of the Spatial Connectome Model which produces the results "SCM_Predictions.csv", in addition to the updated signed wiring rule networks "ovec" and "O~". Running this code requires the files "X.csv", "Y.csv", "Known Network.csv", "Connectome.csv", and "ConnectomeEdgelist.csv".


“Network Based Prediction Cross Validation.py” - Performs a 10-fold cross validation procedure using the SCM, SL3, SPA, and SL2 methods, and produces a precision rank plot like that of Figure 3, panel E in this work. Running this code requires the files "X.csv", "Y.csv", "Known Network.csv", "Connectome.csv", and "ConnectomeEdgelist.csv".
