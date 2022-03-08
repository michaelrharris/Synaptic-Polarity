“Known Network.csv" - Network of known positive and negative connections


“Connectome.csv” - Contact adjacency matrix for the entire verified space of connections, known and unknown


“ConnectomeEdgelist.csv” - The contact adjacency matrix in edge-list format


“X.csv” - Each of the 295 neurons of interest and their accompanying neurotransmitters


“Y.csv” - Each of the 295 neurons of interest and their accompanying receptor genes


“O.csv” - Signed wiring rule network of interactions between neurotransmitters and receptor genes


"ovec.csv" - Vector format of the updated signed wiring rule network produced by the SCM


"O~.csv" - Matrix format of the updated signed wiring rule network produced by the SCM


"Kprime.csv" - Kronecker product of X and Y truncated to only include connections with known polarity 


“SCM.py” - Implementation of the Spatial Connectome Model. 


“Network Based Prediction Cross Validation.py” - performs a 10-fold cross validation procedure using the SCM, SL3, SPA, and SL2 methods, and produces a precision rank plot like that of Figure 3, panel E in this work.


"Source Data - Fenyves et al..xlsx" - Source data used in this work


"Neuron Tags.csv" - Labels of the 295 neurons present in "X", "Y", the known network, and connectome constructions in ascending order.


"Neurotransmitters.csv" - Labels of the 3 neurotransmitters present in "X", "O", and "O~" in ascending order


"Receptor Genes.csv" - Labels of the 42 receptor genes present in "Y", "O", and "O~" in ascending order
