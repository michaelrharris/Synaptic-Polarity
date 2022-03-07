import numpy as np
import csv

reader=csv.reader(open("Connectome.csv"), delimiter=",")
c=np.array(list(reader)).astype('float')
reader=csv.reader(open("X.csv"), delimiter=",")
X=np.array(list(reader)).astype("float")
reader=csv.reader(open("Y.csv"), delimiter=",")
Y=np.array(list(reader)).astype("float")
reader=csv.reader(open("Known Network.csv"), delimiter=",")
a=np.array(list(reader)).astype("float")
reader=csv.reader(open("ConnectomeEdgelist.csv"), delimiter=",")
connedge=np.array(list(reader)).astype("str")

for x in range(len(connedge)):
    if connedge[x][2] == "no pred":
        c[int(connedge[x][0])][int(connedge[x][1])] = 0
    if connedge[x][2] == "complex":
        c[int(connedge[x][0])][int(connedge[x][1])] = 0
                     
k = np.kron(X,Y)
cf=c.flatten()
kp=np.delete(k,np.nonzero(1-cf),0)
af=a.flatten()
ap=af[cf==1]
alpha = 31.92
o = np.linalg.inv(kp.T@kp + alpha*(np.identity(len(kp.T@kp))))@kp.T@ap
O = np.reshape(o,(len(X.T),len(Y.T)))

reader=csv.reader(open("Connectome.csv"), delimiter=",")
c=np.array(list(reader)).astype('float')

A = X@O@Y.T
for x in range(len(A)):
    for y in range(len(A)):
        if c[x][y] == 0:
            A[x][y] = 0
            
