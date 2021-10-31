import numpy as np
import csv


reader=csv.reader(open("C:/Users/youve/OneDrive/Desktop/Data/Connectome.csv"), delimiter=",")
ta=np.array(list(reader)).astype('float')
reader=csv.reader(open("C:/Users/youve/OneDrive/Desktop/Data/X.csv"), delimiter=",")
tx=np.array(list(reader)).astype("float")
reader=csv.reader(open("C:/Users/youve/OneDrive/Desktop/Data/Y.csv"), delimiter=",")
ty=np.array(list(reader)).astype("float")
reader=csv.reader(open("C:/Users/youve/OneDrive/Desktop/Data/known.csv"), delimiter=",")
tb=np.array(list(reader)).astype("float")
reader=csv.reader(open("C:/Users/youve/OneDrive/Desktop/Data/ConnectomeEdgelist.csv"), delimiter=",")
connedge=np.array(list(reader)).astype("str")


for x in range(len(connedge)):
    if connedge[x][2] == "no pred":
        ta[int(connedge[x][0])][int(connedge[x][1])] = 0
    if connedge[x][2] == "complex":
        ta[int(connedge[x][0])][int(connedge[x][1])] = 0


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
ta=np.array(list(reader)).astype("float")

B = tx@tof@ty.T
for x in range(len(B)):
    for y in range(len(B)):
        if ta[x][y] == 0:
            B[x][y] = 0

            
