import math
import numpy as np
from matplotlib import pyplot as plt
X=[40,40]

def waterVel(X):
    x,y=X[0],X[1]
    R=math.sqrt(x**2+y**2)
    T=24*60*60 #setter perioden til 1 døgn i sekunder
    factor=2*math.pi/T
    Y=[x,y]
    Y[0]+=factor*(-y)
    Y[1]+=factor*x
    Y=np.array(Y)
    return Y,R

def eulerForEq2(X,t): #Beregner w(i+1) med Eulers metode
    x=X[0]
    y=X[1]
    Y=[x,y]
    T=24*60*60
    factor=2*math.pi/T
    Y=[Y[0]+t*factor*(-y),Y[1]+t*factor*(x)]
    x=Y[0]
    y=Y[1]
    R=math.sqrt(x**2+y**2)
    return Y,R

def ETM(X,t): #Beregner w(i+1) med Eulers metode
    x=X[0]
    y=X[1]
    R=math.sqrt(x**2+y**2)
    T=24*60*60 #setter perioden til 1 døgn i sekunder
    factor=2*math.pi/T
    Y=[x,y]
    W=[0,0]
    W[0]+=factor*(-y)
    W[1]+=factor*x
    #W[0]+=t*factor*-(W[1])
    #W[1]+=t*factor*(W[0])
    Y[0]+=t/2*(factor*(-X[1])+factor*(-(X[1]+t*factor*(-X[1]))))
    Y[1]+=t/2*(factor*X[0]+factor*(X[0]+t*factor*X[0]))
    x=Y[0]
    y=Y[1]
    R=math.sqrt(x**2+y**2)
    return Y,R

vector=[0,100]
vector2=[0,100]
tid=60*60

for i in range(24*60*60//tid):
    vector,R=ETM(vector,tid)
for i in range(24*60*60//tid):
    vector2,R2=eulerForEq2(vector2,tid)
print(vector, R)
print(vector2,R2)

def plotThing(X,h):
    #print(X)
    R=math.sqrt(X[0]**2+X[1]**2)
    R2=0
    x,y=X[0],X[1]
    Y=[x,y]
    for i in range(2*24*60*60//(h)):
        Y,R2=eulerForEq2(Y,h)
    error=abs(R-R2)
    return error
times=[30,60,4*60,8*60,15*60,30*60,3600,2*3600,4*3600,8*3600]
errorArr=[]
X=[0,100]
#print("X:",X)
for h in times:
    errorArr.append(plotThing(X,h))
#print("X:",X)
h0=1 #sekund
thiserror = plotThing(X,h0)
while thiserror<10:
    ##print("error:",thiserror,"må være sant",thiserror<10)
    h0+=1
    thiserror = plotThing(X,h0)
#print("for å få en error større enn 10 meter må tidssteget være",h0,"sekunder")
plt.loglog(times,errorArr)
plt.xlabel("tid")
plt.ylabel("error")
plt.show()