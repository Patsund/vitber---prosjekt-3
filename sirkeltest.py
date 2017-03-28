import math
import numpy as np
from matplotlib import pyplot as plt
X=[40,40]
def eq1(X):
    x=X[0]
    y=X[1]
    R=math.sqrt(x**2+y**2)
    T=24*60*60 #setter perioden til 1 døgn i sekunder
    factor=2*math.pi/T
    Y=[x,y]
    Y[0]+=factor*(-y)
    Y[1]+=factor*x
    return Y
def f(foo,X,t): #Beregner w(i+1) med Eulers metode
    x=X[0]
    y=X[1]
    Z=eq1(X)
    Y=[X[0]+t*Z[0],X[1]+t*Z[1]]
    x=Y[0]
    y=Y[1]
    R=math.sqrt(x**2+y**2)
    return Y,R
def plotThing(X,h):
    R=math.sqrt(X[0]**2+X[1]**2)
    R2=0
    x,y=X[0],X[1]
    Y=[x,y]
    for i in range(2*24*60*60//(h)):
        Y,R2=f(Y,h)
    error=abs(R-R2)
    return error
times=[30,60,4*60,8*60,15*60,30*60,3600,2*3600,4*3600,8*3600]
errorArr=[]
print("X:",X)
for h in times:
    thiserror = plotThing(X,h)
    print(thiserror)
    errorArr.append(thiserror)
print("X:",X)
h0=1 #sekund
thiserror = plotThing(X,h0)
while thiserror<10:
    print("error:",thiserror,"må være sant",thiserror<10)
    X=[40,40] #Hvorfor?
    h0+=1
    thiserror = plotThing(X,h0)
print("for å få en error større enn 10 meter må tidssteget være",h0,"sekunder")
plt.loglog(times,errorArr)
plt.xlabel("tid")
plt.ylabel("error")
plt.show()