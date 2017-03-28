__author__ = 'Patrik'
# -*- coding: utf-8 -*-
import math
import numpy as np
from matplotlib import pyplot as plt
X=[40,40]
def f(X,t):
    x=X[0]
    y=X[1]
    R=math.sqrt(x**2+y**2)
    T=24*60*60
    factor=2*math.pi/T
    Y=[x,y]
    Y[0]+=t*factor*(-y)
    Y[1]+=t*factor*x
    x=Y[0]
    y=Y[1]
    R=math.sqrt(x**2+y**2)
    return Y,R
def findError\
                (X,h):
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
#print("X:",X)
for h in times:
    print(findError
          (X,h))
    errorArr.append(findError
                    (X,h))
#print("X:",X)
X=[0,100]
h0=1
#print("X:",X)
while findError\
            (X,h0)<10:
    #print("X:",X)
    print("error:",findError
    (X,h0),"må være sant",findError
    (X,h0)<10)
    h0+=1
print("for å få en error større enn 10 meter må tidssteget være",h0,"sekunder")
plt.loglog(times,errorArr)
plt.show()