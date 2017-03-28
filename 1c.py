import math
import numpy as np
from matplotlib import pyplot as plt

def Vwater(t,X):
    x,y=X[0],X[1]
    T=24*60*60
    factor=2*math.pi/T
    Y=[factor*(-X[1]),factor*X[0]]
    return Y

def eulerForEq1(X,Vwater,alpha,m,Xdot,h):
    x,y=X[0],X[1]
    VwaterVector=Vwater(X)
    print("Vw",VwaterVector,"Xdot",Xdot)
    x+=h*Xdot[0]
    y+=h*Xdot[1]
    newX=[x,y]
    speedx,speedy=Xdot[0],Xdot[1]
    speedx+=h*(alpha/m)*(VwaterVector[0]-speedx)
    speedy+=h*(alpha/m)*(VwaterVector[1]-speedy)
    newXdot=[speedx,speedy]
    return newX,newXdot

def ETMforEq1(X,Vwater,alpha,m,Xdot,h,t):
    x,y=X[0],X[1]
    waterVelocity=Vwater(t,X)
    R=math.sqrt(X[0]**2+X[1]**2)
    print("R",R,"Vw",waterVelocity,"Xdot",Xdot)
    XdotTilde=[Xdot[0]+h*(alpha/m)*(waterVelocity[0]-Xdot[0]),Xdot[1]+h*(alpha/m)*(waterVelocity[1]-Xdot[1])]
    Xtilde=[X[0]+h*Xdot[0],X[1]+h*Xdot[1]]
    waterVelocityTilde=Vwater(t+h,Xtilde)
    newXdot=[Xdot[0]+h/2*(alpha/m*(waterVelocity[0]-Xdot[0])+alpha/m*(waterVelocityTilde[0]-(Xdot[0]+h*alpha/m*(waterVelocity[0]-Xdot[0])))),Xdot[1]+h/2*(alpha/m*(waterVelocity[1]-Xdot[1])+alpha/m*(waterVelocityTilde[1]-(Xdot[1]+h*alpha/m*(waterVelocity[1]-Xdot[1]))))]
    newX=[X[0]+h/2*(Xdot[0]+newXdot[0]),X[1]+h/2*(Xdot[1]+newXdot[1])]
    return newX,newXdot

X=[100,0]
Xdot=[0,0]
step=30
alpha=5e-5
m=1e-2
print(X)
for i in range(2*24*60*60//step):
    X,Xdot=ETMforEq1(X,Vwater,alpha,m,Xdot,step,step)
print(X)