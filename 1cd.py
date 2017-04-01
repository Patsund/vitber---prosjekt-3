import math
import numpy as np
from matplotlib import pyplot as plt

T = 24 * 60 * 60
factor = 2 * math.pi / T

def analyticSolution(L,alpha,m):
    k = alpha/m
    w = 2*np.pi/(24*60*60)*k
    A = np.array([
    [ 0, 0, 1, 0],
    [ 0, 0, 0, 1],
    [ 0,-w,-k, 0],
    [ w, 0, 0,-k]
    ])
    lams , V = np.linalg.eig(A)
    X0 = np.array([L,0.,0.,0.])
    C = np.linalg.solve(V,X0)
    t = 2*24*3600
    X = V.dot(C*np.exp(lams*t)) # pun
    X = [X[0].real, X[1].real]
    #print("Analytic position after 48 hours = ", X[0].real, X[1].real)
    return X

def Vwater(t,X):
    Y = [factor*(-X[1]), factor*X[0]]
    return Y

def eulerForEq1(X,Vwater,alpha,m,Xdot,h):
    x,y=X[0],X[1]
    VwaterVector=Vwater(h,X)
    x+=h*Xdot[0]
    y+=h*Xdot[1]
    newX=[x,y]
    speedx,speedy=Xdot[0],Xdot[1]
    speedx+=h*(alpha/m)*(VwaterVector[0]-speedx)
    speedy+=h*(alpha/m)*(VwaterVector[1]-speedy)
    newXdot=[speedx,speedy]
    return newX,newXdot

def ETMforEq1(X,Vwater,alpha,m,Xdot,h,t):
    waterVelocity=np.array(Vwater(t,X))
    XdotTilde=Xdot+h*(alpha/m)*(waterVelocity-Xdot)
    Xtilde=X+h*Xdot
    waterVelocityTilde = np.array(Vwater(t + h, Xtilde))
    newXdot=Xdot+h/2*(alpha/m*(waterVelocity-Xdot)+alpha/m*(waterVelocityTilde-(XdotTilde)))
    newX=X+h/2*(Xdot+newXdot)
    return newX,newXdot

def task1c():
    ##Constants
    L = 1.0E+02
    alpha = 5e-5
    m = 1e-2

    numberOfTimesteps = 10
    analyticEndpoint = np.array(analyticSolution(L, alpha, m))
    print("Analytisk ende", analyticEndpoint)
    timestepArray = np.linspace(10, 300, num=numberOfTimesteps)
    errorArray = np.zeros(numberOfTimesteps)
    for i in range(numberOfTimesteps): #numberOfTimesteps
        h = timestepArray[i]
        plt.figure(i)
        plt.title("tidssteg " + str(h))
        coordinateArray = np.array([[0, 0]])
        X = np.array([L, 0])
        Xdot = np.array([0, 0.000001])
        for j in range(int(2 * 24 * 60 * 60 / h) + 1):
            X, Xdot = ETMforEq1(X, Vwater, alpha, m, Xdot, h, h * j)
            dummyarray = np.array([X])
            # print(dummyarray)
            coordinateArray = np.concatenate((coordinateArray, dummyarray), axis=0)
        h = (2 * 24 * 60 * 60) % h
        X, Xdot = ETMforEq1(X, Vwater, alpha, m, Xdot, h, 2 * 24 * 60 * 60 - h)
        print("thisEndpoint:", X)

        # Her sorteres arrayet så det kan plottes
        coordinateArray = np.delete(coordinateArray, 0, 0)  # Fjerner 0,0 som ble brukt for å initialisere arrayet
        coordinateArray = sorted(coordinateArray, key=lambda e: e[0])
        xValueArray = [k[0] for k in coordinateArray]
        yValueArray = [k[1] for k in coordinateArray]
        plt.plot(xValueArray, yValueArray, 'ro', markersize=1)
        plt.plot([analyticEndpoint[0]], [analyticEndpoint[1]], 'ko', markersize=4)
        errorArray[i] = np.linalg.norm(np.array(X) - analyticEndpoint)
        print(errorArray)

    plt.figure()
    plt.plot(timestepArray,errorArray)
    plt.xlabel("timestep / sekunder")
    plt.ylabel("error / meter")
    plt.title("Avvik fra analytisk løsning")
    plt.show()

def task1d():
    ##Constants
    L = 1.0E+02
    alpha = 5e-5
    m = 1e-2
    timeNow = 0
    totalTime = 48* 3600
    deviationLimit = 0.1
    X = np.array([L, 0])
    Xdot = np.array([0, 0.000001])
    timeStepArray = []
    deviationArray = []
    coordinateArray = []
    analyticEndpoint = analyticSolution(L, alpha, m)
    print("Analytisk ende", analyticEndpoint)
    while (timeNow < totalTime):
        #print("starter løkke")
        if timeNow==0:
            h=100
            print("initializing h from",h,"seconds")
            Xeul, XdotEul = eulerForEq1(X,Vwater,alpha,m,Xdot,h)
            Xtrap, XdotTrap = ETMforEq1(X, Vwater, alpha, m, Xdot, h, timeNow)
            thisDeviation = np.linalg.norm(Xtrap - Xeul)
            while thisDeviation>deviationLimit:
                h-=1
                Xeul, XdotEul = eulerForEq1(X,Vwater,alpha,m,Xdot,h)
                Xtrap, XdotTrap = ETMforEq1(X, Vwater, alpha, m, Xdot, h, timeNow)
                thisDeviation = np.linalg.norm(Xtrap - Xeul)
            print("h initialized to",h,"seconds")
        h = min(h, totalTime - timeNow)
        Xeul, XdotEul = eulerForEq1(X,Vwater,alpha,m,Xdot,h)
        Xtrap, XdotTrap = ETMforEq1(X, Vwater, alpha, m, Xdot, h, timeNow)
        thisDeviation = np.linalg.norm(Xtrap - Xeul)
        while (thisDeviation > deviationLimit):
            print(timeNow,"###\nNY LØKKE\n###")
            h /= 2
            Xeul, XdotEul = eulerForEq1(X, Vwater, alpha, m, Xdot, h)
            Xtrap, XdotTrap = ETMforEq1(X, Vwater, alpha, m, Xdot, h, timeNow)
            thisDeviation = np.linalg.norm(Xtrap - Xeul)
        #Legger til timestep og deviation idet de blir godkjent
        timeStepArray.append([timeNow, h])
        deviationArray.append(thisDeviation)
        timeNow+=h
        if (thisDeviation < deviationLimit/10):
            print(timeNow,"###\nØKER H\n###")
            h *= 2
        coordinateArray.append(Xtrap)
        X, Xdot = Xtrap, XdotTrap
    xValueArray = [k[0] for k in coordinateArray]
    yValueArray = [k[1] for k in coordinateArray]
    timeStepArray_xValues = [k[0] for k in timeStepArray]
    timeStepArray_yValues = [k[1] for k in timeStepArray]
    print(timeStepArray)
    plt.figure()
    plt.title("Tidssteg")
    plt.plot(timeStepArray_xValues, timeStepArray_yValues, 'go', markersize=1)
    plt.xlabel("tid / sekunder")
    plt.ylabel("tidssteg / sekunder")
    plt.figure()
    plt.title("Avvik mellom Euler og Trapes")
    plt.plot(timeStepArray_xValues, deviationArray, 'ko', markersize=1)
    plt.xlabel("tid / sekunder")
    plt.ylabel("avvik / meter")
    plt.figure()
    plt.title("Bane med variabelt tidssteg")
    plt.plot(xValueArray, yValueArray, 'ro', markersize=0.4)
    print("Final error:", np.linalg.norm(X - analyticEndpoint))
    plt.show()

#task1c()
task1d()