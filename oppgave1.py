import numpy as np
from matplotlib import pyplot as plt

alpha = 5e-5
m = 1e-2

def Vwater(t,X):
    T=24*3600
    factor = 2 * np.pi / T
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


def fForEq2(X,t):
    dX=np.array(Vwater(t,X))
    return dX

def fForEq1(X,t):
    positions=[X[0],X[1]]
    Vw=Vwater(t,positions)
    dx=X[2]
    dy=X[3]
    dvx=alpha/m*(Vw[0]-X[2])
    dvy=alpha/m*(Vw[1]-X[3])
    return [dx,dy,dvx,dvy]


def rk2(X,Vwater,f,h,t): #X er en array som kan inneholde enten x,y eller x,y,vx,vy
    K1=f(X,t)
    K1=np.array(K1)
    K2=f(X+h*K1,t+h)
    X_=X+h/2*(K1+K2)
    return X_

def errorForRK2(X,h):
    newX=np.array(X)
    xList,yList=[],[]
    xList.append(newX[0])
    yList.append(newX[1])
    timeFinal,timeNow=48*3600,0
    n=int(timeFinal/h)
    for i in range(n+1):
        h=min(h,timeFinal-timeNow)
        newX=rk2(newX,Vwater,fForEq2,h,timeNow)
        xList.append(newX[0])
        yList.append(newX[1])
        timeNow+=h
    X=np.array(X)
    newX=np.array(newX)
    vecError=X-newX
    plt.plot(xList,yList)
    plt.show()
    error=np.sqrt(vecError[0]**2+vecError[1]**2)
    return error

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

def task1b():
    X=[0,100]
    errorForRK2(X,10000)
def task1c():
    ##Constants
    L = 1.0E+02
    alpha = 5e-5
    m = 1e-2

    numberOfTimesteps = 10
    analyticEndpoint = np.array(analyticSolution(L, alpha, m))
    print("Analytisk ende", analyticEndpoint)
    timestepArray = np.linspace(390, 402, num=numberOfTimesteps)
    errorArray = np.zeros(numberOfTimesteps)
    for i in range(numberOfTimesteps): #numberOfTimesteps
        h = timestepArray[i]
        plt.figure(i)
        plt.title("tidssteg " + str(h))
        coordinateArray = np.array([[0, 0]])
        X = np.array([L, 0, 0, 0])
        timeNow=0
        timeFinal=48*3600
        for j in range(int(2 * 24 * 60 * 60 / h) + 1):
            h=min(h,timeFinal-timeNow)
            X=rk2(X,Vwater,fForEq1,h,timeNow)
            timeNow+=h
            dummyarray=[X[0],X[1]]
            dummyarray = np.array([dummyarray])
            coordinateArray = np.concatenate((coordinateArray, dummyarray), axis=0)

        # Her sorteres arrayet så det kan plottes
        coordinateArray = np.delete(coordinateArray, 0, 0)  # Fjerner 0,0 som ble brukt for å initialisere arrayet
        coordinateArray = sorted(coordinateArray, key=lambda e: e[0])
        xValueArray = [k[0] for k in coordinateArray]
        yValueArray = [k[1] for k in coordinateArray]
        plt.plot(xValueArray, yValueArray, 'ro', markersize=1)
        plt.plot([analyticEndpoint[0]], [analyticEndpoint[1]], 'ko', markersize=4)
        Xcoord=[X[0],X[1]]
        errorArray[i] = np.linalg.norm(np.array(Xcoord) - analyticEndpoint)

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
    Xdot = np.array([0, 0])
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
            Xtrap = rk2(X,Vwater,fForEq1,h,timeNow)
            thisDeviation = np.linalg.norm(Xtrap - Xeul)
            while thisDeviation>deviationLimit:
                h-=1
                Xeul, XdotEul = eulerForEq1(X,Vwater,alpha,m,Xdot,h)
                Xtrap
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


#task1b()
task1c()