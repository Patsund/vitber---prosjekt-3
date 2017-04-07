import numpy as np
from matplotlib import pyplot as plt

alpha = 5e-5
m = 1e-2

def Vwater(t,X):
    T=24*3600
    factor = 2 * np.pi / T
    newX = [factor*(-X[1]), factor*X[0]]
    return newX

def eulerForEq1(X,t,Vwater,Xdot,h):
    x,y=X[0],X[1]
    VwaterVector=Vwater(t,X)
    x+=h*Xdot[0]
    y+=h*Xdot[1]
    newX=[x,y]
    dx,dy=Xdot[0],Xdot[1]
    dx+=h*(alpha/m)*(VwaterVector[0]-dx)
    dy+=h*(alpha/m)*(VwaterVector[1]-dy)
    newXdot=[dx,dy]
    return newX,newXdot

def eulerForEq2(X,Vwater,Xdot,h):
    x=X[0]
    y=X[1]
    newX=[x,y]#ny matrise for å ikke endre den gamle
    T=24*60*60
    factor=2*np.pi/T
    newX=[newX[0]+h*factor*(-y),newX[1]+h*factor*(x)]
    return newX

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


def rk2(X,Vwater,f,h,t): #X er en array som kan inneholde enten x,y eller x,y,vx,vy avhening av om hhv eq 2 eller eq 1 skal brukes
    K1=f(X,t)
    K1=np.array(K1)
    K2=f(X+h*K1,t+h)
    X_=X+h/2*(K1+K2)
    return X_

def errorAndTrajectoryForRK2(X,h,figname):
    newX=np.array(X)
    xList,yList=[
        X=[]
                ],[]
    xList.append(newX[0])
    X=[]
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
    if figname=="dummy":
        X=[]
        plt.figure(str(h)+"Trap")
        plt.title("Trajectory for particle using trapezoid, h = "+str(h))
        plt.plot(xList,yList)
    elif figname=="errorOnly":
        pass
    else:
        plt.figure(figname)
        plt.title("Tracectory for particle using trapezoid, h = "+str(h))
        plt.plot(xList,yList)
        plt.savefig(figname)
        plt.close()
    error=np.sqrt(vecError[0]**2+vecError[1]**2)
    return error

X=[]
def errorAndTrajectoryForEuler(X,h,figname):
    x,y=X[0],X[1]
    xList,yList=[],[]
    newX=[x,y] #oppretter en ny matrise for å unngå å endre på X. Må skje siden vi endrer på newX ila neste for-løkke
    xList.append(newX[0])
    yList.append(newX[1])
    timeFinal=2*24*60*60
    timeNow=0
    for i in range(int(2*24*60*60/h)+1):
        h=min(h,timeFinal-timeNow)
        timeNow+=h
        newX=eulerForEq2(newX,h)
        xList.append(newX[0])
        yList.append(newX[1])
    X=np.array(X)#gjør om til numpyarrays for å kunne gjøre vektoraritmetikk. Kunne i utgangspunktet tatt inn nparrays,
    newX=np.array(newX)#men funksjonen er nå kompatibel med både vanlige arrays og nparrays.
    vecError=X-newX #finner vektoren fra X til newX
    if figname=="dummy":
        plt.figure(str(h)+"Eul")
        plt.title("Trajectory for particle using Euler, h = "+str(h))
        plt.plot(xList,yList)
    elif figname=="errorOnly":
        pass
    else:
        plt.figure(figname)
        plt.title("Trajectory for particle using Euler, h = "+str(h))
        plt.plot(xList,yList)
        plt.savefig(figname)
        plt.close()
    error=np.sqrt(vecError[0]**2+vecError[1]**2)#finner avstanden mellom startpunktet og sluttpunktet
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

def task1a(bool=False):
    X=[100,0]
    h0=60*60 #sekund
    thiserror = errorAndTrajectoryForEuler(X,h0,"errorOnly")
    while thiserror>10:#finner høyeste tidssteg som gir feil mindre enn 10 meter
        h0-=1
        thiserror = errorAndTrajectoryForEuler(X,h0,"errorOnly")
    timestepEuler=h0
    print("tidssteget som kreves for å få en global feil mindre enn 10 m er h =",timestepEuler)
    if bool:
        X=[100,0]
        h=100
        errorAndTrajectoryForEuler(X,h,"1ahlik100.pdf")
        times=np.logspace(2,4,10)
        errorArr=[]
        for h in times:
            X=[100,0]
            errorArr.append(errorAndTrajectoryForEuler(X,h,"errorOnly"))
        plt.figure("1aerror")
        plt.title("Error for Euler")
        plt.loglog(times,errorArr,"ro")
        plt.axvline(x=timestepEuler)
        plt.savefig("1aerror.pdf")
    else:
        X=[100,0]
        h=100
        errorAndTrajectoryForEuler(X,h,"dummy")
        times=np.logspace(2,4,10)
        errorArr=[]
        for h in times:
            X=[100,0]
            errorArr.append(errorAndTrajectoryForEuler(X,h,"errorOnly"))
        plt.figure("1aerror")
        plt.title("Error for Euler")
        plt.loglog(times,errorArr,"ro")
        plt.axvline(x=timestepEuler)
        plt.show()

def task1b(bool=False):
    X=[100,0]
    h0=60*60 #sekund
    thiserror = errorAndTrajectoryForEuler(X,h0,"errorOnly")
    while thiserror>10:#finner høyeste tidssteg som gir feil mindre enn 10 meter
        h0-=1
        thiserror = errorAndTrajectoryForEuler(X,h0,"errorOnly")
    timestepEuler=h0
    h0=60*60 #sekund
    thiserror = errorAndTrajectoryForRK2(X,h0,"errorOnly")
    while thiserror>10:
        h0-=1
        thiserror = errorAndTrajectoryForRK2(X,h0,"errorOnly")
    timestepTrap=h0
    print("tidssteget som kreves for å få en global feil mindre enn 10 m med trapesmetoden er h =",timestepTrap)
    if bool:
        X=[100,0]
        h=100
        errorAndTrajectoryForRK2(X,h,"1ahlik100.pdf")
        times=np.logspace(2,4,10)
        errorArr1=[]
        errorArr2=[]
        for h in times:
            X=[100,0]
            errorArr1.append(errorAndTrajectoryForRK2(X,h,"errorOnly"))
            errorArr2.append(errorAndTrajectoryForEuler(X,h,"errorOnly"))
        plt.figure("1aerror")
        plt.title("Error for Euler and trapezoid with ")
        plt.loglog(times,errorArr1,"bo",label="Trapezoid")
        plt.loglog(times,errorArr2,"ro",label="Euler")
        plt.axvline(x=timestepEuler, color="r")
        plt.axvline(x=timestepTrap, color="b")
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig("1berror.pdf")
    else:
        X=[100,0]
        h=100
        errorAndTrajectoryForRK2(X,h,"dummy")
        times=np.logspace(2,4,10)
        errorArr1=[]
        errorArr2=[]
        for h in times:
            X=[100,0]
            errorArr1.append(errorAndTrajectoryForRK2(X,h,"errorOnly"))
            errorArr2.append(errorAndTrajectoryForEuler(X,h,"errorOnly"))
        plt.figure("1aerror")
        plt.title("Error for Euler and Trapezoid")
        plt.loglog(times,errorArr2,"ro",label="Trapezoid")
        plt.loglog(times,errorArr1,"bo",label="Euler")
        plt.axvline(x=timestepEuler, color=)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.show()

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