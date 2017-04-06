import math
import numpy as np
from matplotlib import pyplot as plt
import time



def Vwater(t,X):
    x,y=X[0],X[1]
    T=24*60*60 #periode
    factor=2*math.pi/T
    Y=[factor*(-X[1]),factor*X[0]] #returnerer hastigheten basert på koordinater i henhold til Vw fra oppgave 1
    Y=np.array(Y)
    return Y

def eulerForEq2(X,t): #Beregner w(i+1) med Eulers metode
    x=X[0]
    y=X[1]
    Y=[x,y]#ny matrise for å ikke endre den gamle
    T=24*60*60
    factor=2*math.pi/T
    Y=[Y[0]+t*factor*(-y),Y[1]+t*factor*(x)]
    x=Y[0]
    y=Y[1]
    return Y

def ETMforEq2(X,t,h,slopeFunction):
    Xvel=slopeFunction(t,X) #finner hastighet i nåværende punkt
    Xvel=np.array(Xvel)
    Xnext=X+h*Xvel #finner approksimert neste punkt
    XvelNext=slopeFunction(t+h,Xnext) #finner hastighet i approksimert neste punkt
    XvelNext=np.array(XvelNext)
    Xfinal=X+h/2*(Xvel+XvelNext)
    #print("Xfinal",Xfinal)
    return Xfinal


def genLogSpace( array_size, num ):
    lspace = np.around(np.logspace(1,np.log10(array_size),num)).astype(np.uint64)
    return np.array(sorted(set(lspace.tolist())))-1

#vector=[0,100]
#vector2=[0,100]
#tid=60
### DETTE MÅ FIKSES I HENHOLD TIL NÅVÆRENDE FUNKSJONER ###
# for i in range(24*60*60//tid):
#     vector=ETM(vector,tid)
# for i in range(24*60*60//tid):
#     vector2=eulerForEq2(vector2,tid)
# print(vector)
# print(vector2)

def errorForEuler(X,h,figname):
    #print(X) #debugging
    x,y=X[0],X[1]
    xList,yList=[],[]
    Y=[x,y] #oppretter en ny matrise for å unngå å endre på X. Må skje siden vi endrer på Y ila neste for-løkke
    xList.append(Y[0])
    yList.append(Y[1])
    timeFinal=2*24*60*60
    timeNow=0
    for i in range(int(2*24*60*60/h)+1):
        h=min(h,timeFinal-timeNow)
        timeNow+=h
        Y=eulerForEq2(Y,h)
        xList.append(Y[0])
        yList.append(Y[1])
    X=np.array(X)#gjør om til numpyarrays for å kunne gjøre vektoraritmetikk. Kunne i utgangspunktet tatt inn nparrays,
    Y=np.array(Y)#men funksjonen er nå kompatibel med både vanlige arrays og nparrays.
    vecError=X-Y #finner vektoren fra X til Y
    if figname!="dummy.pdf":
        plt.figure(figname)
        plt.plot(xList,yList)
        plt.savefig(figname)
        plt.close()
    error=np.sqrt(vecError[0]**2+vecError[1]**2)#finner avstanden mellom startpunktet og sluttpunktet
    return error

def errorForETM(X,h,figname):
    x,y=X[0],X[1]
    xList,yList=[],[]
    timeNow=0
    finalTime=2*24*60*60
    n=int(finalTime/h)
    newX=np.array([x,y])
    xList.append(newX[0])
    yList.append(newX[1])
    for i in range(n+1):
        h=min(h,finalTime-timeNow) #setter h til å være gjenværende tid dersom gjenværende tid er mindre enn h
        newX=ETMforEq2(newX,timeNow,h,Vwater)
        xList.append(newX[0])
        yList.append(newX[1])
        timeNow+=h
    X=np.array(X)
    newX=np.array(newX)
    vecError=X-newX
    #print(vecError)
    if figname!="dummy.pdf":
        plt.figure(figname)
        plt.plot(xList,yList)
        plt.savefig(figname)
        plt.close()
    error=np.sqrt(vecError[0]**2+vecError[1]**2)
    return error

def plotTrajectory(X,figname):
    plt.figure(figname)
    plt.plot(X)
times=np.logspace(2,4,10)
errorArr1=[]
errorArr2=[]
X=[0,100]
for h in times: #lager en array med avstand mellom start- og sluttpunkt i henhold til økende tidssteg
    newname="pdfer/figure"
    newname+=str(h)
    eulname=newname+"eul.pdf"
    etmname=newname+"etm.pdf"
    errorArr2.append(errorForETM(X,h,etmname))
    errorArr1.append(errorForEuler(X,h,eulname))


h0=24*60*60 #sekund
print("Høyeste error Euler",errorArr1[-4])
print("Høyeste error ETM",errorArr2[-4])
thiserror = errorForEuler(X,h0,"dummy.pdf")
while thiserror>1:#finner høyeste tidssteg som gir feil mindre enn 10 meter
    h0-=1
    thiserror = errorForEuler(X,h0,"dummy.pdf")
    #print(thiserror)
tidsstegEuler=h0
print("for å akkurat få en error akkurat enn 10 meter med preisjon på 1 sekund må tidssteget være",h0,"sekunder")
h0=60*60
thiserror = errorForETM(X,h0,"dummy.pdf")
print("første verdi",thiserror)
while thiserror>1:
    h0-=1
    thiserror = errorForETM(X,h0,"dummy.pdf")
    #print(thiserror)
tidsstegETM=h0
print("for å akkurat få en error mindre enn 10 meter med presisjon på 1 sekund med ETM må tidssteget være",h0,"sekunder")
plt.figure()
plt.title("Global error")
plt.loglog(times,errorArr1,"ro")
plt.loglog(times,errorArr2,"bo")
plt.xlabel("tid / sekunder")
plt.ylabel("error / meter")
plt.savefig("pdfer/globalError.pdf")
X=(100,0)
n=30
endTimeEuler,endTimeETM=0,0
for i in range(n):
    startTimeEuler=time.time()
    errorEuler=errorForEuler(X,tidsstegEuler,"dummy.pdf")
    endTimeEuler+=time.time()-startTimeEuler
    startTimeETM=time.time()
    errorETM=errorForETM(X,tidsstegETM,"dummy.pdf")
    endTimeETM+=time.time()-startTimeEuler
endTimeEuler/=n
endTimeETM/=n
print("ETMtid/Eulertid",endTimeETM/endTimeEuler)
print("Error for Euler",errorEuler,"brukte",endTimeEuler,"sekunder \nError for ETM",errorETM,"brukte",endTimeETM,"sekunder.")