import math
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import xarray as xr
from scipy.interpolate import RectBivariateSpline
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pyproj
import random
import time
plt.style.use('bmh')
startTime=time.time()

#Denne funksjonen gir oss vannhastigheten med bakgrunn i dataene
class Interpolator():
    def __init__(self, dataset):
        self.dataset = dataset

    def get_interpolators(self, X, it):
        # Add a buffer of cells around the extent of the particle butt
        buf  = 3
        # Find extent of particle butt in terms of indices
        imax = np.searchsorted(self.dataset.X, np.amax(X[0,:])) + buf
        imin = np.searchsorted(self.dataset.X, np.amin(X[0,:])) - buf
        jmax = np.searchsorted(self.dataset.Y, np.amax(X[1,:])) + buf
        jmin = np.searchsorted(self.dataset.Y, np.amin(X[1,:])) - buf
        # Take out subset of array, to pass to
        # interpolation object
        # Fill NaN values (land cells) with 0, otherwise
        # interpolation won't work
        u    = self.dataset.u[it, 0, jmin:jmax, imin:imax].fillna(0.0)
        v    = self.dataset.v[it, 0, jmin:jmax, imin:imax].fillna(0.0)
        # RectBivariateSpline essentially returns a function,
        # which can be called to get value at arbitrary position
        # kx and ky sets order of spline interpolation along either direction (must be 1 <= kx <= 5)
        # transpose arrays to switch order of coordinates
        fu   = RectBivariateSpline(self.dataset.X[imin:imax], self.dataset.Y[jmin:jmax], u.T)#, kx = 3, ky = 3)
        fv   = RectBivariateSpline(self.dataset.X[imin:imax], self.dataset.Y[jmin:jmax], v.T)#, kx = 3, ky = 3)
        return fu, fv

    def get_time_index(self, t):
        # Get index of largest timestamp smaller than (or equal to) t
        return np.searchsorted(self.dataset.time, t, side='right') - 1

    def __call__(self, X, t):
        # get index of current time in dataset
        it = self.get_time_index(t)
        # get interpolating functions,
        # covering the extent of the particle
        fu, fv = self.get_interpolators(X, it)
        # Evaluate velocity at position(x[:], y[:])
        dx = fu(X[0,:], X[1,:], grid = False)
        dy = fv(X[0,:], X[1,:], grid = False)
        return np.array([dx, dy])

datapath = 'C:/Users/Patrik/Downloads/NorKyst-800m.nc' #her må egen datapath skrives inn
d  = xr.open_dataset(datapath)
f  = Interpolator(dataset = d)
t = np.datetime64('2017-02-01T12:00:00')
X = np.array([-3000000, -1200000]).reshape(2, 1)  #reshape (2,Np)

t0   = np.datetime64('2017-02-01T12:00:00')
h  = np.timedelta64(3600, 's')
#Step forward
i = 2
t = t0 + i*h
#Get number of seconds in h
h_seconds = h / np.timedelta64(1, 's')

def Vwater(t,X):
    return f(X,t)

def slopeFunctionForEq2(X,t):
    dX=np.array(Vwater(t,X))
    return dX


def rk2(X,Vwater,slopeFunction,h,t): #X er en array som kan inneholde enten x,y eller x,y,vx,vy avhening av om hhv eq 2 eller eq 1 skal brukes
    K1=slopeFunction(X,t)
    K1=np.array(K1)
    K2=slopeFunction(X+h*K1,t+h)
    X_=X+h/2*(K1+K2)
    return X_

def particleTrajectory(X0, time_final, h, time_initial, velocityField, integrator):
    numberOfTimeSteps = int((time_final - time_initial) / h)
    X = np.zeros((numberOfTimeSteps + 1, *X0.shape))
    # X lagrer alle partiklenes posisjon gjennom hele tidsforløpet
    X[0, :] = X0
    time_now = time_initial
    for step in range(numberOfTimeSteps):  # Tok bort +1 pga hele timer
        X[step + 1, :] = integrator(X[step, :],velocityField, slopeFunctionForEq2, h, time_now )
        time_now += h
    return X

def randomX0Array(lowendY,highendY,lowendX,highendX,numberOfParticles):
    finalArray=[]
    for i in range(numberOfParticles):
        finalArray.append(random.randint(lowendY,highendY))
    for i in range(numberOfParticles):
        finalArray.append(random.randint(lowendX,highendX))
    return finalArray
def task2a():
    numberOfParticles = 1
    initArray=randomX0Array(-1.21e6,-1.19e6,-3.01e6,-2.99e6,numberOfParticles)
    print(len(initArray))
    X0 = np.array(initArray).reshape(2, numberOfParticles)  # reshape (2,Np)
    t0 = np.datetime64('2017-02-01T12:00:00')
    tEnd = np.datetime64('2017-02-11T12:00:00')
    h = np.timedelta64(3600, 's')
    trajectories = particleTrajectory(X0, tEnd, h, t0, Vwater, rk2)
    trajectories = np.hsplit(trajectories, len(trajectories[0]))
    xArray = trajectories[0]
    yArray = trajectories[1]
    xArrayParticleSplit = np.array([ [xArray[i][0][0] for i in range(len(xArray))] ])
    yArrayParticleSplit = np.array([ [yArray[i][0][0] for i in range(len(yArray))] ])
    print("xArrayParticleSplit\n", xArrayParticleSplit[0][:5])
    for particle in range(1,numberOfParticles): #Append smeller alt i samme brackets. Prøver vstack
        xArrayParticleSplit = np.vstack((xArrayParticleSplit, np.array([ [xArray[i][0][particle] for i in range(len(xArray))] ])))
        yArrayParticleSplit = np.vstack((yArrayParticleSplit, np.array([ [yArray[i][0][particle] for i in range(len(yArray))] ])))
    plt.figure()
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    land_10m = cfeature.NaturalEarthFeature('physical','land','10m',color='#00aa00')
    ax.add_feature(land_10m)
    ax.coastlines(resolution='10m')
    p1 = pyproj.Proj(d.projection_stere.proj4)
    p2 = pyproj.Proj(proj='latlong')
    plt.title("Partikkelens bane")
    ax.set_extent((-4, 15, 57, 67))
    for index in range(numberOfParticles):
        lons, lats = pyproj.transform(p1, p2, xArrayParticleSplit[index], yArrayParticleSplit[index])
        ax.plot(lons, lats,".", transform=ccrs.PlateCarree(), zorder=2)
    print("Viser plott nå")
    endTime=time.time()
    print("tid brukt",endTime-startTime)
    plt.show()