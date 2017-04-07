import math
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr
from scipy.interpolate import RectBivariateSpline
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pyproj
import random
import time
plt.style.use('bmh')

#Denne funksjonen gir oss vannhastigheten med bakgrunn i dataene
class Interpolator():
    def __init__(self, dataset):
        self.dataset = dataset

    def get_interpolators(self, X, it):
        # Add a buffer of cells around the extent of the particle cloud
        buf  = 3
        # Find extent of particle cloud in terms of indices
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

datapath = 'C:/Users/Even/Downloads/Norkyst-800m.nc'
d  = xr.open_dataset(datapath)
f  = Interpolator(dataset = d)
t = np.datetime64('2017-02-01T12:00:00')
X = np.array([-3000000, -1200000]).reshape(2, 1)  #reshape (2,Np)
#print("Hastighet i dette punktet ved denne tiden:\n",f(X, t))

t0   = np.datetime64('2017-02-01T12:00:00')
h  = np.timedelta64(3600, 's')
#Step forward
i = 2
t = t0 + i*h
#Get number of seconds in h
h_seconds = h / np.timedelta64(1, 's')

def Vwater(t,X):
    return f(X,t)

def ETMforEq2(X,t,h,Vwater):
    h_seconds = h / np.timedelta64(1, 's')
    Xvel=np.array(Vwater(t,X)).reshape(X.shape) #finner hastighet i nåværende punkt
    Xnext=X+h_seconds*Xvel #finner approksimert neste punkt
    XvelNext=np.array(Vwater(t+h,Xnext)).reshape(X.shape) #finner hastighet i approksimert neste punkt
    Xfinal=X+h_seconds/2*(Xvel+XvelNext)
    #print("Xfinal",Xfinal)
    return Xfinal

def particleTrajectory(X0, time_final, h, time_initial, velocityField, integrator):
    numberOfTimeSteps = int((time_final - time_initial) / h)
    X = np.zeros((numberOfTimeSteps + 1, *X0.shape))
    # X lagrer alle partiklenes posisjon gjennom hele tidsforløpet
    X[0, :] = X0
    time_now = time_initial
    for step in range(numberOfTimeSteps):  # Tok bort +1 pga hele timer
        time_now += h  # numpy håndterer time64-opplegg
        X[step + 1, :] = integrator(X[step, :], time_now, h, velocityField)
        # Denne bør virkelig returnere koordinater på formen (2,1) - og det gjør den
    return X

def task2a():
    clockStart = time.time()
    numberOfParticles = 1
    X0 = np.array([-3e6, -1.2e6]).reshape(2, numberOfParticles)  # reshape (2,Np)
    #Funker ikke med alt for store verdier
    print("X0: \n", X0)
    t0 = np.datetime64('2017-02-05T12:00:00')
    tEnd = np.datetime64('2017-02-15T12:00:00')
    h = np.timedelta64(3600, 's')
    trajectories = particleTrajectory(X0, tEnd, h, t0, Vwater, ETMforEq2)
    #print("trajectories før splitting:", trajectories[:5])
    trajectories = np.hsplit(trajectories, len(trajectories[0])) #Splitter i x og y
    xArray = trajectories[0]
    yArray = trajectories[1]
    #print("xArray in all its glory:\n", xArray)
    #print("Arrays hver for seg")
    #print(xArray[:5])
    #print(yArray[:5])
    #Her er det nye i denne fila:
    #print(xArray[0][0][0])
    #print(xArray[0][0][1])
    xArrayParticleSplit = np.array([ [xArray[i][0][0] for i in range(len(xArray))] ])
    yArrayParticleSplit = np.array([ [yArray[i][0][0] for i in range(len(yArray))] ])
    print("xArrayParticleSplit\n", xArrayParticleSplit[0][:5]) #Denne funker fint
    for particle in range(1,numberOfParticles): #Append smeller alt i samme brackets. Derfor vstack
        xArrayParticleSplit = np.vstack((xArrayParticleSplit, np.array([ [xArray[i][0][particle] for i in range(len(xArray))] ])))
        yArrayParticleSplit = np.vstack((yArrayParticleSplit, np.array([ [yArray[i][0][particle] for i in range(len(yArray))] ])))
    print(xArrayParticleSplit)
    plt.figure()
    plt.title("05. - 15. februar")
    for index in range(numberOfParticles): #endret fra len(xArrayParticleSplit)
        plt.plot(xArrayParticleSplit[index], yArrayParticleSplit[index])
    clockEnd = time.time()
    print("Tid brukt:", clockEnd - clockStart)
    print("Viser plott nå")
    plt.show()

def randomX0Array(lowendY,highendY,lowendX,highendX,numberOfParticles):
    finalArray=np.zeros(numberOfParticles*2)
    for i in range(numberOfParticles):
        finalArray[i] = random.randint(lowendY,highendY)
    for i in range(numberOfParticles):
        finalArray[numberOfParticles + i] = random.randint(lowendX,highendX)
    return finalArray

def frequencyCounter(X,Nx,Ny,numberOfParticles):
    cellSize = 800
    x = np.array((X[0, :]+3.01e6)/cellSize, dtype= int)
    y = np.array((X[1, :]+1.3e6) / cellSize, dtype= int)
    frequencyGrid = np.zeros((Ny, Nx))
    for i in range(numberOfParticles):
        frequencyGrid[y[i], x[i]] += 1
    return frequencyGrid

#task2a()
def task3b():
    Nx, Ny = 800, 300
    x = -3010000 + 800 * np.arange(Nx)
    y = -1300000 + 800 * np.arange(Ny)
    print("heyheyheyyy")
    x, y = np.meshgrid(x, y)
    print("ferdig med meshgrid")
    # Randomize startposisjoner og kalkuler sanns.
    # Lager en grid. Skal koordinatene i disse endres eller hva? Hvordan flytter vi partiklene?
    fig = plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', color='#dddddd')
    ax.add_feature(land_10m)
    ax.coastlines(resolution='10m')
    p1 = pyproj.Proj(d.projection_stere.proj4)
    p2 = pyproj.Proj(proj='latlong')
    lons, lats = pyproj.transform(p1, p2, x, y)
    print("ferdig med å preppe ax")
    numberOfParticles=1000
    X0 = randomX0Array(-3.001e6, -3e6, -1.201e6, -1.2e6, numberOfParticles).reshape(2, numberOfParticles)
    print("laget X0")
    t0 = np.datetime64('2017-02-01T12:00:00')
    tEnd = np.datetime64('2017-02-11T12:00:00')
    h = np.timedelta64(3600, 's')
    trajectories = particleTrajectory(X0, tEnd, h, t0, Vwater, ETMforEq2)
    print("halvveis ferdig med trajectories")
    trajectories = np.hsplit(trajectories, len(trajectories[0]))  # Splitter i x og y
    print("trajectories er fiksa")
    xArray = trajectories[0]
    yArray = trajectories[1]
    for day in range(0,6):
        print("dag",2*day)
        X0=np.array(np.concatenate((xArray[day*2*24],yArray[day*2*24]))).reshape(2,numberOfParticles)
        concentration = frequencyCounter(X0,Nx,Ny,numberOfParticles)
        concentration = np.ma.masked_array(concentration, mask=concentration == 0)
        # Her brukes C
        # C = (x + 2950000) ** 2 + (y + 1150000) ** 2
        # ax.pcolormesh(lons, lats, C, transform=ccrs.PlateCarree(), zorder=2)
        ax.pcolormesh(lons, lats, concentration, transform=ccrs.PlateCarree(), zorder=2)
    ax.set_extent((-5, 15, 57, 67))
    print("begynner å save")
    plt.savefig("figurgitter.pdf")
    print("ferdig")
    plt.show()


def task3b_new():
    Nx, Ny = 600, 300
    # Her settes ...
    x = -3010000 + 800 * np.arange(Nx)
    y = -1300000 + 800 * np.arange(Ny)
    x, y = np.meshgrid(x, y)
    fig = plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', color='#dddddd')
    ax.add_feature(land_10m)
    ax.coastlines(resolution='10m')
    p1 = pyproj.Proj(d.projection_stere.proj4)
    p2 = pyproj.Proj(proj='latlong')
    lons, lats = pyproj.transform(p1, p2, x, y)
    #print(lons)
    #print(lats)
    numberOfParticles=100000
    X0 = randomX0Array(-3.11e6, -3.1e6, -1.31e6, -1.3e6, numberOfParticles).reshape(2, numberOfParticles)
    print("laget X0")
    t0 = np.datetime64('2017-02-01T12:00:00')
    tEnd = np.datetime64('2017-02-11T12:00:00')
    h = np.timedelta64(3600, 's')
    trajectories = particleTrajectory(X0, tEnd, h, t0, Vwater, ETMforEq2)
    print("trajectories er fiksa")
    xArray = trajectories[:, 0, :]
    yArray = trajectories[:, 1, :]
    dataGrid_X = d.X.values
    dataGrid_Y = d.Y.values
    colormap = ['Reds', 'Oranges', 'Greens', 'Blues', 'Purples', 'PuRd'] #['autumn', 'cool', 'copper' , 'summer', 'winter', 'autumn']
    for day in range(0, 6):
        print("dag", 2*day)
        plt.figure(day)
        concentration, coordinates_X, coordinates_Y = np.histogram2d(xArray[day*2*24], yArray[day*2*24], bins=(dataGrid_X, dataGrid_Y))
        concentration = np.ma.masked_array(concentration, mask=concentration == 0)
        coordinates_X, coordinates_Y = np.meshgrid(coordinates_X, coordinates_Y)
        lons, lats = pyproj.transform(p1, p2, coordinates_X, coordinates_Y)
        ax.pcolormesh(lons, lats, concentration.T, transform=ccrs.PlateCarree(), zorder=2, cmap='gist_heat_r')
        plt.savefig("oppgave3b"+str(2*day)+"pdf")
    ax.set_extent((0, 7, 58, 64))
    #print("begynner å save")
    #plt.savefig("figurgitter.pdf")
    print("ferdig")
    plt.show()

startTime=time.time()
task3b_new()
