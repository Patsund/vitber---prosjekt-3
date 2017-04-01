import math
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr
from scipy.interpolate import RectBivariateSpline
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pyproj
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

datapath = 'C:/Users/Even/Downloads/NorKyst-800m.nc'
d  = xr.open_dataset(datapath)
f  = Interpolator(dataset = d)
t = np.datetime64('2017-02-01T12:00:00')
X = np.array([-3000000, -1300000]).reshape(2,1) #reshape (2,Np)
print("Hastighet i dette punktet ved denne tiden:",f(X, t))

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
    Xvel=Vwater(t,X) #finner hastighet i nåværende punkt
    Xvel=np.array(Xvel)
    Xnext=X+h_seconds*Xvel #finner approksimert neste punkt
    XvelNext=Vwater(t+h,Xnext) #finner hastighet i approksimert neste punkt
    XvelNext=np.array(XvelNext)
    Xfinal=X+h_seconds/2*(Xvel+XvelNext)
    #print("Xfinal",Xfinal)
    return Xfinal

def task2a():
    #POSITION
    X = np.array([-3000000, -1200000])
    coordinateArray = np.array([[0, 0]])
    coordinateArray = np.delete(coordinateArray, 0, 0)  # Fjerner 0,0 som ble brukt for å initialisere arrayet
    #TIME
    t0 = np.datetime64('2017-02-01T12:00:00')
    h = np.timedelta64(3600, 's')
    for i in range(24*10):
        X = ETMforEq2(X, t0 + h*i, h, Vwater)
        dummyarray = np.array([X])
        coordinateArray = np.concatenate((coordinateArray, dummyarray), axis=0)
    xValueArray = [k[0] for k in coordinateArray]
    yValueArray = [k[1] for k in coordinateArray]
    plt.figure()
    plt.title("Trajectory")
    plt.plot(xValueArray, yValueArray, 'ro', markersize=1)
    plt.show()

task2a()
########################################
#### Plotting trajectories on a map ####
########################################
