
import numpy as np
import timeit
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d
import os.path
from os import path
from scipy import interpolate,signal
import random 
from matplotlib.colors import LogNorm

fname='dt0p00001/'
x1=np.loadtxt(fname+'continua_x.out', delimiter = ',')
y1=np.loadtxt(fname+'continua_y.out', delimiter = ',')
z1=np.loadtxt(fname+'continua_z.out', delimiter = ',')
fname='dt0p00005/'
x2=np.loadtxt(fname+'continua_x.out', delimiter = ',')
y2=np.loadtxt(fname+'continua_y.out', delimiter = ',')
z2=np.loadtxt(fname+'continua_z.out', delimiter = ',')
fname='dt0p0001/'
x3=np.loadtxt(fname+'continua_x.out', delimiter = ',')
y3=np.loadtxt(fname+'continua_y.out', delimiter = ',')
z3=np.loadtxt(fname+'continua_z.out', delimiter = ',')
fname='dt0p0005/'
x4=np.loadtxt(fname+'continua_x.out', delimiter = ',')
y4=np.loadtxt(fname+'continua_y.out', delimiter = ',')
z4=np.loadtxt(fname+'continua_z.out', delimiter = ',')
fname='dt0p001/'
x5=np.loadtxt(fname+'continua_x.out', delimiter = ',')
y5=np.loadtxt(fname+'continua_y.out', delimiter = ',')
z5=np.loadtxt(fname+'continua_z.out', delimiter = ',')

def myerr(x1,y1,z1,x2,y2,z2):
    dx2=x2-x1
    dy2=y2-y1
    dz2=z2-z1
    err=np.sqrt(dx2*dx2+dy2*dy2+dz2*dz2)
    return sum(sum(err))

print(myerr(x1,y1,z1,x2,y2,z2))
print(myerr(x1,y1,z1,x3,y3,z3))
print(myerr(x1,y1,z1,x4,y4,z4))
print(myerr(x1,y1,z1,x5,y5,z5))
