
import numpy as np
import globals
from io_continua import *
from useful_scripts import *
import os.path
from os import path

def mydistspace(start,stop,n,lmin):
    lmax=(stop-start)/2.0
    lmin=lmin/2.0
    x1=2.0            # use log curve from x1 to x2
    x2=5.0
    logx1=np.log10(x1)
    logx2=np.log10(x2)
    c1=(lmax-lmin)/(x2-x1)
    c2=(lmax*x1-lmin*x2)/(x1-x2)
    xf=c1*np.logspace(logx1, logx2, num=int(n/2), endpoint=True)+c2
    xb=-xf[::-1]    # reverse
    x=np.concatenate((xb,xf), axis=0)+lmax+start
    return x

def myinterpolation(x):
    len1,len2=np.shape(x)                               # setup original indices
    lin1=np.linspace(0, 1, len1)
    lin2=np.linspace(0, 1, len2)
    x=interpolate.interp2d(lin2, lin1, x, kind='cubic') # interpolation function | forgot why indices are swapped
    lin1=np.linspace(0, 1, N1+2)                        # setup new indices
    lin2=np.linspace(0, 1, N2+2)
    x=x(lin2,lin1)                                      # interpolate | forgot why indices are swapped
    return  x

def read_previous():
    N1=globals.N1
    N2=globals.N2
    x,y,z,l1,l2,l3,l4=read_continua()
    if np.shape(x)!=(N1+2,N2+2):
        len1,len2=np.shape(x)
        print('Interpolating continua data from %ix%i to %ix%i' % (len1,len2,N1+2,N2+2))
        x=myinterpolation(x)
        y=myinterpolation(y)
        z=myinterpolation(z)
        l1=myinterpolation(l1)
        l2=myinterpolation(l2)
        l3=myinterpolation(l3)
        l4=myinterpolation(l4)
    else:
        print('Continua data have the same shape')
    return x,y,z,l1,l2,l3,l4

def create_positions():
    N1=globals.N1
    N2=globals.N2
    netY=globals.netY
    netZ=globals.netZ
    print('Using default initial conditions...')
    if globals.net_geometry==1: # plane net
        x=np.zeros((N1+2,N2+2))
        #lin1=(netY/N1)*mydistspace(0, N1+1, N1+2,0.005)
        #lin2=(netZ/N2)*mydistspace(0, N2+1, N2+2,0.001)
        lin1=(netY/N1)*np.linspace(0, N1+1, N1+2)
        lin2=(netZ/N2)*np.linspace(0, N2+1, N2+2)
        z,y=np.meshgrid(lin2,lin1)
    if globals.net_geometry==2: # net cage
        lin1=(netY/N1)*np.linspace(0, N1+1, N1+2)
        lin2=np.linspace(0, N2+1, N2+2)/(N2+1)
        linx=globals.radius*np.sin(2.0*np.pi*lin2)
        linz=globals.radius*(1.0-np.cos(2.0*np.pi*lin2))
        x,na=np.meshgrid(linx,lin1)
        y,na=np.meshgrid(lin1,linx,indexing='ij')
        z,na=np.meshgrid(linz,lin1)
    #ly,lz=y[1:,:]-y[:-1,:],z[:,1:]-z[:,:-1]
    xyz1,xyz2,xyz3,xyz4=neighbor_info(x,y,z)    # xyz3,xyz4 are transposed for convenience
    dx1,dy1,dz1=relative_var(xyz1)
    dx2,dy2,dz2=relative_var(xyz2)
    dx3,dy3,dz3=relative_var(xyz3)
    dx4,dy4,dz4=relative_var(xyz4)
    l1=np.sqrt(dx1*dx1+dy1*dy1+dz1*dz1)         # relaxed distance for spring 1
    l2=np.sqrt(dx2*dx2+dy2*dy2+dz2*dz2)         # relaxed distance for spring 2
    l3=np.sqrt(dx3*dx3+dy3*dy3+dz3*dz3)         # relaxed distance for spring 3
    l4=np.sqrt(dx4*dx4+dy4*dy4+dz4*dz4)         # relaxed distance for spring 4
    return x,y,z,l1,l2,l3,l4

def set_positions():
    tmp=path.exists('continua_x.out')
    if tmp: x,y,z,l1,l2,l3,l4=read_previous()
    else:   x,y,z,l1,l2,l3,l4=create_positions()
    return x,y,z,l1,l2,l3,l4

def set_velocities():
    N1=globals.N1
    N2=globals.N2
    vx=np.zeros((N1+2,N2+2))
    vy=np.zeros((N1+2,N2+2))
    vz=np.zeros((N1+2,N2+2))
    middle1=int(N1/2.0)
    #middle2=int(N2/2.0)
    #vx[middle1:middle+2,middle2:middle2+2]=5.0
    #vx[0,:]=5.0
    return vx,vy,vz
