
import matplotlib.animation as animation
import numpy as np
from matplotlib import pyplot as plt
import params

def setup_animation(): 
    FFMpegWriter=animation.writers['ffmpeg']
    writer=FFMpegWriter(fps=1.0/params.dt/100)
    return writer

def mylims():
    myxlim=[-1.75,1.75]#][-0.25,2.0]
    myylim=[-0.20,1.70]#[-0.25,0.5]
    myzlim=[-0.20,1.95]#[-0.25,1.5]
    return myxlim,myylim,myzlim

def figAspectRatios(myxlim,myylim,myzlim):
    return np.ptp(myxlim), np.ptp(myylim), np.ptp(myzlim)

def setup_mainFigure(x,y,z):
    fig1=plt.figure(1)
    axs1=fig1.add_subplot(projection='3d')
    axs1.view_init(elev=0, azim=0, vertical_axis='y')
    axs1.plot_wireframe(x, y, z, rstride=1, cstride=1, color='gray')
    axs1.set_xlabel('x')
    axs1.set_ylabel('y')
    axs1.set_zlabel('z')
    update_animation(axs1,x,y,z)
    return fig1,axs1

def update_animation(axs1,x,y,z):
    axs1.cla()
    # update net
    myxlim,myylim,myzlim=mylims()
    axs1.set_xlim(myxlim)
    axs1.set_ylim(myylim)
    axs1.set_zlim(myzlim)
    rx,ry,rz=figAspectRatios(myxlim,myylim,myzlim)
    axs1.set_box_aspect((rz,rx,ry))                                     # aspect ratio = 1:1:1 in data space
    axs1.plot_wireframe(x, y, z, rstride=1, cstride=1, color='gray')
    # highlight heavier nodes
    mmax=np.max(params.m)
    axs1.plot(x[params.m==mmax], y[params.m==mmax], z[params.m==mmax], 'bo', ms=4)
