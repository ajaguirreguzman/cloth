
import matplotlib.animation as animation
import numpy as np
from matplotlib import pyplot as plt
import globals

def setup_animation(): 
    FFMpegWriter=animation.writers['ffmpeg']
    writer=FFMpegWriter(fps=0.5/globals.dt)
    return writer

def mylims():
    myxlim=[-.2,1.2]#][-0.25,2.0]
    myylim=[-.2,1.2]#[-0.25,0.5]
    myzlim=[-.2,1.2]#[-0.25,1.5]
    return myxlim,myylim,myzlim

def figAspectRatios(myxlim,myylim,myzlim):
    return np.ptp(myxlim), np.ptp(myylim), np.ptp(myzlim)

def setup_mainFigure(x,y,z):
    fig1=plt.figure(1)
    axs1=fig1.add_subplot(projection='3d')
    axs1.view_init(elev=10, azim=-70, vertical_axis='y')
    #net,=axs1.plot(x.flatten(), y.flatten(), z.flatten(), 'ro', markersize=1)
    axs1.plot_wireframe(x, y, z, rstride=1, cstride=1, color='gray')
    axs1.set_xlabel('x')
    axs1.set_ylabel('y')
    axs1.set_zlabel('z')
    myxlim,myylim,myzlim=mylims()
    axs1.set_xlim(myxlim)
    axs1.set_ylim(myylim)
    axs1.set_zlim(myzlim)
    rx,ry,rz=figAspectRatios(myxlim,myylim,myzlim)
    axs1.set_box_aspect((rz,rx,ry))  # aspect ratio is 1:1:1 in data space
    #plt.show()
    return fig1,axs1
