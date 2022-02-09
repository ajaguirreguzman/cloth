
# Included:
# - Last state (x,y,z,ly,lz) is saved to "continua" files. If they exist, these files are read
# - If N1,N2 are different, x,y,z,ly,lz are interpolated.
#
# Observations:
# - For validation (static analysis): consider a single line of springs subject to gravity and compare the resulting shape and (static) tension with the catenary solution. See Hou et al. 2018 PMR. Also, replicate the dynamic analysis therein.

# strain = deformation
# stress = applied forces

import numpy as np
import timeit
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d
from scipy import interpolate,signal
import random 
from matplotlib.colors import LogNorm

# my functions
import globals
from init_conditions import *
from iterative_methods import *
from energy import *
from plots import *
from useful_scripts import *
from io_continua import *


start = timeit.default_timer()

if __name__ == "__main__": 
    globals.initialize() 
    N1=globals.N1
    N2=globals.N2
    tmax=globals.tmax
    dt=globals.dt
    see_elasticity=globals.see_elasticity
    see_net=globals.see_net
    see_hmap=globals.see_hmap
    method=globals.method

    # initial conditions
    x,y,z,globals.ly,globals.lz=set_positions()
    vx,vy,vz=set_velocities()
    x0,y0,z0=x,y,z

    # visualize elasticity
    if see_elasticity: plot_elasticity(globals.k1)

    # movie settings
    nframes=int(tmax/dt)

    # movie setup
    FFMpegWriter=animation.writers['ffmpeg']
    writer=FFMpegWriter(fps=10)
    # animation
    fig1=plt.figure(1)
    axs1=fig1.add_subplot(projection='3d')
    axs1.view_init(elev=10, azim=-70, vertical_axis='y')
    #net,=axs1.plot(x.flatten(), y.flatten(), z.flatten(), 'ro', markersize=1)
    axs1.plot_wireframe(x, y, z, rstride=1, cstride=1, color='gray')
    axs1.set_xlabel('x')
    axs1.set_ylabel('y')
    axs1.set_zlabel('z')
    myxlim=[-.2,1.2]#][-0.25,2.0]
    myylim=[-.2,1.2]#[-0.25,0.5]
    myzlim=[-.2,1.2]#[-0.25,1.5]
    axs1.set_xlim(myxlim)
    axs1.set_ylim(myylim)
    axs1.set_zlim(myzlim)
    rx,ry,rz=np.ptp(myxlim), np.ptp(myylim), np.ptp(myzlim)
    axs1.set_box_aspect((rz,rx,ry))  # aspect ratio is 1:1:1 in data space
    #plt.show()

    if see_net:
        plt.show()
        exit()

    # visualize heatmap
    if see_hmap: hmap=plot_hmap(reference_area(x,y,z)[1])

    # initial energy
    K0,V0=calc_energy(x,y,z,vx,vy,vz)

    # file to write
    file=open('energy.out','a')

    # initialize
    t=0.0
    timet=0.0
    dmax=0.0
    xm=0.0
    ym=0.0
    zm=0.0

    with writer.saving(fig1, "movie.mp4", dpi=200):
        # main loop
        print('Time loop... ', end='')
        for i in range(nframes):

            time1 = timeit.default_timer()

            if method==1:
                if i==0: print('using forward (explicit) Euler method')
                xp,yp,zp,vx,vy,vz=euler_method(x,y,z,vx,vy,vz)
            if method==2:
                if i==0: print('using semi-implicit Euler method')
                xp,yp,zp,vx,vy,vz=semiImplicitEuler_method(x,y,z,vx,vy,vz)
            if method==3:
                if i==0: print('using backward (implicit) Euler method')
                xp,yp,zp,vx,vy,vz=implicitEuler_method(x,y,z,vx,vy,vz)
            if method==4:
                if i==0: print('using midpoint method')
                xp,yp,zp,vx,vy,vz=midpoint_method(x,y,z,vx,vy,vz)
            if method==5:
                if i==0: print('using RK2 method')
                xp,yp,zp,vx,vy,vz=rk2_method(x,y,z,vx,vy,vz)
            if method==6:
                if i==0: print('using RK4 method')
                xp,yp,zp,vx,vy,vz=rk4_method(x,y,z,vx,vy,vz)
            if method==7:
                if i==0: print('using Verlet method')
                xp,yp,zp,vx,vy,vz=verlet_method(x,y,z,xm,ym,zm,vx,vy,vz,i)
            if method==8:
                if i==0: print('using velocity Verlet method')
                xp,yp,zp,vx,vy,vz=velocityVerlet_method(x,y,z,vx,vy,vz)
            if method==9:
                if i==0: print('using Yoshida method')
                xp,yp,zp,vx,vy,vz=yoshida_method(x,y,z,vx,vy,vz)

            # iteration info
            time2 =timeit.default_timer()
            timet+=(time2-time1)
            timea =timet/(i+1)

            # update (no need to store two levels of velocity) 
            xm,ym,zm=x,y,z  # old
            x,y,z=xp,yp,zp  # current
            t+=dt

            # stop if diverge
            if max(np.sqrt(x*x+y*y+z*z).flatten())>100:
                print('Solution diverged!')
                break

            # update movie every 200 frames
            if i % 10 == 0:
                #net.set_data_3d(x.flatten(), y.flatten(), z.flatten())
                axs1.cla()
                axs1.set_xlim(myxlim)
                axs1.set_ylim(myylim)
                axs1.set_zlim(myzlim)
                axs1.set_box_aspect((rz,rx,ry))  # aspect ratio is 1:1:1 in data space
                axs1.plot_wireframe(x, y, z, rstride=1, cstride=1, color='gray')
                if see_hmap: hmap.set_data(reference_area(x,y,z)[1])
                plt.pause(0.01)
                if globals.mk_movie: writer.grab_frame()

            maxChange=calc_posChange(x,y,z,x0,y0,z0)    # change in position
            K,V=calc_energy(x,y,z,vx,vy,vz)             # kinetic and potential energy
            file.write('%f\t%f\t%f\n' % (K+V,K,V))

            # print iteration countdown, physical time, iter time, average iter time
            #print('%i \t %.6fs \t %.6fs \t %.6fs \t %.2f%%' % (nframes-i,t,time2-time1,timea,maxChange))

    file.close()

    write_continua(x,y,z,globals.ly,globals.lz)

    damping_ratio=globals.c/2.0/np.sqrt(globals.k1[0,0])    # damping ratio
    Ei=K0+V0                                # initial energy
    Ef=K+V                                  # final energy
    Ee=100.0*(Ef-Ei)/Ei                     # percentage error of final energy to initial energy
    stop = timeit.default_timer()           # final time

    # info
    print('Damping ratio     = %.4f   ' % damping_ratio)
    print('Error (Ef-Ei)/Ei  = %.6f %%' % Ee)
    print('Average iter time = %.6f s'  % timea)
    print('Total run time    = %.2f s'  % (stop - start))

    plt.show()
