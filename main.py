
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
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d
#from scipy import interpolate,signal
#import random 
from matplotlib.colors import LogNorm

# my functions
import globals
from init_conditions import *
from iterative_methods import *
from energy import *
from plots import *
from useful_scripts import *
from io_continua import *
from animation_setup import *
from print_info import *


start = timeit.default_timer()

if __name__ == "__main__": 
    globals.initialize() 
    N1=globals.N1
    N2=globals.N2
    tmax=globals.tmax
    dt=globals.dt
    see_net=globals.see_net
    see_hmap=globals.see_hmap
    method=globals.method

    # print some info
    print_preInfo()

    # initial conditions
    x,y,z,globals.l1,globals.l2,globals.l3,globals.l4=set_positions()
    vx,vy,vz=set_velocities()
    x0,y0,z0=x,y,z

    # visualize elasticity
    #if globals.see_elasticity: plot_elasticity(globals.k1)
    if globals.see_elasticity: plot_elasticity(globals.m)

    # movie settings
    nframes=int(tmax/dt)

    # movie
    writer=setup_animation()
    
    # animation
    myxlim,myylim,myzlim=mylims()
    rx,ry,rz=figAspectRatios(myxlim,myylim,myzlim)
    fig1,axs1=setup_mainFigure(x,y,z)

    if see_net:
        plt.show()
        exit()

    # visualize heatmap
    if see_hmap: hmap=plot_hmap(reference_area(x,y,z)[1])

    # initial energy
    K0,V0=calc_energy(x,y,z,vx,vy,vz)

    # file to write
    file1=open('position.out','a')
    file2=open('energy.out','a')

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

            file1.write('%f,%f\n' % (t,z[0,1]-z0[0,1]))

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
                if i==0: print('using predictor-corrector method')
                xp,yp,zp,vx,vy,vz=predictorCorrector_method(x,y,z,vx,vy,vz)
            if method==5:
                if i==0: print('using midpoint method')
                xp,yp,zp,vx,vy,vz=midpoint_method(x,y,z,vx,vy,vz)
            if method==6:
                if i==0: print('using RK2 method')
                xp,yp,zp,vx,vy,vz=rk2_method(x,y,z,vx,vy,vz)
            if method==7:
                if i==0: print('using RK4 method')
                xp,yp,zp,vx,vy,vz=rk4_method(x,y,z,vx,vy,vz)
            if method==8:
                if i==0: print('using Verlet method')
                xp,yp,zp,vx,vy,vz=verlet_method(x,y,z,xm,ym,zm,vx,vy,vz,i)
            if method==9:
                if i==0: print('using velocity Verlet method')
                xp,yp,zp,vx,vy,vz=velocityVerlet_method(x,y,z,vx,vy,vz)
            if method==10:
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
            if max(np.sqrt(x*x+y*y+z*z).flatten())>200:
                print('Solution diverged!')
                #exit()
                break

            # update movie every 200 frames
            if i % 100 == 0:
                #net.set_data_3d(x.flatten(), y.flatten(), z.flatten())
                update_animation(axs1,x,y,z)
                if see_hmap: hmap.set_data(reference_area(x,y,z)[1])
                plt.pause(0.001)
                if globals.mk_movie: writer.grab_frame()

            maxChange=calc_posChange(x,y,z,x0,y0,z0)    # change in position
            K,V=calc_energy(x,y,z,vx,vy,vz)             # kinetic and potential energy
            file2.write('%f\t%f\t%f\n' % (K+V,K,V))

            if (i+1) % (0.1*nframes) == 0: print('%.0f%%' % (100*i/nframes))
            #print('%i \t %.6fs \t %.6fs \t %.6fs \t %.2f%%' % (nframes-i,t,time2-time1,timea,maxChange))

    file1.close()
    file2.close()

    write_continua(x,y,z,globals.l1,globals.l2,globals.l3,globals.l4)

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
