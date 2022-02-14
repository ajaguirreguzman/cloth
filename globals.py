
# Observations:
# - Mass can be set per each particle
# - Modulus of elasticity can be set per each spring-like connection

import numpy as np

def initialize(): 
    global net_geometry,radius
    global see_elasticity,see_net,see_hmap,method
    global mk_movie
    global add_hydroForce,add_gravityForce
    global dt,tmax,N1,N2,netY,netZ,model,dw,theta
    global rho,ux,uy,uz
    global ly,lz,L,dW,c,k1,k2,m
    net_geometry=2
    radius=0.5
    see_elasticity=0
    see_net=0
    see_hmap=0
    mk_movie=0
    add_hydroForce=0
    add_gravityForce=0
    method=4                # 1 forward Euler, 2 semi-implicit Euler, 3 backward (implicit) Euler,
                            # 4 predictor-corrector, 5 midpoint, 6 RK2, 7 RK4, 8 Verlet,
                            # 9 velocity Verlet, 10 Yoshida
    dt=1e-4                 # [s] time step
    tmax=1.0                # [s] maximum integration time
    c =0.0                  # viscous damping coefficient
    N1=20                   # number of particles along y (edges not included) | USE EVEN NUMBER
    N2=20                   # number of particles along z (edges not included) | USE EVEN NUMBER
    netY=0.5                # net size in y
    netZ=2.0*np.pi          # net size in z
    model='S1'              # model for drag and lift coefficients = S1, S4, none (set both to 1)
    dW=0.01                 # [m] twine diameter
    theta=0.0               # = pi/2 - alpha, angle between plane net normal vector and flow
    rho=999.06              # [kg/m3] fluid density
    ux,uy,uz=3.0,0.0,0.0    # [m/s] fluid velocity
    ly=netY/N1              # [m] half mesh size = L_twine + D_twine = interparticle distance along y
    lz=netZ/N2              # [m] half mesh size = L_twine + D_twine = interparticle distance along z
    L=np.mean((ly,lz))      # [m] half mesh size used to calc Sn in S1 model -> CHECK THIS
    k1=np.ones((N1+2,N2))*1e5   # [N/m2] spring constant along z
    #k1[-1,:]=1e6
    k2=np.ones((N1,N2+2))*1e5   # [N/m2] spring constant along y
    k2=k2.transpose()
    m=np.ones((N1+2,N2+2))      # mass
