
# Observations:
# - Mass can be set per each particle
# - Modulus of elasticity can be set per each spring-like connection

import numpy as np

def initialize(): 
    global net_geometry,radius
    global see_elasticity,see_net,see_hmap,method
    global mk_movie
    global add_hydroForce,add_gravityForce
    global dt,tmax,N1,N2,netY,netZ,model,dw
    global rho,ux,uy,uz
    global l1,l2,l3,l4,L,dW,c,k1,k2,m
    net_geometry='cage'          # 'panel', 'cage'
    radius=0.5
    see_elasticity=0
    see_net=0
    see_hmap=0
    mk_movie=0
    add_hydroForce=1
    add_gravityForce=1
    method=7                # 1 forward Euler, 2 semi-implicit Euler, 3 backward (implicit) Euler,
                            # 4 predictor-corrector, 5 midpoint, 6 RK2, 7 RK4, 8 Verlet,
                            # 9 velocity Verlet, 10 Yoshida
    dt=5e-4                 # [s] time step
    tmax=10.0               # [s] maximum integration time
    c =0.0                  # viscous damping coefficient
    N1=12                   # number of particles along y (edges not included) | USE EVEN NUMBER
    N2=40                   # number of particles along z (edges not included) | USE EVEN NUMBER
    netY=1.0                # net size in y
    netZ=2.0*np.pi          # net size in z
    model='S1'              # model for drag and lift coefficients = S1, S4, none (set both to 1)
    dW=0.01                 # [m] twine diameter
    rho=999.06              # [kg/m3] fluid density
    ux,uy,uz=3.0,0.0,0.0    # [m/s] fluid velocity
    l1=netY/N1              # [m] half mesh size = L_twine + D_twine = interparticle distance along y
    l3=netZ/N2              # [m] half mesh size = L_twine + D_twine = interparticle distance along z
    L=np.mean((l1,l3))      # [m] half mesh size used to calc Sn in S1 model -> CHECK THIS
    k1=np.ones((N1+2,N2))*1e5   # [N/m2] spring constant along z
    k2=np.ones((N1,N2+2))*1e5   # [N/m2] spring constant along y
    k2=k2.transpose()
    m=np.ones((N1+2,N2+2))      # mass
