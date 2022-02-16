
# Included:
# - Forward (explicit) Euler (1st order)
# - Semi-implicit Euler (1st order)
# - Backward (implicit) Euler (1st order)
# - Midpoint method (2nd order)
# - RK2 (2nd order)
# - RK4 (4th order)
# - Verlet (2nd order)
# - velocity Verlet / leapfrog (2nd order, though a(t+dt) uses an Euler-advanced velocity which might affect order)
# - Yoshida integrator (4th order)
# To be implemented:
# - (two-step) Adams–Bashforth method
# - Runge-Kutta-Verner fifth-order and sixth-order method.
# Observations:
# - Vel. Verlet performs worse than Euler, possibly because a(x,v) instead of a(x). Here, a(t+dt) uses v(t+dt/2)
# Resources:
# - Various relevant Wikipedia links
# - https://www.phys.hawaii.edu/~varner/PHYS305-Spr12/P305_lab8.html
# - https://www.compadre.org/PICUP/resources/Numerical-Integration/

import params
from forces import a
from useful_scripts import newton_raphson
import numpy as np

def euler_step(x,y,z,fx,fy,fz,dt):
    xp=x+dt*fx
    yp=y+dt*fy
    zp=z+dt*fz
    return xp,yp,zp

def euler_method(x,y,z,vx,vy,vz):                       # standard forward (explicit) Euler method
    dt=params.dt
    ax ,ay ,az =a(x,y,z,vx,vy,vz)                       # acceleration at t_n
    vxp,vyp,vzp=euler_step(vx,vy,vz,ax,ay,az,dt)        # velocity at t_(n+1)
    xp ,yp ,zp =euler_step(x ,y ,z ,vx,vy,vz,dt)        # position at t_(n+1)
    return xp,yp,zp,vxp,vyp,vzp

def semiImplicitEuler_method(x,y,z,vx,vy,vz):
    dt=params.dt
    ax ,ay ,az =a(x,y,z,vx,vy,vz)                       # acceleration at t_n
    vxp,vyp,vzp=euler_step(vx,vy,vz,ax ,ay ,az ,dt)     # velocity at t_(n+1)
    xp ,yp ,zp =euler_step(x ,y ,z ,vxp,vyp,vzp,dt)     # notice that v(t+dt) is used instead of v(t)
    return xp,yp,zp,vxp,vyp,vzp

def newtonRaphson_func(xp,yp,zp,vx,vy,vz,vxg,vyg,vzg):  # function for Newton-Raphson algorithm
    dt=params.dt
    axg,ayg,azg=a(xp,yp,zp,vxg,vyg,vzg)                 # x(t+dt) and velocity guess
    fxg=vxg-vx-axg*dt
    fyg=vyg-vy-ayg*dt
    fzg=vzg-vz-azg*dt
    return fxg,fyg,fzg

def implicitEuler_method(x,y,z,vx,vy,vz):               # or backward Euler method
    dt=params.dt
    N1=params.N1
    N2=params.N2
    xp,yp,zp=euler_step(x,y,z,vx,vy,vz,dt)
    # velocity initial guess
    vxg=np.ones((N1+2,N2+2))
    vyg=np.ones((N1+2,N2+2))
    vzg=np.ones((N1+2,N2+2))
    # initialize
    dvg=1.0                                             # solution might be independent of dvg
    tol=1e-8
    i=1
    # Newton-Raphson loop
    while i<1000:
        # function(guess)
        fxg,fyg,fzg=newtonRaphson_func(xp,yp,zp,vx,vy,vz,vxg,vyg,vzg)
        fx1,fy1,fz1=newtonRaphson_func(xp,yp,zp,vx,vy,vz,vxg-dvg,vyg-dvg,vzg-dvg)
        fx2,fy2,fz2=newtonRaphson_func(xp,yp,zp,vx,vy,vz,vxg+dvg,vyg+dvg,vzg+dvg)
        # derivative at guess
        dfxg,dfyg,dfzg=(fx2-fx1)/2.0/dvg,(fy2-fy1)/2.0/dvg,(fz2-fz1)/2.0/dvg
        # tolerance
        if np.max(abs(fxg/dfxg))<tol: break
        if np.max(abs(fyg/dfyg))<tol: break
        if np.max(abs(fzg/dfzg))<tol: break
        # iterate
        vxg=vxg-fxg/dfxg
        vyg=vyg-fyg/dfyg
        vzg=vzg-fzg/dfzg
        i+=1
    vxp,vyp,vzp=vxg,vyg,vzg
    return xp,yp,zp,vxp,vyp,vzp

def predictorCorrector_method(x,y,z,vx,vy,vz):
    # Heun's method predicts using (explicit) Euler method and corrects using (implicit) trapezoidal rule
    dt=params.dt
    ncorrectors=1
    # predictor step
    ax ,ay ,az =a(x,y,z,vx,vy,vz)
    vxg,vyg,vzg=euler_step(vx,vy,vz,ax,ay,az,dt)
    xg ,yg ,zg =euler_step(x ,y ,z ,vx,vy,vz,dt)
    # corrector step: implicit but with predicted value
    for i in range(ncorrectors):
        #print('Corrector step #%i: x[12,12]=%.8f' % (i+1,x[12,12]))
        axg,ayg,azg=a(xg,yg,zg,vxg,vyg,vzg)
        axm=0.5*(axg+ax)
        aym=0.5*(ayg+ay)
        azm=0.5*(azg+az)
        vxm=0.5*(vxg+vx)
        vym=0.5*(vyg+vy)
        vzm=0.5*(vzg+vz)
        vxg,vyg,vzg=euler_step(vx,vy,vz,axm,aym,azm,dt)
        xg ,yg ,zg =euler_step(x ,y ,z ,vxm,vym,vzm,dt)
    return xg,yg,zg,vxg,vyg,vzg

def midpoint_method(x,y,z,vx,vy,vz):                    # order-2 for position, 1 for velocity
    dt=params.dt
    ax ,ay ,az =a(x,y,z,vx,vy,vz)                       # acceleration at t_n
    vxp,vyp,vzp=euler_step(vx,vy,vz,ax ,ay ,az ,dt)     # velocity at t_(n+1)
    vxm=0.5*(vxp+vx)
    vym=0.5*(vyp+vy)
    vzm=0.5*(vzp+vz)
    xp ,yp ,zp =euler_step(x ,y ,z ,vxm,vym,vzm,dt)     # ( v(t+dt) + v(t) ) / 2 is used!!
    return xp,yp,zp,vxp,vyp,vzp

def rk2_method(x,y,z,vx,vy,vz):
    dt=params.dt
    ax ,ay ,az =a(x,y,z,vx,vy,vz)                       # acceleration at t_n
    vxh,vyh,vzh=euler_step(vx,vy,vz,ax,ay,az,dt/2.0)    # velocity at t_(n+1/2)
    xh ,yh ,zh =euler_step(x ,y ,z ,vx,vy,vz,dt/2.0)    # position at t_(n+1/2)
    axh,ayh,azh=a(xh,yh,zh,vxh,vyh,vzh)                 # acceleration at t_(n+1/2)
    vxp,vyp,vzp=euler_step(vx,vy,vz,axh,ayh,azh,dt)     # velocity at t_(n+1)
    xp ,yp ,zp =euler_step(x ,y ,z ,vxh,vyh,vzh,dt)     # position at t_(n+1)
    return xp,yp,zp,vxp,vyp,vzp

def rk4_method(x1,y1,z1,vx1,vy1,vz1):
    dt=params.dt
    # k1
    ax1,ay1,az1=a(x1,y1,z1,vx1,vy1,vz1)
    # k2
    x2 ,y2 ,z2 =euler_step(x1 ,y1 ,z1 ,vx1,vy1,vz1,dt/2.0)
    vx2,vy2,vz2=euler_step(vx1,vy1,vz1,ax1,ay1,az1,dt/2.0)
    ax2,ay2,az2=a(x2,y2,z2,vx2,vy2,vz2)
    # k3
    x3 ,y3 ,z3 =euler_step(x1 ,y1 ,z1 ,vx2,vy2,vz2,dt/2.0)
    vx3,vy3,vz3=euler_step(vx1,vy1,vz1,ax2,ay2,az2,dt/2.0)
    ax3,ay3,az3=a(x3,y3,z3,vx3,vy3,vz3)
    # k4
    x4 ,y4 ,z4 =euler_step(x1 ,y1 ,z1 ,vx3,vy3,vz3,dt)
    vx4,vy4,vz4=euler_step(vx1,vy1,vz1,ax3,ay3,az3,dt)
    ax4,ay4,az4=a(x4,y4,z4,vx4,vy4,vz4)
    # all together
    vxp=vx1+(1.0/6.0)*dt*(ax1+2.0*ax2+2.0*ax3+ax4)
    vyp=vy1+(1.0/6.0)*dt*(ay1+2.0*ay2+2.0*ay3+ay4)
    vzp=vz1+(1.0/6.0)*dt*(az1+2.0*az2+2.0*az3+az4)
    xp =x1 +(1.0/6.0)*dt*(vx1+2.0*vx2+2.0*vx3+vx4)
    yp =y1 +(1.0/6.0)*dt*(vy1+2.0*vy2+2.0*vy3+vy4)
    zp =z1 +(1.0/6.0)*dt*(vz1+2.0*vz2+2.0*vz3+vz4)
    return xp,yp,zp,vxp,vyp,vzp

def verlet_method(x,y,z,xm,ym,zm,vx,vy,vz,i):
    # xm, ym and zm are old positions
    # first step is treated differently
    dt=params.dt
    ax,ay,az=a(x,y,z,vx,vy,vz)
    if i==0:
        xp=x+vx*dt+0.5*ax*dt*dt  # 2nd-degree polynomial
        yp=y+vy*dt+0.5*ay*dt*dt
        zp=z+vz*dt+0.5*az*dt*dt
        vxp=(xp-x)/dt            # backward difference
        vyp=(yp-y)/dt
        vzp=(zp-z)/dt
    else:
        xp=2.0*x-xm+ax*dt*dt
        yp=2.0*y-ym+ay*dt*dt
        zp=2.0*z-zm+az*dt*dt
        vxp=(xp-xm)/(2.0*dt)     # central difference; velocity is one step behind!
        vyp=(yp-ym)/(2.0*dt)
        vzp=(zp-zm)/(2.0*dt)
    return xp,yp,zp,vxp,vyp,vzp

def velocityVerlet_method(x,y,z,vx,vy,vz):          # aka leapfrog
    dt=params.dt
    # step 1
    ax,ay,az=a(x,y,z,vx,vy,vz)
    xp=x+vx*dt+0.5*ax*dt*dt
    yp=y+vx*dt+0.5*ay*dt*dt
    zp=z+vx*dt+0.5*az*dt*dt
    # step 2
    vxp,vyp,vzp=euler_step(vx,vy,vz,ax,ay,az,dt)    # use Euler to get v(t+dt), because a(x,v), not just a(x)
    axp,ayp,azp=a(xp,yp,zp,vxp,vyp,vzp)             # using an Euler-advanced velocity is a great venture
    # step 3
    axm=0.5*(axp+ax)
    aym=0.5*(ayp+ay)
    azm=0.5*(azp+az)
    vxp,vyp,vzp=euler_step(vx,vy,vz,axm,aym,azm,dt)
    return xp,yp,zp,vxp,vyp,vzp

def yoshida_method(x,y,z,vx,vy,vz):
    dt=params.dt
    # preliminaries
    cbrt2=2.0**(3.0/2.0)
    w0=-cbrt2/(2.0-cbrt2)
    w1=1.0/(2.0-cbrt2)
    c1=w1/2.0
    c2=(w0+w1)/2.0
    c3=c2
    c4=c1
    d1=w1
    d2=w0
    d3=d1
    # step 1
    x1 ,y1 ,z1 =euler_step(x ,y ,z ,vx,vy,vz,c1*dt)         # x(t+c1*dt)
    ax ,ay ,az =a(x ,y ,z ,vx ,vy ,vz )
    vx1,vy1,vz1=euler_step(vx,vy,vz,ax,ay,az,c1*dt)         # vx(t+c1*dt)
    ax1,ay1,az1=a(x1,y1,z1,vx1,vy1,vz1)                     # notice that ax(t+c1*dt) uses vx(t+c1*dt) !!!
    vx1,vy1,vz1=euler_step(vx,vy,vz,ax1,ay1,az1,d1*dt)      # vx(t+d1*dt)
    # step 2
    x2 ,y2 ,z2 =euler_step(x1,y1,z1,vx1,vy1,vz1,c2*dt)      # x(t+c2*dt)
    ax1,ay1,az1=a(x1,y1,z1,vx1,vy1,vz1)
    vx2,vy2,vz2=euler_step(vx1,vy1,vz1,ax1,ay1,az1,c2*dt)   # vx(t+c2*dt)
    ax2,ay2,az2=a(x2,y2,z2,vx2,vy2,vz2)                     # notice that ax(t+c2*dt) uses vx(t+c2*dt) !!!
    vx2,vy2,vz2=euler_step(vx1,vy1,vz1,ax2,ay2,az2,d2*dt)   # vx(t+d2*dt)
    # step 3
    x3 ,y3 ,z3 =euler_step(x2,y2,z2,vx2,vy2,vz2,c3*dt)      # x(t+c3*dt)
    ax2,ay2,az2=a(x2,y2,z2,vx2,vy2,vz2)
    vx3,vy3,vz3=euler_step(vx2,vy2,vz2,ax2,ay2,az2,c3*dt)   # vx(t+c3*dt)
    ax3,ay3,az3=a(x3,y3,z3,vx3,vy3,vz3)                     # notice that ax(t+c3*dt) uses vx(t+c3*dt) !!!
    vx3,vy3,vz3=euler_step(vx2,vy2,vz2,ax3,ay3,az3,d3*dt)   # vx(t+d3*dt)
    # step 4
    x4 ,y4 ,z4 =euler_step(x3,y3,z3,vx3,vy3,vz3,c4*dt)      # x(t+c4*dt)
    vx4,vy4,vz4=vx3,vy3,vz3
    return x4,y4,z4,vx4,vy4,vz4

