
# Testing explicit and implicit methods
# https://en.wikipedia.org/wiki/Explicit_and_implicit_methods
# https://beltoforion.de/en/harmonic_oscillator/
# https://galileo.phys.virginia.edu/classes/152.mf1i.spring02/Oscillations4.htm

import numpy as np
import matplotlib.pyplot as plt
import timeit

def euler_step(x,v):
    return x+v*dt

def homogeneous_solution(c,k):
    delta=c/2.0
    omega0=np.sqrt(k)
    omega=np.sqrt(omega0**2.0-delta**2.0)
    bigSqrt=np.sqrt(2.0*delta*v0*x0+v0**2.0+x0**2.0*(delta**2.0+omega**2.0))
    A=-bigSqrt/2.0/omega
    phi=2.0*np.pi-2.0*np.arctan( (v0+delta*x0) / (x0*omega-bigSqrt) )
    return np.exp(-delta*t)*(2.0*A*np.cos(phi+omega*t))

def inhomogeneous_solution(m,b,k):
    omega0=np.sqrt(k/m)
    phi=np.arctan((b/m)*omegaF/(omegaF**2.0-omega0**2.0))
    x0=(F0/m)/np.sqrt((b/m)**2.0*omegaF**2 + (omega0**2-omegaF**2)**2.0)
    return x0*np.sin(omegaF*t+phi)

def a(x,v,t):
    F0=0.3
    omega=0.5
    return -k*x-c*v+F0*np.sin(omega*t)

def fxANDdfx(x,f,xg):
    x1=x[x<xg][-1]
    x2=x[x>xg][0]
    f1=f[x<xg][-1]
    f2=f[x>xg][0]
    m=(f2-f1)/(x2-x1)
    b=f1-m*x1
    return m*xg+b,m

def newton_raphson(x,f,xg): # xg is initial guess
    tol=1e-4
    i=1
    dx=8.0
    while i<10000:
        print(i)
        #fxg,dfxg=fxANDdfx(x,f,xg)
        fxg,na=fxANDdfx(x,f,xg)
        fx1,na=fxANDdfx(x,f,xg-dx)
        fx2,na=fxANDdfx(x,f,xg+dx)
        dfxg=(fx2-fx1)/2.0/dx
        #print('%i \t %.8f \t %.8f' % (i,xg,fxg))
        #if abs(fxg)<tol: break
        if abs(fxg/dfxg)<tol: break
        xg=xg-fxg/dfxg
        i+=1
    return xg


def semiImplicitEuler_method(x,v,t):
    file = open('data_semiImplicit.out', 'a')
    start_semiImplicit = timeit.default_timer()
    for i in range(int(tmax/dt)):
        file.write('%f\t%f\n' % (t,x))
        vp=euler_step(v,a(x,v,t))
        xp=euler_step(x,vp)
        x=xp
        v=vp
        t+=dt
    stop_semiImplicit = timeit.default_timer()
    file.close()
    print('Tot time semi-implicit  = %f' % (stop_semiImplicit-start_semiImplicit))

def implicitEuler_method(x,v,t):
    file = open('data_implicit.out', 'a')
    start_implicit = timeit.default_timer()
    for i in range(int(tmax/dt)):
        file.write('%f\t%f\n' % (t,x))
        xp=euler_step(x,v)
        # initialize
        tol=1e-4
        i=1
        vg=1.0
        dv=0.01
        # Newton-Raphson loop
        while i<10000:
            print(i)
            fvg=vg-v-a(xp,vg,t+dt)*dt
            v1=vg-dv
            v2=vg+dv
            fv1=v1-v-a(xp,v1,t+dt)*dt
            fv2=v2-v-a(xp,v2,t+dt)*dt
            # fv2-fv1 = ...
            # v2-v1-a(xp,v2)*dt+a(xp,v1)*dt
            # 2*dv-(-k*xp-c*v2)*dt+(-k*xp-c*v1)*dt
            # 2*dv + k*xp*dt + c*v2*dt - k*xp*dt - c*v1*dt
            # 2*dv + c*v2*dt - c*v1*dt
            # 2*dv + c*dt*(v2-v1)
            # 2*dv + c*dt*2*dv
            # (1+c*dt)*2*dv
            # (fv2-fv1)/(2*dv) = ...
            # 1+c*dt
            dfvg=(fv2-fv1)/2.0/dv
            if abs(fvg/dfvg)<tol: break
            vg=vg-fvg/dfvg
            i+=1
        vp=vg
        x=xp
        v=vp
        t+=dt
    stop_implicit = timeit.default_timer()
    file.close()
    print('Tot time implicit       = %f' % (stop_implicit-start_implicit))

def crankNicolson_method(x,v,t):
    file = open('data_crankNicolson.out', 'a')
    start_crankNicolson = timeit.default_timer()
    for i in range(int(tmax/dt)):
        file.write('%f\t%f\n' % (t,x))
        # vp-v-0.5*(-k*x-c*vp)*dt-0.5*(-k*x-c*v)*dt=0
        f=y-v-0.5*( a(x,y,t) + a(x,v,t) )*dt
        vp=newton_raphson(y,f,1.0)
        xp=euler_step(x,vp)
        x=xp
        v=vp
        t+=dt
    stop_crankNicolson = timeit.default_timer()
    file.close()
    print('Tot time Crank-Nicolson = %f' % (stop_crankNicolson-start_crankNicolson))

# global params
tmax=100.0
dt=0.5
k=1.0
c=0.05

# damping-to-spring-force ratio
zeta=c/2.0/np.sqrt(k)
print('zeta=%f' % (c/2.0/np.sqrt(k)))
if zeta>1: print('overdamped')
if zeta==1: print('critically damped')
if zeta<1: print('underdamped')

# initialize
x=0.0
v=1.0
t=0.0
y=np.linspace(0,10000,10001)*20.0/10000-10.0

semiImplicitEuler_method(x,v,t)
implicitEuler_method(x,v,t)
crankNicolson_method(x,v,t)

# read
datae=np.loadtxt('data_semiImplicit.out')
datai=np.loadtxt('data_implicit.out')
datac=np.loadtxt('data_crankNicolson.out')
t=datae[:,0]
xe=datae[:,1]
xi=datai[:,1]
xc=datac[:,1]

# analytical solution
x0=0.0
v0=1.0

xa=homogeneous_solution(c,k) # assumes m=1

F0=0.3
omegaF=0.5
xa+=inhomogeneous_solution(1.0,c,k) # m,b,k


plt.figure()
plt.plot(t,xa, label='Analytical solution', linewidth=15, color='gray')
plt.plot(t,xe, label='Semi-implicit Euler', linewidth=10)
plt.plot(t,xi, label='Implicit Euler', linewidth=5)
plt.plot(t,xc, label='Crank-Nicolson Euler', color='black')
plt.legend()
plt.xlabel('Time')
plt.ylabel('Solution')
plt.ylim((-1.1,1.1))
plt.show()
