
# Included:
# - Linear-elastic interparticle, frictional, hydrodynamic (drag and lift) and gravitational forces
# To be included:
# - Buoyancy forces
# - Constrained dynamics: currently constraints on positions (artificially) advance them one time step ahead, while velocities stay behind. This is incorrect. Implement constraints as a force instead.
# - Wilson 1967's expression for the linear-elastic response of the net. This is used by Bi et al. 2014a,b.

import globals
from useful_scripts import *
from hydro_models import *

def spring_force(k,l,xyz):
    dx,dy,dz=relative_var(xyz)
    r=np.sqrt(dx*dx+dy*dy+dz*dz)            # relative distance
    F=-k*(r-l)                              # restoring force
    V=0.5*abs(k)*(r-l)**2.0
    return F,V

def frictional_force(c,xyz,vxyz):
    dx,dy,dz=relative_var(xyz)
    vx,vy,vz=relative_var(vxyz)
    r=np.sqrt(dx*dx+dy*dy+dz*dz)            # relative distance
    v=np.sqrt(vx*vx+vy*vy+vz*vz)            # relative velocity
    xn,yn,zn=dx/r,dy/r,dz/r                 # unit vector components
    F=-c*(xn*vx+yn*vy+zn*vz)                # frictional force
    return F

def hydrodynamic_force(x,y,z,vx,vy,vz):
    ux=globals.ux
    uy=globals.uy
    uz=globals.uz
    rho=globals.rho
    Cd,Cl=hydrodynamic_coeffs()
    A,ex,ey,ez=reference_area(x,y,z)                # area of net panel element (screen?) and normal vector components
    urx,ury,urz=ux-vx,uy-vy,uz-vz                   # relative velocity components
    ur=np.sqrt(urx*urx+ury*ury+urz*urz)             # relative velocity magnitude
    iDx,iDy,iDz=urx/ur,ury/ur,urz/ur                # unit vector for drag force
    uv_col=set_crossProduct_array(urx,ury,urz)      # setup column vector for cross product
    e_col=set_crossProduct_array(ex,ey,ez)          # setup column vector for cross product
    P_col=np.cross(uv_col,np.cross(e_col,uv_col))   # P = (u-v) x e_n x (u-v)
    Px,Py,Pz=unset_crossProduct_array(P_col)        # back to 2D array 
    P=np.sqrt(Px*Px+Py*Py+Pz*Pz)                    # norm 
    iLx,iLy,iLz=safe_div(Px,P),safe_div(Py,P),safe_div(Pz,P)    # unit vector for lift force
    Fdx=0.5*Cd*rho*A*ur*ur*iDx                      # drag force
    Fdy=0.5*Cd*rho*A*ur*ur*iDy
    Fdz=0.5*Cd*rho*A*ur*ur*iDz
    Flx=0.5*Cl*rho*A*ur*ur*iLx                      # lift force
    Fly=0.5*Cl*rho*A*ur*ur*iLy
    Flz=0.5*Cl*rho*A*ur*ur*iLz
    return Fdx,Fdy,Fdz,Flx,Fly,Flz

def constraint_force():
    lz=globals.lz
    d1=z[:,1:]-z[:,:-1]
    d2=abs(d1)
    d3=(d2-lz)/d2
    z[:,:-1]=z[:,:-1]+0.5*d1*d3
    z[:,1:] =z[:,1:] -0.5*d1*d3
    Fx=0.0
    Fy=0.0
    Fz=0.0
    return Fx,Fy,Fz

def add_anchors(ax,ay,az):
    """
    ax[-1,0]=0.0
    ay[-1,0]=0.0
    az[-1,0]=0.0
    ax[-1,-1]=0.0
    ay[-1,-1]=0.0
    az[-1,-1]=0.0
    ax[0,:]=0.0
    ay[0,:]=0.0
    az[0,:]=0.0
    ax[-1,:]=0.0
    ay[-1,:]=0.0
    az[-1,:]=0.0
    """
    ax[:,0]=0.0
    ay[:,0]=0.0
    az[:,0]=0.0
    ax[:,-1]=0.0
    ay[:,-1]=0.0
    az[:,-1]=0.0
    return ax,ay,az

def force_components(F,xyz):
    dx,dy,dz=relative_var(xyz)
    r=np.sqrt(dx*dx+dy*dy+dz*dz)            # relative distance
    Fx=F*(dx/r)
    Fy=F*(dy/r)
    Fz=F*(dz/r)
    return Fx,Fy,Fz

def add_forces(F1,F2,F3,F4):
    N1=globals.N1
    N2=globals.N2
    # extract components
    F1x,F1y,F1z=F1
    F2x,F2y,F2z=F2
    F3x,F3y,F3z=F3
    F4x,F4y,F4z=F4
    # initialize
    ax=np.zeros((N1+2,N2+2))
    ay=np.zeros((N1+2,N2+2))
    az=np.zeros((N1+2,N2+2))
    # particles with 4 neighbors
    ax[1:-1,1:-1]=(F1x+F2x)[1:-1,:]+(F3x+F4x)[:,1:-1]
    ay[1:-1,1:-1]=(F1y+F2y)[1:-1,:]+(F3y+F4y)[:,1:-1]
    az[1:-1,1:-1]=(F1z+F2z)[1:-1,:]+(F3z+F4z)[:,1:-1]
    # left
    ax[1:-1,0]=-F1x[1:-1,0]+(F3x+F4x)[:,0]
    ay[1:-1,0]=-F1y[1:-1,0]+(F3y+F4y)[:,0]
    az[1:-1,0]=-F1z[1:-1,0]+(F3z+F4z)[:,0]
    # top
    ax[0,1:-1]=(F1x+F2x)[0,:]-F3x[0,1:-1]
    ay[0,1:-1]=(F1y+F2y)[0,:]-F3y[0,1:-1]
    az[0,1:-1]=(F1z+F2z)[0,:]-F3z[0,1:-1]
    # right
    ax[1:-1,-1]=-F2x[1:-1,-1]+(F3x+F4x)[:,-1]
    ay[1:-1,-1]=-F2y[1:-1,-1]+(F3y+F4y)[:,-1]
    az[1:-1,-1]=-F2z[1:-1,-1]+(F3z+F4z)[:,-1]
    # bottom
    ax[-1,1:-1]=(F1x+F2x)[-1,:]-F4x[-1,1:-1]
    ay[-1,1:-1]=(F1y+F2y)[-1,:]-F4y[-1,1:-1]
    az[-1,1:-1]=(F1z+F2z)[-1,:]-F4z[-1,1:-1]
    # top left corner
    ax[0,0]=-F1x[0,0]-F3x[0,0]
    ay[0,0]=-F1y[0,0]-F3y[0,0]
    az[0,0]=-F1z[0,0]-F3z[0,0]
    # top right corner
    ax[0,-1]=-F2x[0,-1]-F3x[0,-1]
    ay[0,-1]=-F2y[0,-1]-F3y[0,-1]
    az[0,-1]=-F2z[0,-1]-F3z[0,-1]
    # bottom left corner
    ax[-1,0]=-F1x[-1,0]-F4x[-1,0]
    ay[-1,0]=-F1y[-1,0]-F4y[-1,0]
    az[-1,0]=-F1z[-1,0]-F4z[-1,0]
    # bottom right corner
    ax[-1,-1]=-F2x[-1,-1]-F4x[-1,-1]
    ay[-1,-1]=-F2y[-1,-1]-F4y[-1,-1]
    az[-1,-1]=-F2z[-1,-1]-F4z[-1,-1]
    return ax,ay,az

def a(x,y,z,vx,vy,vz):
    k1=globals.k1
    k2=globals.k2
    ly=globals.ly
    lz=globals.lz
    c=globals.c
    m=globals.m
    # preliminaries
    xyz1,xyz2,xyz3,xyz4=neighbor_info(x,y,z)            # xyz3 and xyz4 are transposed for convenience
    vxyz1,vxyz2,vxyz3,vxyz4=neighbor_info(vx,vy,vz)     # vxyz3 and vxyz4 are transposed for convenience
    lz1,lz2,ly3,ly4=relaxed_distance(lz,ly)             # ly3 and ly4 are transposed for convenience
    # spring force from neighbors
    Fs1,na=spring_force( k1,lz1,xyz1)           # k1 is modulus of elasticity along z; V1 = potential energy in spring 1
    Fs2,na=spring_force(-k1,lz2,xyz2)           # notice the -k
    Fs3,na=spring_force( k2,ly3,xyz3)           # since all arguments are transposed, then Fs3 is also transposed
    Fs4,na=spring_force(-k2,ly4,xyz4)
    # frictional force
    Ff1=frictional_force( c,xyz1,vxyz1)
    Ff2=frictional_force(-c,xyz2,vxyz2)
    Ff3=frictional_force( c,xyz3,vxyz3)         # same here, Ff3 is transposed
    Ff4=frictional_force(-c,xyz4,vxyz4)
    # total force components per neighbor
    F1=force_components(Fs1+Ff1,xyz1)
    F2=force_components(Fs2+Ff2,xyz2)
    F3=mytr(force_components(Fs3+Ff3,xyz3))     # now we transpose back
    F4=mytr(force_components(Fs4+Ff4,xyz4))
    # acceleration due to interaction with neighbors
    ax,ay,az=add_forces(F1,F2,F3,F4)
    # add drag and lift forces
    if globals.add_hydroForce:
        Fdx,Fdy,Fdz,Flx,Fly,Flz=hydrodynamic_force(x,y,z,vx,vy,vz)
        ax+=Fdx+Flx
        ay+=Fdy+Fly
        az+=Fdz+Flz
    # divide by mass
    ax=ax/m
    ay=ay/m
    az=az/m
    # add gravitational acceleration
    if globals.add_gravityForce:
        ay+=-9.8
    # constraints
    # ADD CONSTRAINTS HERE
    # https://en.wikipedia.org/wiki/Verlet_integration#Constraints
    # https://en.wikipedia.org/wiki/Constraint_(computational_chemistry)
    # add anchors
    ax,ay,az=add_anchors(ax,ay,az)
    return ax,ay,az
