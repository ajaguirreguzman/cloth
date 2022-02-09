
import numpy as np
import globals
from useful_scripts import neighbor_info,relaxed_distance
from forces import spring_force

def calc_energy(x,y,z,vx,vy,vz):
    k1=globals.k1
    k2=globals.k2
    ly=globals.ly
    lz=globals.lz
    xyz1,xyz2,xyz3,xyz4=neighbor_info(x,y,z)    # xyz3 and xyz4 are transposed for convenience
    lz1,lz2,ly3,ly4=relaxed_distance(lz,ly)     # ly3 and ly4 are transposed for convenience
    na,V1=spring_force( k1,lz1,xyz1)            # k1 is modulus of elasticity along z; V1 = potential energy in spring 1
    na,V2=spring_force(-k1,lz2,xyz2)            # notice the -k
    na,V3=spring_force( k2,ly3,xyz3)            # since all arguments are transposed, then Fs3 is also transposed
    na,V4=spring_force(-k2,ly4,xyz4)
    V3=V3.transpose()
    V4=V4.transpose()
    V=np.sum(V1)+np.sum(V2[:,-1])+np.sum(V3)+np.sum(V4[-1,:])
    K=0.5*np.sum(vx*vx+vy*vy+vz*vz)
    return K,V

def calc_posChange(x,y,z,x0,y0,z0):
    changex=x-x0
    changey=y-y0
    changez=z-z0
    change=np.sqrt(changex*changex+changey*changey+changez*changez)
    return np.max(change)