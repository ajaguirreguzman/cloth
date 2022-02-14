
import numpy as np
import globals

# For definition of theta, see Fig.5 in Cheng et al. 2020.

def hydrodynamic_coeffs(ex,ey,ez):
    if globals.model=='S1':
        ux=globals.ux
        uy=globals.uy
        uz=globals.uz
        Sn=2.0*globals.dW/globals.L
        ex=abs(ex)
        ey=abs(ey)
        ez=abs(ez)
        udote=ux*ex+uy*ey+uz*ez
        normu=np.sqrt(ux*ux+uy*uy+uz*uz)
        norme=np.sqrt(ex*ex+ey*ey+ez*ez)
        theta=np.arccos(udote/normu/norme)
        Cd=0.04+(-0.04+Sn-1.24*Sn*Sn+13.7*Sn*Sn*Sn)*np.cos(theta)
        Cl=(0.57*Sn-3.54*Sn*Sn+10.1*Sn*Sn*Sn)*np.sin(2.0*theta)
    else:
        Cd=1.0
        Cl=1.0
    return Cd,Cl
