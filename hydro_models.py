
import numpy as np
import globals

def hydrodynamic_coeffs():
    model=globals.model
    L=globals.L
    dW=globals.dW
    theta=globals.theta
    if model=='S1':
        Sn=2.0*dW/L
        Cd=0.04+(-0.04+Sn-1.24*Sn*Sn+13.7*Sn*Sn*Sn)*np.cos(theta)
        Cl=(0.57*Sn-3.54*Sn*Sn+10.1*Sn*Sn*Sn)*np.sin(2.0*theta)
    else:
        Cd=1.0
        Cl=1.0
    return Cd,Cl
