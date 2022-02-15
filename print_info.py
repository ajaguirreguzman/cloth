
import globals
import numpy as np


def print_preInfo():
    Sn=2.0*globals.twineD/globals.twineL
    omega=np.sqrt(globals.k1[-1,-1]/globals.m[-1,-1])
    print('Solidity ratio    = %.4f   ' % Sn)
    print('SHO frequency     = %.4f   ' % omega)
    print('Number of nodes   = %i     ' % (globals.N1*globals.N2))
