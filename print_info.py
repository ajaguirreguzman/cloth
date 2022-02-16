
import params
import numpy as np


def print_preInfo():
    Sn=2.0*params.twineD/params.twineL
    omega=np.sqrt(params.k1[-1,-1]/params.m[-1,-1])
    print('Solidity ratio    = %.4f   ' % Sn)
    print('SHO frequency     = %.4f   ' % omega)
    print('Number of nodes   = %i     ' % (params.N1*params.N2))
