
import numpy as np

def write_continua(x,y,z,l1,l2,l3,l4):
    print('Saving state...')
    np.savetxt('continua_x.out', x, delimiter = ',')
    np.savetxt('continua_y.out', y, delimiter = ',')
    np.savetxt('continua_z.out', z, delimiter = ',')
    np.savetxt('continua_l1.out', l1, delimiter = ',')
    np.savetxt('continua_l2.out', l2, delimiter = ',')
    np.savetxt('continua_l3.out', l3, delimiter = ',')
    np.savetxt('continua_l4.out', l4, delimiter = ',')

def read_continua():
    print('Reading continua files...')
    x=np.loadtxt('continua_x.out', delimiter = ',')
    y=np.loadtxt('continua_y.out', delimiter = ',')
    z=np.loadtxt('continua_z.out', delimiter = ',')
    l1=np.loadtxt('continua_l1.out', delimiter = ',')
    l2=np.loadtxt('continua_l2.out', delimiter = ',')
    l3=np.loadtxt('continua_l3.out', delimiter = ',')
    l4=np.loadtxt('continua_l4.out', delimiter = ',')
    return x,y,z,l1,l2,l3,l4
