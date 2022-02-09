
import numpy as np

def write_continua(x,y,z,ly,lz):
    print('Saving state...')
    np.savetxt('continua_x.out', x, delimiter = ',')
    np.savetxt('continua_y.out', y, delimiter = ',')
    np.savetxt('continua_z.out', z, delimiter = ',')
    np.savetxt('continua_ly.out', ly, delimiter = ',')
    np.savetxt('continua_lz.out', lz, delimiter = ',')

def read_continua():
    print('Reading continua files...')
    x=np.loadtxt('continua_x.out', delimiter = ',')
    y=np.loadtxt('continua_y.out', delimiter = ',')
    z=np.loadtxt('continua_z.out', delimiter = ',')
    ly=np.loadtxt('continua_ly.out', delimiter = ',')
    lz=np.loadtxt('continua_lz.out', delimiter = ',')
    return x,y,z,ly,lz
