
import numpy as np
import params

def relative_var(xyz):
    x,y,z=xyz
    dx=x[:,1:]-x[:,0:-1]
    dy=y[:,1:]-y[:,0:-1]
    dz=z[:,1:]-z[:,0:-1]
    return dx,dy,dz

def mytr(xyz):
    x,y,z=xyz
    return x.transpose(),y.transpose(),z.transpose()

def convert2col(f):
    f_col=f.reshape((len(f),1))
    return f_col
    
def set_crossProduct_array(x,y,z):
    x_flat=x.flatten()
    y_flat=y.flatten()
    z_flat=z.flatten()
    x_col=convert2col(x_flat)
    y_col=convert2col(y_flat)
    z_col=convert2col(z_flat)
    myarray=np.concatenate((x_col, y_col, z_col), axis=1)
    return myarray

def unset_crossProduct_array(f):
    fx=f[:,0]
    fy=f[:,1]
    fz=f[:,2]
    fx=fx.reshape((params.N1+2,params.N2+2))
    fy=fy.reshape((params.N1+2,params.N2+2))
    fz=fz.reshape((params.N1+2,params.N2+2))
    return fx,fy,fz

# P1___P2
#  |  /|    Area assigned to
#  | / |    P1 is that of the 
#  |/__|    triangle P1P2P4
# P4   P3
def reference_area(x,y,z):
    A=np.zeros(np.shape(x))
    P1P2x,P1P2y,P1P2z=relative_var((x,y,z))
    P1P4x,P1P4y,P1P4z=mytr(relative_var(mytr((x,y,z))))
    P1P2x=np.append(P1P2x,convert2col(P1P2x[:,-1]),axis=1)      # duplicate last column
    P1P2y=np.append(P1P2y,convert2col(P1P2y[:,-1]),axis=1)
    P1P2z=np.append(P1P2z,convert2col(P1P2z[:,-1]),axis=1)
    P1P4x=np.append(P1P4x,[P1P4x[-1,:]],axis=0)
    P1P4y=np.append(P1P4y,[P1P4y[-1,:]],axis=0)
    P1P4z=np.append(P1P4z,[P1P4z[-1,:]],axis=0)
    P1P2_col=set_crossProduct_array(P1P2x,P1P2y,P1P2z)
    P1P4_col=set_crossProduct_array(P1P4x,P1P4y,P1P4z)
    P1P2xP1P4_col=np.cross(P1P2_col,P1P4_col)                   # cross product
    P1P2xP1P4x,P1P2xP1P4y,P1P2xP1P4z=unset_crossProduct_array(P1P2xP1P4_col)
    P1P2xP1P4=np.sqrt(P1P2xP1P4x*P1P2xP1P4x+P1P2xP1P4y*P1P2xP1P4y+P1P2xP1P4z*P1P2xP1P4z)
    A=0.5*P1P2xP1P4
    # unit normal vector of the panel
    ex,ey,ez=P1P2xP1P4x/P1P2xP1P4,P1P2xP1P4y/P1P2xP1P4,P1P2xP1P4z/P1P2xP1P4
    return A,ex,ey,ez

def safe_div(a,b):
    c = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
    return c

# Neighbors
# reference
# numbering:
#     3
#     |
# 1 --o-- 2
#     |
#     4
def neighbor_info(x,y,z):
    xyz1=x[:, :-1],y[:, :-1],z[:, :-1]
    xyz2=x[:,1:  ],y[:,1:  ],z[:,1:  ]
    xyz3=x[ :-1,:],y[ :-1,:],z[ :-1,:]
    xyz4=x[1:  ,:],y[1:  ,:],z[1:  ,:]
    xyz3=mytr(xyz3)
    xyz4=mytr(xyz4)
    return xyz1,xyz2,xyz3,xyz4

def relaxed_distance(lz,ly):    # relaxed interparticle distance matrix
    lz1=lz[:,:-1]
    lz2=lz[:,1:]
    ly3=np.transpose(ly[:-1,:])
    ly4=np.transpose(ly[1:,:])
    return lz1,lz2,ly3,ly4

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
    while i<10000:
        fxg,dfxg=fxANDdfx(x,f,xg)
        #print('%i \t %.8f \t %.8f' % (i,xg,fxg))
        #if abs(fxg)<tol: break
        if abs(fxg/dfxg)<tol: break
        xg=xg-fxg/dfxg
        i+=1
    return xg

