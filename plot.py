
import matplotlib.pyplot as plt
import numpy as np

def peaks(t,y):
    yp=np.ones(len(y))
    tp=np.ones(len(y))
    k=0
    for i in range(len(y)-1):
        j=i+1
        if y[i]>y[i-1] and y[i]>y[i+1]:
            tp[k]=t[i]
            yp[k]=y[i]
            k+=1
    tp=tp[:k]
    yp=yp[:k]
    return tp,yp

data1=np.loadtxt('position_semiImplicitEuler.out', delimiter = ',')
data2=np.loadtxt('position_predictorCorrector.out', delimiter = ',')
data3=np.loadtxt('position_RK4.out', delimiter = ',')
data4=np.loadtxt('position.out', delimiter = ',')
t=data1[:,0]
z1=data1[:,1]
z2=data2[:,1]
z3=data3[:,1]
t4=data4[:,0]
z4=data4[:,1]

omega=np.sqrt(2e5)
za=(5.0/omega)*np.sin(omega*t)
za4=(5.0/omega)*np.sin(omega*t4)
err_trend=(1.0/2.0)*(5e-4)**2.0*(-5.0*omega*np.sin(omega*t))

tp1,zp1=peaks(t,za-z1)
tp2,zp2=peaks(t,za-z2)
tp3,zp3=peaks(t,za-z3)
tp4,zp4=peaks(t4,za4-z4)

fig=plt.figure()
axs=plt.gca()
#axs.plot(t,z)
#axs.plot(t,za)
axs.plot(tp1,zp1, label='Semi-implicit Euler')
axs.plot(tp2,zp2, label='PC, dt=5e-4')
axs.plot(tp3,zp3, label='RK4')
axs.plot(tp4,zp4, label='PC, dt=1e-4')
#axs.plot(t,err_trend)
#axs.set_ylim([-.1,.1])
plt.xlabel('Time (s)')
plt.ylabel('peaks( y_{analytical} - y_{numerical} )')
plt.legend()
plt.show()
