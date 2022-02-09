
import numpy as np
import matplotlib.pyplot as plt
from root_finding import *

x=np.linspace(0,10000,10001)*1.0/10000-0.0
f=np.cos(np.pi*x)
#f=(x-1.0/3.0)**3.0

plt.figure()
plt.plot(x,f,'o')
plt.show()

r=newton_raphson(x,f,0.75)
print(r)
