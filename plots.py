
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

def plot_elasticity(k1):
    fig=plt.figure()
    axs=plt.gca()
    h=axs.imshow(k1, cmap='Greys', norm=LogNorm(vmin=k1.min(), vmax=k1.max()),
                 interpolation='nearest', origin='lower')
    fig.colorbar(h, ax=axs)
    plt.show()
    exit()

def plot_hmap(A):
    fig=plt.figure()
    axs=plt.gca()
    hmap=axs.imshow(A, cmap='RdBu', interpolation='nearest', vmin=-1, vmax=1, origin='lower')
    cbar=fig.colorbar(hmap, ax=axs)
    return hmap
