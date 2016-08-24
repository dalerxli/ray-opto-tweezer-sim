import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

import scipy.interpolate as interpol

values = np.loadtxt("2d-na068-linear-a1.tsv")

# Coordinates
xx = values[:,0].flatten()
zz = values[:,2].flatten()

# Forces
ux = values[:,3].flatten()
uz = values[:,5].flatten()

# Add the negative x coordinates (for a symmetric plot)
xx = np.hstack([xx, -xx])
zz = np.hstack([zz, zz])
ux = np.hstack([ux, -ux])
uz = np.hstack([uz, uz])

#norms = np.sqrt(uz**2 + ux**2)

#ux = 15*ux/norms
#uz = 15*uz/norms

### Only take every nth row
#n = 4

#plt.quiver(zz[::n], xx[::n], uz[::n], ux[::n], norms[::n], cmap=cm.seismic, scale=400)
#plt.colorbar()
#plt.show()

#coords = np.vstack([zz,xx]).transpose()
#vals = np.vstack([uz,ux]).transpose()

zi = np.linspace(zz.min(), zz.max(), 100)
xi = np.linspace(xx.min(), xx.max(), 100)

ipts_z, ipts_x = np.meshgrid(zi, xi)

ux_int = interpol.griddata((zz, xx), ux, (ipts_z, ipts_x), method='cubic')
uz_int = interpol.griddata((zz, xx), uz, (ipts_z, ipts_x), method='cubic')

speed = np.sqrt(ux_int**2 + uz_int**2)

plt.streamplot(ipts_z, ipts_x, uz_int, ux_int,          # data
               color=speed,         # array that determines the colour
               cmap=cm.cool,        # colour map
               linewidth=2,         # line thickness
               arrowstyle='->',     # arrow style
               arrowsize=1.5)       # arrow size

#plt.colorbar()                      # add colour bar on the right

plt.title('Stream Plot, Dynamic Colour')
plt.show()