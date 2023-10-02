import yt

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.ion()

basename='edr'

base = basename
base = "%snew" % basename
base = "%s_aux_" % basename
# Load the dataset.

xprobmin1     = -2.5e0
xprobmax1     = 2.5e0
xprobmin2     = -2.5e0
xprobmax2     = 2.5e0
xprobmin3     = 0.7e0
xprobmax3     = 12e0

mx=64
my=64
mz=64

from os.path import join
unitFile=join(".","units.dat")
unitLength,unitTime, unitDens, unitVel,unitPres,unitMag,unitTemp = np.loadtxt(unitFile, usecols=(0,1,2,3,4,5,6),unpack=True)
unitJ = unitDens * unitVel/(unitMag * unitTime)
print("UNITJ=",unitJ)


zheight=1.0
zheight=1.4
zz=np.linspace(xprobmin3,xprobmax3,mz)
zindex=np.argmin(np.absolute(zz-zheight))

center = [0.5*(xprobmin1+xprobmax1),0.5*(xprobmin2+xprobmax2),xprobmax3]
width=[xprobmax1-xprobmin1,xprobmax2-xprobmin2,0]

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


# Create a 3D matplotlib figure for visualizing the surface
fig = plt.figure()
ax = fig.add_subplot(projection="3d")


listFiles=range(87,10000)
for kk in listFiles:
    ax.cla()
    ax.invert_yaxis()  
    ds = yt.load("%s%04d.dat" % (base,kk))
    ds.periodicity = (True, True, True)
    print(ds.field_list)
    print(ds['time'])
    #ds.add_field("amrvac", "")
    # Create a slice plot for the dataset.  With no additional arguments,
    # the width will be the size of the domain and the center will be the
    # center of the simulation box

    # Identify the isodensity surface in this sphere with density = 1e-24 g/cm^3
    rr = ds.r[xprobmin1:xprobmax1,xprobmin2:xprobmax2,xprobmin3:xprobmax3]
    rr = ds.r[-2:2,-2:2,xprobmin3:xprobmax3]
    rr = ds.r[-2:2,-2:2,0.75:11.95]
    rr = ds.all_data()
    print("MINMAX",np.min(rr["jz"]), np.max(rr["jz"]))

    vals=[-4e-3,4e-3]
    #vals=[-2e-4,2e-4]

    for vv in vals:
      surface = ds.surface(rr,"jz", vv)
      
      # Color this isodensity surface according to the log of the temperature field
      #colors = yt.apply_colormap(np.log10(surface["jz"]), cmap_name="hot")
      colors = yt.apply_colormap(surface["jz"], cmap_name="seismic")
      
      p3dc = Poly3DCollection(surface.triangles, linewidth=0.0)
      
      # Set the surface colors in the right scaling [0,1]
      p3dc.set_facecolors(colors[0, :, :] / 255.0)
      ax.add_collection(p3dc)
    
    # Let's keep the axis ratio fixed in all directions by taking the maximum
    # extent in one dimension and make it the bounds in all dimensions
    max_extent = (surface.vertices.max(axis=1) - surface.vertices.min(axis=1)).max()
    centers = (surface.vertices.max(axis=1) + surface.vertices.min(axis=1)) / 2
    bounds = np.zeros([3, 2])
    bounds[:, 0] = centers[:] - max_extent / 2
    bounds[:, 1] = centers[:] + max_extent / 2
    print("BONDS",bounds)
    ax.auto_scale_xyz(bounds[0, :], bounds[1, :], bounds[2, :])
    
    # Save the figure
    plt.savefig("%d_Surface.png" % kk)

