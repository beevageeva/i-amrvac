import yt

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
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


zheight=6.4
zheight=1.4
#zheight=1.0
zz=np.linspace(xprobmin3,xprobmax3,mz)
xx=np.linspace(xprobmin1,xprobmax1,mx)
yy=np.linspace(xprobmin2,xprobmax2,my)
zindex=np.argmin(np.absolute(zz-zheight))

center = [0.5*(xprobmin1+xprobmax1),0.5*(xprobmin2+xprobmax2),xprobmax3]
width=[xprobmax1-xprobmin1,xprobmax2-xprobmin2,0]


fig,ax=plt.subplots(nrows=1,ncols=1)

cb=None

listFiles=range(0,10000)
for kk in listFiles:
    ax.cla()
    ds = yt.load("%s%04d.dat" % (base,kk))
    print(ds.field_list)
    print(ds['time'])
    #ds.add_field("amrvac", "")
    # Create a slice plot for the dataset.  With no additional arguments,
    # the width will be the size of the domain and the center will be the
    # center of the simulation box
    data = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
    #arr = data['b3'].to_ndarray()[:,:,mz-1]
    #arr = data['b3'].to_ndarray()[mx//2,:,:]
    #arr = data['m_c2'].to_ndarray()[mx//2,:,:]/data['rho_c'].to_ndarray()[mx//2,:,:]
    #arr = data['m_c2'].to_ndarray()[mx//2,:,:]/data['rho_c'].to_ndarray()[mx//2,:,:]
    #arr = data['jz'].to_ndarray()[mx//2,:,:]
    #arr = data['jz'].to_ndarray()[:,:,mz-1]
    arr = data['jz'].to_ndarray()[:,:,zindex]*unitJ*1e9
    #arr = data['jz'].to_ndarray()[:,:,zindex]
    #arr = data['jz'].to_ndarray()[:,my//2,:]*unitJ*1e9
    #arr = data['jz'].to_ndarray()[mx//3,:,:]*unitJ*1e9
    #arr = data['rho_c'].to_ndarray()[:,:,zindex]
    print("MINMAX2 ", np.min(arr), np.max(arr))
#    print("MINMAX ", np.min(data['b3'].to_ndarray()), np.max(data['b3'].to_ndarray()))
#    print("MINMAX mc1 ", np.min(data['m_c1'].to_ndarray()), np.max(data['m_c1'].to_ndarray()))
#    print("MINMAX mc2 ", np.min(data['m_c2'].to_ndarray()), np.max(data['m_c2'].to_ndarray()))
#    print("MINMAX mc3 ", np.min(data['m_c3'].to_ndarray()), np.max(data['m_c3'].to_ndarray()))
    print(type(arr))
    mm=max(abs(np.min(arr)),abs(np.max(arr)) )

    arr=arr.transpose()
    im=ax.imshow(arr,origin='lower',norm=matplotlib.colors.Normalize(vmin=-mm,vmax=mm),cmap='seismic',extent=(xprobmin1,xprobmax1,xprobmin2,xprobmax2))    
    if(not cb is None):
      cb.remove()
    cb=plt.colorbar(im)
    cb.set_label(r"$J_z$ [nA/m$^2$]")
    levels=[-2,-1,1,2]
    levels=[-2,-1,-0.5,-0.25,0.25,0.5,1,2]
    cp = ax.contour(xx,yy,arr,levels,cmap="PiYG",linewidths=1.5)
    ax.set_title("z=140 km")
    ax.text(0.92, 0.9, ("%.1f s"% (float(ds["time"])*unitTime) ),size=10, color='k',ha="center", va="center",transform=ax.transAxes)
    fig.canvas.draw()    
    plt.draw()
    plt.show()
    plt.savefig("B3%04d.png" % kk)
