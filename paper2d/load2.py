import yt

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ion()

basename='edr'


xprobmin1     = -2.5e0
xprobmax1     = 2.5e0
xprobmin2     = 0.9e0
xprobmax2     = 12.5e0


mx=128
my=256

xx=np.linspace(xprobmin1,xprobmax1,mx)
yy=np.linspace(xprobmin2,xprobmax2,my)

xmin=-1.5
xmax=1.5
ymin=0.9
ymax=3.05


xminindex=np.argmin(np.absolute(xx-xmin))
xmaxindex=np.argmin(np.absolute(xx-xmax))
yminindex=np.argmin(np.absolute(yy-ymin))
ymaxindex=np.argmin(np.absolute(yy-ymax))
dir1="."

from os.path import join
unitFile=join(dir1,"units.dat")
unitLength,unitTime, unitDens, unitVel,unitPres,unitMag,unitTemp = np.loadtxt(unitFile, usecols=(0,1,2,3,4,5,6),unpack=True)
unitJ = unitDens * unitVel/(unitMag * unitTime)
unitVel=unitLength/unitTime

newUnitJ=unitJ*1e6

print("UNITJ=",unitJ)


xx=np.linspace(xprobmin1,xprobmax1,mx)
yy=np.linspace(xprobmin2,xprobmax2,my)



fig,ax=plt.subplots(nrows=1,ncols=1)


dataset="vel"
dataset="b2"
dataset="jy"
dataset="rho_c1"
dataset="b3"
dataset="vc2"
dataset="jx"
dataset="jz"
dataset="vc3"
dataset="vc1"


labels={}
labels["rho_c1"]=r"$\rho_{\rm c1}$ [kg/m$^3$]"
labels["vc3"]=r"$v_z$ [m/s]"
labels["vc2"]=r"$v_y$ [m/s]"
labels["vc1"]=r"$v_x$ [m/s]"
labels["jz"]=r"$J_z$ [$\mu$ A/m$^2$]"
labels["jy"]=r"$J_y$ [$\mu$ A/m$^2$]"
labels["jx"]=r"$J_x$ [$\mu$ A/m$^2$]"
labels["b3"]=r"$B_{\rm z1}$ [T]"
labels["b2"]=r"$B_{\rm y1}$ [T]"
labels["b1"]=r"$B_{\rm x1}$ [T]"

#https://matplotlib.org/users/colormapnorms.html
class MidpointNormalize(matplotlib.colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))



cb=None

listFiles=range(0,10000)
#listFiles=[0,1]
for kk in listFiles:
    ax.cla()
    if(dataset in ["jz", "jy", "jx"]):
      base = "%s_aux_" % basename
    elif(dataset in ["vc2", "vc3", "vc1", "vel", "rho_c"]):
      base = "%snew" % basename
    else:
      base = basename
    # Load the dataset.
    ds = yt.load(join(dir1,("%s%04d.dat" % (base,kk))))
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
    #arr = data['jz'].to_ndarray()[mx//2,:,:]
    #arr = data['jz'].to_ndarray()[:,:,mz-1]
    #arr = data['jy'].to_ndarray()[:,:,0]*unitJ*1e9
    if(dataset in ["jz","jy","jx"]):
      arr = data[dataset].to_ndarray()[:,:,0]*newUnitJ
    elif(dataset in ["b1","b2","b3"]):
      arr = data[dataset].to_ndarray()[:,:,0]*unitMag
    elif dataset=="vc1":  
      arr = (data['m_c1'].to_ndarray()[:,:,0]/data['rho_c'].to_ndarray()[:,:,0])*unitVel
    elif dataset=="vc2":  
      arr = (data['m_c2'].to_ndarray()[:,:,0]/data['rho_c'].to_ndarray()[:,:,0])*unitVel
    elif dataset=="vc3":  
      arr = (data['m_c3'].to_ndarray()[:,:,0]/data['rho_c'].to_ndarray()[:,:,0])*unitVel
    elif dataset=="vel":  
      arr = (np.sqrt( 
             (data['m_c1'].to_ndarray()[:,:,0])**2 + 
             (data['m_c2'].to_ndarray()[:,:,0])**2 +
             (data['m_c3'].to_ndarray()[:,:,0])**2 ) /data['rho_c'].to_ndarray()[:,:,0])*unitVel
    
    elif dataset=="rho_c1":
      arr=data["rho_c"].to_ndarray()[:,:,0] * unitDens  
    else:
      arr=data[dataset].to_ndarray()[:,:,0] 

    #arr = data['jz'].to_ndarray()[:,:,zindex]
    #arr = data['jz'].to_ndarray()[:,my//2,:]*unitJ*1e9
    #arr = data['jz'].to_ndarray()[mx//3,:,:]*unitJ*1e9
    #arr = data['rho_c'].to_ndarray()[:,:,zindex]
    arr=arr[xminindex:xmaxindex+1,yminindex:ymaxindex+1]
    print("MINMAX2 ", np.min(arr), np.max(arr))
#    print("MINMAX ", np.min(data['b3'].to_ndarray()), np.max(data['b3'].to_ndarray()))
#    print("MINMAX mc1 ", np.min(data['m_c1'].to_ndarray()), np.max(data['m_c1'].to_ndarray()))
#    print("MINMAX mc2 ", np.min(data['m_c2'].to_ndarray()), np.max(data['m_c2'].to_ndarray()))
#    print("MINMAX mc3 ", np.min(data['m_c3'].to_ndarray()), np.max(data['m_c3'].to_ndarray()))
    print(type(arr))
    mm=max(abs(np.min(arr)),abs(np.max(arr)) )
    if(dataset=="vc1"):
      mm=800
    norm=matplotlib.colors.Normalize(vmin=-mm,vmax=mm)
    cmap="seismic"
    if(dataset=="jz"):
      norm=MidpointNormalize(vmin=-1.1,vmax=0.2,midpoint=0.0)
    elif(dataset=="vel"):
      norm=matplotlib.colors.Normalize()
      cmap="binary"

    arr=arr.transpose()
    im=ax.imshow(arr,origin='lower',norm=norm,cmap=cmap,extent=(xx[xminindex],xx[xmaxindex],yy[yminindex],yy[ymaxindex]))    
    if(not cb is None):
      cb.remove()
    cb=plt.colorbar(im)
    cb.set_label(labels[dataset])
    levels=[-2,-1,1,2]
    levels=[-2,-1,-0.5,-0.25,0.25,0.5,1,2]
    #cp = ax.contour(xx,yy,arr,levels,cmap="PiYG",linewidths=1.5)
    ax.set_title("z=140 km")
    ax.text(0.92, 0.9, ("%.1f s"% (float(ds["time"])*unitTime) ),size=10, color='k',ha="center", va="center",transform=ax.transAxes)
    fig.canvas.draw()    
    plt.draw()
    plt.show()
    fn=join(dir1,("%s%04d.png" % (dataset, kk)))
    plt.savefig(fn)
