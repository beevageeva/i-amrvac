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

from os.path import join
unitFile=join(".","units.dat")
unitLength,unitTime, unitDens, unitVel,unitPres,unitMag,unitTemp = np.loadtxt(unitFile, usecols=(0,1,2,3,4,5,6),unpack=True)
unitJ = unitDens * unitVel/(unitMag * unitTime)
print("UNITJ=",unitJ)


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


center = [0.5*(xprobmin1+xprobmax1),0.5*(xprobmin2+xprobmax2),0]
width=[xprobmax1-xprobmin1,xprobmax2-xprobmin2,0]


fig,ax=plt.subplots(nrows=1,ncols=1)



cb=None

listFiles=range(0,10000)


#listFiles=[15]
for kk in listFiles:
    ax.cla()
    base = "%s_aux_" % basename
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
    #arr = data['jz'].to_ndarray()[mx//2,:,:]
    #arr = data['jz'].to_ndarray()[:,:,mz-1]
    #arr = data['jy'].to_ndarray()[:,:,0]*unitJ*1e9
    jy = data['jy'].to_ndarray()[:,:,0]*unitJ*1e9
    jx = data['jx'].to_ndarray()[:,:,0]*unitJ*1e9

    #arr = data['jz'].to_ndarray()[:,:,zindex]
    #arr = data['jz'].to_ndarray()[:,my//2,:]*unitJ*1e9
    #arr = data['jz'].to_ndarray()[mx//3,:,:]*unitJ*1e9
    #arr = data['rho_c'].to_ndarray()[:,:,zindex]
#    print("MINMAX ", np.min(data['b3'].to_ndarray()), np.max(data['b3'].to_ndarray()))
#    print("MINMAX mc1 ", np.min(data['m_c1'].to_ndarray()), np.max(data['m_c1'].to_ndarray()))
#    print("MINMAX mc2 ", np.min(data['m_c2'].to_ndarray()), np.max(data['m_c2'].to_ndarray()))
#    print("MINMAX mc3 ", np.min(data['m_c3'].to_ndarray()), np.max(data['m_c3'].to_ndarray()))

    print("LEN xx ", len(xx))
    print("LEN yy ", len(yy))
    print("JXSHAPE ", jx.shape)
    print("JYSHAPE ", jy.shape)

    xplot=xx[xminindex:xmaxindex+1]
    yplot=yy[yminindex:ymaxindex+1]


    jx=jx[xminindex:xmaxindex+1,yminindex:ymaxindex+1]
    jy=jy[xminindex:xmaxindex+1,yminindex:ymaxindex+1]

    #ax.quiver(xx,yy,jx.T,jy.T,scale=1e4, units='xy',
    #              color="orange", width=1.0, headwidth=0.8)
    ax.quiver(xplot,yplot,jx.T,jy.T)
    #ax.streamplot(xplot,yplot,jx.T,jy.T)

    ax.text(0.92, 0.9, ("%.1f s"% (float(ds["time"])*unitTime) ),size=10, color='k',ha="center", va="center",transform=ax.transAxes)
    fig.canvas.draw()    
    plt.draw()
    plt.show()
    plt.savefig("J%04d.png" % kk)
