import yt

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.ion()

base = 'implicit_coupled'
base = 'implicit_couplednew'
base = 'implicit_coupled_aux_'
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


center = [0.5*(xprobmin1+xprobmax1),0.5*(xprobmin2+xprobmax2),xprobmax3]
width=[xprobmax1-xprobmin1,xprobmax2-xprobmin2,0]


fig,ax=plt.subplots(nrows=1,ncols=1)


listFiles=range(10000)
for kk in listFiles:
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
    arr = data['jz'].to_ndarray()[:,:,mz-1]
    print("MINMAX2 ", np.min(arr), np.max(arr))
#    print("MINMAX ", np.min(data['b3'].to_ndarray()), np.max(data['b3'].to_ndarray()))
#    print("MINMAX mc1 ", np.min(data['m_c1'].to_ndarray()), np.max(data['m_c1'].to_ndarray()))
#    print("MINMAX mc2 ", np.min(data['m_c2'].to_ndarray()), np.max(data['m_c2'].to_ndarray()))
#    print("MINMAX mc3 ", np.min(data['m_c3'].to_ndarray()), np.max(data['m_c3'].to_ndarray()))
    print(type(arr))
    im=ax.imshow(arr,origin='lower',norm=matplotlib.colors.Normalize(vmin=-1e-8,vmax=1e-8),cmap='seismic')    
    fig.canvas.draw()    
    plt.draw()
    plt.show()
    plt.savefig("B3%04d.png" % kk)
