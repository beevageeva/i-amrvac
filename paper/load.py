import yt

import numpy as np

base = 'C2/edrnew'
# Load the dataset.

xprobmin1     = -2.5e0
xprobmax1     = 2.5e0
xprobmin2     = -2.5e0
xprobmax2     = 2.5e0
xprobmin3     = 0.7e0
xprobmax3     = 12e0

center = [0.5*(xprobmin1+xprobmax1),0.5*(xprobmin2+xprobmax2),xprobmax3]
width=[xprobmax1-xprobmin1,xprobmax2-xprobmin2,0]


listFiles=range(40,41)
for kk in listFiles:
    ds = yt.load("%s%04d.dat" % (base,kk))
    print(ds.field_list)
    #ds.add_field("amrvac", "")
    # Create a slice plot for the dataset.  With no additional arguments,
    # the width will be the size of the domain and the center will be the
    # center of the simulation box
    slc = yt.SlicePlot(ds, "z", ("amrvac", "b3"), center=center, width=width)
    data = ds.all_data()
    slc.set_log(("amrvac", "b3"), False)
    print("MINMAX ", np.min(data["b3"]), np.max(data["b3"]))


    slc.save("B3%04d.png" % kk)
