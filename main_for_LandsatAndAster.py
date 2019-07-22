import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import ndimage
import csv
import sys
import os, re
from datetime import datetime
import matplotlib
import matplotlib.patches as patches

import plume_abstraction_clustering as L
import ASTER_plume_detection as A



# in both scripts it is done for HINKLEY
#change inside those scripts for other locations

mask = '/home/users/mp877190/CODE/plume_detection_clustering/Hinkley_landmask.nc'

plume_cube_aster = A.ASTER_plume_main()
plume_cube_landsat = L.Landsat_plume_main()
	
#plume_cube_all = np.concatenate((plume_cube_aster, plume_cube_landsat), axis=0)
plume_cube_all = plume_cube_landsat
#for i in plume_cube_all:
	#plt.imshow(i)
	#plt.show()
end_result = np.nansum(plume_cube_all, axis=0)
end_result = (end_result / np.amax(end_result)) * 100

print mask

landmask = np.flipud(np.array(nc.Dataset(mask).variables['Stat_landmask']))
a = np.ma.masked_equal(landmask, 0)
aMask = np.ma.getmaskarray(a)
result_masked = np.ma.array(end_result, mask=aMask, fill_value=np.nan)

## PLOT THE PROBABILITY DENSITY MAP 
my_cmap = matplotlib.cm.hot_r 
my_cmap.set_under(color='paleturquoise') #, alpha=0.5)
my_cmap.set_bad(color='sienna')
plt.imshow(result_masked, cmap= my_cmap, vmin=1)
# Create a Rectangle patch
#rect = patches.Rectangle((250,0),150,150,linewidth=1,edgecolor='k',facecolor='none')
#ax.add_patch(rect)
plt.colorbar()
plt.title('HINKLEY - Probability of plume extent from Landsat and ASTER \n based on 1.5 deg above ambient') #3 st dev, no mask extension')
plt.tight_layout()
plt.show()
