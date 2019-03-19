import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys, os
import netCDF4 as nc

import Landsat_ncfiles
import water_detection_mask as wdm

def heysham_bt_to_sst(bt_array):

	bt_sim = bt_array * 0.89207433 
	bt_sim = bt_sim + 30.02847216 
	sst_array = bt_sim * 1.08778792
	sst_array = sst_array - 22.97856999
	return sst_array

def hartlepool_bt_to_sst(bt_array):

	bt_sim = bt_array * 0.7079483 
	bt_sim = bt_sim + 82.10687098 
	sst_array = bt_sim * 1.0202827
	sst_array = sst_array - 4.01934945
	return sst_array
	
def hinkley_bt_to_sst(bt_array):

	bt_sim = bt_array *  0.86914992 
	bt_sim = bt_sim + 36.40805221 
	sst_array = bt_sim * 1.13316041
	sst_array = sst_array - 35.73586162
	return sst_array
	
def hunterston_bt_to_sst(bt_array):

	bt_sim = bt_array * 0.7427912 
	bt_sim = bt_sim + 71.99667606 
	sst_array = bt_sim * 0.99927464 
	sst_array = sst_array + 2.21099584
	return sst_array
	
def dungeness_bt_to_sst(bt_array):

	bt_sim = bt_array * 1.15388349 
	bt_sim = bt_sim - 45.08048914
	sst_array = bt_sim * 1.11363911
	sst_array = sst_array - 30.29539824
	return sst_array
	
def sizewell_bt_to_sst(bt_array):

	bt_sim = bt_array * 0.8643254
	bt_sim = bt_sim + 38.06518136
	sst_array = bt_sim *  1.07763576
	sst_array = sst_array - 20.1082804
	return sst_array
	
def torness_bt_to_sst(bt_array):

	bt_sim = bt_array * 0.84219701
	bt_sim = bt_sim + 44.27977808
	sst_array = bt_sim * 1.00370779
	sst_array = sst_array + 0.49992172
	return sst_array

#--------------------------------------------------------------------

# read in the netcdf file, file given as argument in the command line	

"""
filename = sys.argv[1]	
date = filename[21:29]


#dataset = nc.Dataset(filename)
nc_file = nc.Dataset(filename, 'r+' , format="NETCDF4")

BT10 = np.array(nc_file.variables['BT_band10'])
BT11 = np.array(nc_file.variables['BT_band11'])
#pick your variables for the mask and for BT
refl_band3 = np.array(nc_file.variables['reflectance_band3'])
refl_band6 = np.array(nc_file.variables['reflectance_band6'])
# Get the MNDWI
aMNDWIMask = wdm.landmask(refl_band3, refl_band6, BT10)
# Apply the mask to the TIR data.
BT10_masked = np.ma.array(BT10, mask=aMNDWIMask, fill_value=np.nan)
BT11_masked = np.ma.array(BT11, mask=aMNDWIMask, fill_value=np.nan)

# Convert DN to BT to SST.
SST = heysham_bt_to_sst(BT_band10_masked)
SST_Celsius = SST - 273.15
SST_Celsius[SST_Celsius > 350] = np.nan
my_cmap = mpl.cm.coolwarm
my_cmap.set_under(color='white')
plt.imshow(SST_Celsius, cmap=my_cmap)
plt.colorbar()
plt.clim(5, 15)
plt.title('Heysham ' + str(date))
plt.show()
#plt.savefig('/home/users/mp877190/getting_netcdf/SST_plots/Sizewell/Sizewell__'+str(date)+'.png')
"""

