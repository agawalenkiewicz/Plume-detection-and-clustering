import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
import os, re
from datetime import datetime
import matplotlib

import Landsat_ncfiles
import water_detection_mask as wdm
import make_RGB
import bt_to_sst_double_regression as sst_regres




def geo_idx(dd, dd_array):
   """
     search for nearest decimal degree in an array of decimal degrees and return the index.
     np.argmin returns the indices of minium value along an axis.
     so subtract dd from all values in dd_array, take absolute value and find index of minium.
    """
   geo_idx = (np.abs(dd_array - np.float(dd))).argmin()
   return geo_idx

def stack_water_detection(filename, path):
	"""
	Input:
	filename - list of netcdf files with MNDWI mask
	path - absolute path to where the netcdf files are stored
	---
	loops through all files in the given list, finds the MNDWI masked layer
	for each layer gives value of 1 if land, 0 is water
	stacks the masks together as a 3D array
	---
	Output: 3D array of 0-1 masks
	"""
	#result_array = np.zeros((len(filename), 401, 401))
	result_array = np.empty((len(filename), 401, 401))
	dates_array = []
	for i, element in enumerate(filename):
		#date_of_the_satellite_obs = element[21:29]
		element_path = os.path.join(path, element)
		nc_file = nc.Dataset(element_path)
		rgb_layer = make_RGB.file_to_rgb(nc_file)
		BT_masked = masked_bt(nc_file)
		my_cmap = matplotlib.cm.coolwarm #BuPu #seismic #magma #coolwarm
		my_cmap.set_over(color='white') #only to get the RGB landmask
		my_cmap.set_under(color='white') #only to get the RGB landmask
		plt.imshow(np.rot90(rgb_layer))
		plt.imshow(np.rot90(BT_masked), cmap=my_cmap, vmax=20, vmin=19)
		#plt.show()

		#BT_masked = np.array(nc_file.variables['BT10_masked'])
		BT_masked[BT_masked > 350] = np.nan 
		BT_masked[BT_masked < 273] = np.nan
		BT_masked[BT_masked.mask == True] = np.nan
		result_array[i,:,:] = BT_masked
		print "results array" , result_array[i,:,:]
		date = re.findall(r"\d{8}", str(element)[50:61]) #[46:54]) 
		#print date
		#print str(element)[50:61] #[46:54] 
		value = str(date[0].strip().strip("'"))
		#print value
		dates_array.append(value)
		
		#plt.imshow(np.rot90(result_array[i,:,:]))
		#plt.colorbar()
		#plt.clim()
		#plt.title('Sizewell  ' + value)
		#plt.show()
		
	#print np.shape(result_array)
	return result_array, dates_array

def masked_bt(nc_file):	
	BT10 = np.array(nc_file.variables['BT_band10'])
	#pick your variables for the mask and for BT
	refl_band3 = np.array(nc_file.variables['reflectance_band3'])
	refl_band6 = np.array(nc_file.variables['reflectance_band6'])
	# Get the MNDWI
	aMNDWIMask = wdm.landmask(refl_band3, refl_band6, BT10)
	# Apply the mask to the TIR data.
	BT_masked = np.ma.array(BT10, mask=aMNDWIMask, fill_value=np.nan)
	#np.ma.set_fill_value(BT_masked, np.nan)
	return BT_masked

def stat_landmask(meanBT):
	"""
	Input: 3D array of 0-1 masks
	---
	sums the value of each pixel in the vertical (through all stacked layers)
	created a 2D layer with summed pixel values
	to normalize the output, it calculates the 95th percentile from the sum
	95th percentile is now the threshols value
	where pixel value above thershold - land
	where pixel values below threshold - water
	masks all pixels that were treated as land
	---
	Output: landmask based on 95th percentile of summed pixel values
	"""
	#meanBT[np.isnan(meanBT)] = np.inf
	#sum = meanBT.sum(axis=0)
	#sumPercentail = np.percentile(sum, 99.)
	#sum[sum >= sumPercentail] = 1
	#sum[sum == np.inf] = 1
	#sum[sum > 999] = 1
	meanBT[np.isnan(meanBT)] = 1
	# Apply the mask to remove land pixels from the image.
	sum = np.ma.masked_where(meanBT == 1, meanBT)
	Landmask = np.ma.getmaskarray(sum)
	return Landmask

def avg_BT(result_array):	
	"""
	kkk
	"""
	meanBT = np.nanmean(result_array, axis=0) #sum(axis=0)
	#sumBT /= np.float(len(filename))
	return meanBT
	

def plot_stat_mask(filepath, Landmask):
	"""
	Input: path to netcdf file we want to process, 95th percentile landmask
	----
	choose the variable from the netcdf file (SST, BT)
	apply a landmask over the chosen variable
	plot a 2D image with the landmask active
	---
	Output: shwoing the 2D masked plot
	"""
	data_path = filepath
	data = nc.Dataset(data_path)
	#choose BT at 10.8 micron
	variable = np.array(data.variables['BT_band10'])
	variable = np.ma.array(variable, mask=Landmask)
	plt.imshow(np.rot90(variable))
	plt.colorbar()
	plt.clim(275,300)
	#plt.title('Heysham 01/10/2015 from ASTER \n with Landsat 8 landmask')
	#plt.show()

def centered_average(nums):
	#print np.nansum(nums)
	#print np.nanmax(nums)
	#print np.nanmin(nums)
	#print len(nums) - 2
	return (np.nansum(nums) - np.nanmax(nums) - np.nanmin(nums)) / (np.count_nonzero(~np.isnan(nums)) - 2) 

######################################################
######################################################
######################################################
"""
# if you want different place, change it after the dot
filename1 = Landsat_ncfiles.heysham_10_50
path1 = Landsat_ncfiles.heysham_10_50_path
masked_stack1, dates_stack1 = stack_water_detection(filename1, path1)

filename2 = Landsat_ncfiles.heysham
path2 = Landsat_ncfiles.heysham_path
masked_stack2, dates_stack2 = stack_water_detection(filename2, path2)

masked_stack = np.concatenate((masked_stack1, masked_stack2),axis=0)

max = np.empty(len(masked_stack))
cent_avg = np.empty(len(masked_stack))
for i, layer in enumerate(masked_stack):
	#print layer
	#print np.count_nonzero(~np.isnan(layer))
	layer_max = np.nanmax(layer)
	layer_avg = centered_average(layer)
	max[i] = layer_max
	cent_avg[i] = layer_avg
	#print "inflow" , layer[191, 194]
	#print "outflow" , layer[194, 197]
print "max" , max
print "centered average" , cent_avg


#averaged_BT = np.ma.array(avg_BT(masked_stack), mask=landmask)
averaged_BT = avg_BT(masked_stack)
averaged_BT[averaged_BT <= 280 ] = np.nan
my_cmap = matplotlib.cm.BuPu
my_cmap.set_bad(color='khaki')
#plt.imshow(np.rot90(averaged_BT), cmap=my_cmap)
#plt.title('Sizewell averaged BT for 10.8 micron \n for 10-50% cloud cover for all scenes')
#plt.colorbar()
#plt.clim(284,289)
#plt.show()

landmask = stat_landmask(averaged_BT)
print landmask


for i, layer in enumerate(masked_stack): 
	meanBT = np.nanmean(layer)
	print meanBT
	layer_amplitude = layer - meanBT
	layer_amplitude[np.isnan(layer_amplitude)] = np.nan #-999 #this is a new added line
	print layer_amplitude
	layer_amplitude[landmask] = np.nan #999
	my_cmap = matplotlib.cm.hot_r #coolwarm #BuPu #seismic #magma #coolwarm
	my_cmap.set_bad(color='tan', alpha=1.0) #olive #black #khaki #silver
	#my_cmap.set_over(color='khaki') #black #olive
	my_cmap.set_under(color='khaki') #khaki #grey

	#levels = [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0] #[2.0] #[0.0, 1.0, 2.0, 3.0]
	#contours = plt.contour(np.rot90(layer_amplitude), levels, colors='black')
	#plt.clabel(contours, inline=True, fmt='%1.1f', fontsize=10)
	#contour_filled = plt.contourf(np.rot90(layer_amplitude), levels, extend='both', cmap='hot_r', alpha=0.5) #cmap='coolwarm'
	#plt.colorbar(contour_filled)
	calculated_sst = sst_regres.heysham_bt_to_sst(layer_amplitude)
	SST_Celsius = calculated_sst - 273.15
	SST_Celsius = np.rot90(SST_Celsius)
	max_t = np.nanmax(SST_Celsius)
	min_t = np.nanmin(SST_Celsius)
	plt.imshow(SST_Celsius, cmap= my_cmap, vmax=max_t, vmin=min_t)
	#plt.imshow(np.rot90(SST_Celsius), cmap= my_cmap, vmax=max_t, vmin=min_t) #alpha=0.5
 	##plt.imshow(np.rot90(layer_amplitude),cmap= my_cmap, vmax = 3, vmin=-2)
	##plt.contour(np.rot90(layer_amplitude), 5, cmap='RdGy')
	plt.colorbar()
	plt.clim()
	#plt.show()
	#plt.savefig('/glusterfs/surft/users/mp877190/scenes/dungeness_plume_' +str(i) + '.png')
	plt.clf()
	##crash

#exampleFilePath = os.path.join(path, filename[-1])
exampleFilePath = '/glusterfs/surft/users/mp877190/data/datastore/EE/ASTER_L1T/heysham/2194252766_tir/AST_L1T_00309302015215121_20151001103139_9227.nc'
#plot_stat_mask(exampleFilePath, landmask)
"""