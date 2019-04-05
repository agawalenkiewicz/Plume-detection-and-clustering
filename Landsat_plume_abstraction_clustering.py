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

import landmask_95percentile as lm95
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

def subimage(image_as_array, step):
	"""
	Get a sub-image from original image.
	Sub-image and original image share the centre.
	The outer bounds of sub-image are specified by the step_count_from_centre argument.
	step - step size from the centre expressed as number of pixels
	"""
	subimage_2d_array = image_as_array[200-int(step):200+int(step)]
	return subimage_2d_array
	
def find_max(subimage):
	"""
	Find the maximum value of the smaller subimage
	concentrated on the power ststaion location.
	"""
	max_val_subimage = np.nanmax(subimage)
	return max_val_subimage
	
def morphological_dilation(masked_image, n): #n=3
	"""
	Extending the landmask.
	Should extend the landmask over the image for 3 pixels (0.0015 degrees)
	----------
	from stackoverflow:
	def dilation(a, n):
	m = np.isnan(a)
	s = np.full(n, True, bool)
	return ndimage.binary_dilation(m, structure=s, origin=-(n//2))
	-----------
	For sparse initial masks and small n this one is also pretty fast:
	def index_expansion(a, n):
	mask = np.isnan(a)
	idx = np.flatnonzero(mask)
	expanded_idx = idx[:,None] + np.arange(1, n)
	np.put(mask, expanded_idx, True, 'clip')
	return mask
	"""
	mask = np.isnan(masked_image)
	s = ndimage.morphology.generate_binary_structure(2, 1)
	extended_mask = ndimage.binary_dilation(mask, structure=s, iterations=3).astype(mask.dtype)
	return extended_mask
	#mask = np.isnan(masked_image)
	#idx = np.flatnonzero(mask)
	#expanded_idx = idx[:,None] + np.arange(1, n)
	#np.put(mask, expanded_idx, True, 'clip')
	#return mask

def detect_peaks(image):
    """
    Takes an image and detect the peaks usingthe local maximum filter.
    Returns a boolean mask of the peaks (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """

    # define an 8-connected neighborhood
    neighborhood = ndimage.morphology.generate_binary_structure(2,2)

    #apply the local maximum filter; all pixel of maximal value 
    #in their neighborhood are set to 1
    local_max = ndimage.filters.maximum_filter(image, footprint=neighborhood)==image
    #local_max is a mask that contains the peaks we are 
    #looking for, but also the background.
    #In order to isolate the peaks we must remove the background from the mask.

    #we create the mask of the background
    background = (image==0)

    #a little technicality: we must erode the background in order to 
    #successfully subtract it form local_max, otherwise a line will 
    #appear along the background border (artifact of the local maximum filter)
    eroded_background = ndimage.morphology.binary_erosion(background, structure=neighborhood, border_value=1)

    #we obtain the final mask, containing only peaks, 
    #by removing the background from the local_max mask (xor operation)
    detected_peaks = local_max ^ eroded_background

    return detected_peaks

	
def centered_average(nums):
	#print np.nansum(nums)
	#print np.nanmax(nums)
	#print np.nanmin(nums)
	#print len(nums) - 2
	return (np.nansum(nums) - np.nanmax(nums) - np.nanmin(nums)) / (np.count_nonzero(~np.isnan(nums)) - 2)	


def choose_plume(image_thresholded):
	
	#now find the objects
	where_are_NaNs = np.isnan(image_thresholded)
	image_thresholded[where_are_NaNs] = 0
	
	labeled_image, numobjects = ndimage.label(image_thresholded)
	
	#plt.imshow(labeled_image)
	#plt.show()
	
	object_areas = np.bincount(labeled_image.ravel())[:]
	#to exclude the first object which is background , index from 1
	object_idx = [i for i in range(1, numobjects) if object_areas[i] > 2]
	print('object area' , object_areas)
	print('object idx' , object_idx)
	# Remove small white regions
	#labeled_image = ndimage.binary_opening(labeled_image)
	# Remove small black hole
	#labeled_image = ndimage.binary_closing(labeled_image)
	chosen_object = [0,50]
	for object in object_idx: #range(0,numobjects):
		#object = object + 1
		#print('object' , object)
		iy, ix = np.where(labeled_image == object)
		centridx_y = 200
		centridx_x = 190 #200
		min_dist = np.min(np.sqrt((centridx_y - iy)**2 + (centridx_x - ix)**2))
		if min_dist < chosen_object[1]:
			chosen_object = [object, min_dist]
		#print(object , min_dist)
	#print('Chosen object' , chosen_object)
	#chosen_plume = np.where((labeled_image == chosen_object[0]), chosen_object[0], 0)
	if chosen_object[1] == 50:
		chosen_plume = np.zeros_like(labeled_image)
	else:
		chosen_plume = np.where((labeled_image == chosen_object[0]), 1, 0)
	area = sum(sum(i == True for i in chosen_plume))
	print "Detected plume area (number of pixels):" , area
	return chosen_plume
	
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------	
def Landsat_plume_main():
	filename1 = Landsat_ncfiles.hinkley
	path1 = Landsat_ncfiles.hinkley_path
	masked_stack1, dates_stack1 = lm95.stack_water_detection(filename1, path1)

	filename2 = Landsat_ncfiles.hinkley_10_50
	path2 = Landsat_ncfiles.hinkley_10_50_path
	masked_stack2, dates_stack2 = lm95.stack_water_detection(filename2, path2)

	masked_stack = np.concatenate((masked_stack1, masked_stack2),axis=0)
	dates_stack = dates_stack1 + dates_stack2

	print dates_stack

	averaged_BT = lm95.avg_BT(masked_stack)
	landmask = lm95.stat_landmask(averaged_BT)

	
	plume_array = np.zeros_like(masked_stack)
	for i, layer in enumerate(masked_stack): 
		date_in_image_name = dates_stack[i]
		# make sure that the values with statistical landmask are set to NaN
		layer[landmask] = np.nan
		"""
		plt.imshow(np.rot90(layer))
		plt.title('SIZEWELL \n Brightness temperature')
		plt.colorbar()
		plt.show()
		"""
		# expand land mask by morphological dilation
		#extended_landmask = morphological_dilation(layer, 20)
		#layer[extended_landmask] = np.nan

		# depending on location convert BT to SST
		calculated_sst = sst_regres.heysham_bt_to_sst(layer)
		SST_Celsius = calculated_sst - 273.15
		SST_Celsius = np.rot90(SST_Celsius)
		meanSST = np.nanmean(SST_Celsius)
		#print meanSST
		SST_amplitude = SST_Celsius - meanSST
		"""
		my_cmap = matplotlib.cm.BuPu #hot_r #coolwarm #BuPu #seismic #magma #coolwarm
		my_cmap.set_bad(color='khaki', alpha=1.0) #olive #black #khaki #silver
		plt.imshow(SST_amplitude, cmap= my_cmap, vmax=3, vmin=-3)
		plt.title('Hartlepool SST during ebb'+' '+str(date_in_image_name))
		plt.colorbar()
		plt.show()
		#plt.savefig('/home/users/mp877190/getting_netcdf/SST_plots/Hartlepool/Hartlepool__'+str(date_in_image_name)+'.png')
		plt.clf()
		"""
		
		# Create a Rectangle patch
		#rect = patches.Rectangle((250,0),150,150,linewidth=1,edgecolor='k',facecolor='none')
		#ax.add_patch(rect)
		#plt.show()

		# get a central part of the image
		sub_layer = subimage(SST_Celsius, 25)
		max_val = find_max(sub_layer)
		#print('max value' , max_val)
		ambient = centered_average(SST_Celsius[150:300, 250:400]) #[0:150, 250:400])
		accepted_thresh = np.float(ambient) + 1.0 #np.float(max_val) -  np.float(ambient)
		#print('accepted threshold' , accepted_thresh)
		threshold   = accepted_thresh
		image_thresh = np.copy(SST_Celsius)
		image_thresh[image_thresh<threshold] = np.nan

		#plt.imshow(image_thresh)
		#plt.colorbar()
		#plt.title('Image thresholded')
		#plt.show()
		
		plume = choose_plume(image_thresh)
		
		#plt.imshow(plume)
		#plt.title('Detected plume' + str(date_in_image_name))
		#plt.show()
		
		print "Date of the Landsat image" , date_in_image_name
		
		plume_array[i,:,:] = plume

		end_result = np.nansum(plume_array, axis=0)
		end_result = (end_result / np.amax(end_result)) * 100
	return plume_array

	

"""
## PLOT THE PROBABILITY DENSITY MAP 
my_cmap = matplotlib.cm.hot_r 
my_cmap.set_under(color='paleturquoise') #, alpha=0.5)
plt.imshow(end_result, cmap= my_cmap, vmin=1)
plt.colorbar()
plt.title('HEYSHAM \n Probability of plume extent during ebb \n with waters warmer than ambient by 1.0 degree')
plt.tight_layout()
plt.show()
"""




"""
#-- Plotting... -------------------------------------
fig, ax = plt.subplots()
ax.imshow(image_thresh)
ax.set_title('Original Data')

fig, ax = plt.subplots()
ax.imshow(labeled_image)
ax.set_title('Labeled objects')

fig, axes = plt.subplots(ncols=numobjects)
for ax, sli in zip(axes.flat, slices):
	ax.imshow(labeled_image[sli], vmin=0, vmax=numobjects)
	tpl = 'BBox:\nymin:{0.start}, ymax:{0.stop}\nxmin:{1.start}, xmax:{1.stop}'
	ax.set_title(tpl.format(*sli))
fig.suptitle('Individual Objects')

plt.show()
"""
