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

import coral_reef_files


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
		print('object' , object)
		iy, ix = np.where(labeled_image == object)
		centridx_y = 200
		centridx_x = 200
		min_dist = np.min(np.sqrt((np.abs(centridx_y - iy))**2 + (np.abs(centridx_x - ix))**2))
		if min_dist < chosen_object[1]:
			chosen_object = [object, min_dist]
		print(object , min_dist)
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


	filename1 = Landsat_ncfiles.hinkley_ebb
	path1 = Landsat_ncfiles.hinkley_path
	masked_stack1 = lm95.stack_water_detection(filename1, path1)

	filename2 = Landsat_ncfiles.hinkley_10_50_ebb
	path2 = Landsat_ncfiles.hinkley_10_50_path
	masked_stack2 = lm95.stack_water_detection(filename2, path2)

	filename3 = Landsat_ncfiles.hinkley_2018_ebb
	path3 = Landsat_ncfiles.hinkley_2018_path
	masked_stack3 = lm95.stack_water_detection(filename3, path3)

	masked_stack = np.concatenate((masked_stack1, masked_stack2, masked_stack3),axis=0)
	#dates_stack = dates_stack1 + dates_stack2 + dates_stack3
	
	"""
	# Plot a scatter of differences between sst from oe and observed bt10
	# based on matched scenes
	# from all three directories concatenated
	diff10_array = diff10_array1 + diff10_array2 #+ diff10_array3
	diff11_array = diff11_array1 + diff11_array2 #+ diff11_array3
	tcwv_array = tcwv_array1 + tcwv_array2 #+ tcwv_array3
	plt.scatter(diff10_array, tcwv_array, c='r', label='sst - bt10.8')
	plt.scatter(diff11_array, tcwv_array, c='b', label='sst - bt12')
	plt.title('Differences between SST and BT for different SST \n Torness')
	plt.xlabel('Difference [K]')
	plt.ylabel('SST prior [K]') #('TCWV for that scene [kg/m2]')
	plt.legend(loc='upper left')
	plt.show()
	"""
	#print dates_stack

	averaged_BT = lm95.avg_BT(masked_stack)
	landmask = lm95.stat_landmask(averaged_BT)

	
	plume_array = np.zeros_like(masked_stack)
	for i, layer in enumerate(masked_stack): 
		#date_in_image_name = dates_stack[i]
		# make sure that the values with statistical landmask are set to NaN
		layer[landmask] = np.nan
		# expand land mask by morphological dilation
		extended_landmask = morphological_dilation(layer, 20)
		layer[extended_landmask] = np.nan
		
		# depending on location convert BT to SST
		calculated_sst = sst_regres.hinkley_bt_to_sst(layer)
		SST_Celsius = calculated_sst - 273.15
		SST_Celsius = np.rot90(SST_Celsius)
		meanSST = np.nanmean(SST_Celsius)
		#print meanSST
		
		# Create a Rectangle patch
		#rect = patches.Rectangle((250,0),150,150,linewidth=1,edgecolor='k',facecolor='none')
		#ax.add_patch(rect)
		#plt.show()

		# get a central part of the image
		sub_layer = subimage(SST_Celsius, 25)
		max_val = find_max(sub_layer)
		#print('max value' , max_val)
		
		
                # indexing into ambient - first row (so y or vertical) then column (so x or horizontal)
		ambient = centered_average(SST_Celsius[0:150, 250:400]) #[0:150, 250:400])
		accepted_thresh = np.float(ambient) + 1.5 #np.float(max_val) -  np.float(ambient)
		
		#accepted_thresh = np.float(centered_average(SST_Celsius)) + (3.0 * np.nanstd(SST_Celsius))

		threshold   = accepted_thresh
		image_thresh = np.copy(SST_Celsius)
		image_thresh[image_thresh<threshold] = np.nan

		#plt.imshow(image_thresh)
		#plt.colorbar()
		#plt.title('Image thresholded')
		#plt.show()
		
		plume = choose_plume(image_thresh)
		SST_amplitude = SST_Celsius - ambient
		temp_val = SST_amplitude * plume
		
		#plt.imshow(plume)
		#plt.title(i)
		#plt.show()

		plume_array[i,:,:] = plume #temp_val

		end_result = np.nansum(plume_array, axis=0)
		end_result = (end_result / np.amax(end_result)) * 100 #(end_result / len(end_result))
	#plt.imshow(end_result)
	#plt.colorbar()
	#plt.clim()
	#plt.title('Delta temp warmer than amabient')
	#plt.show()
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
