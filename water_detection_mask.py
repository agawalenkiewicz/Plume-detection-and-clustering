import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
import os



def landmask(aB3, aB6, aB11):
    """
    Return an array that can be used to mask land pixels.
    """

    # Covert the DN into reflectance values.
    aB3Reflectance = aB3
    aB6Reflectance = aB6

    # Apply the mask from band 11 to remove black pixels from the image.
    aB11Mask = np.ma.getmaskarray(aB11)
    aB3Reflectance = np.ma.array(aB3Reflectance, mask=aB11Mask)
    aB6Reflectance = np.ma.array(aB6Reflectance, mask=aB11Mask)

    # Calculate the MNDWI to find the water.
    aLandmask = mndwi(aB3Reflectance, aB6Reflectance, fCutoff=0.8)
    return aLandmask

def mndwi(aB3, aB6, fCutoff=0.8):
    """
    Modification of Normalized Difference Water Index
    The B3 and B6 arrays need to be reflectance values.
    High values (over cutoff) are water, low values are land.
    """
	
    aMNDWI = (aB3 - aB6) / (aB3 + aB6)

    max = np.amax(aMNDWI)
    min = np.amin(aMNDWI)
    # Normalise so ranges from 0-1
    aMNDWI -= min
    aMNDWI /= max - min

    # Remove the points that lie outside of the range.
    aMNDWIMask = np.ma.masked_outside(aMNDWI, fCutoff , 1.0)
    aMNDWIMask = np.ma.getmaskarray(aMNDWIMask)
    return aMNDWIMask

"""
# read in the netcdf file, file given as argument in the command line	
filename = sys.argv[1]	
#nc_file = nc.Dataset(filename)
dataset = nc.Dataset(filename, 'r+' , format="NETCDF4")


#pick your variables for the mask and for BT
refl_band3 = dataset.variables['reflectance_band3']
refl_band6 = dataset.variables['reflectance_band6']

BT_band10 = dataset.variables['BT_band10']
BT_band11 = dataset.variables['BT_band11']

# Get the land mask.
aLandmask = landmask(refl_band3, refl_band6, BT_band10)

# Apply the mask to the TIR data.
BT_band10_masked = np.ma.array(BT_band10, mask=aLandmask)
BT_band11_masked = np.ma.array(BT_band11, mask=aLandmask)


plt.imshow(BT_band10_masked)
plt.show()

############
#Get the masked BT 10 and BT 11 as new variables into existing netcdf
############
lat = dataset.dimensions['lat'] 
lon = dataset.dimensions['lon']
# Create the actual 4-d variable
BT10_masked = dataset.createVariable('BT10_masked', np.float32, ('lat','lon'))
BT11_masked = dataset.createVariable('BT11_masked', np.float32, ('lat','lon'))
#bands.units = 'micro_meters'
BT10_masked.units = 'K'
BT11_masked.units = 'K'
#Assign values to the variables
BT10_masked[:,:] = BT_band10_masked
BT11_masked[:,:] = BT_band11_masked
#Close dataset
dataset.close()
"""