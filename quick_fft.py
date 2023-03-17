# Extracting from FITS and putting it back to it later

"""
Here I will define three modules here to:
1-] extracts ['SCI'] image data from a FITS file
2-] makes an FFT with a "threshold_radius" parameter
3-] put processed image (FFT-filtered) data back into ['SCI']
in a new file with _fft.fits suffix.
"""

#####################################
#####################################
#                                   #
#                                   #
#      1. Extracting SCI image      #
#                                   #
#                                   #
#####################################
#####################################

# Loading necessary libraries

import os # operating system package
import pathlib # custom listing of files in a folder
from astropy.io import fits # to open FITS file
from astropy.nddata.utils import Cutout2D # cutting parts from FOV from numpy arrays

# Below, files are considered to be in the work directory
# if not, modify "work_dir" address

# Custom listing of the files to retrieve list of FITS
work_dir = pathlib.Path(os.getcwd())
i2d_files = list(work_dir.glob("*_f1000w_i2d.fits")) # resampled i2d files in F1000W filter
i2d_files = [os.path.basename(i2d) for i2d in i2d_files] # getting base file name

# Retrieving Header Data Unit of FITS
with fits.open(i2d_files[0]) as hdu_fits: 
    image_data = hdu_fits['SCI'].data # this is a numpy nd array

print("### *** \\\ EXTRACTION TASK COMPLETE //// **** ###")

# Now cutting out the critical FOV parts we were looking for
# in two different 512 x 512 numpy.ndarrays for FFT

ready_for_fft_top = Cutout2D(image_data, (685,753), (512,512))
ready_for_fft_bottom = Cutout2D(image_data, (685,314), (512,512))

#####################################
#####################################
#                                   #
#      Fast Fourier Transform       #
#           to the image            #
#                                   #
#####################################
#####################################

# Loading necessary libraries
import numpy as np # Numerical Python, fft and other numpy array operations

# Apply FFT to the image data
img_fft = np.fft.fft2(ready_for_fft_bottom.data) # Cutout requires .data to point data

# Shift the zero-frequency component to the center of the spectrum
img_fft = np.fft.fftshift(img_fft)

# Set a threshold radius to filter high frequency components
threshold_radius = 1

# Create a mask to keep low frequency components
mask = np.zeros(ready_for_fft_bottom.shape)
mask_inv = np.ones(ready_for_fft_bottom.shape)
center = (mask.shape[0] // 2, mask.shape[1] // 2)

# np.ogrid creates a mesh for indexing
y, x = np.ogrid[:mask.shape[0], :mask.shape[1]] # y and x are arrays for indices

"""
Now we obtain mask after array operation
More precisely, "mask_area" is a True False matrix, 
it will show whether a specific indexed data will
be included in the final masked numpy array.
"""
mask_area = (x - center[1])**2 + (y - center[0])**2 <= threshold_radius**2
# The following assignment put 1 to "True" locations in "mask_area" mask

mask[mask_area] = 1
mask_inv[mask_area] = 0

# Apply the mask to the FFT output
img_fft_filtered = img_fft * mask
img_fft_filtered_inv = img_fft * mask_inv

# Shift the zero-frequency component back to the top-left corner of the spectrum
img_fft_filtered = np.fft.ifftshift(img_fft_filtered)
img_fft_filtered_inv = np.fft.ifftshift(img_fft_filtered_inv)

# Apply inverse FFT to get the filtered image
img_filtered = np.abs(np.fft.ifft2(img_fft_filtered)) # np.abs calculates absolute values
img_filtered_inv = np.abs(np.fft.ifft2(img_fft_filtered_inv)) # np.abs calculates absolute values

print("### *** \\\ FFT TASK COMPLETE //// **** ###")

#####################################
#####################################
#                                   #
#                                   #
#     3. Reingesting edited SCI     #
#                                   #
#                                   #
#####################################
#####################################

with fits.open(i2d_files[0]) as hdu_fits:
    hdu_fits['SCI'].data = img_filtered
    hdu_fits.writeto(i2d_files[0][:-5]+"_cutfft2.fits", overwrite=True)


with fits.open(i2d_files[0]) as hdu_fits2:
    hdu_fits2['SCI'].data = img_filtered_inv
    hdu_fits2.writeto(i2d_files[0][:-5]+"_cutfftinv2.fits", overwrite=True)
    
print("### *** \\\ REINGESTION TASK COMPLETE //// **** ###")