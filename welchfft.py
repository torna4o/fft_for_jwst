# FFT filtering of recognized pattern from JWST MIRI calib_3 level images


import os # operating system package
import matplotlib.pyplot as plt # plotting functionality
import pathlib # managing and filtering files in a directory
import numpy as np # Numerical Python

from astropy.io import fits # plotting fits images
from astropy.nddata.utils import Cutout2D # Cuts a part of the image from FITS

from scipy.ndimage import convolve1d # Required for the functions below
from scipy.signal import firwin, welch # firwin filter creation and Welch's method

def estimate_distortion_freq(image, min_frequency=1/25):
  """
  Estimates distortion frequency as spectral peak in vertical dimension.
  The idea is summing along one axis and dealing with as if dealing a time-series etc.
  The same will be done for horizontal axis below
  """
  f, pxx = welch(np.reshape(image, (len(image), -1), 'C').sum(axis=1))
  
  """
  Welch's method documentation
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html#r34b375daf612-1
  Welch's method scientific paper
  https://ieeexplore.ieee.org/document/1161901
  
  """
  pxx[f < min_frequency] = 0.0
  return f[pxx.argmax()]

def remove_lines(image, distortion_freq=None, num_taps=65, eps=0.025):
  """Removes horizontal line artifacts from scanned image.
  Args:
    image: 2D or 3D array.
    distortion_freq: Float, distortion frequency in cycles/pixel, or
      `None` to estimate from spectrum.
    num_taps: Integer, number of filter taps to use in each dimension.
    eps: Small positive param to adjust filters cutoffs (cycles/pixel).
  Returns:
    Denoised image.
    
  Source of the code:
  https://stackoverflow.com/questions/65480162/how-to-remove-repititve-pattern-from-an-image-using-fft/65482958#65482958
  
  """
  image = np.asarray(image, float) # makes the integer units float type
  if distortion_freq is None:
    distortion_freq = estimate_distortion_freq(image)

  # firwin() is a choice of filter creation for the following operations
  hpf = firwin(num_taps, distortion_freq - eps,
               pass_zero='highpass', fs=1)
  lpf = firwin(num_taps, eps, pass_zero='lowpass', fs=1)
  aaaa = convolve1d(convolve1d(image, hpf, axis=0), lpf, axis=1)
  """
  bb = np.fft.fft2(aaaa)
  bb = np.fft.fftshift(bb)
  plt.imshow(np.abs(bb), cmap='gray', vmin=0, vmax= 1)
  plt.show()
  """
  return image - convolve1d(convolve1d(image, hpf, axis=0), lpf, axis=1)


def remove_lines_h(image, distortion_freq=None, num_taps=65, eps=0.025):
  # Horizontal dimension counterpart of the function above
  image = np.asarray(image, float)
  if distortion_freq is None:
    distortion_freq = estimate_distortion_freq(image)

  hpf = firwin(num_taps, distortion_freq - eps,
               pass_zero='highpass', fs=1)
  lpf = firwin(num_taps, eps, pass_zero='lowpass', fs=1)
  return image - convolve1d(convolve1d(image, hpf, axis=1), lpf, axis=0)

def estimate_distortion_freq_h(image, min_frequency=1/25):
  # Horizontal dimension counterpart of the function above
  f, pxx = welch(np.reshape(image, (len(image), -1), 'C').sum(axis=0))
  pxx[f < min_frequency] = 0.0
  return f[pxx.argmax()]
  


# Locating .fits files
work_dir = pathlib.Path("####Your working directory####") # or os.getcwd()
i2d_files = list(work_dir.glob("*f770w_i2d.fits")) # this retrieve list of F770W, modify at your will
i2d_files = [os.path.basename(i2d) for i2d in i2d_files] # obtain the names
i2d_files[0] = "####Your working directory####" + i2d_files[0]

with fits.open(i2d_files[0]) as hdu_fits: 
    image_data = hdu_fits['SCI'].data # this is a numpy nd array of SCIENCE data in FITS

"""
The following lines create a 1024x1024 zero-padded science array
For the following Fourier-like methods.

According to the focused regions, you may modify them to other sizes
or dimensions.
"""
zeroo = np.zeros((1024,1024))

# Cutout demans the data, then a point location (here 725, 630), 
# and sizes of the margins of the rectangular, which centers the point location

dat = Cutout2D(image_data, (725,630), (780,560))
dat = dat.data # Cutout2D's SCIENCE data part.

zeroo[50:(50 + dat.shape[0]), 51:(51 + dat.shape[1])] = dat # To replicate study's analysis


aaa = estimate_distortion_freq(zeroo) # Remove vertical lines
bbbb = remove_lines(zeroo, distortion_freq=aaa, num_taps=65, eps=0.025)
eee = estimate_distortion_freq_h(zeroo) # Remove horizontal lines
ffff = remove_lines_h(zeroo, distortion_freq=eee, num_taps=65, eps=0.025)
ccc = estimate_distortion_freq_h(bbbb) # Remove first vertical, then horizontal lines
dddd = remove_lines_h(bbbb, distortion_freq=ccc, num_taps=65, eps=0.025)

# Plotting the results, in this case, for F770W filter of JWST MIRI

plt.subplot(221)
plt.imshow(zeroo[50:(50 + 780), 50:(50 + 560)], origin='lower', cmap='afmhot', vmin=5.95, vmax=6.3)
plt.title('Original Image F770W')
plt.subplot(222)
plt.imshow(bbbb[50:(50 + 780), 50:(50 + 560)], origin='lower', cmap='afmhot', vmin=5.3, vmax=6)
plt.title('Horizontal Strips Removed F770W')
plt.subplot(223)
plt.imshow(ffff[50:(50 + 780), 50:(50 + 560)], origin='lower', cmap='afmhot', vmin=5.3, vmax=6)
plt.title('Vertical Strips Removed F770W')
plt.subplot(224)
plt.imshow(dddd[50:(50 + 780), 50:(50 + 560)], origin='lower', cmap='afmhot', vmin=4.85, vmax=5.3)
plt.title('Horizontal, then Vertical Strips Removed F770W')
plt.show()