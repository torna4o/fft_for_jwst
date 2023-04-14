# FFT for JWST
Processing JWST images with FFT - 

This repository is intended for putting all FFT-related convenience functions in my research regarding JWST imaging preprocessing and source finder.

You will need numpy, os, pathlib, and astropy packages. These can work even in Windows 10 environment, and also in Google Colab.

"fft_astro.py"

This mainly stores the two FFT function codes I used in the research note I published as a preprint in arXiv reporistory:
https://arxiv.org/abs/2304.00728
You may just keep the functions and remove the rest, as the lines below the function is just for the people who would like to replicate the results.

"gaia_query.py"

A code to easily query the region around the center of the FITS image one feeds to the function. In addition to GAIA, it has functions for querying VizieR as well, and then plotting the results as a scatter plot and overlaying them with the image at the beginning.

"miripsf.py"

This code quickly prints a Basic Photometry results with default and available methods without background subtraction. Only photutils and astropy packages are sufficient.

"quick_fft.py"

Quickly applying a standard FFT filter to an image provided. It applies a 512,512 cutout around a specified point in the image, then with a threshold radius circle, separates the data into two according to its location in the frequency domain. 

"welchfft.py"

After several trial-error and look-up, Welch's method was by far the most successful in finding and removing the vertical and horizontal lines in JWST MIRI images. Hence, a code is provided here to be able to quickly apply this method. One example application to JWST MIRI images is here:
https://arxiv.org/abs/2304.00728

