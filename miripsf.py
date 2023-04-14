from photutils.psf import BasicPSFPhotometry
from photutils.detection import IRAFStarFinder

from photutils.psf import IntegratedGaussianPRF, DAOGroup
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.modeling.fitting import LevMarLSQFitter
from photutils.background import MADStdBackgroundRMS

def defpsf(img, sigma_psf = 2.0):

    """
    This is a default Basic Photometry to quickly
    try new images PSF results.
    """
    
    bkgrms = MADStdBackgroundRMS() # Background RMS using Median Absolute Deviation
    std = bkgrms(img)
    print("Median Absolute Deviation: ", std)
    daogroup = DAOGroup(2.0 * sigma_psf * gaussian_sigma_to_fwhm)
    my_finder = IRAFStarFinder(threshold=3.5*std, fwhm=sigma_psf * gaussian_sigma_to_fwhm)
    my_psf_model = IntegratedGaussianPRF(sigma=sigma_psf)
    my_fitter = LevMarLSQFitter()

    my_photometry = BasicPSFPhotometry(finder=my_finder,
                                   group_maker = daogroup,
                                   bkg_estimator = None,
                                   psf_model = my_psf_model,
                                   fitter = my_fitter,
                                   fitshape=(7,7))
    return my_photometry(image=img), my_photometry.get_residual_image()
