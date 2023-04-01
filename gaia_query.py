##################################################
##     ##     ##    ## ### ##     ##    ###     ##
#### #### ### ## ## ##  ## ## ### ## ### ## ### ##
#### #### ### ##   ### # # ## ### ## ### ## ### ##
#### #### ### ## ## ## ##  ##     ## ### ## ### ##
#### ####     ## ## ## ### ## ### ##    ###     ##
##################################################

#####################################
#####################################
#                                   #
#                                   #
#    Loading Required Libraries     #
#                                   #
#                                   #
#####################################
#####################################



import os # Operating System package
import pathlib # Easier and automatic file querying and retrieval

from astroquery.gaia import Gaia # To Query Gaia data releases
from astroquery.vizier import Vizier # To query catalogs in VizieR
from astropy.wcs import WCS # World Coordinate System wrapper and manipulator
from astropy.io import fits # FITS file opener, knowledge extractor and manipulator

from astropy.coordinates import SkyCoord # Celestial coordinate repr. interface
from astropy.coordinates import ICRS # A common choice of coordinate ref. system
import astropy.units as u # ensuring standard unit to represent custom lengths

import matplotlib.pyplot as plt # plotting

#####################################
#####################################
#                                   #
#                                   #
#           Data Loading            #
#                                   #
#                                   #
#####################################
#####################################

work_dir = pathlib.Path(os.getcwd()) # This work_dir object facilitates the later operations
i2d_files = list(work_dir.glob("*_f1000w_i2d.fits")) # To query JWST MIRI F1000W images
i2d_files = [os.path.basename(i2d) for i2d in i2d_files]

# Opening FITS file via "with" to ensure its closure after extracting data and info

with fits.open(i2d_files[0]) as hdu_fits: 
    image_data = hdu_fits['SCI'].data # this is a numpy nd array
    wcs_helix = WCS(hdu_fits[1].header) # WCS object from FITS

# crval[0] and [1] gives the center of the FITS file.

coord = SkyCoord(ra=wcs_helix.wcs.crval[0], dec=wcs_helix.wcs.crval[1], unit="deg", frame="icrs")
exp_coeff = 0.7 # A multiplier of width and height for querying
width = u.Quantity(0.1*exp_coeff, u.deg)
height = u.Quantity(0.1*exp_coeff, u.deg)

#####################################
#####################################
#                                   #
#                                   #
#         Querying Catalogs         #
#                                   #
#                                   #
#####################################
#####################################

"""
Catalog search may result in tremendous amount of rows and columns.
Thus, there are "max_lines" and "max_width" parameters to tweak the output.
"""

print("|||| ---- |||| Querying GAIA Data Release III |||| ---- ||||")

# Default GAIA data source, Gaia.MAIN_GAIA_TABLE, is Gaia DR III.

r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
r.pprint(max_lines=12, max_width=130)

Gaia.MAIN_GAIA_TABLE = "gaiadr2.gaia_source"  # Select Data Release 2

print("|||| ---- |||| Querying GAIA Data Release II |||| ---- ||||")

r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
r.pprint(max_lines=12, max_width=130)

catalog_k=["VII/233"]

print("|||| ---- |||| Querying VizieR |||| ---- ||||")

"""
List of catalog IDs can be found in the following webpage:
https://cdsarc.cds.unistra.fr/cgi-bin/cat
The following are some examples
"""

catalog_k = 'II/363/unwise'
catalog_2mass = 'II/246/out'
catalog_most = 'B/gemini/obscore'
cat_ssts ='II/368/sstsl2'

Vizier.ROW_LIMIT = -1 # alters the 50 rows default limit

# Unspecified "catalog" parameter in Vizier.query_region will list from all catalogs.

r2 = Vizier.query_region(coord, width=width, catalog = [catalog_k, catalog_2mass, cat_ssts])
print(r2)

"""
The following parameters retrieves the names and coordinates from query result dict to a list.

RAJ2000 is Right Ascension for J2000 in degree.decimal
DEJ2000 is Declination for J2000 in degree.decimal
"""

ralist = r2[cat_ssts]['RAJ2000'].tolist()
declist = r2[cat_ssts]['DEJ2000'].tolist()
names = r2[cat_ssts]['Name'].tolist()

#####################################
#####################################
#                                   #
#                                   #
#        Overlaying Catalogs        #
#                                   #
#                                   #
#####################################
#####################################

ax = plt.subplot(projection=wcs_helix) # subplot with specific projection

ax.imshow(hdu_sci, origin='lower', cmap='afmhot', alpha=1, vmin=17.5, vmax=18.9)

overlay = ax.get_coords_overlay('fk5') #Overlaying ICRS coordinates on counter axes
overlay.grid(color='red', ls='dotted')
overlay[0].set_axislabel('Right Ascension (J2000)')
overlay[1].set_axislabel('Declination (J2000)')

# The following line puts objects in the selected catalogs on image
ax.scatter(ralist,declist,transform=ax.get_transform('fk5'), s=15, edgecolor='red')
"""
The following line put name labels on catalog objects drawn on the figure
    va is vertical alignment of the label, and ha is horizontal
    ax.get_transform('fk5') permits coinciding the coordinate systems of the Figure
    and scatter plots
"""
for (xi, yi, ni) in zip(ralist, declist, names):
    plt.text(xi, yi, ni, va='bottom', ha='center', fontsize=5, color= 'yellow', transform=ax.get_transform('fk5'))
plt.show()
