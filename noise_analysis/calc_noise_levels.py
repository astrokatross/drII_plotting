#!/usr/bin/python

import numpy as np
import os
import re
import glob
import astropy.io.fits as fits
from astropy.table import Table, Column
from astropy import wcs
from astropy import units as u
from astropy import coordinates
from astropy.coordinates import SkyCoord
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib import rc
from astropy.visualization.wcsaxes import Quadrangle
from matplotlib.colors import LogNorm
import healpy as hp 
from matplotlib.ticker import LogFormatter 
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import pandas as pd
import cartopy.crs as ccrs
from astropy.nddata import Cutout2D

aspect = 20
pad_fraction = 0.5
subplot_cols = 1
subplot_rows = 1
save_path = "/data/gleam_x/drII_plotting/plots"
data_path = "/data/gleam_x/processing/GLEAMX-DRII/"
vmin = -0.003e3
vmax = 0.01e3


CHANSN = ['072-080MHz', '080-088MHz', '088-095MHz', '095-103MHz', '103-111MHz', '111-118MHz', '118-126MHz',
 '126-134MHz', '139-147MHz', '147-154MHz', '154-162MHz', '162-170MHz', '170-177MHz', 
 '177-185MHz', '185-193MHz', '193-200MHz', '200-208MHz', '208-216MHz', '216-223MHz', '223-231MHz']

WIDECHANS = ["072-103MHz", '103-134MHz', "139-170MHz", "170-200MHz", "200-231MHz"]

for chan in WIDECHANS:
    im_rms_dr2 = fits.open(f"{data_path}/GLEAMX_DRII_{chan}_rms.fits")
    wcs_rms_dr3 = wcs.WCS(im_rms_dr2[0].header)
    centre = SkyCoord("00:30:00 -30:00:00", frame="fk5", unit=(u.hour,u.deg))
    size = u.Quantity((120,140),u.deg)


    dr2region = Cutout2D(im_rms_dr2[0].data*1000, centre, size, wcs = wcs_rms_dr3)

    q1 = np.nanpercentile(dr2region.data, 25)
    q2 = np.nanpercentile(dr2region.data, 50)
    q3 = np.nanpercentile(dr2region.data, 75)

    print(f"freq: {chan} = {q2} plus {q3-q2} minus {q2-q1}")
    