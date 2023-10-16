#!/usr/bin/env python

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

import sys

# TODO: make into argparse
cat = sys.argv[1]
out = sys.argv[2]

# Nights in IDR1   RA start    RA end
# 2018-02-04         50           220
# 2018-02-09         40           215
# 2018-02-20         60           240
# 2018-03-05         85           270
# So at maximum 50 to 240
# Also need to subtract ~10 deg from each end because of PB cutoff
# 60 to 215
# Visually, image quality is much lower past 190, probably because of Cen A
# So let's go with 60 to 195 = 4h to 13h


# Nights in Dec-13 teststrip   RA start    RA end
# 2020-10-02                    338          97
# 2018-02-09                    302          
# 2018-02-20                    60           240
# 2018-03-05                    85           270
# So at maximum 50 to 240
# Also need to subtract ~10 deg from each end because of PB cutoff
# 60 to 215
# Visually, image quality is much lower past 190, probably because of Cen A
# So let's go with 60 to 195 = 4h to 13h

ramin = 95
ramax = 310
decmax = 30

hdu = fits.open(cat)
idx = np.where(hdu[1].data["ref_dec"] < decmax)
idx1 = np.intersect1d(idx, np.where(hdu[1].data["ref_ra"] <= ramin))
# idx = np.intersect1d(idx, np.where(hdu[1].data["dec"] > decmin))
idx2 = np.intersect1d(idx, np.where(hdu[1].data["ref_ra"] >= ramax ))

idx = np.union1d(idx1, idx2)
hdu[1].data = hdu[1].data[idx]

hdu.writeto(out,overwrite=True)
