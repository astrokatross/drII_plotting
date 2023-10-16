#!/usr/bin/python

import numpy as np
import os
import re
import glob
import astropy.io.fits as fits
from astropy import wcs
from astropy import units as u
from astropy import coordinates
from astropy.coordinates import SkyCoord
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib import rc
from astropy.visualization.wcsaxes import Quadrangle
# rc('text', usetex=True)
# rc('font',**{'family':'serif','serif':['serif']})

"""
Script to plot nice coverage of the final mosaic and (hopefully) overlay the region to be considered during the data release 
"""

plt.rcParams["axes.axisbelow"] = False
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 10,
    "axes.titlesize": 10
})
cm = 1/2.54
sbplt_pad_left  = 0.125  # the left side of the subplots of the figure
sbplt_pad_right = 0.9    # the right side of the subplots of the figure
sbplt_pad_bottom = 0.1   # the bottom of the subplots of the figure
sbplt_pad_top = 0.9      # the top of the subplots of the figure
sbplt_pad_wspace = 0.1   # the amount of width reserved for blank space between subplots
sbplt_pad_hspace = 0.5   # the amount of height reserved for white space between subplots

subplot_cols = 1
subplot_rows = 1
save_path = "/data/gleam_x/drII_plotting/plots"
data_path = "/data/gleam_x/drII_plotting/data"
vmin = -0.003e3
vmax = 0.01e3
def redo_axis_labels(ticklist):
    new_label_list = []
    for label in ticklist:
        string = label.get_position()[0]
##        # Get numeric equivalent, change to hours:
        tick_num = np.degrees(float(string))
        if tick_num <= 0.:
            tick_num = tick_num + 360.
        tick_num = 360. - tick_num
        tick_num /= 15.
    # Note that this rounds to integer, might not be what you want:
        tick_string = "{:2.0f}".format(tick_num)
    # Add a raised "h" to denote hours:
        new_label = "$" + tick_string + "^h$"
#        new_label = "$" + str(tick_num) + "^h$"
        new_label_list.append(new_label)
    return new_label_list

def hr2rad(x):
    ''' convert hours to radians for the molleweide projection '''
    if x > 12.:
        x -= 24.
    return -np.radians(x*15.)

def unwrap(x):
    if x > np.radians(180.):
        x -= 2*np.pi
    return -x
vunwrap = np.vectorize(unwrap)

im_dr3 = fits.open(f"{data_path}/GLEAMX_DRII_170-231MHz.fits")
wcs_dr3 = wcs.WCS(im_dr3[0].header)


fig = plt.figure(figsize=(25*cm,20*cm))
ax = fig.add_subplot(1,1,1,projection=wcs_dr3)
im = ax.imshow(im_dr3[0].data*1000,origin="lower",vmin = vmin, vmax = vmax,cmap="gnuplot2")
overlay = ax.get_coords_overlay('icrs')
overlay.grid(color='black', ls='dotted')
r = Quadrangle((100, -85)*u.deg, -150*u.deg, 115*u.deg,
               edgecolor='white', facecolor='none',
               transform=ax.get_transform('fk5'))

ax.add_patch(r)
lon = ax.coords['dec']
lon.set_axislabel("Declination")
lat = ax.coords['ra']
lat.set_axislabel("Right Ascension")
lat.set_ticks(number=3)


fig.savefig(f"{save_path}/mosaic_coverage.png",bbox_inches="tight")
plt.clf()


plotnum = 1

fig = plt.figure(figsize=(20*cm,15*cm))

ax = fig.add_subplot(111, projection="mollweide")
ax.grid(True)

# Define already-observed region
tmp = [0., 0., 12, 12]
x_obs_1 = np.array([hr2rad(x) for x in tmp])
y_obs_1 = np.radians([-90, 30., 30, -90.])

tmp = [12.000001, 12.000001, 23.9999, 23.9999]
x_obs_2 = np.array([hr2rad(x) for x in tmp])
y_obs_2 = np.radians([-90., 30., 30., -90.])

# Define DR2
tmp = [20., 20., 6.5, 6.5]
x_obs_3 = np.array([hr2rad(x) for x in tmp])
y_obs_3 = np.radians([-90, 30., 30, -90.])

# Define DR1
tmp = [4., 4., 12., 12.]
x_es = np.array([hr2rad(x) for x in tmp])
y_es = np.radians([-32.7, -20.7, -20.7, -32.7])

# Define DR3
tmp = [6., 6., 12., 12.]
x_obs_4 = np.array([hr2rad(x) for x in tmp])
y_obs_4 = np.radians([-15, 30., 30, -15])

# Define DR3
tmp = [12.00001, 12.00001, 13., 13.]
x_obs_5 = np.array([hr2rad(x) for x in tmp])
y_obs_5 = np.radians([-32.7, -20.7, -20.7, -32.7])

#width = np.radians(30.)
#height = np.radians(30.)
#ax.add_patch(Ellipse((hr2rad(18.),np.radians(-30.)),
#                         width=width, height=height, angle = 0.0,
#                        ))
# First one is just to get it plotted so I can get the ticklabels
# fig.savefig("dummy.png")

# Observed region
ax.plot(x_obs_1, y_obs_1, color="C3", alpha=0.1, zorder = -10)
ax.fill(x_obs_1, y_obs_1, color="C3", alpha=0.1, zorder = -10)
ax.plot(x_obs_2, y_obs_2, color="C3", alpha=0.1, zorder = -10)
ax.fill(x_obs_2, y_obs_2, color="C3", alpha=0.1, zorder = -10, label="")
# DR1
ax.plot(x_es, y_es, color="C4", alpha=0.3)
ax.fill(x_es, y_es, color="C4", alpha=0.3, label="DR1", zorder=-10)

# DR2
ax.plot(x_obs_3, y_obs_3, color="C6", alpha=0.5, zorder = -5)
ax.fill(x_obs_3, y_obs_3, color="C6", alpha=0.5, zorder = -5, label="This release")

# DR3
# ax.plot(x_obs_4, y_obs_4, color="C1", alpha=0.3, zorder = -5)
# ax.fill(x_obs_4, y_obs_4, color="C1", alpha=0.3, zorder = -5, label="Processing")
ax.plot(x_obs_5, y_obs_5, color="C4", alpha=0.3, zorder = -5)
ax.fill(x_obs_5, y_obs_5, color="C4", alpha=0.3, zorder = -5)

# Plot the Galactic plane
l = np.arange(0, 360, 3)
b1 = -10.*np.ones(len(l))
b2 = 10.*np.ones(len(l))

gal_1 = SkyCoord(l, b1, unit=(u.deg, u.deg), frame="galactic")
gal_2 = SkyCoord(l, b2, unit=(u.deg, u.deg), frame="galactic")

ax.scatter(vunwrap(np.radians(gal_1.fk5.ra.value)), np.radians(gal_1.fk5.dec.value), color="black", zorder=100, marker=".", s=0.5)
ax.scatter(vunwrap(np.radians(gal_2.fk5.ra.value)), np.radians(gal_2.fk5.dec.value), color="black", zorder=100, marker=".", s=0.5)

gc = SkyCoord(0.0, 0.0, unit=(u.deg, u.deg), frame="galactic")
gsp = SkyCoord(0.0, -90, unit=(u.deg, u.deg), frame="galactic")

hyda = coordinates.get_icrs_coordinates("Hydra A")
vira = coordinates.get_icrs_coordinates("Virgo A")
crab = coordinates.get_icrs_coordinates("Crab")
cena = coordinates.get_icrs_coordinates("Centaurus A")
pica = coordinates.get_icrs_coordinates("Pictor A")
cyga = coordinates.get_icrs_coordinates("Cygnus A")

ax.scatter(unwrap(np.radians(gc.fk5.ra.value)), np.radians(gc.fk5.dec.value), marker="*",color="black")

# ax.scatter(unwrap(np.radians(gsp.fk5.ra.value)), np.radians(gsp.fk5.dec.value), marker="*")

# for a in hyda, vira, crab, cena, pica, cyga:
#    ax.scatter(unwrap(np.radians(a.fk5.ra.value)), np.radians(a.fk5.dec.value), marker="*", color="orange")

ticklist = ax.get_xmajorticklabels()
new_label_list = redo_axis_labels(ticklist)
new = ax.set_xticklabels(new_label_list)

ax.legend()


fig.savefig(f"{save_path}/gleamx_coverage_dr2.pdf", dpi=200, pad_inches=0.0, bbox_inches='tight')



