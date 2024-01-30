#!/usr/bin/python

import numpy
#tables and votables
from astropy.io.votable import parse_single_table
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.visualization.wcsaxes import SphericalCircle
import matplotlib.pyplot as plt
from reproject import reproject_interp
import matplotlib as mpl
from matplotlib import rc
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 8})

cm = 1/2.54

# Not used at the moment
mwa_bmaj=0.0208 # degrees
mwa_bmin=0.0165 # degrees
mwa_bpa=-31.88 # degrees

datadir="/data/gleam_x/drII_plotting/data/"
savedir = "/data/gleam_x/drII_plotting/plots/"
mosaic="GLEAMX_DRII_170-231MHz_cut.fits"
splodge_mosaic="GLEAMX_DRII_200-231MHz.fits"
bkg="GLEAMX_DRII_170-231MHz_bkg.fits"
rms="GLEAMX_DRII_170-231MHz_rms.fits"
glm="GLEAM_rg_cut.fits"

vmin_mos = -0.003e3
vmax_mos = 0.01e3
vmin_glm = -0.02e3
vmax_glm = 0.06e3
vmin_bkg = -0.0002e3
vmax_bkg = 0.0002e3
vmin_rms = 0.0006e3
vmax_rms = 0.001e3

c = SkyCoord("01:40:00 -30:00:00", frame="fk5", unit=(u.hour, u.deg))
framesize = u.Quantity((3.5, 3.5), u.deg)

cmap="gnuplot2"
#cmap="Greys_r"

# table = parse_single_table("XG_mosaic.vot")
# my_data = table.array

# Set up the giant wrapper fits figure
fig = plt.figure(figsize=(19*cm,8.5*cm))

hdu_mos = fits.open(f"{datadir}{mosaic}")[0]
w_mos = wcs.WCS(hdu_mos.header)
hdu_glm = fits.open(f"{datadir}{glm}")[0]
w_glm = wcs.WCS(hdu_glm.header)

glm_im, glm_footprint = reproject_interp(hdu_glm, hdu_mos.header)


# cutout_mos = Cutout2D(1000*hdu_mos[0].data, c, framesize, w_mos)
# Left-most panel: GLEAM
# cutout_glm = Cutout2D(1000*hdu_glm[0].data, c, framesize, w_glm)




ax_glm = fig.add_axes([0.085, 0.085, 0.3, 0.85], projection = w_mos)
im_glm = ax_glm.imshow(glm_im*1000, origin="lower", cmap=cmap, vmin = vmin_glm, vmax = vmax_glm)
cax_glm = fig.add_axes([0.095, 0.86, 0.3-0.015, 0.015])
cb_glm = plt.colorbar(im_glm, cax = cax_glm, orientation="horizontal")
cb_glm.ax.xaxis.set_ticks_position('top')
cb_glm.ax.xaxis.set_label_position('top')

# Next panel: attractive slice of GLEAM-X

ax_mos = fig.add_axes([0.395,0.085,0.3,0.85], projection = w_mos)
im_mos = ax_mos.imshow(hdu_mos.data*1000, origin="lower", cmap=cmap, vmin = vmin_mos, vmax = vmax_mos)
cax_mos = fig.add_axes([0.4, 0.86, 0.3-0.015, 0.015])
cb_mos = plt.colorbar(im_mos, cax = cax_mos, orientation="horizontal")
cb_mos.ax.xaxis.set_ticks_position('top')
cb_mos.ax.xaxis.set_label_position('top')

## Top right panel: background
hdu_bkg = fits.open(f"{datadir}{bkg}")
w_bkg = wcs.WCS(hdu_bkg[0].header)
cutout_bkg = Cutout2D(1000*hdu_bkg[0].data, c, framesize, w_bkg)
ax_bkg = fig.add_axes([0.74,0.6,0.2,0.3], projection = cutout_bkg.wcs)
im_bkg = ax_bkg.imshow(cutout_bkg.data, origin="lower", cmap=cmap, vmin = vmin_bkg, vmax = vmax_bkg)
cax_bkg = fig.add_axes([0.92, 0.6, 0.008, 0.3])
cb_bkg = plt.colorbar(im_bkg, cax = cax_bkg)
#
## Bottom right panel: RMS
hdu_rms = fits.open(f"{datadir}{rms}")
w_rms = wcs.WCS(hdu_rms[0].header)
cutout_rms = Cutout2D(1000*hdu_rms[0].data, c, framesize, w_rms)
ax_rms = fig.add_axes([0.74,0.15,0.2,0.3], projection = cutout_rms.wcs)
im_rms = ax_rms.imshow(cutout_rms.data, origin="lower", cmap=cmap, vmin = vmin_rms, vmax = vmax_rms)
cax_rms = fig.add_axes([0.92, 0.15, 0.008, 0.3])
cb_rms = plt.colorbar(im_rms, cax = cax_rms)

#
lon = ax_glm.coords['dec']
lon.set_axislabel("Declination")
for ax in ax_glm, ax_mos:
    lat = ax.coords['ra']
    lat.set_axislabel("Right Ascension")
lon = ax_mos.coords['dec']
lon.set_axislabel(" ")
#lon.set_ticks_visible(False)
lon.set_ticklabel_visible(False)

for ax in ax_bkg, ax_rms:
    lon = ax.coords['dec']
    lon.set_axislabel("Dec")
#    lon.set_major_formatter('dd')
    lat = ax.coords['ra']
    lat.set_axislabel("RA")
# RA ticks overlap if left alone
    lat.set_ticks(number=3)


for cb in cb_glm, cb_mos, cb_bkg, cb_rms:
    cb.set_label("Brightness / mJy beam$^{-1}$")#, labelpad=-2)
#for cb in cb_bkg, cb_rms:
#    cb.set_label("Flux density / mJy beam$^{-1}$")

fig.savefig(f"{savedir}XG_mosaic.pdf", dpi=900, bbox_inches='tight')
fig.savefig(f"{savedir}XG_mosaic.png", dpi=900, bbox_inches='tight')

# This is just another plot to have a cutout of the splode phantom nuisance 

# c = SkyCoord("00:00:00 +12:13:00", frame="fk5", unit=(u.hour, u.deg))
# framesize = (4*u.deg,10*u.deg)

# fig = plt.figure(figsize=(19*cm,8.5*cm))


# hdu_mos = fits.open(f"{datadir}{splodge_mosaic}")
# w_mos = wcs.WCS(hdu_mos[0].header)
# cutout_mos = Cutout2D(1000*hdu_mos[0].data, c, framesize, w_mos)
# ax_mos = fig.add_subplot(111, projection = cutout_mos.wcs)
# im_mos = ax_mos.imshow(cutout_mos.data, origin="lower", cmap=cmap,  vmin = vmin_mos, vmax = vmax_mos)
# # cax_mos = fig.add_axes([0.4, 0.86, 0.3-0.015, 0.015])
# # cax_rms = fig.add_axes([0.92, 0.15, 0.008, 0.3])
# cb_mos = plt.colorbar(im_mos, ax = ax_mos,fraction=0.018,pad=0.05)

# cb_mos.ax.xaxis.set_ticks_position('top')
# cb_mos.ax.xaxis.set_label_position('top')
# cb_mos.set_label("Brightness / mJy beam$^{-1}$")

# phantom1 = SkyCoord("00:09:03.9 +12:07:25.3", unit=(u.hourangle, u.deg))
# phantom2 = SkyCoord("00:03:59.4 +12:09:31.3", unit=(u.hourangle, u.deg))
# phantom3 = SkyCoord("23:58:57.3 +12:09:23.7", unit=(u.hourangle, u.deg))
# phantom4 = SkyCoord("23:53:59.5 +12:07:52.6", unit=(u.hourangle, u.deg))
# phantom5 = SkyCoord("23:48:51.8 +12:08:32.3", unit=(u.hourangle, u.deg))
# phantom6 = SkyCoord("23:43:46.2 +12:05:38.0", unit=(u.hourangle, u.deg))
# phantom7 = SkyCoord("00:14:07.1 +12:05:47.8", unit=(u.hourangle, u.deg))
# phantom8 = SkyCoord("00:19:13.3 +12:01:59.0", unit=(u.hourangle, u.deg))

# for phan in [phantom1, phantom2, phantom3, phantom4, phantom5, phantom6,phantom7,phantom8]:
#     r = SphericalCircle(phan, 0.2*u.deg, edgecolor="white", lw=1, facecolor="none", transform=ax_mos.get_transform('fk5'))
#     ax_mos.add_patch(r)

# lon = ax_mos.coords['dec']
# lon.set_axislabel("Declination")
# # for ax in ax_mos:
# lat = ax_mos.coords['ra']
# lat.set_axislabel("Right Ascension")
# # lon = ax_mos.coords['dec']
# # lon.set_axislabel(" ")
# #lon.set_ticks_visible(False)
# # lon.set_ticklabel_visible(False)

# fig.savefig(f"{savedir}phantom_cutout.pdf", dpi=900, bbox_inches='tight')
# fig.savefig(f"{savedir}phantom_cutout.png", dpi=900, bbox_inches='tight')


# c = SkyCoord("23:33:55.2 -23:43:40", frame="fk5", unit=(u.hour, u.deg))
# framesize = 2*u.deg

# fig = plt.figure(figsize=(8.5*cm,8.5*cm))


# hdu_mos = fits.open(f"{datadir}{mosaic}")
# w_mos = wcs.WCS(hdu_mos[0].header)
# cutout_mos = Cutout2D(1000*hdu_mos[0].data, c, framesize, w_mos)
# ax_mos = fig.add_subplot(111, projection = cutout_mos.wcs)
# im_mos = ax_mos.imshow(cutout_mos.data, origin="lower", cmap=cmap,  vmin = vmin_mos, vmax = vmax_mos)
# # cax_mos = fig.add_axes([0.4, 0.86, 0.3-0.015, 0.015])
# # cax_rms = fig.add_axes([0.92, 0.15, 0.008, 0.3])
# cb_mos = plt.colorbar(im_mos, ax = ax_mos,fraction=0.018,pad=0.05)

# cb_mos.ax.xaxis.set_ticks_position('top')
# cb_mos.ax.xaxis.set_label_position('top')
# cb_mos.set_label("Brightness / mJy beam$^{-1}$")

# lon = ax_mos.coords['dec']
# lon.set_axislabel("Declination")
# # for ax in ax_mos:
# lat = ax_mos.coords['ra']
# lat.set_axislabel("Right Ascension")
# # lon = ax_mos.coords['dec']
# # lon.set_axislabel(" ")
# #lon.set_ticks_visible(False)
# # lon.set_ticklabel_visible(False)

# fig.savefig(f"{savedir}GLEAMX_J233355-234340_cutout.pdf", dpi=900, bbox_inches='tight')
# fig.savefig(f"{savedir}GLEAMX_J233355-234340_cutout.png", dpi=900, bbox_inches='tight')