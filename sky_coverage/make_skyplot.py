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
# rc('text', usetex=True)
# rc('font',**{'family':'serif','serif':['serif']})

"""
Script to plot nice coverage of the final mosaic and (hopefully) overlay the region to be considered during the data release 
"""

plt.rcParams["axes.axisbelow"] = False
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 9,
    "axes.titlesize": 9
})
cm = 1/2.54
sbplt_pad_left  = 0.125  # the left side of the subplots of the figure
sbplt_pad_right = 0.9    # the right side of the subplots of the figure
sbplt_pad_bottom = 0.1   # the bottom of the subplots of the figure
sbplt_pad_top = 0.9      # the top of the subplots of the figure
sbplt_pad_wspace = 0.1   # the amount of width reserved for blank space between subplots
sbplt_pad_hspace = 0.5   # the amount of height reserved for white space between subplots

aspect = 20
pad_fraction = 0.5
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

def deg_unwrap(x):
    if x > 180.:
        x -= 360
    return -x

cmap=plt.cm.gnuplot2

vunwrap = np.vectorize(unwrap)
deg_vunwrap = np.vectorize(deg_unwrap)


# im_dr3 = fits.open(f"{data_path}/GLEAMX_DRII_170-231MHz.fits")
# wcs_dr2 = wcs.WCS(im_dr3[0].header)

im_rms_dr2 = fits.open(f"{data_path}/GLEAMX_DRII_170-231MHz_rms.fits")
wcs_rms_dr3 = wcs.WCS(im_rms_dr2[0].header)
centre = SkyCoord("01:30:00 -30:00:00", frame="fk5", unit=(u.hour,u.deg))
size = u.Quantity((115,140),u.deg)
# rms_cutout = Cutout2D(im_rms_dr2[0].data, centre, size, wcs=wcs_rms_dr3)

print(f"mean and std rms: {np.nanmean(im_rms_dr2[0].data*1000)} +/- {np.nanstd(im_rms_dr2[0].data*1000)}")

print(f"min: {np.nanmin(im_rms_dr2[0].data*1000)}")
print(f"max: {np.nanmax(im_rms_dr2[0].data*1000)} ")
print(f"min (minus in plusminus): {np.nanmean(im_rms_dr2[0].data*1000)- np.nanmin(im_rms_dr2[0].data*1000)}")
print(f"min (plus in plusminus): {np.nanmax(im_rms_dr2[0].data*1000)-np.nanmean(im_rms_dr2[0].data*1000)}")



fig = plt.figure(figsize=(8.9*cm,8*cm))
ax = fig.add_subplot(1,1,1,projection=wcs_rms_dr3)
# im = ax.imshow(im_rms_dr2[0].data*1000,origin="lower",cmap="gnuplot2", norm=LogNorm())
im = ax.imshow(im_rms_dr2[0].data*1000, origin="lower", cmap="gnuplot2", norm=LogNorm(vmin=0.5,vmax=13))
# overlay = ax.get_coords_overlay('icrs')
# overlay.grid(color='black', ls='dotted')
# r = Quadrangle((100, 30)*u.deg, -140*u.deg, -115*u.deg,
#                edgecolor='white', facecolor='none',
#                transform=ax.get_transform('fk5'))

r = Quadrangle((96, 0)*u.deg, -140*u.deg, -40*u.deg,
               edgecolor='white', facecolor='none',
               transform=ax.get_transform('fk5'))
ax.add_patch(r)

formatter = LogFormatter(10, labelOnlyBase=False) 
cb = fig.colorbar(im, ax=ax,format=formatter, fraction=0.0345, pad=0.05)
cb.set_label("RMS (mJy/beam)")
# overlay = ax.get_coords_overlay('fk5')
# overlay.grid(color='black', ls='dotted')

# # Hide the right and top spines
# ax.coords["ra"].set_visible(False)
# overlay.set_visible(False)


# # Only show ticks on the left and bottom spines

# ax[0].set_axislabel(' ')
# ax[1].set_axislabel(' ')

ax.update({'xlabel': "Right Ascension (deg)", 'ylabel': "Declination (deg)"})
ax.grid(color="black", ls="dotted")
# ax.update({"labelbottom": False,  "labeltop": True})

# overlay = ax.get_coords_overlay('icrs')
# overlay[0].set_axislabel(' ')
# overlay[0].set_axislabel_position("lt")
# overlay[1].set_axislabel(' ')
# overlay[1].set_axislabel_position("lt")
# overlay.grid(color='black', ls='dotted')
# ax[0].set_axislabel(' ')
# ax[0].set_axislabel_position("lt")
# ax[1].set_axislabel(' ')
# ax[1].set_axislabel_position("lt")

lon = ax.coords['dec']
lon.set_ticks(number=5)
# lon.set_axislabel("Declination")
# lon.set_ticklabel('white')
lat = ax.coords['ra']
# lat.set_axislabel("Right Ascension")
lat.set_axislabel_position("t")
lat.set_ticklabel_position("t")
lat.set_format_unit(u.deg)
lat.set_ticks(number=7)
# lat.set_ticklabel("white")

fig.savefig(f"{save_path}/mosaic_coverage.pdf",bbox_inches="tight")
plt.clf()


# just a quick addition for some noise analysis 


im_rms_dr2 = fits.open(f"{data_path}/GLEAMX_DRII_170-231MHz_rms.fits")
wcs_rms_dr3 = wcs.WCS(im_rms_dr2[0].header)
centre = SkyCoord("00:30:00 -30:00:00", frame="fk5", unit=(u.hour,u.deg))
size = u.Quantity((120,140),u.deg)


dr2region = Cutout2D(im_rms_dr2[0].data*1000, centre, size, wcs = wcs_rms_dr3)

q1 = np.nanpercentile(dr2region.data, 25)
q2 = np.nanpercentile(dr2region.data, 50)
q3 = np.nanpercentile(dr2region.data, 75)

print(f"q2: {q2} plus {q3-q2} minus {q2-q1}")
print(f"Other vals: {q1,q2,q3}")
print(f"{np.nanmean(dr2region.data)}")
print(f"{np.nanmean(im_rms_dr2[0].data)}")
plotnum = 1

# fig = plt.figure(figsize=(20*cm,15*cm))

# ax = fig.add_subplot(111, projection="mollweide")
# ax.grid(True)

# # Define already-observed region
# tmp = [0., 0., 12, 12]
# x_obs_1 = np.array([hr2rad(x) for x in tmp])
# y_obs_1 = np.radians([-90, 30., 30, -90.])

# tmp = [12.000001, 12.000001, 23.9999, 23.9999]
# x_obs_2 = np.array([hr2rad(x) for x in tmp])
# y_obs_2 = np.radians([-90., 30., 30., -90.])

# # Define DR2
# tmp = [20.4, 20.4, 6.4, 6.4]
# x_obs_3 = np.array([hr2rad(x) for x in tmp])
# y_obs_3 = np.radians([-90, 30., 30, -90.])

# # Define Brandon region 
# tmp = [21.07, 21.07, 6.4, 6.4]
# x_obs_b = np.array([hr2rad(x) for x in tmp])
# y_obs_b = np.radians([-40, 0., 0, -40.])

# # Define DR1
# tmp = [4., 4., 12., 12.]
# x_es = np.array([hr2rad(x) for x in tmp])
# y_es = np.radians([-32.7, -20.7, -20.7, -32.7])

# # Define DR3
# tmp = [6., 6., 12., 12.]
# x_obs_4 = np.array([hr2rad(x) for x in tmp])
# y_obs_4 = np.radians([-15, 30., 30, -15])

# # Define DR3
# tmp = [12.00001, 12.00001, 13., 13.]
# x_obs_5 = np.array([hr2rad(x) for x in tmp])
# y_obs_5 = np.radians([-32.7, -20.7, -20.7, -32.7])

# # Define cursed region GLEAM
# tmp = [22., 22., 23.99999, 23.9999]
# x_obs_curse = np.array([hr2rad(x) for x in tmp])
# y_obs_curse = np.radians([30., 0., 0, 30])

# width = np.radians(30.)
# height = np.radians(30.)
# # ax.add_patch(Ellipse((hr2rad(18.),np.radians(-30.)),
#                         # width=width, height=height, angle = 0.0, 
#                     #    ))
# # First one is just to get it plotted so I can get the ticklabels
# fig.savefig(f"{save_path}/dummy.png")

# # Observed region
# ax.plot(x_obs_1, y_obs_1, color="C3", alpha=0.1, zorder = -10)
# ax.fill(x_obs_1, y_obs_1, color="C3", alpha=0.1, zorder = -10)
# ax.plot(x_obs_2, y_obs_2, color="C3", alpha=0.1, zorder = -10)
# ax.fill(x_obs_2, y_obs_2, color="C3", alpha=0.1, zorder = -10, label="")
# # # DR1
# ax.plot(x_es, y_es, color="C4", alpha=0.3)
# ax.fill(x_es, y_es, color="C4", alpha=0.3, label="DR1", zorder=-10)

# # DR2
# ax.plot(x_obs_3, y_obs_3, color="C6", alpha=0.4, zorder = -5)
# ax.fill(x_obs_3, y_obs_3, color="C6", alpha=0.4, zorder = -5, label="DR2")

# # Cursed
# # ax.plot(x_obs_curse, y_obs_curse, color="C9", alpha=0.4, zorder = -6)
# # ax.fill(x_obs_curse, y_obs_curse, color="C9", alpha=0.4, zorder = -6, label="Cursed region")

# # # Brandon region
# # ax.plot(x_obs_b, y_obs_b, color="C6", alpha=0.5, zorder = -5)
# # ax.fill(x_obs_b, y_obs_b, color="C6", alpha=0.5, zorder = -5, label="GXDS")

# # DR3
# # ax.plot(x_obs_4, y_obs_4, color="C1", alpha=0.3, zorder = -5)
# # ax.fill(x_obs_4, y_obs_4, color="C1", alpha=0.3, zorder = -5, label="Processing")
# ax.plot(x_obs_5, y_obs_5, color="C4", alpha=0.3, zorder = -5)
# ax.fill(x_obs_5, y_obs_5, color="C4", alpha=0.3, zorder = -5)

# # Plot the Galactic plane
# l = np.arange(0, 360, 3)
# b1 = -10.*np.ones(len(l))
# b2 = 10.*np.ones(len(l))

# gal_1 = SkyCoord(l, b1, unit=(u.deg, u.deg), frame="galactic")
# gal_2 = SkyCoord(l, b2, unit=(u.deg, u.deg), frame="galactic")

# ax.scatter(vunwrap(np.radians(gal_1.fk5.ra.value)), np.radians(gal_1.fk5.dec.value), color="black", zorder=100, marker=".", s=0.5)
# ax.scatter(vunwrap(np.radians(gal_2.fk5.ra.value)), np.radians(gal_2.fk5.dec.value), color="black", zorder=100, marker=".", s=0.5)

# gc = SkyCoord(0.0, 0.0, unit=(u.deg, u.deg), frame="galactic")
# gsp = SkyCoord(0.0, -90, unit=(u.deg, u.deg), frame="galactic")

# hyda = coordinates.get_icrs_coordinates("Hydra A")
# vira = coordinates.get_icrs_coordinates("Virgo A")
# crab = coordinates.get_icrs_coordinates("Crab")
# cena = coordinates.get_icrs_coordinates("Centaurus A")
# pica = coordinates.get_icrs_coordinates("Pictor A")
# cyga = coordinates.get_icrs_coordinates("Cygnus A")

# ax.scatter(unwrap(np.radians(gc.fk5.ra.value)), np.radians(gc.fk5.dec.value), marker="*",color="black")

# # phantom1 = SkyCoord("00:09:03.9 +12:07:25.3", unit=(u.hourangle, u.deg))
# # phantom2 = SkyCoord("00:03:59.4 +12:09:31.3", unit=(u.hourangle, u.deg))
# # phantom3 = SkyCoord("23:58:57.3 +12:09:23.7", unit=(u.hourangle, u.deg))
# # phantom4 = SkyCoord("23:53:59.5 +12:07:52.6", unit=(u.hourangle, u.deg))
# # phantom5 = SkyCoord("23:48:51.8 +12:08:32.3", unit=(u.hourangle, u.deg))
# # phantom6 = SkyCoord("23:43:46.2 +12:05:38.0", unit=(u.hourangle, u.deg))
# # phantom7 = SkyCoord("00:14:07.1 +12:05:47.8", unit=(u.hourangle, u.deg))
# # phantom8 = SkyCoord("00:19:13.3 +12:01:59.0", unit=(u.hourangle, u.deg))
# # for phan in [phantom1, phantom2, phantom3, phantom4, phantom5, phantom6,phantom7]:
# #     ax.scatter(unwrap(np.radians(phan.fk5.ra.value)), np.radians(phan.fk5.dec.value), marker=".", color=cmap(0.5))
# # ax.scatter(unwrap(np.radians(phan.fk5.ra.value)), np.radians(phan.fk5.dec.value), marker=".", color=cmap(0.5), label="Phantom")
# # ax.scatter(unwrap(np.radians(gsp.fk5.ra.value)), np.radians(gsp.fk5.dec.value), marker="*")

# # for a in hyda, vira, crab, cena, pica, cyga:
# #    ax.scatter(unwrap(np.radians(a.fk5.ra.value)), np.radians(a.fk5.dec.value), marker="*", color="orange")

# ticklist = ax.get_xmajorticklabels()
# new_label_list = redo_axis_labels(ticklist)
# new = ax.set_xticklabels(new_label_list)

# ax.legend()


# fig.savefig(f"{save_path}/gleamx_coverage.pdf", dpi=300, pad_inches=0.0, bbox_inches='tight')
# fig.savefig(f"{save_path}/gleamx_coverage.png", dpi=300, pad_inches=0.0, bbox_inches='tight')
# plt.clf()



# negative_srcs = Table.read(f"{data_path}/GLEAMX_DRII_170-231MHz_psf_comp_negative_comp.fits")

# ras = negative_srcs["ra"]
# decs = negative_srcs["dec"]

# ramin = 90
# ramax = 320
# decmax = 30
# decmin = -85

# idx = np.where(decs < decmax)
# # idx = np.intersect1d(idx, np.where(cat_neg["dec"] > decmin))
# idx1 = np.intersect1d(idx, np.where(ras <= ramin))
# idx2 = np.intersect1d(idx, np.where(ras >= ramax ))
# idx = np.union1d(idx1, idx2)

# ras = ras[idx]
# decs = decs[idx]

# # hpmap = cat2hpx(ras, decs, 12, radec=True)


# pk_flux = np.array(negative_srcs["peak_flux"])
# pk_flux = pk_flux[idx]
# rms = np.array(negative_srcs["local_rms"])
# rms = rms[idx]

# coords = SkyCoord(ra=ras, dec=decs, frame='icrs', unit='deg')

# ra=np.radians(vunwrap(coords.ra.value))
# dec = np.radians(decs)


# weights = -pk_flux/rms 

# weighted_counts = np.histogram2d(ra,dec, weights=np.array(-negative_srcs["peak_flux"])/np.array(negative_srcs["local_rms"]), bins=[30,10])[0]


# theta = np.radians(90. - b)
# phi = np.radians(l)

# nsources = len(negative_srcs)
# nside=12
# npix=hp.nside2npix(nside)
# hpx_map = np.zeros(npix, dtype=int)
# fig = plt.figure(figsize=(25*cm,20*cm))
# fignum = plt.gcf()


# # ipix = hp.ang2pix(nside, theta, phi)
# ipix = hp.ang2pix(nside, (90-decs)/180*np.pi, ras/180*np.pi, nest=False)
# idx, counts = np.unique(ipix, return_counts=True)
# hpx_map = np.zeros(npix, dtype=int)
# hpx_map[idx] = counts

# density_map = hp.mollview(np.log10(hpx_map+1), cmap="gnuplot2_r", fig=fignum, title="Source Density Map",norm="hist", return_projected_map=True)
# hp.graticule()

# plt.clf()


# fig = plt.figure(figsize=(25*cm,20*cm))
# ax = fig.add_subplot(1,1,1,projection="mollweide")

# im = ax.imshow(hpx_map,cmap="gnuplot2")
# # overlay = ax.get_coords_overlay('icrs')
# overlay.grid(color='black', ls='dotted')
# r = Quadrangle((96, -40)*u.deg, -140*u.deg, 40*u.deg,
#                edgecolor='white', facecolor='none',
#                transform=ax.get_transform('fk5'))
# ax.add_patch(r)

# formatter = LogFormatter(10, labelOnlyBase=False) 
# cb = fig.colorbar(im, ax=ax,format=formatter, fraction=0.0345, pad=0.05)
# cb.set_label("RMS (mJy/beam)")
# overlay = ax.get_coords_overlay('fk5')
# overlay.grid(color='white', ls='dotted')

# # Hide the right and top spines
# ax["ra"].set_visible(False)
# overlay["dec"].set_visible(False)


# # Only show ticks on the left and bottom spines
# overlay[0].set_axislabel(' ')
# overlay[1].set_axislabel(' ')


# ax.update({'xlabel': "Right Ascension", 'ylabel': "Declination"})
# lon = ax.coords['dec']
# lon.set_axislabel("Declination")
# lat = ax.coords['ra']
# lat.set_axislabel("Right Ascension")
# lat.set_ticks(number=3)
# lat.set_ticklabel("white")

# # fig.savefig(f"{save_path}/mosaic_coverage.pdf",bbox_inches="tight")
# # plt.clf()


# fig.savefig(f"{save_path}/testing.png")
# plt.clf()


# fig = plt.figure(figsize=(20*cm,15*cm))

# ax = fig.add_subplot(111,projection='mollweide')
# ax.grid(True)
# width = np.radians(30.)
# height = np.radians(30.)



# overlay = ax.get_coords_overlay('icrs')
# overlay.grid(color='black', ls='dotted')
# r = Quadrangle((96, -40)*u.deg, -140*u.deg, 40*u.deg,
#                edgecolor='white', facecolor='none',
#                transform=ax.get_transform('fk5'))


# First one is just to get it plotted so I can get the ticklabels
# fig.savefig(f"{save_path}/dummy.png")


# # ax.plot(ra,dec,'o',markersize=5,color="C6")
# # hex_map = ax.hist2d(ra,dec,bins=45,cmap='gnuplot_r')

# hex_map = ax.hexbin(ra,dec,cmap="gnuplot2_r",C=weights, gridsize=50,bins='log',vmax=5)
# # ax.grid(True)
# # ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
# # r = Quadrangle((np.radians(90), np.radians(30))*u.rad, -(np.radians(140))*u.rad, -np.radians(115)*u.rad,edgecolor="black", facecolor='none')

# # ax.add_patch(Quadrangle((np.radians(90), np.radians(30))*u.rad, -(np.radians(140))*u.rad, -np.radians(115)*u.rad,edgecolor="black", facecolor='none'))
# # ax.plot(x_obs_3, y_obs_3, color="C3", alpha=0.1, zorder = -5)
# # ax.fill(x_obs_3, y_obs_3, color="C3", alpha=0.1, zorder = -5, label="GLEAM-X DR2")
# # cb = plt.colorbar(hex_map, spacing='uniform', ax=ax,fraction=0.025, pad=0.05,format=formatter)
# # cb = fig.colorbar(hex_map)
# # cb.set_label("RMS (mJy/beam)")
# ticklist = ax.get_xmajorticklabels()
# new_label_list = redo_axis_labels(ticklist)
# new = ax.set_xticklabels(new_label_list)
# plt.xlabel('RA')
# plt.ylabel('Declination')

# fig.savefig(f"{save_path}/testing_2.png")