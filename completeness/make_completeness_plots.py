from argparse import ArgumentParser

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as maxes
from scipy.interpolate import interp1d
from astropy.io import fits
from astropy.wcs import WCS
#from astropy.visualization.wcsaxes import Quadrangle
from astropy.coordinates import SkyCoord
import astropy.units as u
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.mpl_axes import Axes
from matplotlib.ticker import FormatStrFormatter,MultipleLocator,FixedFormatter, FixedLocator,MaxNLocator


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 9})

cm = 1/2.54

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

def redo_axis_labels(ticklist):
    new_label_list = []
    for label in ticklist:
        # string = label.get_position()[0]
##        # Get numeric equivalent, change to hours:
        tick_num = label
        # tick_num = np.degrees(float(string))+50
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

LIMITS1a = [SkyCoord(i, unit=(u.hourangle, u.deg), frame='icrs') for i in ['06:00:00 +30:00:00','00:00:00 00:00:00']]
LIMITS1b = [SkyCoord(i, unit=(u.hourangle, u.deg), frame='icrs') for i in ['23:59:59 +30:00:00','21:00:00 00:00:00']]

LIMITS2a = [SkyCoord(i, unit=(u.hourangle, u.deg), frame='icrs') for i in ['06:00:00 00:00:00','00:00:00 -80:00:00']]
LIMITS2b = [SkyCoord(i, unit=(u.hourangle, u.deg), frame='icrs') for i in ['23:59:59 00:00:00','21:00:00 -80:00:00']]

LIMITS3a = [SkyCoord(i, unit=(u.hourangle, u.deg), frame='icrs') for i in ['06:00:00 -60:00:00','00:00:00 -90:00:00']]
LIMITS3b = [SkyCoord(i, unit=(u.hourangle, u.deg), frame='icrs') for i in ['23:59:59 -60:00:00','21:00:00 -90:00:00']]

LIMITS4a = [SkyCoord(i, unit=(u.hourangle, u.deg), frame='icrs') for i in ['06:00:00 -20:42:00','00:00:00 -32:42:00']]
LIMITS4b = [SkyCoord(i, unit=(u.hourangle, u.deg), frame='icrs') for i in ['23:59:59 -20:42:00','21:00:00 -32:42:00']]


LIMITSalla = [SkyCoord(i, unit=(u.hourangle, u.deg), frame='icrs') for i in ['06:00:00 +30:00:00','00:00:00 -90:00:00']]

LIMITSallb = [SkyCoord(i, unit=(u.hourangle, u.deg), frame='icrs') for i in ['23:59:59 +30:00:00','21:00:00 -90:00:00']]


LIMITS = [SkyCoord(i, unit=(u.hourangle, u.deg), frame='icrs') for i in ['13:00:00 -20:42:00','04:00:00 -32:42:00']]

def make_curve_plot(stats, flux_levels, base_out):
    curve = interp1d(
                stats[1],
                flux_levels,
    )

    fig, ax = plt.subplots(1,1, figsize=(8.9*cm,6.8*cm))

    ax.errorbar(
        flux_levels,
        stats[1],
        ms=2.0, 
        color="k",
        marker='o',
        linewidth=0.7,
        yerr=[stats[1]-stats[0], stats[2]-stats[1]]
    )

    ax.axhline(
        100,
        ls='--'
    )
    ax.grid(
        'major'
    )

    ax.set(
        xscale='log',
        xlabel='Flux density (mJy)',
        ylabel='Completeness (\%)',
    )
# Remove unnecessary scientific notation
    ax.xaxis.set_major_formatter(FormatStrFormatter('%3.0f'))

    fig.tight_layout()
    fig.savefig(f'{base_out}_curve.pdf', bbox_inches="tight")
    fig.savefig(f'{base_out}_curve.png', bbox_inches="tight")


def overlay_box(ax, text, x=0.02, y=0.125):
    ax.text(
        x, 
        y,
        text, 
        transform=ax.transAxes,
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='black', boxstyle='round,pad=0.25')
    )


def make_spatial_plot(comp_cube, flux_levels, w, base_out, cmap='inferno'):

# TODO these shouldn't be defined twice (here and at the top)
    limits = SkyCoord(
        np.array((359,0))*u.deg,
        np.array((-89, 30))*u.deg
        )
    # limits_ra = limits.ra.wrap_at(180*u.deg)
    # pix_limits = w.all_world2pix(limits_ra-180*u.deg, limits.dec, 0)
    pix_limits = w.all_world2pix(limits.ra, limits.dec, 0)
    x_lim = pix_limits[0][::-1]
    y_lim = pix_limits[1]
    print(limits)
    print(pix_limits)
    print(x_lim)
    print(y_lim)
    # x_lim = [90, 230]
    start_x, start_y = 0.15, 0.1
    delta_x, delta_y = 0.7, 0.2
    pad_y = 0.01
    offset_y = lambda x: x * (delta_y + pad_y) 

    #    ax4_loc = [start_x, start_y, delta_x, delta_y]

    ax3_loc = [start_x, start_y, delta_x, delta_y]
    ax2_loc = [start_x, start_y+offset_y(1), delta_x, delta_y]
    ax1_loc = [start_x, start_y+offset_y(2), delta_x, delta_y]

    #    for loc in (ax1_loc, ax2_loc, ax3_loc): #ax4_loc
    #        print(loc)

    fig = plt.figure(figsize=(17*cm, 7.3*cm))

    cax = fig.add_axes([0.865, 0.1, 0.0085, offset_y(2) + delta_y])
    ax1 = fig.add_axes(ax1_loc)
    # 
    ax1.imshow(
        np.roll(comp_cube[0].data[6],-180), 
        # comp_cube[0].data[6],
        vmin=0,
        vmax=100, cmap=cmap,
        aspect='auto',
        extent=[-180,180,90,-90]
    )
    ax1.set(
        xlim=[-90,50],
        ylim=[-90,30],
        ylabel='Dec'
    )
    overlay_box(ax1, f"{flux_levels[6]:.0f} mJy")
    overlay_box(ax1, "(a)", y=0.75)

    # lon = ax1.coords[0]
    # lon.set_format_unit(u.deg)
    # # lon.set_coord_type("longitude", coord_wrap=-180)
    # lon.set_ticklabel_visible(False)
    # lon.set_axislabel('')


    ax2 = fig.add_axes(ax2_loc)
    # ax2 = divider.append_axes('bottom', pad='3%', size='100%')
    ax2.imshow(
        np.roll(comp_cube[0].data[9],-180),
        # comp_cube[0].data[9],
        vmin=0,
        vmax=100, 
        cmap=cmap,
        aspect='auto',
        extent=[-180,180,90,-90]
    )
    ax2.set(
        xlim=[-90,50],
        ylim=[-90,30],
        ylabel='Dec'
    )
    # ax2.coords
    # lon = ax2.coords[0]

    # lon.set_format_unit(u.deg)
    # # lon.set_coord_type("longitude", coord_wrap=-180)
    # lon.set_ticklabel_visible(False)
    # lon.set_axislabel('')
    overlay_box(ax2, f"{flux_levels[9]:.0f} mJy")
    overlay_box(ax2, "(b)", y=0.75)


    ax3 = fig.add_axes(ax3_loc)

    # ax3 = fig.add_subplot(3,1,3, projection=w)
    # # ax3 = divider.append_axes('bottom', pad='3%', size='100%')
    cim = ax3.imshow(
        np.roll(comp_cube[0].data[12],-180),
        # comp_cube[0].data[12],
        vmin=0,
        vmax=100, 
        cmap=cmap,
        aspect='auto',
        extent=[-180,180,90,-90]
    )
    ax3.set(
        xlim=[-90,50],
        ylim=[-90,30],
        xlabel='RA',
        ylabel='Dec',
        # xticklabels=["6hr", "4hr", "2hr", "0hr", "22hr"],
    )

    xticks = ax3.get_xticks().tolist()

    print(xticks)
    new_label_list = redo_axis_labels(xticks)
    new = ax3.set_xticklabels(new_label_list)

    overlay_box(ax3, f"{flux_levels[12]:.0f} mJy")
    overlay_box(ax3, "(c)", y=0.75)


    cbar = fig.colorbar(cim, cax=cax, label='Completeness (\%)')

    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')


    # fig.tight_layout()
    fig.savefig(f"{base_out}_spatial.png", bbox_inches="tight")
    fig.savefig(f"{base_out}_spatial.pdf", bbox_inches="tight")




def make_completeness_plots(comp_cube, base_out='Completeness', s_min=-3, s_max=-0.5, s_step=0.1):
    comp_cube_fits = fits.open(comp_cube)

    w = WCS(comp_cube_fits[0].header).celestial

    x, y = np.indices(comp_cube_fits[0].data[0].shape)
    coords = w.wcs_pix2world(y, x, 0)

    mask1a = (LIMITS1a[1].ra.deg < coords[0]) & (coords[0] <= LIMITS1a[0].ra.deg) &\
            (LIMITS1a[1].dec.deg < coords[1]) & (coords[1] <= LIMITS1a[0].dec.deg)
    mask1b = (LIMITS1b[1].ra.deg < coords[0]) & (coords[0] <= LIMITS1b[0].ra.deg) &\
            (LIMITS1b[1].dec.deg < coords[1]) & (coords[1] <= LIMITS1b[0].dec.deg)
    mask2a = (LIMITS2a[1].ra.deg < coords[0]) & (coords[0] <= LIMITS2a[0].ra.deg) &\
            (LIMITS2a[1].dec.deg < coords[1]) & (coords[1] <= LIMITS2a[0].dec.deg)
    mask2b = (LIMITS2b[1].ra.deg < coords[0]) & (coords[0] <= LIMITS2b[0].ra.deg) &\
            (LIMITS2b[1].dec.deg < coords[1]) & (coords[1] <= LIMITS2b[0].dec.deg)
    mask3a = (LIMITS3a[1].ra.deg < coords[0]) & (coords[0] <= LIMITS3a[0].ra.deg) &\
            (LIMITS3a[1].dec.deg < coords[1]) & (coords[1] <= LIMITS3a[0].dec.deg)
    mask3b = (LIMITS3b[1].ra.deg < coords[0]) & (coords[0] <= LIMITS3b[0].ra.deg) &\
            (LIMITS3b[1].dec.deg < coords[1]) & (coords[1] <= LIMITS3b[0].dec.deg)

    mask4a = (LIMITS4a[1].ra.deg < coords[0]) & (coords[0] <= LIMITS4a[0].ra.deg) &\
            (LIMITS4a[1].dec.deg < coords[1]) & (coords[1] <= LIMITS4a[0].dec.deg)
    mask4b = (LIMITS4b[1].ra.deg < coords[0]) & (coords[0] <= LIMITS4b[0].ra.deg) &\
            (LIMITS4b[1].dec.deg < coords[1]) & (coords[1] <= LIMITS4b[0].dec.deg)

    maskalla = (LIMITSalla[1].ra.deg < coords[0]) & (coords[0] <= LIMITSalla[0].ra.deg) &\
            (LIMITSalla[1].dec.deg < coords[1]) & (coords[1] <= LIMITSalla[0].dec.deg)
    maskallb = (LIMITSallb[1].ra.deg < coords[0]) & (coords[0] <= LIMITSallb[0].ra.deg) &\
            (LIMITSallb[1].dec.deg < coords[1]) & (coords[1] <= LIMITSallb[0].dec.deg)

    mask = (LIMITS[1].ra.deg < coords[0]) & (coords[0] <= LIMITS[0].ra.deg) &\
            (LIMITS[1].dec.deg < coords[1]) & (coords[1] <= LIMITS[0].dec.deg)

    mask1 = np.logical_xor(mask1a,mask1b)
    mask2 = np.logical_xor(mask2a,mask2b)
    mask3 = np.logical_xor(mask3a,mask3b)
    mask4 = np.logical_xor(mask4a,mask4b)
    maskall = np.logical_xor(maskalla, maskallb)
    dec_regions = [mask, mask1, mask2, mask3, mask4, maskall]
    names = ["dr1", "highdec", "middec","lowdec", "dr1dec", "dr2"]
    i=0
    # convert to mJy in linear space
    s = np.arange(s_min, s_max + 0.0001, s_step)
    sdim = len(s)
    slin = 10 ** (s + 3.0)
    for mk in dec_regions:
        stats = np.nanpercentile(
                comp_cube_fits[0].data[:, mk], 
                [16, 50, 84], 
                axis=1
        )
        print(stats)
        print(np.sum(np.isfinite(comp_cube_fits[0].data[:, mk])))

        savename = f"{base_out}completeness_{names[i]}"
        make_curve_plot(stats, slin, savename)
        i+=1
    

    make_spatial_plot(comp_cube_fits, slin, w, f"{base_out}completeness")


if __name__ == '__main__':
    parser = ArgumentParser(description='Make some completeness figures up')
    parser.add_argument('completeness_cube', type=str, help='Path to the FITS cube produced by the completeness simulations')
    parser.add_argument('-b','--base-out', type=str, default='Completeness', help='Basename to use when making up output files')

    args = parser.parse_args()
    comp_cube = fits.open(args.completeness_cube)
    w = WCS(comp_cube[0].header).celestial
    make_completeness_plots(
        args.completeness_cube,
        base_out=args.base_out
    )
