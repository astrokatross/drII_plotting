from argparse import ArgumentParser
import numpy as np
from astropy.io import fits 
from astropy.nddata import Cutout2D
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
import copy
from matplotlib import rc
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 8,
    "axes.titlesize": 8
})
cm = 1/2.54
fit_g = fitting.LevMarLSQFitter()


# things we need for plotting: 
# bkgsub_sig, masked_bkgsub_sigma, resid_bkgsub_sigma, img, rms, bkg

# note: stdev=1, nsig=5 --> sigma
# stdev=0.025/0.002, nsig=1 --> raw

def fit_hists(data, stdev = 1.0, nsig=5):
    bins=100
    xmin = -nsig*stdev
    xmax = nsig*stdev

    n, binss = np.histogram(data,bins = np.linspace(xmin,xmax,bins))
    n = np.array(n,dtype="float32")/len(n)
    bin_centers = (binss+((binss[1]-binss[0])/2))[:-1]
    g_init = models.Gaussian1D(amplitude=1., mean=0., stddev=stdev)
    gu = fit_g(g_init, bin_centers, n)

    return n, bin_centers, binss, gu 

def plot_distribution(imstring, stdev=1.0, nsig=5., mean_rms=1.0, mean_bkg=0., savenm=""):

# Image names I'll need: 
# _bkgsub_sigma.fits, _bkgsub_masked_sigma.fits, _bkgsub_resid_sigma.fits
# _priorsub_bkgsub_sigma.fits, _priorsub_bkgsub_masked_sigma.fits, _priorsub_bkgsub_resid_sigma.fits 

    mean_bkg = 0.
    mean_rms = 1.
    xmin = -nsig*stdev
    xmax = nsig*stdev

    im = fits.getdata(f"{args.datadir}{imstring}_bkgsub_sigma.fits")
    im_masked = fits.getdata(f"{args.datadir}{imstring}_bkgsub_masked_sigma.fits")
    im_resid = fits.getdata(f"{args.datadir}{imstring}_bkgsub_resid_sigma.fits")

    n_im, bin_centers_im, bins_im, gu_im = fit_hists(im)
    n_im_mask, bin_centers_im_mask, bins_im_mask, gu_im_mask = fit_hists(im_masked)
    n_im_resid, bin_centers_im_resid, bins_im_resid, gu_im_resid = fit_hists(im_resid)

    tgu_im = copy.deepcopy(gu_im)
    tgu_im.mean.value = mean_bkg
    tgu_im.stddev.value = mean_rms
    tgu_im_mask = copy.deepcopy(gu_im_mask)
    tgu_im_mask.mean.value = mean_bkg
    tgu_im_mask.stddev.value = mean_rms
    tgu_im_resid = copy.deepcopy(gu_im_resid)
    tgu_im_resid.mean.value = mean_bkg
    tgu_im_resid.stddev.value = mean_rms

    if args.compare is True: 
        im_gd = fits.getdata(f"{args.datadir}{imstring}_priorsub_bkgsub_sigma.fits")
        im_masked_gd = fits.getdata(f"{args.datadir}{imstring}_priorsub_bkgsub_masked_sigma.fits")
        im_resid_gd = fits.getdata(f"{args.datadir}{imstring}_priorsub_bkgsub_resid_sigma.fits")
        n_im_gd, bin_centers_im_gd, bins_im_gd, gu_im_gd = fit_hists(im_gd)
        n_im_mask_gd, bin_centers_im_mask_gd, bins_im_mask_gd, gu_im_mask_gd = fit_hists(im_masked_gd)
        n_im_resid_gd, bin_centers_im_resid_gd, bins_im_resid_gd, gu_im_resid_gd = fit_hists(im_resid_gd)

        tgu_im_gd = copy.deepcopy(gu_im_gd)
        tgu_im_gd.mean.value = mean_bkg
        tgu_im_gd.stddev.value = mean_rms
        tgu_im_mask_gd = copy.deepcopy(gu_im_mask_gd)
        tgu_im_mask_gd.mean.value = mean_bkg
        tgu_im_mask_gd.stddev.value = mean_rms
        tgu_im_resid_gd = copy.deepcopy(gu_im_resid_gd)
        tgu_im_resid_gd.mean.value = mean_bkg
        tgu_im_resid_gd.stddev.value = mean_rms


    print("calculated distributions, plotting now...")

    if args.compare is True: 
        fig = plt.figure(figsize=(19*cm, 12.4*cm))
        gs = fig.add_gridspec(2,3, wspace=0.2,hspace=0.4)
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[0,1])
        ax3 = fig.add_subplot(gs[0,2])
        ax4 = fig.add_subplot(gs[1,0])
        ax5 = fig.add_subplot(gs[1,1])
        ax6 = fig.add_subplot(gs[1,2])
        axes_total = [ax1, ax2, ax3, ax4, ax5, ax6]
    else: 
        fig = plt.figure(figsize=(19*cm, 6.2*cm))
        gs = fig.add_gridspec(1,3, wspace=0.2)
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[0,1])
        ax3 = fig.add_subplot(gs[0,2])
        axes_total = [ax1, ax2, ax3]

    ax1.bar(bin_centers_im, n_im, color = 'C6', edgecolor = "none", width=(bins_im[1]-bins_im[0]),alpha=0.5) 
    ax1.plot(np.linspace(xmin,xmax,1000),gu_im(np.linspace(xmin,xmax,1000)), color='k',lw=0.5, alpha=0.8, label='Fitted')
    ax1.plot(np.linspace(xmin,xmax,1000),tgu_im(np.linspace(xmin,xmax,1000)), color='k',lw=1, linestyle="--", label='BANE')
    ax1.set_ylim([1.0,1.2*max(n_im)])
    ax1.set_title(f"Bkg subtracted, S/N image", fontsize=7)

    ax2.bar(bin_centers_im_mask, n_im_mask, color = 'C6', edgecolor = "none", width=(bins_im_mask[1]-bins_im_mask[0]),alpha=0.5) 
    ax2.plot(np.linspace(xmin,xmax,1000),gu_im_mask(np.linspace(xmin,xmax,1000)), color='k',lw=0.5, alpha=0.8, label='Fitted')
    ax2.set_ylim([1.0,1.2*max(n_im_mask)])
    ax2.set_title(f"Bkg subtracted, $>5\sigma$ sources masked", fontsize=7)
    ax2.plot(np.linspace(xmin,xmax,1000),tgu_im_mask(np.linspace(xmin,xmax,1000)), color='k',lw=1, linestyle="--", label='BANE')

    ax3.bar(bin_centers_im_resid, n_im_resid, color = 'C6', edgecolor = "none", width=(bins_im_resid[1]-bins_im_resid[0]),alpha=0.5) 
    ax3.plot(np.linspace(xmin,xmax,1000),gu_im_resid(np.linspace(xmin,xmax,1000)), color='k',lw=0.5, alpha=0.8, label='Fitted')
    ax3.set_ylim([1.0,1.2*max(n_im_resid)])
    ax3.set_title(f"$>5\sigma$ Sources and bkg subtrackted ", fontsize=7)
    ax3.plot(np.linspace(xmin,xmax,1000),tgu_im_resid(np.linspace(xmin,xmax,1000)), color='k',lw=1, linestyle="--", label='BANE')

    if args.compare is True: 
        ax4.bar(bin_centers_im_gd, n_im_gd, color = 'C6', edgecolor = "none", width=(bins_im_gd[1]-bins_im_gd[0]),alpha=0.5) 
        ax4.plot(np.linspace(xmin,xmax,1000),gu_im_gd(np.linspace(xmin,xmax,1000)), color='k',lw=0.5, alpha=0.8, label='Fitted')
        ax4.set_title(f"New bkg subtracted, S/N image", fontsize=7)
        ax4.set_ylim([1.0,1.2*max(n_im_gd)])
        ax4.plot(np.linspace(xmin,xmax,1000),tgu_im_gd(np.linspace(xmin,xmax,1000)), color='k',lw=1, linestyle="--", label='BANE')

        ax5.bar(bin_centers_im_mask_gd, n_im_mask_gd, color = 'C6', edgecolor = "none", width=(bins_im_mask_gd[1]-bins_im_mask_gd[0]),alpha=0.5) 
        ax5.set_ylim([1.0,1.2*max(n_im_mask_gd)])
        ax5.set_title(f"New bkg subtracted, $>5\sigma$ sources masked", fontsize=7)
        ax5.plot(np.linspace(xmin,xmax,1000),gu_im_mask_gd(np.linspace(xmin,xmax,1000)), color='k',lw=0.5, alpha=0.8, label='Fitted')
        ax5.plot(np.linspace(xmin,xmax,1000),tgu_im_mask_gd(np.linspace(xmin,xmax,1000)), color='k',lw=1, linestyle="--", label='BANE')
        
        ax6.bar(bin_centers_im_resid_gd, n_im_resid_gd, color = 'C6', edgecolor = "none", width=(bins_im_resid_gd[1]-bins_im_resid_gd[0]),alpha=0.5) 
        ax6.plot(np.linspace(xmin,xmax,1000),gu_im_resid_gd(np.linspace(xmin,xmax,1000)), color='k',lw=0.5, alpha=0.8, label='Fitted')
        ax6.set_ylim([1.0,1.2*max(n_im_resid_gd)])
        ax6.set_title(f"$>5\sigma$ sources and new bkg subtracted", fontsize=7)
        ax6.plot(np.linspace(xmin,xmax,1000),tgu_im_resid_gd(np.linspace(xmin,xmax,1000)), color='k',lw=1, linestyle="--", label='BANE')

    for ax in axes_total:
        ax.axvline(x=0.0, lw=0.5, color='k', linestyle='-')
        ax.axvline(x=-1, lw=0.5, color='k', linestyle='--')
        ax.axvline(x=-2, lw=0.5, color='k', linestyle='-.')
        ax.axvline(x=-5, lw=0.5, color='k', linestyle=':')
        ax.axvline(x=1, lw=0.5, color='k', linestyle='--')
        ax.axvline(x=2, lw=0.5, color='k', linestyle='-.')
        ax.axvline(x=5, lw=0.5, color='k', linestyle=':')
        ax.set_yscale('log')
        ax.set_xlabel("$\sigma$")
        ax.set_xlim([xmin,xmax])
        # ax.legend()

    plt.savefig(f"{args.savedir}/noise_distribution_{savenm}.pdf",bbox_inches='tight')

    return 


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Script for plotting various histograms for paper"
    )
    parser.add_argument(
        '--imagenm',
        type=str,
        default="XG.fits",
        help="Image name to import, if rms and bkg aren't specified, then will find rms and bkg for image of same name"
    )
    parser.add_argument(
        '--datadir',
        type=str,
        default="data/",
        help="directory for the images to load (default=./data)"
    )
    parser.add_argument(
        '--savedir',
        type=str,
        default="plots/",
        help="directory to save the plots made (default=./plots)"
    )
    parser.add_argument(
        '--compare',
        action='store_true',
        default=False,
        help="Do you want to compare with the accurate background distributions? Default=False"
    )


    args = parser.parse_args()
    imstring = args.imagenm.split(".")[0]
    savenm_str = imstring.split("_")[1:]

# Image names I'll need: 
# _bkgsub_sigma.fits, _bkgsub_masked_sigma.fits, _bkgsub_resid_sigma.fits
# _priorsub_bkgsub_sigma.fits, _priorsub_bkgsub_masked_sigma.fits, _priorsub_bkgsub_resid_sigma.fits 
# _bkg.fits, _rms.fits, _priorsub_resid_bkg.fits, _priorsub_resid_rms.fits



    plot_distribution(
        imstring,
        savenm=imstring)


