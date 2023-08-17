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

def fit_hists(data, stdev = 1.0, nsig=1):
    bins=100
    xmin = -nsig*stdev
    xmax = nsig*stdev

    n, binss = np.histogram(data,bins = np.linspace(xmin,xmax,bins))
    n = np.array(n,dtype="float32")/len(n)
    bin_centers = (binss+((binss[1]-binss[0])/2))[:-1]
    g_init = models.Gaussian1D(amplitude=1., mean=0., stddev=stdev)
    gu = fit_g(g_init, bin_centers, n)

    return n, bin_centers, binss, gu 


# todo: make better savename and title 
def plot_distribution(n, bin_centers, binss, gu, stdev, mean_rms=1.0, mean_bkg=0., title="", sigma=False,savenm=""):


    if sigma is False:
        tgu = copy.deepcopy(gu)
        if args.bkgsub is False: 
            tgu.mean.value = mean_bkg
        else: 
            tgu.mean.value = 0
        tgu.stddev.value = mean_rms
        nsig=1
        xlabel="Jy/beam"
    else:
        tgu = copy.deepcopy(gu)
        tgu.mean.value = 0
        tgu.stddev.value = 1.

        # stdev=1
        nsig=5
        xlabel="$\sigma$"


    xmin = -nsig*stdev
    xmax = nsig*stdev

    fig = plt.figure(figsize=(10*cm, 10*cm))
    ax = fig.add_subplot(111)
    ax.bar(bin_centers, n, color = 'C6', edgecolor = "none", width=(binss[1]-binss[0]),alpha=0.5) 
    ax.set_yscale('log')
    ax.set_xlabel(f"Pixel distribution ({xlabel})")
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([1.0,1.2*max(n)])

    ax.plot(np.linspace(xmin,xmax,1000),gu(np.linspace(xmin,xmax,1000)), color='k',lw=0.5, alpha=0.8, label='Pixel dist fit')

    ax.plot(np.linspace(xmin,xmax,1000),tgu(np.linspace(xmin,xmax,1000)), color='k',lw=1, linestyle="--", label='BANE measured')
    ax.legend()
    ax.set_title(f"Pixel distribution {title}")
    plt.savefig(f"noise_distribution_{savenm}.png")


    return 


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Script for plotting various histograms for paper"
    )
    parser.add_argument(
        '--imagenm',
        type=str,
        default="XG",
        help="input image to plot the pxl distribution to (default=XG), not exclude the .fits to make easier for rms and bkg loading"
    )
    parser.add_argument(
        '--rms',
        help="The prenm for the rms, really only need to define if you want to measure using a different rms and bkg to the one that would be default from BANE for the image default=imagenm_rms.fits"
    )
    parser.add_argument(
        '--bkg',
        help="See rms description, but it's for bkg default=imagenm_bkg.fits"
    )
    parser.add_argument(
        '--sigma',
        default=False,
        help="Are you using the sigma images or just raw pixels"
    )
    parser.add_argument(
        '--stdev',
        default=0.025,
        help="the std for fitting, 0.25 for red, 0.025 for white, 1 for sigma"
    )
    parser.add_argument(
        "--bkgsub",
        default=False,
        help="If using a bkgsub one, then doesn't update the gaus fit to bkgmean"
    )
    parser.add_argument(
        '--savenm',
        type=str,
        default="",
        help="output save name"
    )
    args = parser.parse_args()

    if args.rms:
        rms = fits.getdata(f"{args.rms}")
    else:
        rms = fits.getdata(f"{args.imagenm}_rms.fits")
    if args.bkg:
        bkg = fits.getdata(f"{args.bkg}")
    else:
        bkg = fits.getdata(f"{args.imagenm}_bkg.fits")

    image = fits.getdata(f"{args.imagenm}.fits")
    
    
    mean_rms = np.nanmean(rms)
    mean_bkg = np.nanmean(bkg)

    if args.sigma is False:
        print("Sigma is false, raw pxl dist")
        nsig=1
        stdev=args.stdev
        # stdev=0.25
    else:
        nsig=5
        stdev=1

    n, bin_centers, binss, gu = fit_hists(image,nsig=nsig,stdev=stdev)

    plot_distribution(
        n, 
        bin_centers, 
        binss, 
        gu,
        stdev=stdev,
        mean_rms=mean_rms, 
        mean_bkg=mean_bkg, 
        title=args.imagenm, 
        sigma=args.sigma,
        savenm=args.savenm)


