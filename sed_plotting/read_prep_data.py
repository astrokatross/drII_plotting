#! /usr/bin/env python3
from astropy.table import Table, Column
from astropy.io.votable import from_table, writeto
from astropy.coordinates import SkyCoord, search_around_sky
import numpy as np
import astropy.units as u
import argparse
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import fits
import logging 
import pandas as pd
import matplotlib.pyplot as plt 
import os
from scipy.stats import chi2

__author__ = "Kat Ross"
__date__ = "2023-04-14"

"""
Script to organise the data from various catalogues and make a database style structure for easier handling per source in things like sed fitting, variability analysis etc 
"""

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(module)s:%(levelname)s:%(lineno)d %(message)s")
logger.setLevel(logging.INFO)
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 8,
    "axes.titlesize": 8
})


majorticklength=2.5
minorticklength=2
tickwidth=0.5

REF_NU = 200.
CHANSN = ['072_080MHz', '080_088MHz', '088_095MHz', '095_103MHz', '103_111MHz', '111_118MHz', '118_126MHz',
 '126_134MHz', '139_147MHz', '147_154MHz', '154_162MHz', '162_170MHz', '170_177MHz', 
 '177_185MHz', '185_193MHz', '193_200MHz', '200_208MHz', '208_216MHz', '216_223MHz', '223_231MHz']

GLEAMX_FREQS = np.array([np.mean([float(i) for i in name.replace('MHz','').split('_')[-2:]]) for name in CHANSN])
GLEAMX_INT = [f'int_flux_{c}' for c in CHANSN]
GLEAMX_ERR = [f'err_int_flux_{c}' for c in CHANSN]

class source:
    def __init__(self, *args, **kwargs):
        if len(args) == 1 and len(kwargs) == 0:
            self.pos = args[0]
        else:
            self.pos = SkyCoord(*args, **kwargs)
        
        
    def _make_str(self, tex=False):
        _format = None if tex is False else 'latex'
        
        return (f'GLEAM-X '\
                f'J{self.pos.ra.to_string(unit=u.hourangle, sep="", precision=1, pad=True, format=_format)}' \
                f'{self.pos.dec.to_string(sep="", precision=0, alwayssign=True, pad=True, format=_format)}')

    def __str__(self):
        return self._make_str()
    
    def __repr__(self):
        return str(self)
     
    @property
    def tex(self):
        return self._make_str(tex=True)


def plot_sed(freq, flux, fluxerr, fit_params, coord_src, outpath):
    fig, ax = plt.subplots(1,1)

    nu = np.geomspace(
        np.min(freq),
        np.max(freq),
        100
    )

    ax.errorbar(
        freq,
        flux,
        yerr=fluxerr,
        ls='None',
        marker='.',
        color="C6"
    )

    legend = False

    if fit_params is not None:
        legend = True
        
        if len(fit_params) == 2:

            ax.plot(
                nu,
                power_law(nu, *fit_params),
                ls='--',
                color='C9',
                label=f"Power-law, alpha {fit_params[1]:.3f}"
            )
        elif len(fit_params) == 3: 
            ax.plot(
                nu,
                curved_power_law(nu, *fit_params),
                ls=':',
                color='C4',
                label=f'Curved power-law, q {fit_params[2]:.3f}'
            )



    if legend is True:
        ax.legend()

    ax.loglog()
    ax.set(
        xlabel='Frequency (MHz)',
        ylabel='Integrated Flux (Jy)',
        title=coord_src.tex
    )
    ax.tick_params(
        axis="both", which="major", direction="in", length=majorticklength, width=tickwidth, pad=5
    )
    ax.tick_params(
        axis="both", which="minor", direction="in", length=minorticklength, width=tickwidth
    )
    savenm = str(coord_src)
    savenm = savenm.replace(" ","_")
    fig.tight_layout()
    fig.savefig(f"{outpath}/{savenm}.png")
    plt.close(fig)


def power_law(nu, norm, alpha):
    return norm * (nu / REF_NU) ** alpha


def curved_power_law(nu, norm, alpha, q):
    spec_nu = nu / REF_NU
        
    return norm * spec_nu ** alpha * \
            np.exp(q * np.log(spec_nu)**2)


def create_source_name(row):
    coord_src = source(
        SkyCoord(row.ref_ra*u.deg, row.ref_dec*u.deg)
    )

    return str(coord_src)


# probs don't need unless i want to expand to do something useful.... 
def read_catalogue(cat:str, source:source, colnames, search_radius = 1*u.degree, cat_ext = "fits"):

    try: 
        cats = Table.read(cat, format=cat_ext).to_pandas()
        catalogue = SkyCoord(cats.ref_ra, cats.ref_dec, frame="fk5", unit=u.deg)
    except FileNotFoundError:
        logger.debug("Can't find the catalogue you're trying to search!?")
        return np.nan 
    except: 
        logger.debug("Can't make catalogue?")
        return np.nan 

    return cats

def read_fit_params(row,pl,cpl):

    if np.isfinite(cpl):
        fit_params = [row[i] for i in ["cpl_norm","cpl_alpha", "cpl_q"]]
    elif np.isfinite(pl): 
        fit_params = [row[i] for i in ["pl_norm","pl_alpha"]]
    else: 
        src_name = row["src_name"]
        logger.warning(f"Both params are none?? No model for {src_name}")
        logger.debug(f"Going to try fitting just a pl")
        freq, int_flux, err_flux = get_freq_flux_err(row)
        p0 = (np.median(int_flux), -0.8)
        try:
            fit_res_pl = curve_fit(
                power_law,
                freq,
                int_flux,
                p0=p0,
                sigma=err_flux,
                absolute_sigma=True
            )
            best_p_pl, covar_pl = fit_res_pl
            err_p_pl = np.sqrt(np.diag(covar_pl))
            dof_pl = len(freq) - 2
            chi2_pl = np.sum(
                ((int_flux - power_law(freq, *best_p_pl)) / err_flux)**2
                )
            rchi2_pl = chi2_pl / dof_pl   
        except RuntimeError:
            fit_pl =  None
            rchi2_pl = np.inf 

        try: 
            p0 = (np.median(int_flux), -0.8, 0)
            fit_res_cpl = curve_fit(
                curved_power_law,
                freq,
                int_flux,
                p0=p0,
                sigma=err_flux,
                absolute_sigma=True
            )
            best_p_cpl, covar_cpl = fit_res_cpl
            err_p_cpl = np.sqrt(np.diag(covar_cpl))
            dof_cpl = len(freq) - 3
            chi2_cpl = np.sum(
                ((int_flux - curved_power_law(freq, *best_p_cpl)) / err_flux)**2
                )
            rchi2_cpl = chi2_cpl / dof_cpl            
        except RuntimeError:
            fit_cpl = None
            rchi2_cpl = np.inf  
        
        if rchi2_cpl < rchi2_pl:
            logger.debug(f"Fit both cpl and pl and cpl preferred, using that as model")
            if np.abs(best_p_cpl[2]) < 0.2 or np.abs(best_p_cpl[2])/err_p_cpl[2] < 3:
                logger.debug(f"This fit would've been cut in the original fitting")
            elif chi2_cpl > chi2.ppf(0.99, dof_cpl):
                logger.debug(f"This would've been cut from original fitting because it's not a confident fit")
            fit_params = [best_p_cpl[0], best_p_cpl[1], best_p_cpl[2]]
        elif rchi2_pl < rchi2_cpl:
            logger.debug(f"Fit both cpl and pl and pl preferred, using that as model")
            if chi2_pl > chi2.ppf(0.99, dof_pl):
                logger.debug("This would've been cut in original fitting from poor fit")
            fit_params = [best_p_pl[0], best_p_pl[1]]
        elif fit_cpl is None and fit_pl is None: 
            logger.warning(f"Tried to fit both but failed both? Don't know what you want from me. ")
            fit_params = None 
        elif rchi2_cpl == rchi2_pl:
            logger.warning(f"Unbelievable chances, they're equally fit, preference for pl ")
            fit_params = [best_p_pl[0], best_p_pl[1]]

    return fit_params

def get_freq_flux_err(row, apply_mask = True, internal_scale=0.02):
    
    freq = GLEAMX_FREQS
    int_flux = np.array([row[i] for i in GLEAMX_INT])
    err_flux = np.array([row[i] for i in GLEAMX_ERR])
    
    # Gleam-x internal flux error
    err_flux = np.sqrt( err_flux ** 2 + (int_flux * internal_scale)**2)
        
    if apply_mask:
        mask = np.isfinite(int_flux) & (int_flux > 0) \
            & np.isfinite(err_flux) & (err_flux > 0)
        
        freq = freq[mask]
        int_flux = int_flux[mask]
        err_flux = err_flux[mask]
    
    return np.squeeze(freq), np.squeeze(int_flux), np.squeeze(err_flux)


def search_cat(source: source, catalogue, search_radius=1*u.arcmin):

    target = SkyCoord(ra=source.ra,dec=source.dec)
    try: 
        distances = target.separation(catalogue, search_radius)
        idc = np.where(distances < search_radius)[0]
        if len(idc) > 1: 
            logger.warning("Found multiple matches! Which one do you want?")
            logger.warning(f"Indecies of matches: {idc}: {distances[idc]}")
            raise IndexError
        else: 
            return catalogue[idc]
        
    except: 
        logger.debug(f"Couldn't find {source.name} in catalogue")
        raise Exception 

    # TODO/WHEREIM UP TO: Use the index from above to get all entries from catalogue for hte source and then create array based on the column names (also add that as variable) then you can return the array with values from the column names 

     


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform fitting and other todos on catalogue')

    parser.add_argument('table', type=str, help='Path to table to read')
    parser.add_argument('-p','--plot', default=False, action='store_true', help='Make plots for SEDs with a successful fit')
    parser.add_argument('-v', '--verbose', action='store_true', help='More output information')
    parser.add_argument('-num', default=None, type=int, help="How many plots do you want? If given number, will randomly select n rows to plot")
    parser.add_argument('-d', default="SEDs", help="directory to save the seds, default (SEDs) ")
    parser.add_argument("-src", default=None,type=str, help="src name to search and plot")
    parser.add_argument("-coord", default=None, type=str, help="Coordinate string to use for a crossmatch and plots all srcs within 2arcmin. assumes hh:mm:ss +dd:mm:ss format")

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.debug(f"Reading the table in now, it's big and will take ages!")
    df = Table.read(args.table).to_pandas()
    logger.debug(f"Have successully read in!")

    # TODO: currently assumes table has columns RA and Dec in degs, maybe add option to change that if not? 
    if args.coord: 
        logger.debug(f"Given coordinates, will perform crossmatch")
        catalogue = SkyCoord(df.ref_ra, df.ref_dec, frame="fk5", unit=u.deg)
        if os.path.isfile(args.coord):
            logger.debug(f"Found file of coords, will cross match for all coords ")
            coords = Table.read(args.coord).to_pandas()
            cross_coords = SkyCoord(coords.RA, coords.Dec, unit=u.deg)
            idxc, idxcatalog, d2d, d3d = catalogue.search_around_sky(cross_coords, 2*u.arcmin)
            plt_inds = idxcatalog 
            logger.debug(f"Found {len(plt_inds)} sources in crossmatch of coords ")
        else: 
            logger.debug(f"Given coordinates, will plot all srcs in crossmatch within 2armin of {args.coord}")
            cross_coord = SkyCoord(args.coord, unit=(u.hourangle, u.deg))
            

            d2d = cross_coord.separation(catalogue)
            catmask = d2d <= 2*u.arcmin
            plt_inds = np.where(catmask)[0]
            if len(plt_inds) == 0:
                logger.warning(f"There are no matches for 2arcmin... going to try again with 5... ")
                catmask = d2d <= 5*u.arcmin
                plt_inds = np.where(catmask)[0]
                if len(plt_inds) == 0:
                    logger.warning(f"Still found no sources... too bad")
                    pass 
        mask_fit = df
        logger.debug(f"Found {len(plt_inds)} sources in crossmatch, plotting them all")
    elif args.src is not None and args.num is None:
        src_row = df.loc[df.src_name == args.src.encode()]
        plt_inds = src_row.index.values
        mask_fit = df
        logger.debug(f"Only plotting a single source: {args.src}")
    elif args.num is not None and args.src is None:
        plot_num = args.num
        logger.debug(f"Going to plot only subset: {plot_num} sources")
        mask_fit = df.loc[np.isfinite(df.cpl_chi2) | np.isfinite(df.pl_rchi2)]
        max_num = mask_fit.shape[0]
        logger.debug(f"Total with a fit: {max_num}")
        plt_inds = np.random.randint(0,max_num,plot_num)
    else:
    # if args.numum is None and args.src is None:
        plot_num = len(df)
        logger.debug(f"Going to plot everything: {plot_num} sources")


    for i in plt_inds:
        src_row = mask_fit.iloc[i]
        coord_src = source(SkyCoord(src_row.ref_ra*u.deg, src_row.ref_dec*u.deg))
        coordnm = str(coord_src).replace(" ","_")
        logger.debug(f"Plotting SED for {coordnm}")
        
        freq, int_flux, err_flux = get_freq_flux_err(src_row, apply_mask = False, internal_scale=0.02)

        fit_params = read_fit_params(src_row,src_row.pl_chi2,src_row.cpl_chi2)

        plot_sed(freq, int_flux, err_flux, fit_params, coord_src, args.d)