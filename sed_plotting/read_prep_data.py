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
from matplotlib.ticker import ScalarFormatter

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

x_ticks=[80, 90, 100, 120, 130, 150, 200, 250]
majorticklength=2.5
minorticklength=2
tickwidth=0.5
cm = 1/2.54


REF_NU = 200.
CHANSN = ['072_080MHz', '080_088MHz', '088_095MHz', '095_103MHz', '103_111MHz', '111_118MHz', '118_126MHz',
 '126_134MHz', '139_147MHz', '147_154MHz', '154_162MHz', '162_170MHz', '170_177MHz', 
 '177_185MHz', '185_193MHz', '193_200MHz', '200_208MHz', '208_216MHz', '216_223MHz', '223_231MHz']

VARCHANSN = [
    "107",
    "115",
    "123",
    "130",
    "143",
    "150",
    "158",
    "166",
    "174",
    "181",
    "189",
    "197",
    "204",
    "212",
    "220",
    "227",
    ]
CHANSN_FINAL = [
    "076",
    "084",
    "092",
    "099",
    "107",
    "115",
    "122",
    "130",
    "143",
    "151",
    "158",
    "166",
    "174",
    "181",
    "189",
    "197",
    "204",
    "212",
    "220",
    "227",
    ]



GLEAMX_FREQS = np.array([np.mean([float(i) for i in name.replace('MHz','').split('_')[-2:]]) for name in CHANSN])
GLEAMX_INT = [f"int_flux_{c}" for c in CHANSN]
GLEAMX_ERR = [f"err_int_flux_{c}" for c in CHANSN]
GLEAMX_RMS = [f"local_rms_{c}" for c in CHANSN]

GLEAMX_FREQS_FINAL = np.array([float(i) for i in CHANSN_FINAL])
GLEAMX_INT_FINAL = [f"int_flux_{c}" for c in CHANSN_FINAL]
GLEAMX_ERR_FINAL = [f"err_int_flux_{c}" for c in CHANSN_FINAL]
GLEAMX_RMS_FINAL = [f"local_rms_{c}" for c in CHANSN_FINAL]

VAR_FREQS = np.array([float(c) for c in VARCHANSN])
VARY1_INT = [f"S_{c}_yr1" for c in VARCHANSN]
VARY2_INT = [f"S_{c}_yr2" for c in VARCHANSN]
VARY1_RMS = [f"local_rms_{c}_yr1" for c in VARCHANSN]
VARY2_RMS = [f"local_rms_{c}_yr2" for c in VARCHANSN]



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


def plot_sed(freq, flux, fluxerr, fit_params, src_row, outpath):
    # fig = plt.figure(figsize=(8*cm, 8*cm))
    fig, ax = plt.subplots(1,1, figsize=(8*cm,5.5*cm))

    # ax = fig.add_axes(111)



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
        color="hotpink"
    )

    legend = False

    if fit_params is not None:
        # legend = True
        if len(fit_params) == 2:
            
            fit_p0 = fit_params

            fit_res = curve_fit(power_law, freq, flux, fit_p0,sigma=fluxerr,absolute_sigma=True)

            no_samps = 1000
            samps = np.random.multivariate_normal(fit_res[0], fit_res[1],size=no_samps).swapaxes(0,1)

            models = power_law(nu[:,None], *samps)
            q16, q50, q84 = np.percentile(models, [16, 50, 84], axis=1)

            ax.plot(
                nu,
                q50,
                lw=0.5,
                color='hotpink',
                # label=f"Power-law, alpha {fit_params[1]:.3f}"
            )
            ax.fill_between(nu,q16,q84,alpha=0.3, color="hotpink")
            ax.text(0.04, 0.125, "Power Law", transform=ax.transAxes,bbox=dict(facecolor="white",alpha=0.8, edgecolor="black",boxstyle="round,pad=0.25"))


        elif len(fit_params) == 3: 
            fit_p0 = fit_params

            fit_res = curve_fit(curved_power_law, freq, flux, fit_p0,sigma=fluxerr,absolute_sigma=True)

            no_samps = 1000
            samps = np.random.multivariate_normal(fit_res[0], fit_res[1],size=no_samps).swapaxes(0,1)

            models = curved_power_law(nu[:,None], *samps)
            q16, q50, q84 = np.percentile(models, [16, 50, 84], axis=1)

            ax.plot(
                nu,
                q50,
                lw=0.5,
                color='hotpink',
                # label=f"Power-law, alpha {fit_params[1]:.3f}"
            )
            ax.fill_between(nu,q16,q84,alpha=0.3, color="hotpink")
            # ax.text(0.315, 0.125, "Curved Power Law", transform=ax.transAxes,bbox=dict(facecolor="white",alpha=0.8, edgecolor="black",boxstyle="round,pad=0.25"))
            ax.text(0.315, 0.825, "Curved Power Law", transform=ax.transAxes,bbox=dict(facecolor="white",alpha=0.8, edgecolor="black",boxstyle="round,pad=0.25"))

    if args.var is not None: 
        src_var_row = search_cat(src_row, args.varcat, search_radius=1.5*u.arcmin)
        if src_var_row.size == 0: 
            logger.warning("There was no crossmatch data for this source, not plotting")
        else: 
            try:
                freq1, flux1, errflux1 = get_freq_flux_err(src_var_row, freqs=VAR_FREQS, fluxes=VARY1_INT, errs = VARY1_RMS, internal_scale=0.03, apply_mask=False)
                freq2, flux2, errflux2 = get_freq_flux_err(src_var_row, freqs=VAR_FREQS, fluxes=VARY2_INT, errs = VARY2_RMS, internal_scale=0.03, apply_mask=False)
            except: 
                freq_var, int_flux_var, err_flux_var = get_freq_flux_err(src_var_row, apply_mask = False, internal_scale=0.02)
            
            if args.var==1:
                logger.debug("GOING TO PLOT FOR 2013")
                ax.errorbar(
                freq1,
                flux1,
                yerr=errflux1,
                ls='None',
                marker='.',
                color="darkviolet",
                alpha=0.5,
                label="2013"
                )
            elif args.var==2:
                logger.debug("GOING TO PLOT FOR 2014")
                ax.errorbar(
                freq2,
                flux2,
                yerr=errflux2,
                ls='None',
                marker='.',
                color="mediumblue",
                alpha=0.5,
                label="2014"
                )  
            elif args.var==0:
                logger.debug("GOING TO PLOT 2013,2014")
                ax.errorbar(
                freq1,
                flux1,
                yerr=errflux1,
                ls='None',
                marker='.',
                color="darkviolet",
                alpha=0.5,
                label="2013"
                )
                ax.errorbar(
                freq2,
                flux2,
                yerr=errflux2,
                ls='None',
                marker='.',
                color="mediumblue",
                alpha=0.5,
                label="2014"
                )
            elif args.var == 3: 
                ax.errorbar(
                freq_var,
                int_flux_var,
                yerr=err_flux_var,
                ls='None',
                marker='.',
                color="darkviolet",
                alpha=0.5,
                )                


    if legend is True:
        ax.legend()


    
    coord_src = source(SkyCoord(src_row.ref_ra*u.deg, src_row.ref_dec*u.deg))
    ax.loglog()
    ax.set(
        xlabel='Frequency (MHz)',
        ylabel='Integrated Flux (Jy)',
        title=coord_src.tex
    )
    ax.tick_params(
        axis="both", which="major", direction="in", length=majorticklength, width=tickwidth, pad=5, labelsize=7
    )
    ax.tick_params(
        axis="both", which="minor", direction="in", length=minorticklength, width=tickwidth, pad=5, labelsize=7
    )
    ax.tick_params(axis="x", which="both", pad=5, labelsize=7)
    ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.xaxis.set_minor_formatter(ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.yaxis.set_minor_formatter(ScalarFormatter(useMathText=True))

    savenm = str(coord_src)
    savenm = savenm.replace(" ","_")
    fig.tight_layout()
    try: 
        fig.savefig(f"{outpath}/{savenm}.png")
    except FileNotFoundError:
        logger.debug(f"Couldnt save the sed because there's no folder: {outpath}, will just save here...")
        fig.savefig(f"./{savenm}.pdf")

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

def try_ransac_fit(freq, flux,err_flux,p0):
    chi2 = np.zeros_like(freq)
    for i in range(len(freq)):
        idx2 = np.arange(0, len(freq))
        idx2 = np.delete(idx2, i)
        temp_freq = freq[idx2]

        temp_model, covar_temp_model = curve_fit(power_law,temp_freq, flux[idx2],p0=p0,sigma=err_flux[idx2],absolute_sigma=True)
        chi2[i] = np.sum(
            ((flux[idx2] - power_law(temp_freq, *temp_model)) / err_flux[idx2])**2
            )

    return chi2


def fit_pl(freq, flux, err_flux):
    p0 = (np.median(flux), -0.8)
    try: 
        fit_res_pl = curve_fit(
            power_law, 
            freq, 
            flux, 
            p0=p0, 
            sigma=err_flux,
            absolute_sigma=True
        )
        best_p_pl, covar_pl = fit_res_pl

        err_p_pl = np.sqrt(np.diag(covar_pl))
        dof_pl = len(freq) - 2
        chi2_pl = np.sum(
            ((flux - power_law(freq, *best_p_pl)) / err_flux)**2
            )
        rchi2_pl = chi2_pl / dof_pl
        ransac_chi2 = try_ransac_fit(freq,flux,err_flux,p0)
        logger.debug(f"The ransac chi2 is: {ransac_chi2}")
        if chi2_pl > chi2.ppf(0.999, dof_pl):
            logger.debug(f"PL fit would've been set to none before any checks: {chi2_pl} > {chi2.ppf(0.999, dof_pl)}")
            ransac_chi2 = try_ransac_fit(freq,flux,err_flux,p0)
            # logger.debug(f"The ransac chi2 is: {ransac_chi2}")
            fit_pl = None
            rchi2_pl = None
            chi2_pl = None 
        else: 
            fit_pl = best_p_pl
    except RuntimeError:
        fit_pl =  None
        rchi2_pl = np.inf 

    return best_p_pl, chi2_pl, rchi2_pl

def fit_cpl(freq, flux, err_flux):
    try: 
        p0 = (np.median(flux), -0.8, 0)
        fit_res_cpl = curve_fit(
            curved_power_law,
            freq,
            flux,
            p0=p0,
            sigma=err_flux,
            absolute_sigma=True
        )
        best_p_cpl, covar_cpl = fit_res_cpl
        err_p_cpl = np.sqrt(np.diag(covar_cpl))
        dof_cpl = len(freq) - 3
        chi2_cpl = np.sum(
            ((flux - curved_power_law(freq, *best_p_cpl)) / err_flux)**2
            )
        rchi2_cpl = chi2_cpl / dof_cpl  
        if chi2_cpl > chi2.ppf(0.99, dof_cpl):
            logger.debug(f"CPL would've been set to none before any checks ")
            fit_cpl = None
            chi2_cpl = None
            rchi2_cpl = None 
        if np.abs(best_p_cpl[2]) < 0.2 or np.abs(best_p_cpl[2])/err_p_cpl[2] <3:
            logger.debug(f"CPL would've been cut because the curvature was poor or flipped")
            logger.debug(f"{round(np.abs(best_p_cpl[2]),3)}<0.2 or {np.abs(best_p_cpl[2])/err_p_cpl[2]} < 3")
            fit_cpl = None 
            rchi2_cpl = None
            chi2_cpl = None
    except RuntimeError:
        best_p_cpl = None
        rchi2_cpl = np.inf
        chi2_cpl = np.inf
    except:
        print("other issue!?!? ")
        best_p_cpl = None
        rchi2_cpl = np.inf
        chi2_cpl = np.inf
    return best_p_cpl, chi2_cpl, rchi2_cpl

def read_fit_params(row,pl,cpl,try_fit=False):

    if np.isfinite(cpl):
        fit_params = [row[i] for i in ["cpl_norm","cpl_alpha", "cpl_q"]]
    elif np.isfinite(pl): 
        fit_params = [row[i] for i in ["pl_norm","pl_alpha"]]
    else: 
        logger.debug(f"No fitting params found for {row['src_name']}")
        fit_params = None 
    if try_fit is True: 
        src_name = row["src_name"]
        logger.warning(f"Both params are none?? No model for {src_name}")
        logger.debug(f"Going to try fitting just a pl")
        freq, int_flux, err_flux = get_freq_flux_err(row)
        if len(freq) < 15:
            logger.debug(f"Would've been cut from few freqs to fit: len(freq) = {len(freq)}")
        best_p_pl, chi2_pl, rchi2_pl = fit_pl(freq, int_flux, err_flux)
        best_p_cpl, chi2_cpl, rchi2_cpl = fit_cpl(freq, int_flux, err_flux)   

        
        logger.debug(f"Comparing the fits: rchi2_pl={rchi2_pl} vs rchi2_cpl={rchi2_cpl}")
        if rchi2_cpl is None and rchi2_pl is not None: 
            logger.debug(f"Only fit PL so using this as preferred!")
            fit_params = [best_p_pl[0], best_p_pl[1]]
        elif rchi2_pl is None and rchi2_cpl is not None: 
            logger.debug(f"Only fit CPL so using this as preferred!")
            fit_params = [best_p_cpl[0], best_p_cpl[1], best_p_cpl[2]]
        elif rchi2_cpl is None and rchi2_pl is None: 
            logger.warning(f"Tried to fit both but failed both? Don't know what you want from me. ")
            fit_params = None 
        elif rchi2_cpl < rchi2_pl:
            logger.debug(f"Fit both cpl and pl and cpl preferred, using that as model")
            fit_params = [best_p_cpl[0], best_p_cpl[1], best_p_cpl[2]]
        elif rchi2_pl < rchi2_cpl:
            logger.debug(f"Fit both cpl and pl and pl preferred, using that as model")
            fit_params = [best_p_pl[0], best_p_pl[1]]
        elif rchi2_cpl == rchi2_pl:
            logger.warning(f"Unbelievable chances, they're equally fit, preference for pl ")
            fit_params = [best_p_pl[0], best_p_pl[1]]

    return fit_params

def get_freq_flux_err(row, freqs=GLEAMX_FREQS, fluxes=GLEAMX_INT, errs = GLEAMX_ERR, apply_mask = True, internal_scale=0.02):
    
    freq = freqs
    int_flux = np.array([row[i] for i in fluxes])
    err_flux = np.array([row[i] for i in errs])
    
    # Gleam-x internal flux error
    err_flux = np.sqrt( err_flux ** 2 + (int_flux * internal_scale)**2)
    if np.nanmean(err_flux) >= 1:
        logger.debug("This is poorly fit, using rms for errs but watch out.")
        rms = np.array([row[i] for i in GLEAMX_RMS])
        err_flux = np.sqrt(rms **2 +(int_flux * internal_scale)**2)

        
    if apply_mask:
        if errs == GLEAMX_ERR:
            if (src_row.ref_ra < 45) & (src_row.ref_ra > 30) & (src_row.ref_dec < -5) & (src_row.ref_dec > -7):
                logger.debug(f"This is in a bad region, flagging freq 99")
                mask_reg1 = (freq == 99.)
                freq = freq[~mask_reg1]
                int_flux = int_flux[~mask_reg1]
                err_flux = err_flux[~mask_reg1]
            if (src_row.ref_ra < 55) & (src_row.ref_ra > 45 )& (src_row.ref_dec < -13) & (src_row.ref_dec > -15):
                logger.debug(f"This is in a bad region, flagging freq 130")
                mask_reg2 = (freq == 130.)
                freq = freq[~mask_reg2]
                int_flux = int_flux[~mask_reg2]
                err_flux = err_flux[~mask_reg2]
            if (src_row.ref_ra < 72) & (src_row.ref_ra >66) & (src_row.ref_dec <2) & (src_row.ref_dec > -2):
                logger.debug(f"This is in a bad region, flagging freq 130")
                mask_reg3 = (freq == 130.)
                freq = freq[~mask_reg3]
                int_flux = int_flux[~mask_reg3]
                err_flux = err_flux[~mask_reg3]                       
            if (src_row.ref_ra < 76) & (src_row.ref_ra >70) & (src_row.ref_dec <9) & (src_row.ref_dec > 5):
                logger.debug(f"This is in a bad region, flagging freq 99")
                mask_reg4 = (freq == 99.)
                freq = freq[~mask_reg4]
                int_flux = int_flux[~mask_reg4]
                err_flux = err_flux[~mask_reg4]   
            if (src_row.ref_ra < 8 )& (src_row.ref_ra > 5) & (src_row.ref_dec < 8 )& (src_row.ref_dec > 6):
                logger.debug(f"This is in a bad region, flagging freq 99")
                mask_reg5 = (freq == 99.)
                freq = freq[~mask_reg5]
                int_flux = int_flux[~mask_reg5]
                err_flux = err_flux[~mask_reg5]  

        mask = np.isfinite(int_flux) & (int_flux > 0) \
            & np.isfinite(err_flux) & (err_flux > 0)
        
        freq = freq[mask]
        int_flux = int_flux[mask]
        err_flux = err_flux[mask]
    

    
    return np.squeeze(freq), np.squeeze(int_flux), np.squeeze(err_flux)


def search_cat(source, catalogue, search_radius=1*u.arcmin):

    target = SkyCoord(src_row.ref_ra*u.deg, src_row.ref_dec*u.deg, frame="icrs")
    try: 
        catalogue = Table.read(catalogue)
        catalogue["RA"] = catalogue["RA"].astype(float)
        catalogue["Dec"] = catalogue["Dec"].astype(float)
        catalogue = catalogue.to_pandas()
        cats = SkyCoord(catalogue.RA, catalogue.Dec, frame="fk5", unit=u.deg)
    except:
        catalogue = Table.read(catalogue).to_pandas()
        cats = SkyCoord(catalogue.ref_ra, catalogue.ref_dec,frame="fk5", unit=u.deg)
    
    try: 
        distances = target.separation(cats)
        idc = distances <= search_radius
        src_inds = np.where(idc)[0]
        if len(src_inds) > 1: 
            logger.warning("Found multiple matches! Which one do you want?")
            logger.warning(f"Indecies of matches: {src_inds}: {distances[idc]}")
            logger.debug(f"Going to choose nearer source, but may be wrong/resolved source")
            minid = src_inds[np.argmin(distances[idc])]
            return catalogue.iloc[minid]
        else: 
            return catalogue.iloc[idc]
        
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
    parser.add_argument("-var", default=None,type=int, help="Do you also want to plot 2013, 2014 from Ross et al. 2020? 0 = both, 1 = 2013, 2=2014, will also try fitting just pl and cpl for each yr" )
    parser.add_argument("-varcat", type=str, default="/data/MWA/master_pop_extended.fits", help="If plotting the var, this is the dest of var catalogue.")

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
            
            try: 
                coords = Table.read(args.coord).to_pandas()
                cross_coords = SkyCoord(coords.RA, coords.Dec, unit=u.deg)
            except: 
                coords = pd.read_csv(args.coord, sep=",", header=None)
                cross_coords = SkyCoord(np.array(coords[0]), np.array(coords[1]), unit=u.deg)
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
        # mask_fit = df.loc[np.isfinite(df.cpl_chi2) | np.isfinite(df.pl_rchi2)]
        max_num = df.shape[0]
        mask_fit = df 
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
        
        freq, int_flux, err_flux = get_freq_flux_err(src_row, apply_mask = True, internal_scale=0.02)
        p0 = (np.median(int_flux), -0.8)
        ransac_chi2 = try_ransac_fit(freq,int_flux,err_flux,p0)
        bad_freq_mask = ransac_chi2 < (np.mean(ransac_chi2)-(5*np.std(ransac_chi2)))
        if any(bad_freq_mask):
            logger.warning(f"THERE IS A BAD FREQ FOR THIS SOURCE: {coordnm}")
            logger.warning(f"You should check mosaic {freq[bad_freq_mask]}")
        fit_params = read_fit_params(src_row,src_row.pl_chi2,src_row.cpl_chi2)
        
        plot_sed(freq, int_flux, err_flux, fit_params, src_row, args.d)
