#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from astropy.io.votable import from_table, writeto
from astropy.coordinates import SkyCoord, search_around_sky
import astropy.units as u
import argparse
from scipy.optimize import curve_fit
from astropy.io import fits
import logging 
import pandas as pd
import os
from scipy.stats import chi2
from matplotlib.ticker import ScalarFormatter, FuncFormatter
import read_prep_data
import pickle 
import sed_models
import sed_fitting
import scipy.special as special
import cmasher as cmr

plt.rcParams["font.family"] = "serif"
plt.rcParams["axes.grid"] = False

majorticklength=2.5
minorticklength=2
tickwidth=0.5

data_dir = "/data/gleam_x/drII_plotting/data"
save_dir = "/data/gleam_x/drII_plotting/plots/SEDs"

u.hms = u.def_unit("h:m:s")
plot_models = ["curve", "singhomobremss", "doubssa","plbreak"]

CHANS = [
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

GLEAMX_FREQ = np.array([float(i) for i in CHANS])
GLEAMX_INT = [f"int_flux_{c}" for c in CHANS]
GLEAMX_ERR = [f"err_int_flux_{c}" for c in CHANS]
GLEAMX_RMS = [f"local_rms_{c}" for c in CHANS]
x_ticks=[70, 100, 150, 200, 300, 500, 800, 1400, 3000, 5000, 10000]

def get_freq_flux_err(row, freqs=GLEAMX_FREQ, fluxes=GLEAMX_INT, errs = GLEAMX_ERR, internal_scale=0.02):
    
    freq = [i*u.MHz for i in freqs]
    int_flux = np.array([row[i] for i in fluxes])
    int_flux = [i*u.Jy for i in np.squeeze(int_flux)]
    err_flux = np.array([row[i] for i in errs])
    err_flux = [i*u.Jy for i in np.squeeze(err_flux)]

    return freq, int_flux, err_flux


def plot_sed(freq, flux, err_flux, fit_params, source_name, outpath="/data/gleam_x/drII_plotting/plots/", plot_extras = None, survey_legend = False):
    fig, ax = plt.subplots(1,1)

    freq = [i.to(u.GHz) for i in freq]
    fit_freq = np.array([i.value for i in freq])
    fit_flux = np.array([i.value for i in flux])
    fit_err_flux = np.array([i.value for i in err_flux])
    fluxnegind = np.where((fit_flux <= 0))[0]
    
    if len(fluxnegind) > 0:
        print(f"There are some negative fluxes???")
        fluxposind = np.where(fit_flux > 0)[0]
        fit_flux = fit_flux[fluxposind]
        fit_err_flux = fit_err_flux[fluxposind]


    
    nu = np.geomspace(
        np.min(fit_freq),
        20,
        100
    )

    ax.errorbar(
        [i.to(u.GHz).value for i in freq],
        [i.to(u.Jy).value for i in flux],
        yerr=[i.to(u.Jy).value for i in err_flux],
        ls='None',
        marker='.',
        color="hotpink",
        label="GLEAM-X"
    )
    if plot_extras is not None: 
        for nm in plot_extras.keys():
            if nm in ["atca"]:
                extra_dict = plot_extras[nm]
                if nm == 'gleamyr1':
                    color="mediumblue"
                elif nm == "gleamyr2":
                    color="darkviolet"
                else: 
                    color="k"
                # if len(extra_dict['flux']) > 1: 
                plot_freq = [i.to(u.GHz) for i in extra_dict["freq"]]
                plot_freq = np.array([i.value for i in plot_freq])
                plot_flux = [i.to(u.Jy) for i in extra_dict["flux"]]
                plot_flux = np.array([i.value for i in plot_flux])
                plot_errflux = [i.to(u.Jy) for i in extra_dict["err_flux"]]
                plot_errflux = np.array([i.value for i in plot_errflux])
                # else: 
                #     plot_freq = extra_dict['freq'].to(u.GHz).value
                #     plot_flux = extra_dict['flux'].to(u.Jy).value
                #     plot_errflux = extra_dict['err_flux'].to(u.Jy).value
                if survey_legend is True: 
                    ax.errorbar(
                        plot_freq,
                        plot_flux,
                        yerr=plot_errflux,
                        ls='None',
                        marker=extra_dict["marker"],
                        label=extra_dict["label"],
                        color=color,
                        alpha=0.5,
                    )
                else: 
                    ax.errorbar(
                        plot_freq,
                        plot_flux,
                        yerr=plot_errflux,
                        ls='None',
                        marker=extra_dict["marker"],
                        color=color,
                        alpha=0.5,
                    )
                fit_flux = np.append(fit_flux, plot_flux)
                fit_err_flux = np.append(fit_err_flux, plot_errflux)
     
    legend = survey_legend

    if fit_params is not None:
        legend = True
        colourmap = cmr.neon 
        colour_indices = np.linspace(0.2,0.8,len(fit_params.keys()))
        i=0
        for fits in fit_params.keys():
            fit_dict = fit_params[fits]
            if fit_dict['best_params'] is None: 
                pass 
            # if fit_dict['name'] in plot_models:

            else: 
                fit_p0 = fit_dict["best_params"]


                ax.plot(
                    nu,
                    fit_dict["model"](nu, *fit_dict["best_params"]),
                    alpha=0.6,
                    # color=colourmap(colour_indices[i]),
                    color=fit_dict['colour'],
                    # label=fit_dict['name'],
                )
                # pass
            i+=1

    if legend is True:
        ax.legend()
    
    ax.loglog()
    ax.set(
        xlabel='Frequency (GHz)',
        ylabel='Integrated Flux (Jy)',
        title=f"{source_name.format('latex')}",
    )

    min_ind = np.where(fit_flux == np.nanmin(fit_flux))[0][0]
    max_ind = np.where(fit_flux == np.nanmax(fit_flux))[0][0]
    y_min = (fit_flux[min_ind]-fit_err_flux[min_ind])*0.7
    if y_min < 0: 
        print(f"The ymin is less than zero, going to set limit to 0.8*flux")
        y_min = np.abs(fit_flux[min_ind]*0.8) 
    y_max = (fit_flux[max_ind]+fit_err_flux[max_ind])*1.2

    ax.tick_params(
        axis="both", which="major", direction="in", length=majorticklength, width=tickwidth, pad=5
    )
    ax.tick_params(
        axis="both", which="minor", direction="in", length=minorticklength, width=tickwidth
    )

    # ax.set_xticks(x_ticks)
    # if fit_params is None and plot_extras is None: 
    ax.set_xlim([0.07, 6.])
    # ax.xaxis.set_minor_formatter(ScalarFormatter(useMathText=True))

    ax.set_ylim([y_min, y_max])

    formatter = FuncFormatter(lambda y, _: '{:0.16g}'.format(y))

    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.yaxis.set_minor_formatter(ScalarFormatter(useMathText=True))

    savenm = str(source_name)
    savenm = savenm.replace(" ","_")
    fig.tight_layout()
    try: 
        fig.savefig(f"{outpath}/{savenm}.png")
    except FileNotFoundError:
        print(f"Couldnt save the sed because there's no folder: {outpath}, will just save here...")
        fig.savefig(f"./{savenm}.png")

    plt.close(fig)
    return 



main_cat = Table.read(f"{data_dir}/GLEAMX_DRII_filtered_prime_sedfit_paper_ucd.fits").to_pandas()
main_coords = SkyCoord(main_cat["RAJ2000"], main_cat["DEJ2000"],frame="fk5",unit=u.deg)

# src_name = "GLEAM-X J233355.2-234340"
# src_name = "GLEAM-X J052517.9-460020"
src_name = "GLEAM-X J023558.0-442619"
sav_name = src_name.replace(" ", "_")
src_row = main_cat.loc[main_cat.Name == src_name.encode()]
freq, flux, err_flux = get_freq_flux_err(src_row)

try: 
    with open(f"{data_dir}/SED_fits/{sav_name}_fits.pkl", "rb") as fp:
        fit_params = pickle.load(fp)
except: 
    print("Couldnt find the file with fit params, refitting ")
    # fit_params = None
    freq = [i.to(u.GHz) for i in freq]
    fit_freq = np.array([i.value for i in freq])
    fit_flux = np.array([i.value for i in flux])
    fit_err_flux = np.array([i.value for i in err_flux])

    # low_freq = fit_freq[0:5]
    # low_flux = fit_flux[0:5]
    # low_err_flux = fit_err_flux[0:5]

    high_freq = fit_freq[18:]
    high_flux = fit_flux[18:]
    high_err_flux = fit_err_flux[18:]

    # best_pl_low, chi2_pl_low, rchi2_pl_low = read_prep_data.fit_pl(low_freq, low_flux, low_err_flux)
    # best_pl_hig, chi2_pl_hig, rchi2_pl_hig = read_prep_data.fit_pl(high_freq, high_flux, high_err_flux)
    best_cpl, chi2_cpl, rchi2_cpl = read_prep_data.fit_cpl(fit_freq, fit_flux, fit_err_flux)
    fit_params = {
    # "pl_low": {
    #     "name": "PL",
    #     "model": read_prep_data.power_law,
    #     "best_params": best_pl_low,
    #     "chi2": chi2_pl_low,
    #     "rchi2": rchi2_pl_low,
    #     "colour": "hotpink",
    # },
    # "pl_high" : {
    #     "name": "PL",
    #     "model": read_prep_data.power_law,
    #     "best_params": best_pl_hig,
    #     "chi2": chi2_pl_hig,
    #     "rchi2": rchi2_pl_hig,
    #     "colour": "hotpink" 
    # },
    "cpl_mwa": {
        "name": "CPL",
        "model": read_prep_data.curved_power_law,
        "best_params": best_cpl,
        "chi2": chi2_cpl,
        "rchi2": rchi2_cpl,
        "colour": "darkviolet"
    }
    }


    # curve_freq = freq[5:]
    # curve_flux = flux[5:]
    # curve_err_flux = err_flux[5:]

    curve_freq = freq
    curve_flux = flux
    curve_err_flux = err_flux

#     for nm in extra_fluxes.keys():
#         extra_dict = extra_fluxes[nm]
#         if extra_dict['label'] in ["Parkes", "VLASS19", "VLASS22"]:
#             pass
#         curve_flux = np.append(curve_flux, extra_dict['flux'].to(u.Jy))
#         curve_freq = np.append(curve_freq, extra_dict['freq'].to(u.GHz))
#         curve_err_flux = np.append(curve_err_flux, extra_dict['err_flux'].to(u.Jy))

    curve_freq = [i.to(u.GHz) for i in curve_freq]
    curve_freq = np.array([i.value for i in curve_freq])
    curve_flux = np.array([i.value for i in curve_flux])
    curve_err_flux = np.array([i.value for i in curve_err_flux])

#     best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, model=sed_models.singSSA)
#     curved_dict =  {
#             "name": "singSSA",
#             "model": sed_models.singSSA,
#             "best_params": best_p,
#             "chi2": chi2,
#             "rchi2": rchi2,
#             "colour": "mediumblue"
#         }
#     fit_params["ssa"] = curved_dict

#     best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, model=sed_models.singinhomobremss)
#     curved_dict =  {
#             "name": "singinhomobremss",
#             "model": sed_models.singinhomobremss,
#             "best_params": best_p,
#             "chi2": chi2,
#             "rchi2": rchi2,
#             "colour": "mediumblue"
#         }
#     fit_params["inffa"] = curved_dict

#     best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, model=sed_models.curve)
#     curved_dict =  {
#             "name": "curve",
#             "model": sed_models.curve,
#             "best_params": best_p,
#             "chi2": chi2,
#             "rchi2": rchi2,
#             "colour": "mediumblue"
#         }
#     fit_params["curve"] = curved_dict


#     best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, model=sed_models.doubSSA)
#     curved_dict =  {
#             "name": "doubSSA",
#             "model": sed_models.doubSSA,
#             "best_params": best_p,
#             "chi2": chi2,
#             "rchi2": rchi2,
#             "colour": "mediumblue"
#         }
#     fit_params["doubssa"] = curved_dict


    # best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, model=sed_models.singhomobremss)
    # curved_dict =  {
    #         "name": "singhomobremss",
    #         "model": sed_models.singhomobremss,
    #         "best_params": best_p,
    #         "chi2": chi2,
    #         "rchi2": rchi2,
    #         "colour": "mediumblue"
    #     }
    # fit_params["ffa"] = curved_dict

#     best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, model=sed_models.internalbremss)
#     curved_dict =  {
#             "name": "internalbremss",
#             "model": sed_models.internalbremss,
#             "best_params": best_p,
#             "chi2": chi2,
#             "rchi2": rchi2,
#             "colour": "mediumblue"
#         }
#     fit_params["intbremss"] = curved_dict

#     best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, model=sed_models.powlawbreak_nophys)
#     curved_dict =  {
#             "name": "powlawbreak_nophys",
#             "model": sed_models.powlawbreak_nophys,
#             "best_params": best_p,
#             "chi2": chi2,
#             "rchi2": rchi2,
#             "colour": "mediumblue"
#         }
#     fit_params["plbreak"] = curved_dict

#     best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, model=sed_models.doubhomobremss)
#     curved_dict =  {
#             "name": "doubhomobremss",
#             "model": sed_models.doubhomobremss,
#             "best_params": best_p,
#             "chi2": chi2,
#             "rchi2": rchi2,
#             "colour": "mediumblue"
#         }
#     fit_params["doubffa"] = curved_dict

    with open(f"{data_dir}/SED_fits/{sav_name}_fits.pkl", "wb") as fp:
        pickle.dump(fit_params, fp)


try: 
    with open(f"{data_dir}/SED_fits/{sav_name}_extra_fluxes.pkl", "rb") as fp:
        extra_fluxes = pickle.load(fp)
except:
    print("Couldnt load the extra fluxes... setting to None")
    extra_fluxes = None

curve_freq = freq
curve_flux = flux
curve_err_flux = err_flux
low_freq = freq
low_flux = flux
low_err_flux = err_flux

high_freq = freq[18:]
high_flux = flux[18:]
high_err_flux = err_flux[18:]
# gleamyr1 = extra_fluxes["gleamyr1"]
# gleamyr2 = extra_fluxes["gleamyr2"]
# gleam_freq = gleamyr1['freq']

# flux_yr1 = gleamyr1["flux"]
# err_fluxyr1 = gleamyr1['err_flux']
# flux_yr2 = gleamyr2['flux']
# err_fluxyr2 = gleamyr2['err_flux']

high_freq = [i.to(u.GHz) for i in high_freq]
high_freq = np.array([i.value for i in high_freq])
high_flux = np.array([i.value for i in high_flux])
high_err_flux = np.array([i.value for i in high_err_flux])

for nm in extra_fluxes.keys():
    if nm in ['atca']:
        extra_dict = extra_fluxes[nm]
        extra_flux = [i.to(u.Jy) for i in extra_dict['flux']]
        extra_freq = [i.to(u.GHz) for i in extra_dict['freq']]
        extra_errflux = [i.to(u.Jy) for i in extra_dict['err_flux']]
        extra_flux = np.array([i.value for i in extra_flux])
        extra_freq = np.array([i.value for i in extra_freq])
        extra_errflux = np.array([i.value for i in extra_errflux])
        high_flux = np.append(high_flux, extra_flux)
        high_freq = np.append(high_freq, extra_freq)
        high_err_flux = np.append(high_err_flux, extra_errflux)

# curve_freq = [i.to(u.GHz) for i in curve_freq]
# curve_freq = np.array([i.value for i in curve_freq])
# curve_flux = np.array([i.value for i in curve_flux])
# curve_err_flux = np.array([i.value for i in curve_err_flux])

# gleam_freq = [i.to(u.GHz) for i in gleam_freq]
# gleam_freq = np.array([i.value for i in gleam_freq])
# flux_yr1 = np.array([i.value for i in flux_yr1])
# err_fluxyr1 = np.array([i.value for i in err_fluxyr1])
# flux_yr2 = np.array([i.value for i in flux_yr2])
# err_fluxyr2 = np.array([i.value for i in err_fluxyr2])

# low_freq = [i.to(u.GHz) for i in low_freq]
# low_freq = np.array([i.value for i in low_freq])
# low_flux = np.array([i.value for i in low_flux])
# low_err_flux = np.array([i.value for i in low_err_flux])

fit_params = {}

# def kats_model(freq,S_norm,alpha,p,freq_peak, a, alpha_low, break_freq): 

#     return a*((freq/break_freq)**(alpha_low)) + (S_norm*(p+1)*((freq/freq_peak)**(2.1*(p+1)-alpha))*special.gammainc((p+1),((freq/freq_peak)**(-2.1)))*special.gamma(p+1))
# best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, p0=[ 0.17413116,  0.63740553, -0.34191943,  0.34448735, 0.01990462, -0.60918343, 0.1], model=kats_model)
# curved_dict =  {
#         "name": "inFFA+PL",
#         "model": kats_model,
#         "best_params": best_p,
#         "chi2": chi2,
#         "rchi2": rchi2,
#         "colour": "mediumblue"
#     }
# fit_params["kats_model"] = curved_dict


# best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, p0=[1.6, 2.0, 0.3], model=sed_models.internalbremss)
# curved_dict =  {
#         "name": "internalbremss",
#         "model": sed_models.internalbremss,
#         "best_params": best_p,
#         "chi2": chi2,
#         "rchi2": rchi2,
#         "colour": "mediumblue"
#     }
# fit_params["intbremss"] = curved_dict

best_pl_low, chi2_pl_low, rchi2_pl_low = sed_fitting.fit_model(high_freq, high_flux, high_err_flux, p0=[0.02, 1], model=sed_models.powlaw)
curved_dict =  {
        "name": "PL",
        "model": sed_models.powlaw,
        "best_params": best_pl_low,
        "chi2": chi2_pl_low,
        "rchi2": rchi2_pl_low,
        "colour": "hotpink",
    }
fit_params["pl_high"] = curved_dict

# best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, p0=[ 0.1,  0.3, 1 ,  0.2], model=sed_models.singinhomobremss)
# curved_dict =  {
#         "name": "singinhomobremss",
#         "model": sed_models.singinhomobremss,
#         "best_params": best_p,
#         "chi2": chi2,
#         "rchi2": rchi2,
#         "colour": "mediumblue"
#     }
# fit_params["inffa"] = curved_dict

# best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, p0=[0.2, 0.8, 0.2], model=sed_models.singhomobremss)
# curved_dict =  {
#         "name": "singhomobremss",
#         "model": sed_models.singhomobremss,
#         "best_params": best_p,
#         "chi2": chi2,
#         "rchi2": rchi2,
#         "colour": "mediumblue"
#     }
# fit_params["ffa"] = curved_dict

# best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, p0=[0.2,0.2,2.5, -2], model=sed_models.curve)
# curved_dict =  {
#         "name": "curve",
#         "model": sed_models.curve,
#         "best_params": best_p,
#         "chi2": chi2,
#         "rchi2": rchi2,
#         "colour": "hotpink"
#     }
# fit_params["curve"] = curved_dict

# best_p, chi2, rchi2 = sed_fitting.fit_model(gleam_freq, flux_yr1, err_fluxyr1, p0=[0.2,0.2,2.5, -2], model=sed_models.curve)
# curved_dict =  {
#         "name": "curve",
#         "model": sed_models.curve,
#         "best_params": best_p,
#         "chi2": chi2,
#         "rchi2": rchi2,
#         "colour": "mediumblue"
#     }
# fit_params["curve_yr1"] = curved_dict

# best_p, chi2, rchi2 = sed_fitting.fit_model(gleam_freq, flux_yr2, err_fluxyr2, p0=[0.2,0.2,2.5, -2], model=sed_models.curve)
# curved_dict =  {
#         "name": "curve",
#         "model": sed_models.curve,
#         "best_params": best_p,
#         "chi2": chi2,
#         "rchi2": rchi2,
#         "colour": "darkviolet"
#     }
# fit_params["curve_yr2"] = curved_dict

# # best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, p0=[0.1, 0.1, 2, 0.2, 5, 0.2], model=sed_models.doubSSA)
# # curved_dict =  {
# #         "name": "doubSSA",
# #         "model": sed_models.doubSSA,
# #         "best_params": best_p,
# #         "chi2": chi2,
# #         "rchi2": rchi2,
# #         "colour": "mediumblue"
# #     }
# # fit_params["doubssa"] = curved_dict

# best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, p0=[0.1, 0.6, -0.6, 0.3], model=sed_models.powlawbreak_nophys)
# curved_dict =  {
#         "name": "powlawbreak_nophys",
#         "model": sed_models.powlawbreak_nophys,
#         "best_params": best_p,
#         "chi2": chi2,
#         "rchi2": rchi2,
#         "colour": "mediumblue"
#     }
# fit_params["plbreak"] = curved_dict

# # best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, p0=[0.1, 0.2, 0.3], model=sed_models.duffcurve)
# # curved_dict =  {
# #         "name": "duffcurve",
# #         "model": sed_models.duffcurve,
# #         "best_params": best_p,
# #         "chi2": chi2,
# #         "rchi2": rchi2,
# #         "colour": "mediumblue"
# #     }
# # fit_params["duffcurve"] = curved_dict

# best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, p0=[0.2,0.6,-0.3, 0.2,2.0], model=sed_models.singinhomobremssbreak)
# curved_dict =  {
#         "name": "singinhomobremssbreak",
#         "model": sed_models.singinhomobremssbreak,
#         "best_params": best_p,
#         "chi2": chi2,
#         "rchi2": rchi2,
#         "colour": "mediumblue"
#     }
# fit_params["inffab"] = curved_dict
# [ 0.17413116,  0.63740553, -0.34191943,  0.34448735]
# S_norm,alpha,p,freq_peak
# best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, p0=[0.2,0.6,-0.34,0.344, 0.1], model=sed_models.singinhomobremssbreakexp)
# curved_dict =  {
#         "name": "singinhomobremssbreakexp",
#         "model": sed_models.singinhomobremssbreakexp,
#         "best_params": best_p,
#         "chi2": chi2,
#         "rchi2": rchi2,
#         "colour": "mediumblue"
#     }
# fit_params["inffab"] = curved_dict

# best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, p0=[0.2,0.6,-0.3, 0.2,2.0], model=sed_models.singinhomobremssbreakexp)
# curved_dict =  {
#         "name": "singinhomobremssbreakexp",
#         "model": sed_models.singinhomobremssbreakexp,
#         "best_params": best_p,
#         "chi2": chi2,
#         "rchi2": rchi2,
#         "colour": "mediumblue"
#     }
# fit_params["inffaeb"] = curved_dict

# # best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, model=sed_models.singinhomobremsscurve)
# # curved_dict =  {
# #         "name": "singinhomobremsscurve",
# #         "model": sed_models.singinhomobremsscurve,
# #         "best_params": best_p,
# #         "chi2": chi2,
# #         "rchi2": rchi2,
# #         "colour": "mediumblue"
# #     }
# # fit_params["inffaq"] = curved_dict

# # best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, model=sed_models.singhomobremsscurve)
# # curved_dict =  {
# #         "name": "singhomobremsscurve",
# #         "model": sed_models.singhomobremsscurve,
# #         "best_params": best_p,
# #         "chi2": chi2,
# #         "rchi2": rchi2,
# #         "colour": "mediumblue"
# #     }
# # fit_params["inffa"] = curved_dict



# # best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, model=sed_models.doubhomobremss)
# # curved_dict =  {
# #         "name": "doubhomobremss",
# #         "model": sed_models.doubhomobremss,
# #         "best_params": best_p,
# #         "chi2": chi2,
# #         "rchi2": rchi2,
# #         "colour": "mediumblue"
# #     }
# # fit_params["doubhomobremss"] = curved_dict


# best_p, chi2, rchi2 = sed_fitting.fit_model(curve_freq, curve_flux, curve_err_flux, model=sed_models.quad)
# curved_dict =  {
#         "name": "quad",
#         "model": sed_models.quad,
#         "best_params": best_p,
#         "chi2": chi2,
#         "rchi2": rchi2,
#         "colour": "mediumblue"
#     }
# fit_params["quad"] = curved_dict



plot_sed(freq, flux, err_flux, fit_params, src_name, outpath=save_dir, plot_extras = extra_fluxes, survey_legend=False)