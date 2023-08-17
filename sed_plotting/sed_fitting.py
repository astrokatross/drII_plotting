#! /usr/bin/env python3

import astropy
from astropy.table import Table, Column
from astropy.io.votable import from_table, writeto
import numpy as np
import glob
import argparse
import sys
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import fits
import os
import time
import sed_models
import logging 
import read_prep_data


__author__ = "Kat Ross"
__date__ = "2023-04-14"

"""
Script to fit all sources in the GLEAM-X DRII 
"""

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(module)s:%(levelname)s:%(lineno)d %(message)s")
logger.setLevel(logging.INFO)

def read_flux():

    return 

def fit_model(freqs, flux, flux_err, model=read_prep_data.power_law):
    if model == read_prep_data.power_law:
# def fit_pl(freq, flux, fluxerr):
        p0 = (np.median(flux), -0.8)

        try:
            fit_res = curve_fit(
                model,
                freqs,
                flux,
                p0=p0,
                sigma=flux_err,
                absolute_sigma=True
            )
        except RuntimeError:
            return None

        best_p, covar = fit_res
        err_p = np.sqrt(np.diag(covar))
        dof = len(freqs) - 2
        chi2 = np.sum(
            ((flux - model(freqs, *best_p)) / flux_err)**2
            )
        rchi2 = chi2 / dof

        return dict(
            norm=best_p[0], 
            alpha=best_p[1],
            norm_err=err_p[0],
            alpha_err = err_p[1],
            chi2=chi2,
            rchi2=rchi2,
            dof=dof
            )
    elif model == read_prep_data.curved_power_law:
        p0 = (np.median(flux), -0.8, 0)

        try:
            fit_res = curve_fit(
                read_prep_data.curved_power_law,
                freqs,
                flux,
                p0=p0,
                sigma=flux_err,
                absolute_sigma=True
            )
        except RuntimeError:
            return None

        best_p, covar = fit_res
        err_p = np.sqrt(np.diag(covar))
        dof = len(freqs) - 3
        chi2 = np.sum(
            ((flux - read_prep_data.curved_power_law(freqs, *best_p)) / flux_err)**2
            )
        rchi2 = chi2 / dof
        return dict(
            norm=best_p[0], 
            alpha=best_p[1],
            q=best_p[2],
            norm_err=err_p[0],
            alpha_err = err_p[1],
            q_err=err_p[2],
            chi2=chi2,
            rchi2=rchi2,
            dof=dof
            )
    else: 
        logger.warning("Dunno what you want to fit cause I was lazy and only set up too.... plz give more details for model fitting")
        return np.nan 

def calc_confidence_curve(popt, perr, nstd, freqs=np.linspace(70,1400,1000), model=sed_models.powlaw):

    popt_ul = popt + nstd*perr
    popt_ll = popt - nstd*perr 

    fit_ul = model(freqs, *popt_ul)
    fit_ll = model(freqs, *popt_ll)


    return fit_ll, fit_ul

