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

def fit_model(freqs, flux, flux_err, model=sed_models.powlaw):

    try: 
        popt, pcov = curve_fit(model, freqs, flux, sigma=flux_err)
        perr = np.sqrt(np.diag(pcov))

    except: 
        logger.warning("Couldn't fit the sed??")    
        popt = np.nan
        perr = np.nan

    return popt, perr 

def calc_confidence_curve(popt, perr, nstd, freqs=np.linspace(70,1400,1000), model=sed_models.powlaw):

    popt_ul = popt + nstd*perr
    popt_ll = popt - nstd*perr 

    fit_ul = model(freqs, *popt_ul)
    fit_ll = model(freqs, *popt_ll)


    return fit_ll, fit_ul

