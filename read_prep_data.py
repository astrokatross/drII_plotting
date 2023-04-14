#! /usr/bin/env python3

import astropy
from astropy.table import Table, Column
from astropy.io.votable import from_table, writeto
from astropy.coordinates import SkyCoord, search_around_sky
import numpy as np
import astropy.units as u
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
Script to organise the data from various catalogues and make a database style structure for easier handling per source in things like sed fitting, variability analysis etc 
"""

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(module)s:%(levelname)s:%(lineno)d %(message)s")
logger.setLevel(logging.INFO)

class source():
    """
    Custom class for each source that you can add the fluxes, frequencies epochs etc for it to be used in plotting 
    """

    def __init__(self, name: str, ra: int, dec: int):
        self.name = name
        self.ra = ra * u.deg
        self.dec = dec * u.deg


def read_catalogue(cat:str, source:source, colnames, search_radius = 1*u.degree, cat_ext = "fits"):

    try: 
        cats = Table.read(cat, format=cat_ext)
        catalogue = SkyCoord(cats.ref_ra, cats.ref_dec, frame="fk5", unit=u.deg)
    except FileNotFoundError:
        logger.debug("Can't find the catalogue you're trying to search!?")
        return np.nan 
    except: 
        logger.debug("Can't make catalogue?")
        return np.nan 

    return catalogue 


def search_cat(source: source, catalogue, search_radius=1*u.arcmin):

    target = SkyCoord(ra=source.ra,dec=source.dec)
    try: 
        distances = target.separation(catalogue, search_radius)
        idc = np.where(distances < search_radius)[0]
        if len(idc) > 1: 
            logger.warning("Found multiple matches! Which one do you want?")
            logger.warning(f"Indecies of matches: {idc}")
            raise IndexError
    except: 
        logger.debug(f"Couldn't find {source.name} in catalogue")
        raise Exception 

    # TODO/WHEREIM UP TO: Use the index from above to get all entries from catalogue for hte source and then create array based on the column names (also add that as variable) then you can return the array with values from the column names 

    return 
