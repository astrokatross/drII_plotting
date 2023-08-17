from astropy.io import fits 
import numpy as np 
from argparse import ArgumentParser

if __name__ == "__main__":
    parser = ArgumentParser(
        description="Script to assess the quality of images for obsids and return list of obsids that pass quality assurance to be included in moasics. Note: currently only works on obsids given per channel, not per night. "
    )
    parser.add_argument(
        '--image',
        type=str,
        default="XG.fits",
        help="input image to subtract from"
    )
    parser.add_argument(
        '--bkg',
        type=str,
        default="XG_bkg.fits",
        help="input bkg to subtract"
    )
    parser.add_argument(
        '--output',
        type=str,
        default="XG_bkgsubtracted.fits",
        help="output image"
    )
    args = parser.parse_args()

    image = fits.open(args.image)
    bkg = fits.open(args.bkg)

    image[0].data = image[0].data - bkg[0].data 

    fits.writeto(args.output, image[0].data, header= image[0].header,overwrite=True)
