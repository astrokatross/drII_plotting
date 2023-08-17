from astropy.io import fits 
import numpy as np 
from argparse import ArgumentParser

if __name__ == "__main__":
    parser = ArgumentParser(
        description="Script to just divide two images "
    )
    parser.add_argument(
        '--image',
        type=str,
        default="XG.fits",
        help="input image to divide from"
    )
    parser.add_argument(
        '--rms',
        type=str,
        default="XG_rms.fits",
        help="input rms to divide"
    )
    parser.add_argument(
        '--output',
        type=str,
        default="XG_sigma.fits",
        help="output image"
    )
    args = parser.parse_args()

    image = fits.open(args.image)
    rms = fits.open(args.rms)

    image[0].data = image[0].data/rms[0].data 

    fits.writeto(args.output, image[0].data, header= image[0].header,overwrite=True)
