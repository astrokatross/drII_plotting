#!/bin/bash

# First extract small pieces so we don't run out of memory
ra="02:30:00"
dec="-40:00:00"
imgpix=1250
bkgpix=70

# if [[ ! -e "XG.fits" ]]
# then
getfits -o XG.fits GLEAMX_DRII-170-231MHz.fits $ra $dec $imgpix $imgpix
# fi

exts="bkg rms"
for ext in $exts
do
    if [[ ! -e XG_${ext}.fits ]]
    then
# NB: this is a modified version of BANE that is forced to do 10 sigma clip loops
       BANE --cores 1 XG.fits
    fi
done

# Make a small local catalogue
# TODO only define these numbers once, above, and set up argparse for this python script
if [[ ! -e XG_comp.fits ]]
then
# Make my own catalogue
    aegean \
    --seedclip=4 \
    --regroup-eps 5 \
    --maxsummits=5 \
    --priorized=1 \
    --progress \
    --input=XG_priorcat.fits \
    --autoload \
    --table=XG.fits \
    XG.fits
fi
    
    # aegean \
    # --regroup-eps 5 \
    # --cores ${GXNCPUS} \
    # --background "${imname}_ddmod_bkg.fits" \
    # --noise "${imname}_ddmod_rms.fits" \
    # --psf "${imname}_projpsf_psf.fits" \
    # --table "${imname}_ddmod_prior.fits" \
    # --priorized 1 \
    # --input "${prior_cat}" \
    # --progress \
    # "${imname}_ddmod.fits"

# Subtract catalogue from image
AeRes -c XG_comp.fits -f XG.fits -r XG_residual.fits
# Subtract background from image
python subtract.py --image XG.fits --bkg XG_bkg.fits --output XG_bkgsubtracted.fits
# # Subtract catalogue from background-subtracted image
AeRes -c XG_comp.fits --sigma=1 -f XG_bkgsubtracted.fits -r XG_residual_bkgsubtracted.fits
# # Mask instead of subtract
AeRes --mask --sigma=1 -c XG_comp.fits -f XG_bkgsubtracted.fits -r XG_masked_bkgsubtracted.fits
# # Create S/N map
python div.py --image XG.fits --rms XG_rms.fits --output XG_sigma.fits
# # Background-subtracted S/N map
python div.py --image XG_bkgsubtracted.fits --rms XG_rms.fits --output XG_bkgsubtracted_sigma.fits
# # Masked, background-subtracted S/N map
python div.py --image XG_masked_bkgsubtracted.fits --rms XG_rms.fits --output XG_masked_bkgsubtracted_sigma.fits
# # Residuals, background-subtracted S/N map
python div.py --image XG_residual_bkgsubtracted.fits --rms XG_rms.fits --output XG_residual_bkgsubtracted_sigma.fits


python div.py --image XG_residual.fits --rms XG_rms.fits --output XG_residual_sigma.fits
