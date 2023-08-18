#!/bin/bash

# First extract small pieces so we don't run out of memory
ra="02:30:00"
dec="-40:00:00"
imgpix_white=1250
imgpix_red=520
bkgpix=70
freq="_95-103MHz"

if [[ ! -e "XG_white.fits" ]]
then
    getfits -sv -o XG_white.fits GLEAMX_DRII-170-231MHz.fits $ra $dec $imgpix_white $imgpix_white
fi

# TODO: This doesn't work for all frequencies, but mostly doing in prep anyway 
if [[ ! -e XG_red${freq}.fits ]]
then
    getfits -sv -o XG_red.fits GLEAMX_DRII-072-103MHz.fits $ra $dec $imgpix_red $imgpix_red
fi

exts="bkg rms"
for ext in $exts
do
    if [[ ! -e XG_white_${ext}.fits ]]
    then
# NB: this is a modified version of BANE that is forced to do 10 sigma clip loops
       BANE --cores 1 XG_white.fits
    fi
    if [[ ! -e XG_red${freq}_${ext}.fits ]]
    then
# NB: this is a modified version of BANE that is forced to do 10 sigma clip loops
       BANE --cores 1 XG_red${freq}.fits
    fi   
done

# Make a small local catalogue
if [[ ! -e XG_white_comp.fits ]]
then
    aegean \
    --seedclip=4 \
    --maxsummits=5 \
    --progress \
    --psf=XG_white_psfmap.fits \
    --autoload \
    --table=XG_white.fits \
    XG_white.fits
fi

# Make background and rms and comp for red with no prior stuff 
# Make a small local catalogue
if [[ ! -e XG_red${freq}_comp.fits ]]
then
    aegean \
    --seedclip=4 \
    --maxsummits=5 \
    --progress \
    --psf=XG_red${freq}_psfmap.fits \
    --autoload \
    --table=XG_red${freq}.fits \
    XG_red${freq}.fits
fi

# Make background and rms and comp for red with no prior stuff 
# Make a small local catalogue
if [[ ! -e XG_red${freq}_prior_comp.fits ]]
then
    aegean \
    --seedclip=4 \
    --maxsummits=5 \
    --regroup-eps=5 \
    --progress \
    --psf=XG_red${freq}_psfmap.fits \
    --priorized=1 \
    --input=XG_white_comp.fits \
    --autoload \
    --table=XG_red${freq}_prior.fits \
    XG_red${freq}.fits
fi

# Now to subtract the priorized catalogue and then calculate the nosie and rms again 
AeRes -c XG_red${freq}_prior_comp.fits -f XG_red${freq}.fits -r XG_red${freq}_priorsub_resid.fits
BANE --cores 1 XG_red${freq}_priorsub_resid.fits 

# Rerun source finding both prior and no prior with new bkg and rms 
if [[ ! -e XG_red${freq}_priorsub_prior_comp.fits ]]
then
    aegean \
    --seedclip=4 \
    --maxsummits=5 \
    --regroup-eps=5 \
    --progress \
    --psf=XG_red${freq}_psfmap.fits \
    --background=XG_red${freq}_priorsub_resid_bkg.fits \
    --noise=XG_red${freq}_priorsub_resid_rms.fits \
    --priorized=1 \
    --input=XG_white_comp.fits \
    --table=XG_red${freq}_priorsub_prior.fits \
    XG_red${freq}.fits
fi

# Rerun source finding both prior and no prior with new bkg and rms 
if [[ ! -e XG_red${freq}_priorsub_comp.fits ]]
then
    aegean \
    --seedclip=4 \
    --maxsummits=5 \
    --progress \
    --psf=XG_red${freq}_psfmap.fits \
    --background=XG_red${freq}_priorsub_resid_bkg.fits \
    --noise=XG_red${freq}_priorsub_resid_rms.fits \
    --table=XG_red${freq}_priorsub.fits \
    XG_red${freq}.fits
fi



# Subtract catalogue from image
AeRes -c XG_red${freq}_priorsub_comp.fits -f XG_red${freq}.fits -r XG_red${freq}_residual.fits
AeRes -c XG_white_comp.fits -f XG_white.fits -r XG_white_residual.fits

# # Subtract background from image
python subtract.py --image XG_red${freq}.fits --bkg XG_red${freq}_priorsub_resid_bkg.fits --output XG_red${freq}_priorsub_bkgsub.fits
python subtract.py --image XG_red${freq}.fits --bkg XG_red${freq}_bkg.fits --output XG_red${freq}_bkgsub.fits
python subtract.py --image XG_white.fits --bkg XG_white_bkg.fits --output XG_white_bkgsub.fits


# Subtract catalogue from background-subtracted image
AeRes -c XG_red${freq}_priorsub_prior_comp.fits -f XG_red${freq}_priorsub_bkgsub.fits -r XG_red${freq}_priorsub_bkgsub_resid.fits
AeRes -c XG_red${freq}_prior_comp.fits -f XG_red${freq}_bkgsub.fits -r XG_red${freq}_bkgsub_resid.fits
AeRes -c XG_white_comp.fits -f XG_white_bkgsub.fits -r XG_white_bkgsub_resid.fits


# Mask instead of subtract
AeRes --mask --sigma=1 -c XG_red${freq}_priorsub_prior_comp.fits -f XG_red${freq}_priorsub_bkgsub.fits -r XG_red${freq}_priorsub_bkgsub_masked.fits
AeRes --mask --sigma=1 -c XG_red${freq}_prior_comp.fits -f XG_red${freq}_bkgsub.fits -r XG_red${freq}_bkgsub_masked.fits
AeRes --mask --sigma=1 -c XG_white_comp.fits -f XG_white_bkgsub.fits -r XG_white_bkgsub_masked.fits


# Create S/N map
python div.py --image XG_red${freq}.fits --rms XG_red${freq}_priorsub_resid_rms.fits --output XG_red${freq}_priorsub_sigma.fits
python div.py --image XG_red${freq}.fits --rms XG_red${freq}_rms.fits --output XG_red${freq}_sigma.fits
python div.py --image XG_white.fits --rms XG_white_rms.fits --output XG_white_sigma.fits


# Background-subtracted S/N map
python div.py --image XG_red${freq}_priorsub_bkgsub.fits --rms XG_red${freq}_priorsub_resid_rms.fits --output XG_red${freq}_priorsub_bkgsub_sigma.fits
python div.py --image XG_red${freq}_bkgsub.fits --rms XG_red${freq}_rms.fits --output XG_red${freq}_bkgsub_sigma.fits
python div.py --image XG_white_bkgsub.fits --rms XG_white_rms.fits --output XG_white_bkgsub_sigma.fits


# Masked, background-subtracted S/N map
python div.py --image XG_red${freq}_priorsub_bkgsub_masked.fits --rms XG_red${freq}_priorsub_resid_rms.fits --output XG_red${freq}_priorsub_bkgsub_masked_sigma.fits
python div.py --image XG_red${freq}_bkgsub_masked.fits --rms XG_red${freq}_rms.fits --output XG_red${freq}_bkgsub_masked_sigma.fits
python div.py --image XG_white_bkgsub_masked.fits --rms XG_white_rms.fits --output XG_white_bkgsub_masked_sigma.fits


# Residuals, background-subtracted S/N map
python div.py --image XG_red${freq}_priorsub_bkgsub_resid.fits --rms XG_red${freq}_priorsub_resid_rms.fits --output XG_red${freq}_priorsub_bkgsub_resid_sigma.fits
python div.py --image XG_red${freq}_bkgsub_resid.fits --rms XG_red${freq}_rms.fits --output XG_red${freq}_bkgsub_resid_sigma.fits
python div.py --image XG_white_bkgsub_resid.fits --rms XG_white_rms.fits --output XG_white_bkgsub_resid_sigma.fits