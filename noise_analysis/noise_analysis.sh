#!/bin/bash

usage()
{
echo "script to make all the images for noise analysis" 1>&2;
return;
}
freq=
colour=

while getopts ":f:c:" opt; do
    case $opt in
    f)
        freq="${OPTARG}"
        ;;
    c)
        colour="${OPTARG}" 
        ;;
  esac
done    

if [[ -z $freq ]] || [[ -z ${colour} ]]
then
    echo "You didn't define things..."
    return
fi 

echo $freq
# First extract small pieces so we don't run out of memory
ra="02:30:00"
dec="-40:00:00"
imgpix_white=1250
imgpix_red=520
bkgpix=70
datadir="/data/gleam_x/drII_plotting/data"
basedir="/data/gleam_x/drII_plotting/noise_analysis/"
savedir="/data/gleam_x/drII_plotting/plots"
# freq="_103-111MHz"
# colour="green"


if [[ ! -e "${datadir}/XG_white.fits" ]]
then
    echo "Nope no file, going to create"
    cd $datadir
    getfits -sv -o "${datadir}/XG_white.fits ${datadir}/GLEAMX_DRII_170-231MHz.fits" $ra $dec $imgpix_white $imgpix_white
    cd $basedir
else
    echo "have file all g"
fi

# TODO: This doesn't work for all frequencies, but mostly doing in prep anyway 
if [[ ! -e "${datadir}/XG_${colour}${freq}.fits" ]]
then
    cd $datadir
    getfits -sv -o "${datadir}/XG_${colour}${freq}.fits ${datadir}/GLEAMX_DRII${freq}.fits" $ra $dec $imgpix_red $imgpix_red
    cd $basedir
fi

exts="bkg rms"
for ext in $exts
do
    if [[ ! -e "${datadir}/XG_white_${ext}.fits" ]]
    then
# NB: this is a modified version of BANE that is forced to do 10 sigma clip loops
       BANE --cores 1 "${datadir}/XG_white.fits"
    fi
    if [[ ! -e "${datadir}/XG_${colour}${freq}_${ext}.fits" ]]
    then
# NB: this is a modified version of BANE that is forced to do 10 sigma clip loops
       BANE --cores 1 "${datadir}/XG_${colour}${freq}.fits"
    fi   
done

# Make a small local catalogue
if [[ ! -e "${datadir}/XG_white_comp.fits" ]]
then
    aegean \
    --seedclip=4 \
    --maxsummits=5 \
    --progress \
    --psf=${datadir}/XG_white_psfmap.fits \
    --autoload \
    --table="${datadir}/XG_white.fits" \
    "${datadir}/XG_white.fits"
fi

# Make background and rms and comp for ${colour} with no prior stuff 
# Make a small local catalogue
if [[ ! -e "${datadir}/XG_${colour}${freq}_comp.fits" ]]
then
    aegean \
    --seedclip=4 \
    --maxsummits=5 \
    --progress \
    --psf="${datadir}/XG_${colour}${freq}_psfmap.fits" \
    --autoload \
    --table="${datadir}/XG_${colour}${freq}.fits" \
    "${datadir}/XG_${colour}${freq}.fits"
fi

# Make background and rms and comp for ${colour} with no prior stuff 
# Make a small local catalogue
if [[ ! -e "${datadir}/XG_${colour}${freq}_prior_comp.fits" ]]
then
    aegean \
    --seedclip=4 \
    --maxsummits=5 \
    --regroup-eps=5 \
    --progress \
    --psf="${datadir}/XG_${colour}${freq}_psfmap.fits" \
    --priorized=1 \
    --input="${datadir}/XG_white_comp.fits" \
    --autoload \
    --table="${datadir}/XG_${colour}${freq}_prior.fits" \
    "${datadir}/XG_${colour}${freq}.fits"
fi

# Now to subtract the priorized catalogue and then calculate the nosie and rms again 
AeRes -c ${datadir}/XG_${colour}${freq}_prior_comp.fits -f ${datadir}/XG_${colour}${freq}.fits -r ${datadir}/XG_${colour}${freq}_priorsub_resid.fits
BANE --cores 1 ${datadir}/XG_${colour}${freq}_priorsub_resid.fits 

# Rerun source finding both prior and no prior with new bkg and rms 
if [[ ! -e ${datadir}/XG_${colour}${freq}_priorsub_prior_comp.fits ]]
then
    aegean \
    --seedclip=4 \
    --maxsummits=5 \
    --regroup-eps=5 \
    --progress \
    --psf=${datadir}/XG_${colour}${freq}_psfmap.fits \
    --background=${datadir}/XG_${colour}${freq}_priorsub_resid_bkg.fits \
    --noise=${datadir}/XG_${colour}${freq}_priorsub_resid_rms.fits \
    --priorized=1 \
    --input=${datadir}/XG_white_comp.fits \
    --table=${datadir}/XG_${colour}${freq}_priorsub_prior.fits \
    ${datadir}/XG_${colour}${freq}.fits
fi

# Rerun source finding both prior and no prior with new bkg and rms 
if [[ ! -e ${datadir}/XG_${colour}${freq}_priorsub_comp.fits ]]
then
    aegean \
    --seedclip=4 \
    --maxsummits=5 \
    --progress \
    --psf=${datadir}/XG_${colour}${freq}_psfmap.fits \
    --background=${datadir}/XG_${colour}${freq}_priorsub_resid_bkg.fits \
    --noise=${datadir}/XG_${colour}${freq}_priorsub_resid_rms.fits \
    --table=${datadir}/XG_${colour}${freq}_priorsub.fits \
    ${datadir}/XG_${colour}${freq}.fits
fi



# Subtract catalogue from image
AeRes -c ${datadir}/XG_${colour}${freq}_priorsub_comp.fits -f ${datadir}/XG_${colour}${freq}.fits -r ${datadir}/XG_${colour}${freq}_residual.fits
AeRes -c ${datadir}/XG_white_comp.fits -f ${datadir}/XG_white.fits -r ${datadir}/XG_white_residual.fits

# # Subtract background from image
python noise_analysis/subtract.py --image ${datadir}/XG_${colour}${freq}.fits --bkg ${datadir}/XG_${colour}${freq}_priorsub_resid_bkg.fits --output ${datadir}/XG_${colour}${freq}_priorsub_bkgsub.fits
python noise_analysis/subtract.py --image ${datadir}/XG_${colour}${freq}.fits --bkg ${datadir}/XG_${colour}${freq}_bkg.fits --output ${datadir}/XG_${colour}${freq}_bkgsub.fits
python noise_analysis/subtract.py --image ${datadir}/XG_white.fits --bkg ${datadir}/XG_white_bkg.fits --output ${datadir}/XG_white_bkgsub.fits


# Subtract catalogue from background-subtracted image
AeRes -c ${datadir}/XG_${colour}${freq}_priorsub_prior_comp.fits -f ${datadir}/XG_${colour}${freq}_priorsub_bkgsub.fits -r ${datadir}/XG_${colour}${freq}_priorsub_bkgsub_resid.fits
AeRes -c ${datadir}/XG_${colour}${freq}_prior_comp.fits -f ${datadir}/XG_${colour}${freq}_bkgsub.fits -r ${datadir}/XG_${colour}${freq}_bkgsub_resid.fits
AeRes -c ${datadir}/XG_white_comp.fits -f ${datadir}/XG_white_bkgsub.fits -r ${datadir}/XG_white_bkgsub_resid.fits


# Mask instead of subtract
AeRes --mask --sigma=1 -c ${datadir}/XG_${colour}${freq}_priorsub_prior_comp.fits -f ${datadir}/XG_${colour}${freq}_priorsub_bkgsub.fits -r ${datadir}/XG_${colour}${freq}_priorsub_bkgsub_masked.fits
AeRes --mask --sigma=1 -c ${datadir}/XG_${colour}${freq}_prior_comp.fits -f ${datadir}/XG_${colour}${freq}_bkgsub.fits -r ${datadir}/XG_${colour}${freq}_bkgsub_masked.fits
AeRes --mask --sigma=1 -c ${datadir}/XG_white_comp.fits -f ${datadir}/XG_white_bkgsub.fits -r ${datadir}/XG_white_bkgsub_masked.fits


# Create S/N map
python noise_analysis/div.py --image ${datadir}/XG_${colour}${freq}.fits --rms ${datadir}/XG_${colour}${freq}_priorsub_resid_rms.fits --output ${datadir}/XG_${colour}${freq}_priorsub_sigma.fits
python noise_analysis/div.py --image ${datadir}/XG_${colour}${freq}.fits --rms ${datadir}/XG_${colour}${freq}_rms.fits --output ${datadir}/XG_${colour}${freq}_sigma.fits
python noise_analysis/div.py --image ${datadir}/XG_white.fits --rms ${datadir}/XG_white_rms.fits --output ${datadir}/XG_white_sigma.fits


# Background-subtracted S/N map
python noise_analysis/div.py --image ${datadir}/XG_${colour}${freq}_priorsub_bkgsub.fits --rms ${datadir}/XG_${colour}${freq}_priorsub_resid_rms.fits --output ${datadir}/XG_${colour}${freq}_priorsub_bkgsub_sigma.fits
python noise_analysis/div.py --image ${datadir}/XG_${colour}${freq}_bkgsub.fits --rms ${datadir}/XG_${colour}${freq}_rms.fits --output ${datadir}/XG_${colour}${freq}_bkgsub_sigma.fits
python noise_analysis/div.py --image ${datadir}/XG_white_bkgsub.fits --rms ${datadir}/XG_white_rms.fits --output ${datadir}/XG_white_bkgsub_sigma.fits


# Masked, background-subtracted S/N map
python noise_analysis/div.py --image ${datadir}/XG_${colour}${freq}_priorsub_bkgsub_masked.fits --rms ${datadir}/XG_${colour}${freq}_priorsub_resid_rms.fits --output ${datadir}/XG_${colour}${freq}_priorsub_bkgsub_masked_sigma.fits
python noise_analysis/div.py --image ${datadir}/XG_${colour}${freq}_bkgsub_masked.fits --rms ${datadir}/XG_${colour}${freq}_rms.fits --output ${datadir}/XG_${colour}${freq}_bkgsub_masked_sigma.fits
python noise_analysis/div.py --image ${datadir}/XG_white_bkgsub_masked.fits --rms ${datadir}/XG_white_rms.fits --output ${datadir}/XG_white_bkgsub_masked_sigma.fits


# Residuals, background-subtracted S/N map
python noise_analysis/div.py --image ${datadir}/XG_${colour}${freq}_priorsub_bkgsub_resid.fits --rms ${datadir}/XG_${colour}${freq}_priorsub_resid_rms.fits --output ${datadir}/XG_${colour}${freq}_priorsub_bkgsub_resid_sigma.fits
python noise_analysis/div.py --image ${datadir}/XG_${colour}${freq}_bkgsub_resid.fits --rms ${datadir}/XG_${colour}${freq}_rms.fits --output ${datadir}/XG_${colour}${freq}_bkgsub_resid_sigma.fits
python noise_analysis/div.py --image ${datadir}/XG_white_bkgsub_resid.fits --rms ${datadir}/XG_white_rms.fits --output ${datadir}/XG_white_bkgsub_resid_sigma.fits