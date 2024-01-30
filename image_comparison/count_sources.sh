#!/bin/bash

catalogue_gx=../data/GLEAMX_DRII_filtered_prime_sedfit.fits
catalogue_glm=/data/useful_surveys/GLEAM.fits
#157.5-(3.5/2)
#155.75000000000000000000
#157.5+(3.5/2)
#159.25000000000000000000
#-27.5+(3.5/2)
#-25.75000000000000000000
#-27.5-(3.5/2)
#-29.25000000000000000000

# TODO: synchronise this with the plotting code
ra1=23.25
ra2=26.75
dec2=-31.75
dec1=-28.25

stilts tpipe \
        in=$catalogue_glm \
        cmd="select ((RAJ2000<=$ra2)&&(RAJ2000>=$ra1)&&(DEJ2000<=$dec1)&&(DEJ2000>=$dec2)&&(int_flux_wide/local_rms_wide>=5))" \
        cmd='keepcols "RAJ2000 DEJ2000 local_rms_wide"' \
        out='./number_sources_glm.csv'

stilts tpipe \
        in=$catalogue_gx \
        cmd="select ((ref_ra<=$ra2)&&(ref_ra>=$ra1)&&(ref_dec<=$dec1)&&(ref_dec>=$dec2)&&(int_flux/local_rms>=5))" \
        cmd='keepcols "ref_ra ref_dec local_rms"' \
        out='./number_sources.csv'

n=`wc -l number_sources_glm.csv | awk '{print $1}'`
((n-=1))
echo "$n sources in GLEAM region"

n=`wc -l number_sources.csv | awk '{print $1}'`
((n-=1))
echo "$n sources in GLEAM-X region"
