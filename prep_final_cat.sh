#!/usr/bin/env bash

# ------------------------------------------------------------------
# --------------------  CATALOGUE QUALITY CUTS ---------------------
# ------------------------------------------------------------------

# python reliability/negatives.py \
# 	data/GLEAMX_DRII_joined_rescaled_comp.fits \
# 	reliability/GLEAMX_DRII_170-231MHz_psf_comp_negative_comp.fits \
# 	-r ref_ra \
# 	-d ref_dec \
# 	-o data/GLEAMX_DRII_rescaled_filtered.fits

# python catalogue_seds/catalogue_seds.py \
# 	IDR1_subset_filtered.fits \
# 	-c 2 -o IDR1_subset_filtered_SEDs.fits

# python final_cat_prep/catalogue_names.py \
#     data/GLEAMX_DRII_filtered_prime_sedfit.fits

# python add_ucd/add_ucd.py \
#     data/GLEAMX_DRII_filtered_prime_sedfit_paper.fits \
# 	--apply



# ------------------------------------------------------------------
# ------------- NOISE ANALYSIS: PROOF OF CONFUSION -----------------
# ------------------------------------------------------------------

# source noise_analysis/noise_analysis.sh -f "_072-103MHz" -c "red"
# source noise_analysis/noise_analysis.sh -f "_072-080MHz" -c "red"
# source noise_analysis/noise_analysis.sh -f "_080-088MHz" -c "red"
# ./noise_analysis/noise_analysis.sh -f "_088-095MHz" -c "red"
# ./noise_analysis/noise_analysis.sh -f "_095-103MHz" -c "red"

# ------------------------------------------------------------------
# -------------    ASTROMETRY CALCS AND REFINING   -----------------
# ------------------------------------------------------------------

# just for talks 
# python astrometry/astrometry_plot.py data/GLEAMX_DRII_170-231MHz_projpsf_comp.fits /data/gleam_x/GLEAM-X-pipeline/models/NVSS_SUMSS_psfcal.fits  --ra ra --dec dec --plot-output plots/astrometry_for_talks.png --min-snr 50 --flux-col int_flux --flux-err-col err_int_flux

# python3 astrometry/extract_prime_region.py data/GLEAMX_DRII_filtered.fits data/GLEAMX_DRII_filtered_prime.fits

# python3 astrometry/astrometry_plot.py data/GLEAMX_DRII_filtered_prime.fits /data/gleam_x/GLEAM-X-pipeline/models/NVSS_SUMSS_psfcal.fits  --ra ref_ra --dec ref_dec --plot-output plots/astrometry_DRII.pdf --min-snr 50 --flux-col int_flux --flux-err-col err_int_flux --latex

# ------------------------------------------------------------------
# -------------------    NICE PAPER PLOTS   ------------------------
# ------------------------------------------------------------------
# python sky_coverage/make_skyplot.py

# python noise_analysis/plotting_noise.py --imagenm "XG_red_072-103MHz.fits" --compare
# python noise_analysis/plotting_noise.py --imagenm "XG_white.fits"

# python GLEAM_S_alpha_comparison/plot_S_alpha_comparison.py

# python alpha_distribution/plot_alpha.py 

python image_comparison/mosaic_plot.py 