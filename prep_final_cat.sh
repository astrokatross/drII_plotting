#!/usr/bin/env bash

# ------------------------------------------------------------------
# --------------------  CATALOGUE QUALITY CUTS ---------------------
# ------------------------------------------------------------------

# python reliability/negatives.py \
# 	IDR1_subset.fits \
# 	reliability/XG_170-231MHz_psf_comp_negative_dafix_comp.fits \
# 	-r ref_ra \
# 	-d ref_dec \
# 	-o IDR1_subset_filtered.fits

# python catalogue_seds/catalogue_seds.py \
# 	IDR1_subset_filtered.fits \
# 	-c 2 -o IDR1_subset_filtered_SEDs.fits

# python catalogue_names/catalogue_names.py \
#     IDR1_subset_filtered_SEDs.fits

# python add_ucd/add_ucd.py \
#     IDR1_subset_filtered_SEDs_paper.fits \
# 	--apply



# ------------------------------------------------------------------
# ------------- NOISE ANALYSIS: PROOF OF CONFUSION -----------------
# ------------------------------------------------------------------

source noise_analysis/noise_analysis.sh -f "_072-103MHz" -c "red"
# ./noise_analysis/noise_analysis.sh -f "_072-080MHz" -c "red"
# ./noise_analysis/noise_analysis.sh -f "_080-088MHz" -c "red"
# ./noise_analysis/noise_analysis.sh -f "_088-095MHz" -c "red"
# ./noise_analysis/noise_analysis.sh -f "_095-103MHz" -c "red"

# ------------------------------------------------------------------
# -------------    ASTROMETRY CALCS AND REFINING   -----------------
# ------------------------------------------------------------------

# just for talks 
python astrometry/astrometry_plot.py data/GLEAMX_DRII_170-231MHz_projpsf_comp.fits /data/gleam_x/GLEAM-X-pipeline/models/NVSS_SUMSS_psfcal.fits  --ra ra --dec dec --plot-output plots/astrometry_for_talks.png --min-snr 50 --flux-col int_flux --flux-err-col err_int_flux

# python3 astrometry/extract_prime_region.py data/GLEAMX_DRII_combined_v2_sedfit.fits data/GLEAMX_DRII_v2_sedfit_prime.fits

# python3 astrometry/astrometry_plot.py data/GLEAMX_DRII_v2_sedfit_prime.fits /data/gleam_x/GLEAM-X-pipeline/models/NVSS_SUMSS_psfcal.fits  --ra ref_ra --dec ref_dec --plot-output plots/astrometry_DRII.pdf --min-snr 50 --flux-col int_flux --flux-err-col err_int_flux --latex


# ------------------------------------------------------------------
# -------------------    NICE PAPER PLOTS   ------------------------
# ------------------------------------------------------------------
python sky_coverage/make_skyplot.py