export PYTHONPATH=
python3 astrometry_plot.py ../data/GLEAMX_DRII_170-231MHz_projpsf_comp.fits /data/gleam_x/GLEAM-X-pipeline/models/NVSS_SUMSS_psfcal.fits  --ra ra --dec dec --plot-output /data/gleam_x/drII_plotting/plots/astrometry_for_talks.png --min-snr 50 --flux-col int_flux --flux-err-col err_int_flux
# python3 extract_prime_region.py /data/gleam_x/processing/GLEAMX-DRII/GLEAMX_DRII_combined_v2_sedfit.fits /data/gleam_x/processing/GLEAMX-DRII/GLEAMX_DRII_v2_sedfit_prime.fits
# python3 astrometry_plot.py /data/gleam_x/processing/GLEAMX-DRII/GLEAMX_DRII_v2_sedfit_prime.fits /data/gleam_x/GLEAM-X-pipeline/models/NVSS_SUMSS_psfcal.fits  --ra ref_ra --dec ref_dec --plot-output astrometry_DRII.pdf --min-snr 50 --flux-col int_flux --flux-err-col err_int_flux --latex
