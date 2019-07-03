# windpower_GWA
scripts for simulating wind power generation from era5 data with Global Wind Atlas bias correction and validation with wind power time series.


The results of Austria and Brazil were presented at the IEWT 2019 (Internationale Energiewirtschaftstagung 13.-15.Feb 2019, Vienna)
https://iewt2019.eeg.tuwien.ac.at/programme_text
Full Paper: https://iewt2019.eeg.tuwien.ac.at/download/contribution/fullpaper/146/146_fullpaper_20190130_090240.pdf
Version: https://github.com/KatharinaGruber/windpower_GWA/tree/iewt


These results complemented with results from South Africa were presented at EGU 2019 (European Geosciences Union General Assembly 7. - 12. Apr 2019, Vienna).
Link to EGU abstract: https://meetingorganizer.copernicus.org/EGU2019/EGU2019-4752.pdf
Link to presentation: doi.org/10.13140/RG.2.2.31691.85283
Version: https://github.com/KatharinaGruber/windpower_GWA/tree/egu2019

These results complemented with results from the USA were presented at ICEM 2019 (6th International Conference Energy & Meteorology 24. - 27. Jun 2019, Copenhagen).
Link to ICEM abstract: http://icem2019-abstract-submission.p.wemc.currinda.com/days/2019-06-27/abstract/655
Link to presentation: doi.org/10.13140/RG.2.2.14639.79520
Note that there are two issues:
1. Two different resolutions are used in the GWA as it was updated in the meantime: 1km for Austria and Brazil and 250m for South Africa and USA
2. For the simulation in the USA the disph was not taken account of for MERRA-2

Simulation for Austria and Brazil is done in R:

the files "RSkript_aut_3.R" and "RSkript_bra_3.R" contain the main part of the simulation and analysis, including download
the files "ERA5_data.R" and "ERA5_dat_bra.R" contain functions for the download and file conversion of ERA5 wind speed data, the latter is used for Brazil, because the method applied for Austria did not work for Brazil due to larger file sizes
in "functions_aut_2.R" and "functions_bra_2.R" the functions used for simulation and analysis are contained
in "RSkript_bra_country_monthly.R" and "RSkript_bra_states.R" further analysis for ICEM results is conducted

The simulation for South Africa is done in Python:
data for South Africa are downloaded with "download_merra_ZAF.R" and "era5_download_zaf.py"
the files " WP_simulation_ZAF_ERA5.ipynb" and "WP_simulation_ZAF_MERRA2.ipynb" contain the simulation and analysis of wind power generation for South Africa
the wind park information is in "./data/windparks_southafrica.csv"

The simulation for the USA is done in Python:
Reanalysis data for the USA are downloaded with "download_merra2_USA.R" and "download_era5_USA.py"
In "WP_simulation_USA_ERA5.ipynb" and "WP_simulation_USA_MERRA2.ipynb" wind power in the USA is simulated
The analysis is in "WP_simulation_USA_Analysis.ipynb"
"utils" contains functions for the simulation and analysis of wind power generation in the USA

"logging_confic.py" and "config.py" are needed for the download of ERA5 data