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


Simulation for Austria and Brazil is done in R:

the files "RSkript_aut_3.R" and "RSkript_bra_3.R" contain the main part of the simulation and analysis, including download
the files "ERA5_data.R" and "ERA5_dat_bra.R" contain functions for the download and file conversion of ERA5 wind speed data, the latter is used for Brazil, because the method applied for Austria did not work for Brazil due to larger file sizes
in "functions_aut_2.R" and "functions_bra_2.R" the functions used for simulation and analysis are contained


The simulation for South Africa is done in Python:

the file "WP_simulation_ZAF.py" contains the full run of the simulation
the file "era5_download_zaf.py" contains the download function and is called in "WP_simulation_ZAF"
the wind park information is in "./data/windparks_southafrica.csv"
