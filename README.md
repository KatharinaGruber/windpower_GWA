# windpower_GWA
scripts for simulating wind power generation from ERA5 and MERRA-2 reanalysis data with Global Wind Atlas version 2 and 3 bias correction and validation with wind power time series.


The results of Austria (outdated) and Brazil were presented at the IEWT 2019 (Internationale Energiewirtschaftstagung 13.-15.Feb 2019, Vienna)
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
Version: https://github.com/KatharinaGruber/windpower_GWA/tree/icem

Simulation and analysis of wind power generation in Brazil, New Zealand, USA and South Africa is performed in Python.
The scripts
BRA/run_bra_simulation.sh
NZ/run_NZ_simulation.sh
USA/run_usa_simulation.sh
ZAF/run_ZAF_simulation.sh
run the simulation for each country.

Each simulation consists of
- reanalysis data download (ERA5 and MERRA-2)
- wind power simulation for each of the reanalysis with and without GWA
- for the GWA simulation between GWA2 and GWA3 can be chosen
- computation of statistical parameters for simulation

GWA needs to be downloaded manually, as well as wind park information data and validation data
More information can be found in the related publication (Preprint) https://arxiv.org/abs/2012.05648



------
Info for Austria (outdated):

Note that there was an issue: Two different resolutions are used in the GWA as it was updated in the meantime: 1km for Austria and Brazil and 250m for South Africa and USA
More recent simulations are using GWA2 and GWA3, both with a resolution of 250m

Simulation for Austria was done in R:

the files "RSkript_aut_3.R" and "RSkript_bra_3.R" contain the main part of the simulation and analysis, including download
the files "ERA5_data.R" and "ERA5_dat_bra.R" contain functions for the download and file conversion of ERA5 wind speed data, the latter is used for Brazil, because the method applied for Austria did not work for Brazil due to larger file sizes
in "functions_aut_2.R" and "functions_bra_2.R" the functions used for simulation and analysis are contained
in "RSkript_bra_country_monthly.R" and "RSkript_bra_states.R" further analysis for ICEM results is conducted

------
We gratefully acknowledge support from the European Research Council (“reFUEL” ERC2017-STG 758149).