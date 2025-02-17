# Systematic Estimation of Earthquake Source Parameters for Continental Australia: Attenuation and Stress Drop 

This archive provides data used in the production of the manuscript "Systematic Estimation of Earthquake Source Parameters for Continental Australia: Attenuation and Stress Drop" by Trevor I Allen, in review for the Bulletin of the Seismological Society of America.  Files include:

- [brune_source_parameters.csv](https://github.com/GeoscienceAustralia/ga-earthquake-hazard/tree/master/cont_au_source_params_stress_drop/brune_source_parameters.csv): Summary file of earthquake source parameters calulated in this study. The table below provides a description of table attributes.
- [brune_fit](https://github.com/GeoscienceAustralia/ga-earthquake-hazard/tree/master/cont_au_source_params_stress_drop/brune_fit): Folder containing individual [Brune (1970)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/jb075i026p04997) spectral fits to all earthquakes for which source parameters could be resolved in this study
- [within-event_residuals](https://github.com/GeoscienceAustralia/ga-earthquake-hazard/tree/master/cont_au_source_params_stress_drop/within-event_residuals): Folder containing within-event residuals. Residuals are plotted against hypocentral distance for a frequency of 1.0 Hz. Individual stations are labelled.

| ATTRIBUTE | DESCRIPTION |
| --------- | ----------- |
| EVENT | Origin time of earthquake |
| GAID | Geoscience Australia earthquake ID |
| LON | Earthquake longitude (degrees) |
| LAT | Earthquake latitude (degrees) |
| DEP | Earthquake hypocentral depth (km) |
| OMAG | Origial magnitude as per [Geoscience Australia catalogue](https://earthquakes.ga.gov.au/) |
| OMAG_TYPE | Origial magnitude type |
| MB | Body-wave magnitude |
| BRUNE_MAG | Moment magnitude from Brune spectral fitting |
| STRESS_DROP | Stress drop from Brune spectral fitting (MPa) |
| CORN_FREQ | Corner frequency from Brune spectral fitting (Hz) |
| NRECS | Number of recordings used for Brune spectral fitting |
| FMIN | Minimum frequency used to fit Brune spectral model (Hz) |
| FMAX | Maximum frequency used to fit Brune spectral model (Hz) |
