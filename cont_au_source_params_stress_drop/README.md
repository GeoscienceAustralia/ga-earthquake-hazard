# Systematic Estimation of Earthquake Source Parameters for Continental Australia: Attenuation and Stress Drop 

This archive provides data used in the production of the manuscript [Systematic Estimation of Earthquake Source Parameters for Continental Australia: Attenuation and Stress Drop](https://pubs.geoscienceworld.org/ssa/bssa/article/doi/10.1785/0120250054/660495/Systematic-Estimation-of-Earthquake-Source) by Trevor I Allen, published in the Bulletin of the Seismological Society of America.  Files include:

- [brune_source_parameters.csv](https://github.com/GeoscienceAustralia/ga-earthquake-hazard/tree/master/cont_au_source_params_stress_drop/brune_source_parameters.csv): Summary file of earthquake source parameters calculated in this study. The table below provides a description of table attributes.
- [atten_coeffs.csv](https://github.com/GeoscienceAustralia/ga-earthquake-hazard/tree/master/cont_au_source_params_stress_drop/atten_coeffs.csv): Attenuation coefficients used in the evaluation of earthquake source parameters, as applied in Equations 8-11 of the manuscript. The coefficients include the aleatory variability model, such that $`\sigma_T = \sqrt{\sigma_{be}^2 +\sigma_{we}^2}`$, where:
	- $`\sigma_T`$ is the total variability
	- $`\sigma_{be}`$ is the between-event variability
	- $`\sigma_{we}`$ is the within-event variability 
	
- [au_stress_drop_supplemental_material.pdf[(https://github.com/GeoscienceAustralia/ga-earthquake-hazard/blob/master/cont_au_source_params_stress_drop/au_stress_drop_supplemental_material.pdf): additional supplemental material including plots of between- and within-event residuals
- [brune_fit](https://github.com/GeoscienceAustralia/ga-earthquake-hazard/tree/master/cont_au_source_params_stress_drop/brune_fit): Folder containing individual [Brune (1970)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/jb075i026p04997) spectral fits to all earthquakes for which source parameters could be resolved in this study
- [within-event_residuals](https://github.com/GeoscienceAustralia/ga-earthquake-hazard/tree/master/cont_au_source_params_stress_drop/within-event_residuals): Folder containing within-event residuals. Residuals are plotted against hypocentral distance for a frequency of 1.0 Hz. Individual stations are labelled.

| ATTRIBUTE | DESCRIPTION |
| --------- | ----------- |
| EVENT | Origin time of earthquake |
| GAID | Geoscience Australia earthquake ID |
| LON | Earthquake longitude (degrees) |
| LAT | Earthquake latitude (degrees) |
| DEP | Earthquake hypocentral depth (km) |
| OMAG | Original magnitude as per the [Geoscience Australia catalogue](https://earthquakes.ga.gov.au/) |
| OMAG_TYPE | Original magnitude type |
| BRUNE_MAG | Moment magnitude from Brune spectral fitting |
| BRUNE_MAG_STD | Standard deviation of the moment magnitude from Brune spectral fitting |
| STRESS_DROP | Stress drop from Brune spectral fitting (MPa) |
| LOG_STRESS_DROP_STD | log<sub>10</sub> standard deviation of the stress drop from Brune spectral fitting |
| CORN_FREQ | Corner frequency from Brune spectral fitting (Hz) |
| CORN_FREQ_STD | Standard deviation of the corner frequency from Brune spectral fitting (Hz) |
| NRECS | Number of recordings used for Brune spectral fitting |
| FMIN | Minimum frequency used to fit Brune spectral model (Hz) |
| FMAX | Maximum frequency used to fit Brune spectral model (Hz) |
