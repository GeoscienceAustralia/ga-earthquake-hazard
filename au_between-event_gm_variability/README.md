# Systematic Estimation of Earthquake Source Parameters for Continental Australia: Magnitude Conversions and Ground-Motion Variability 

This archive provides data used in the production of the manuscript Systematic Estimation of Earthquake Source Parameters for Continental Australia: Magnitude Conversions and Ground-Motion Variability by Trevor I Allen, published in the Bulletin of the Seismological Society of America. This manuscript is a companion paper to [Systematic Estimation of Earthquake Source Parameters for Continental Australia: Attenuation and Stress Drop](https://pubs.geoscienceworld.org/ssa/bssa/article/doi/10.1785/0120250054/660495/Systematic-Estimation-of-Earthquake-Source). The file for this companion paper is:

- [clustered_brune_source_params.csv](https://github.com/GeoscienceAustralia/ga-earthquake-hazard/tree/master/au_between-event_gm_variability/clustered_brune_source_params.csv): Summary file of earthquake source parameters calculated in this study. The table below provides a description of table attributes.

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
| CLUSTER | Earthquake cluster in which the event is assigned based on <i>k</i>-means unsupervised machine learning |
