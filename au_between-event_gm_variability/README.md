# Systematic Estimation of Earthquake Source Parameters for Continental Australia: Magnitude Conversions and Ground-Motion Variability 

This archive provides data used in the production of the manuscript Systematic Estimation of Earthquake Source Parameters for Continental Australia: Magnitude Conversions and Ground-Motion Variability by Trevor I Allen, submitted for publication in the Bulletin of the Seismological Society of America. This manuscript is a companion paper to [Systematic Estimation of Earthquake Source Parameters for Continental Australia: Attenuation and Stress Drop](https://pubs.geoscienceworld.org/ssa/bssa/article/doi/10.1785/0120250054/660495/Systematic-Estimation-of-Earthquake-Source). The file for this companion paper is:

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
| ML2800 | Local magnitude assuming classical Wood-Anderson instrument calibration as per [Richter (1935)](https://pubs.geoscienceworld.org/ssa/bssa/article/25/1/1/115102/An-instrumental-earthquake-magnitude-scale) |
| ML2080 | Local magnitude assuming revised Wood-Anderson instrument calibration as per [Uhrhammer and Collins (1990)](https://pubs.geoscienceworld.org/ssa/bssa/article/80/3/702/119366/Synthesis-of-Wood-Anderson-seismograms-from) |
| BRUNE_MAG | Moment magnitude from Brune spectral fitting |
| BRUNE_MAG_STD | Standard deviation of the moment magnitude from Brune spectral fitting |
| STRESS_DROP | Stress drop from Brune spectral fitting (MPa) |
| LOG_STRESS_DROP_STD | log<sub>10</sub> standard deviation of the stress drop from Brune spectral fitting |
| CORN_FREQ | Corner frequency from Brune spectral fitting (Hz) |
| CORN_FREQ_STD | Standard deviation of the corner frequency from Brune spectral fitting (Hz) |
| NRECS | Number of recordings used for Brune spectral fitting |
| FMIN | Minimum frequency used to fit Brune spectral model (Hz) |
| FMAX | Maximum frequency used to fit Brune spectral model (Hz) |
| DOMAIN | Neotectonic domain based on [Clark et al. 2011](https://ecat.ga.gov.au/geonetwork/srv/api/records/a05f7892-f6f5-7506-e044-00144fdd4fa6). PRPC = Precambrian cratons and reactivated Proterozoic crust; SOPC = Sprigg Orogen and Phanerozoic accretionary crust; ECPM = extended continental crust and passive margins |
| CLUSTER | Earthquake cluster in which the event is assigned based on <i>k</i>-means unsupervised machine learning |

Orginal magnitude types are defined as:

| ABBREVIATION | DESCRIPTION |
| ------------ | ----------- |
| mb | Body-wave magnitude |
| ML | Local magnitude |
| Mw | Moment magnitude as determined by the [National Earthquake Alerts Centre](https://earthquakes.ga.gov.au/) |
| Mwa | Moment magnitude as determined by [Allen (2012)](https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/74133) |
| Mwg | Moment magnitude as determined by [Ghasemi et al. (2016)](https://aees.org.au/wp-content/uploads/2018/06/341-Ghasemi-et-al.pdf) |
| Mwp | P-phase moment magnitude as determined by the [National Earthquake Alerts Centre](https://earthquakes.ga.gov.au/) |


