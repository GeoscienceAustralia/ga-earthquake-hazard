# A Far-Field Ground-Motion Model for the North Australian Craton from Plate Margin Earthquakes

This archive provides data and codes used in the production of the manuscript "A far-field ground-motion model for the North Australian Craton from plate margin earthquakes" by Trevor I Allen, submitted to the Bulletin of the Seismological Society of America.  Files include:

- **event_table.csv:** List of earthquake parameters and number of recordings used in this study.
- **base_model_flatfile.csv:** File containing the ground-motion intensity measures for recordings used in the development of the base ground-motion model coefficients. The table below provides a description of the flatfile attributes.
- **base_amp_model_flatfile.csv:** File containing the ground-motion intensity measures for recordings used in the development of the amplification model coefficients. The table below provides a description of the flatfile attributes.

| ATTRIBUTE | DESCRIPTION |
| --------- | ----------- |
| ORIGIN_TIME | Origin time of earthquake |
| LON | Earthquake longitude (degrees) |
| LAT | Earthquake latitude (degrees) |
| DEPTH | Earthquake hypocentral depth (km) |
| MW | Moment magnitude |
| STA | Station code |
| NET | FDSN network code |
| INST_TY | Instrument type: BH* = low-sample-rate broadband; HH* = high-sample-rate broadband; SH* = low-sample-rate short period; HN* = high-sample-rate strong motion |
| VS30 | Time averaged shear-wave velocity in the upper 30 m of the site foundation (m/s) |
| RHYP | Hypocentral distance (km) |
| AZIM | Source-to-receiver azimuth (degrees) |
| PGV | Peak ground velocity in units of cm/s |
| PGA | Peak ground velocity in units of g |
| T*X* | Geometric mean of five-percent damped spectral acceleration for the horizontal components for a period of *X* seconds in units of g |

