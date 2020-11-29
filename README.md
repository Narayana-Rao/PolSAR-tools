# SAR Tools
### A python based [QGIS](https://qgis.org/en/site/index.html) plugin

This plugin generates derived SAR parameters (viz. vegetation indices, polarimetric decomposition parameters) from input polarimetric matrix (C3, T3, C2, T2). The input data needs to be in [PolSARpro](https://earth.esa.int/web/polsarpro/home)/[ENVI](https://www.l3harrisgeospatial.com/Software-Technology/ENVI) format (\*.bin and \*.hdr). It requires [numpy](https://numpy.org/), [matplotlib](https://matplotlib.org/) python libraries pre-installed.

### Available functionalities:
	
  **Indices**
  - Radar Vegetation Index (RVI) (Full-pol and dual-pol)
  - Generalized volume Radar Vegetation Index (GRVI)
  - Polarimetric Radar Vegetation Index (PRVI) (Full-pol and dual-pol) 
  - Dual-pol Radar Vegetation Index (DpRVI)
  - Degree of Polarization (DOP) (Full-pol, dual-pol, and compact-pol)
  - Compact-pol Radar Vegetation Index (CpRVI)
  
  **Polarimetric Decompositions**
  - Model free 3-Component decomposition for full-pol data (MF3CF).
  - Model free 3-Component decomposition for compact-pol data (MF3CC) 
  - Improved S-Omega decomposition for compact-pol data (iS-Omega)
	
### References

Chang, J.G., Shoshany, M. and Oh, Y., 2018. Polarimetric Radar Vegetation Index for Biomass Estimation in Desert Fringe Ecosystems. IEEE Transactions on Geoscience and Remote Sensing, 56(12), pp.7102-7108.

Ratha, D., Mandal, D., Kumar, V., McNairn, H., Bhattacharya, A. and Frery, A.C., 2019. A generalized volume scattering model-based vegetation index from polarimetric SAR data. IEEE Geoscience and Remote Sensing Letters, 16(11), pp.1791-1795.

Mandal, D., Kumar, V., Ratha, D., J. M. Lopez-Sanchez, A. Bhattacharya, H. McNairn, Y. S. Rao, and K. V. Ramana, 2020. Assessment of rice growth conditions in a semi-arid region of India using the Generalized Radar Vegetation Index derived from RADARSAT-2 polarimetric SAR data, Remote Sensing of Environment, 237: 111561.

Dey, S., Bhattacharya, A., Ratha, D., Mandal, D. and Frery, A.C., 2020. Target Characterization and Scattering Power Decomposition for Full and Compact Polarimetric SAR Data. IEEE Transactions on Geoscience and Remote Sensing.

Mandal, D., Kumar, V., Ratha, D., Dey, S., Bhattacharya, A., Lopez-Sanchez, J.M., McNairn, H. and Rao, Y.S., 2020. Dual polarimetric radar vegetation index for crop growth monitoring using sentinel-1 SAR data. Remote Sensing of Environment, 247, p.111954.

Mandal, D., Ratha, D., Bhattacharya, A., Kumar, V., McNairn, H., Rao, Y.S. and Frery, A.C., 2020. A Radar Vegetation Index for Crop Monitoring Using Compact Polarimetric SAR Data. IEEE Transactions on Geoscience and Remote Sensing, 58 (9), pp. 6321-6335.

V. Kumar, D. Mandal, A. Bhattacharya, and Y. S. Rao, 2020. Crop Characterization Using an Improved Scattering Power Decomposition Technique for Compact Polarimetric SAR Data. International Journal of Applied Earth Observations and Geoinformation, 88: 102052.
