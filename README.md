# sar_tools
### A python based [QGIS](https://qgis.org/en/site/index.html) plugin

This plugin generates derived SAR parameters (viz. vegetation indices, polarimetric decomposition parameters) from input polarimetric matrix (C3, T3, C2, T2). The input data needs to be in [PolSARpro](https://earth.esa.int/web/polsarpro/home)/[ENVI](https://www.l3harrisgeospatial.com/Software-Technology/ENVI) format (\*.bin and \*.hdr). It requires [numpy](https://numpy.org/), [matplotlib](https://matplotlib.org/) python libraries pre-installed.

### Avialble functionalities:
	
  **Indices**
  - Radar Vegetation Index (RVI)
  - Generalized volume Radar Vegetation Index (GRVI)
  - Polarimetric Radar Vegetation Index (PRVI)  
  - Dual-pol Radar Vegetation Index (DpRVI)

  **Polarimetric Decomposition**
  - Non-Model 3-Component decomposition for full-pol data (NM3CF).
  - Non-Model 3-Component decomposition for compact-pol data (NM3CC) 
	
### References
	
Ratha, D., Mandal, D., Kumar, V., McNairn, H., Bhattacharya, A. and Frery, A.C., 2019. A generalized volume scattering model-based vegetation index from polarimetric SAR data. IEEE Geoscience and Remote Sensing Letters, 16(11), pp.1791-1795.

Dey, S., Bhattacharya, A., Ratha, D., Mandal, D. and Frery, A.C., 2020. Target Characterization and Scattering Power Decomposition for Full and Compact Polarimetric SAR Data. IEEE Transactions on Geoscience and Remote Sensing.

Mandal, D., Kumar, V., Ratha, D., Dey, S., Bhattacharya, A., Lopez-Sanchez, J.M., McNairn, H. and Rao, Y.S., 2020. Dual polarimetric radar vegetation index for crop growth monitoring using sentinel-1 SAR data. Remote Sensing of Environment, 247, p.111954.

Chang, J.G., Shoshany, M. and Oh, Y., 2018. Polarimetric Radar Vegetation Index for Biomass Estimation in Desert Fringe Ecosystems. IEEE Transactions on Geoscience and Remote Sensing, 56(12), pp.7102-7108.
