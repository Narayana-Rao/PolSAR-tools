---
title: 'SAR tools: A QGIS plugin for generating SAR descriptors'
tags:
  - SAR
  - QGIS
  - Vegetation indices
  - Polarimetric decompositions
  
authors:
  - name: Narayanarao Bhogapurapu
    orcid: 0000-0002-6496-7283
    affiliation: 1
  - name: Subhadip Dey
  	orcid: 
  	affliation: 1
  - name: Dipankar Mandal
  	orcid:
  	affliation: 1
  - name: Avik Bhattacharya
  	orcid:
  	affliation: 1
  - name: Rao Y. S.
  	orcid:
  	affliation: 1

affiliations:
 - name: Microwave Remote Sensing Lab, Centre of Studies in Resources Engineering, Indian Institute of Technology Bombay, Mumbai-400076, India
   index: 1
date: 2 December 2020
bibliography: paper.bib
---

# Summary
With increasing number of Synthetic Aperture Radar satellite missions and datasets, the demand for processing tools is also increasing. However, to process Synthetic Aperture Radar data very limited free tools are available ([PolSARpro](https://earth.esa.int/web/polsarpro/home), [SNAP](https://step.esa.int/main/toolboxes/snap/)) with major concentration on pre-processing. In application user point of view there is a neccesity for tools to derive polarimetric descriptors like vegetation indices and decomposition parameters. In addition there are no free tools in a GIS platform, which is very much essential as remote sensing and GIS are highly inter-dependent. So we have developed a plugin which supports data of all the three avaialble polarimetric modes (full, compact and dual).
```SAR tools``` plugin generates polarimetric descriptors (viz. vegetation indices, polarimetric decomposition parameters) from C3/T3/C2/T2 matrices obtained from PolSARpro The input data needs to be in PolSARpro format (```*.bin``` and ```*.hdr```). 	
The plug-in is coded in Python and is dependant of the Quantum GIS framework. More specifically, it makes use of following libraries (bundled with Quantum GIS): [numpy](https://numpy.org/), [gdal](https://gdal.org/) and [QGIS](https://qgis.org/en/site/index.html) core library.

# SAR tools Audience

**SAR tools** is intended for students, researchers and polarimetry experts who would like to derive different SAR descriptors, utilizing the ```QGIS``` and ```python``` ecosystem of diverse tools. Especially for non-domain and application users the plugin interface provides an easy way to process the pre-procesed SAR polarimetric data. 

# SAR tools Functionality

The key functionality of **SAR tools** is organized into three modules:
  - **Full-pol**: 
    - Radar Vegetation Index (RVI) [@Kim_2009]
    - Generalized volume Radar Vegetation Index (GRVI) [@ratha2019generalized]
    - Polarimetric Radar Vegetation Index (PRVI) [@chang2018polarimetric] 
    - Model free 3-Component decomposition for full-pol data (MF3CF) [@dey2020target]
    - Degree of Polarization (DOP) [@barakat1977degree]
  - **Compact-pol**:
    - Model free 3-Component decomposition for compact-pol data (MF3CC) [@dey2020target]
    - Improved S-Omega decomposition for compact-pol data (iS-Omega) [@kumar2020crop]
    - Compact-pol Radar Vegetation Index (CpRVI) [@mandal2020radar]
    - Degree of Polarization (DOP) 
  - **Dual-pol**:
    - Radar Vegetation Index (RVI) [@trudel2012using]
    - Dual-pol Radar Vegetation Index (DpRVI) [@mandal2020dual], 
    - Polarimetric Radar Vegetation Index (PRVI) 
    - Degree of Polarization (DOP) [@barakat1977degree]

# Acknowledgements
The authors would like to thank the developers of [QGIS Plugin Builder](https://github.com/g-sherman/Qgis-Plugin-Builder). 
	
# References
