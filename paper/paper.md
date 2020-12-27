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
  	orcid: 0000-0002-4979-0192
  	affliation: 1
  - name: Dipankar Mandal
  	orcid: 0000-0001-8407-7125
  	affliation: 1
  - name: Avik Bhattacharya
  	orcid: 0000-0001-6720-6108
  	affliation: 1
  - name: Y. S. Rao
  	affliation: 1

affiliations:
 - name: Microwave Remote Sensing Lab, Centre of Studies in Resources Engineering, Indian Institute of Technology Bombay, Mumbai-400076, India
   index: 1
date: 2 December 2020
bibliography: paper.bib
---

# Summary

The demand for processing tools increases with the increasing number of ***Synthetic Aperture Radar (SAR)*** satellite missions and datasets. However, to process SAR data, a minimal number of free tools are available ([PolSARpro](https://earth.esa.int/web/polsarpro/home), [SNAP](https://step.esa.int/main/toolboxes/snap/)), which consolidates all necessary pre-processing steps. Bearing this in mind, there is a need to develop specific tools for the remote sensing user community to derive polarimetric descriptors like the vegetation indices and decomposition parameters. Besides, to the best of our knowledge, there are no such free tools available on the GIS platform, which are quite necessary for SAR remote sensing. 

Hence we have developed a plugin for ```QGIS``` that supports data for all the three available polarimetric modes (i.e., full-, compact, and dual). The ```SAR tools``` plugin generates polarimetric descriptors (viz., vegetation indices, polarimetric decomposition parameters) from the 3x3 (C3/T3) or the 2x2 (C2/T2) covariance (coherency) matrices obtained from the ESA's [PolSARpro](https://earth.esa.int/web/polsarpro/home) software. The input data needs to be in PolSARpro format (```*.bin``` and ```*.hdr```). The plugin is coded in Python and is dependant on the Quantum GIS framework. It uses the following libraries (bundled with Quantum GIS): [numpy](https://numpy.org/), [gdal](https://gdal.org/) and [QGIS](https://qgis.org/en/site/index.html) core library.

# Background
The ***polarimetric decomposition*** techniques which are incorporated in this ```QGIS``` based plugin are model-free, i.e. to compute the decomposition power components no prior assumptions on the volume models is considered. The conventional model-based methods utilize a typical hierarchical process to enumerate power components uses various branching conditions, leading to several limitations. In this regard, these decomposition techniques utilizes some roll-invariant target characterization parameters to decompose the total power into even bounce, odd bounce and diffused power components. The powers obtained from the proposed technique are guaranteed to be non-negative, with the total power being conserved.

***Vegetation indices*** are often used as a proxy to plant growth. While appreciating the potential of vegetation indices derived from optical remote sensing sensors, regional to global products have been supported for operational uses. The Earth Observation (EO) community is relying upon the ```Synthetic Aperture Radar (SAR)``` imaging technology due to its all-weather imaging capability among its numerous advantages. The radar images are presently processed by several downstream users and are more frequently interpreted by non-radar specialists. This shift in paradigm offers the utility of radar-derived vegetation indices quintessential towards the goal of Analysis Ready Data (ARD) products. Recently, we proposed three vegetation indices namely GRVI (Generalized Radar Vegetation Index) [@ratha2019generalized], CpRVI (Compact-pol Radar Vegetation Index) [@mandal2020radar], and Dual-pol Radar Vegetation Index (DpRVI) [@mandal2020dual] for distinct acquisition modes. The vegetation indices have indicated an opportunity to directly estimate biophysical parameters from vegetation index images with fitted models. The retrieval of biophysical parameters from SAR observations is of vital importance for in-season monitoring of crop growth. 


# SAR tools Audience

```SAR tools``` is intended for students, researchers and polarimetry experts who would like to derive different SAR descriptors, utilizing the ```QGIS``` and ```python``` ecosystem of diverse tools. Especially for non-domain and application users the plugin interface provides an easy way to process the pre-procesed SAR polarimetric data. 

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
The authors would like to thank the developers of [QGIS Plugin Builder](https://github.com/g-sherman/Qgis-Plugin-Builder). Authors acknowledge the [GEO-AWS Earth Observation Cloud Credits Program](https://www.earthobservations.org/aws.php), which supported the computation, development, and testing of ```SARtools``` on AWS cloud platform through the project: 'AWS4AgriSAR-Crop inventory mapping from SAR data on cloud computing platform.'
	
# References
