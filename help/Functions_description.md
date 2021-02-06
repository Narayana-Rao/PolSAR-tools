## Functions description

**Full-pol functions**
----------------------

 * Generalized volume based Radar Vegetation Index (GRVI) 
    <center>

    ![grvi](https://latex.codecogs.com/svg.latex?\Large&space;\text{GRVI}=\left(1-\text{GD}_{\text{GV}}\right)\Big(\frac{p}{q}\Big)^{2\,\text{GD}_{\text{GV}}},\quad0\le\text{GRVI}\le1)
    
    </center> 

where, GD<sub>GV</sub> is the geodesic distance between Kennaugh matrices (K) of the observed and the generalized volume scattering model, p,q are minimum and maximum value of distances between K matrices of the observed and elementary targets respectively. A detailed explanation of GRVI is available in [[2]](#2).


````python
    input : input_T3_folder, window_size
    output: GRVI.bin
````

 * Model Free 3-Component decomposition for Full-pol data (MF3CF) 
    
    <center>

    ![mfp](https://latex.codecogs.com/svg.latex?\Large&space;m_{\text{FP}}=\sqrt{1-\frac{27|\mathbf{T3}|}{\big(\mathrm{Trace}(\mathbf{T3})\big)^3}})
    
    ![mfp](https://latex.codecogs.com/svg.latex?\Large&space;\tan\theta_{\text{FP}}=\frac{m_{\text{FP}}{\text{Span}}\left(T_{11}-T_{22}-T_{33}\right)}{T_{11}\left(T_{22}+T_{33}\right)+m_{\text{FP}}^{2}{\text{Span}}^{2}})
    
    </center> 


    <center>

    ![tfp](https://latex.codecogs.com/svg.latex?\Large&space;\noindent\\\P_{d}^{\text{FP}}=\frac{m_{\text{FP}}{\text{Span}}}{2}{\left(1-\sin2\theta_{\text{FP}}\right)}\\\P_{v}^{\text{FP}}={\text{Span}}\left(1-m_{\text{FP}}\right)\\\P_{s}^{\text{FP}}=\frac{m_{\text{FP}}{\text{Span}}}{2}\left(1+\sin2\theta_{\text{FP}}\right))
    
    </center>

    <center>

    ![spanfp](https://latex.codecogs.com/svg.latex?\Large&space;\text{Span}=T_{11}+T_{22}+T_{33})
    
    </center>

````python
    input : input_T3_folder, window_size
    output: Ps.bin,Pd.bin,Pv.bin,Theta_FP.bin
````

 * Radar Vegetation Index (RVI) 
<center>

![grvi](https://latex.codecogs.com/svg.latex?\Large&space;\text{RVI}_{fp}=\frac{8\sigma^\circ_{\text{HV}}}{\sigma^\circ_{\text{HH}}+\sigma^\circ_{\text{VV}}+2\sigma^\circ_{\text{HV}}})
    
</center> 


````python
    input : input_c2_folder, window_size
    output: RVI.bin
````
 * Polarimetric Radar Vegetation Index (PRVI) 
<center>

![grvi](https://latex.codecogs.com/svg.latex?\Large&space;\text{PRVI}_{fp}=(1-\text{DOP}_{fp})\sigma^\circ_{\text{XY}})
    
</center> 


````python
    input : input_c2_folder, window_size
    output: PRVI.bin
````

 * Degree of Polarization (DOP) 
<center>

![grvi](https://latex.codecogs.com/svg.latex?\Large&space;\text{DOP}_{fp}=\sqrt{1-\frac{27\times\text{det([T3])}}{\text{(Trace[T3])}^3}})
    
</center> 

````python
    input : input_c2_folder, window_size
    output: dop_fp.bin
````




**Compact-pol functions**
-------------------------


 * Model Free 3-Component decomposition for Full-pol data (MF3CF) 
    
    <center>

    ![mcp](https://latex.codecogs.com/svg.latex?\Large&space;m_{\text{CP}}=\sqrt{1-\frac{4|\mathbf{C2}|}{\big(\mathrm{Trace}(\mathbf{C2})\big)^2}};\tan\theta_{\text{CP}}=\frac{m_{\text{CP}}{S_0}\left(\text{OC}-\text{SC}\right)}{\text{OC}\times\text{SC}+m_{\text{CP}}^{2}{S_0}^{2}})
    
    </center> 


    <center>

    ![cppowers](https://latex.codecogs.com/svg.latex?\Large&space;P_{d}^{\text{CP}}=\frac{m_{\text{FP}}{S_0}}{2}{\left(1-\sin2\theta_{\text{CP}}\right)};\\\P_{v}^{\text{CP}}={S_0}\left(1-m_{\text{CP}}\right);\\\P_{s}^{\text{CP}}=\frac{m_{\text{CP}}{S_0}}{2}\left(1+\sin2\theta_{\text{CP}}\right))
    
    </center>
    <center>

    ![sparm](https://latex.codecogs.com/svg.latex?\Large&space;S_0=\text{C11+C22};S_1=\text{C11-C22};\\\S_2=\text{C12+C21};S_3=\pm\text{j(C12-C21)})
    
    </center>

    <center>

    ![cparm](https://latex.codecogs.com/svg.latex?\Large&space;\text{SC}=\frac{S_0-S_3}{2};\text{OC}=\frac{S_0+S_3}{2};)
    
    </center>


````python
    input : input_C2_folder, window_size
    output: Ps.bin,Pd.bin,Pv.bin,Theta_FP.bin
````








**Dual-pol**
------------

 * Dual-pol Radar Vegetation Index (DpRVI) 

    <center>

    ![grvi](https://latex.codecogs.com/svg.latex?\Large&space;\text{DpRVI}=1-\Big(\frac{\lambda_1}{\lambda_1+\lambda_2}\Big)\sqrt{1-\frac{4\times\text{det([C2])}}{\text{(Trace[C2])}^2}})
    
    </center> 

where, C2 is co-variance matrix,  and  &lambda;<sub>1</sub> and &lambda;<sub>2</sub> are the eigen values of C2 matrix in descending order.

````python
    input : input_C2_folder, window_size
    output: DpRVI.bin
````

 * Radar Vegetation Index (RVI) 
<center>

![grvi](https://latex.codecogs.com/svg.latex?\Large&space;\text{RVI}_{dp}=\frac{4\sigma^\circ_{\text{XY}}}{\sigma^\circ_{\text{XX}}+\sigma^\circ_{\text{XY}}})
    
</center> 


````python
    input : input_c2_folder, window_size
    output: RVI_dp.bin
````

 * Polarimetric Radar Vegetation Index (PRVI) 
<center>

![grvi](https://latex.codecogs.com/svg.latex?\Large&space;\text{PRVI}_{dp}=(1-\text{DOP}_{dp})\sigma^\circ_{\text{XY}})
    
</center> 


````python
    input : input_c2_folder, window_size
    output: PRVI_dp.bin
````

 * Degree of Polarization (DOP) 
<center>

![grvi](https://latex.codecogs.com/svg.latex?\Large&space;\text{DOP}_{dp}=\sqrt{1-\frac{4\times\text{det([C2])}}{\text{(Trace[C2])}^2}})
    
</center> 

````python
    input : input_c2_folder, window_size
    output: dop_dp.bin
````


## References
-------------
<a id="1">[1]</a> 
Chang, J.G., Shoshany, M. and Oh, Y., 2018. Polarimetric Radar Vegetation Index for Biomass Estimation in Desert Fringe Ecosystems. IEEE Transactions on Geoscience and Remote Sensing, 56(12), pp.7102-7108.

<a id="2">[2]</a> 
Ratha, D., Mandal, D., Kumar, V., McNairn, H., Bhattacharya, A. and Frery, A.C., 2019. A generalized volume scattering model-based vegetation index from polarimetric SAR data. IEEE Geoscience and Remote Sensing Letters, 16(11), pp.1791-1795.

<a id="3">[3]</a> 
Mandal, D., Kumar, V., Ratha, D., J. M. Lopez-Sanchez, A. Bhattacharya, H. McNairn, Y. S. Rao, and K. V. Ramana, 2020. Assessment of rice growth conditions in a semi-arid region of India using the Generalized Radar Vegetation Index derived from RADARSAT-2 polarimetric SAR data, Remote Sensing of Environment, 237: 111561.

<a id="4">[4]</a> 
Dey, S., Bhattacharya, A., Ratha, D., Mandal, D. and Frery, A.C., 2020. Target Characterization and Scattering Power Decomposition for Full and Compact Polarimetric SAR Data. IEEE Transactions on Geoscience and Remote Sensing.

<a id="5">[5]</a> 
Mandal, D., Kumar, V., Ratha, D., Dey, S., Bhattacharya, A., Lopez-Sanchez, J.M., McNairn, H. and Rao, Y.S., 2020. Dual polarimetric radar vegetation index for crop growth monitoring using sentinel-1 SAR data. Remote Sensing of Environment, 247, p.111954.

<a id="6">[6]</a> 
Mandal, D., Ratha, D., Bhattacharya, A., Kumar, V., McNairn, H., Rao, Y.S. and Frery, A.C., 2020. A Radar Vegetation Index for Crop Monitoring Using Compact Polarimetric SAR Data. IEEE Transactions on Geoscience and Remote Sensing, 58 (9), pp. 6321-6335.

<a id="7">[7]</a> 
V. Kumar, D. Mandal, A. Bhattacharya, and Y. S. Rao, 2020. Crop Characterization Using an Improved Scattering Power Decomposition Technique for Compact Polarimetric SAR Data. International Journal of Applied Earth Observations and Geoinformation, 88: 102052.

<a id="8">[8]</a> 
Kim, Y. and van Zyl, J.J., 2009. A time-series approach to estimate soil moisture using polarimetric radar data. IEEE Transactions on Geoscience and Remote Sensing, 47(8), pp.2519-2527.

<a id="9">[9]</a> 
Trudel, M., Charbonneau, F. and Leconte, R., 2012. Using RADARSAT-2 polarimetric and ENVISAT-ASAR dual-polarization data for estimating soil moisture over agricultural fields. Canadian Journal of Remote Sensing, 38(4), pp.514-527.

<a id="10">[10]</a> 
Barakat, R., 1977. Degree of polarization and the principal idempotents of the coherency matrix. Optics Communications, 23(2), pp.147-150.
