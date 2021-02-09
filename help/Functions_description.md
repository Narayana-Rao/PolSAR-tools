## Functions description

**Full-pol functions**
----------------------
Full-pol functionalities require the SAR data in the form of covariance (C3) or coherency matrix (T3). A typical file structures of T3 and C3 matrices are as follows:

<center>

| C3 matrix files |              | T3 matrix files |              |
|-----------------|--------------|----------------:|-------------:|
| C11.bin         | C11.hdr      |         T11.bin |      T11.hdr |
| C12_real.bin    | C12_real.hdr |    T12_real.bin | T12_real.hdr |
| C12_imag.bin    | C12_imag.hdr |    T12_imag.bin | T12_imag.hdr |
| C13_real.bin    | C13_real.hdr |    T13_real.bin | T13_real.hdr |
| C13_imag.bin    | C13_imag.hdr |    T13_imag.bin | T13_imag.hdr |
| C22.bin         | C22.hdr      |         T22.bin |      T22.hdr |
| C23_real.bin    | C23_real.hdr |    T23_real.bin | T23_real.hdr |
| C23_imag.bin    | C23_imag.hdr |    T23_imag.bin | T23_imag.hdr |
| C33.bin         | C33.hdr      |         T33.bin |      T33.hdr |

</center>
<br>

 * ```GRVI``` (Generalized volume based Radar Vegetation Index): This functionality computes the generalized volume based radar vegetation index for full polarimetric SAR data. The required input and the computed output are as follows:

    ````python
        input : input_T3/C3_folder, window_size
        output: GRVI.bin
    ````
    
    The formulation of GRVI is as follows:

    <center>

    ![grvi](https://latex.codecogs.com/svg.latex?\Large&space;\centering\text{GRVI}=\left(1-\text{GD}_{\text{GV}}\right)\Big(\frac{p}{q}\Big)^{2\,\text{GD}_{\text{GV}}},\quad0\le\text{GRVI}\le1)

    </center> 

    where, GD<sub>GV</sub> is the geodesic distance between Kennaugh matrices (K) of the observed and the generalized volume scattering model, p,q are minimum and maximum value of distances between K matrices of the observed and elementary targets respectively. A detailed explanation of GRVI is available in [[2]](#2).


 * ```MF3CF``` (Model Free 3-Component decomposition for Full-pol data): This functionality computes the model free 3 component scattering power decomposition for full polarimetric SAR data. The required input and the computed output are as follows:

    ````python
        input : input_T3/C3_folder, window_size
        output: Ps_FP.bin, Pd_FP.bin, Pv_FP.bin, Theta_FP.bin
    ````
    
    The formulation of the scattering powers (P<sub>s</sub>: Surface, P<sub>d</sub>: Double bounce, P<sub>v</sub>: volume) is as follows:
    <center>

    ![tfp](https://latex.codecogs.com/svg.latex?\Large&space;\noindent\\\P_{d}^{\text{FP}}=\frac{m_{\text{FP}}{\text{Span}}}{2}{\left(1-\sin2\theta_{\text{FP}}\right)}\\\P_{v}^{\text{FP}}={\text{Span}}\left(1-m_{\text{FP}}\right)\\\P_{s}^{\text{FP}}=\frac{m_{\text{FP}}{\text{Span}}}{2}\left(1+\sin2\theta_{\text{FP}}\right))
    
    </center>
    where m<sub>FP</sub> is degree of polarization, &theta;<sub>FP</sub> scattering type parameter, Span is the sum of the diagonal elements os coherence matrix (T3).  The derivation of these parameters in-terms of coherancey matrix (T3) elements is as shown below. Further details can be obtained from [[4]](#4)

    <center>

    ![mfp](https://latex.codecogs.com/svg.latex?\Large&space;m_{\text{FP}}=\sqrt{1-\frac{27|\mathbf{T3}|}{\big(\mathrm{Trace}(\mathbf{T3})\big)^3}};\qquad{}\tan\theta_{\text{FP}}=\frac{m_{\text{FP}}{\text{Span}}\left(T_{11}-T_{22}-T_{33}\right)}{T_{11}\left(T_{22}+T_{33}\right)+m_{\text{FP}}^{2}{\text{Span}}^{2}})
    
    </center> 

    <center>

    ![spanfp](https://latex.codecogs.com/svg.latex?\Large&space;\text{Span}=T_{11}+T_{22}+T_{33})
    
    </center>


 * ```RVI``` (Radar Vegetation Index): This functionality computes the Radar vegetation index for full polarimetric SAR data. The required input and the computed output are as follows:
    ````python
        input : input_T3/C3_folder, window_size
        output: RVI_FP.bin
    ````
    
    The formulation of RVI is as follows:

    <center>

    ![rvifp](https://latex.codecogs.com/svg.latex?\Large&space;\text{RVI}_{fp}=\frac{4\lambda_1}{\lambda_1+\lambda_2+\lambda_3})
        
    </center>
where, &lambda;<sub>1</sub>, &lambda;<sub>2</sub> and &lambda;<sub>3</sub> are the eigen values of coherency matrix (T3) in descending order (&lambda;<sub>1</sub>> &lambda;<sub>2</sub>>&lambda;<sub>3</sub>). Further details can be found in [[8]](#8)

<!-- <center>

![rvifp](https://latex.codecogs.com/svg.latex?\Large&space;\text{RVI}_{fp}=\frac{8\sigma^\circ_{\text{HV}}}{\sigma^\circ_{\text{HH}}+\sigma^\circ_{\text{VV}}+2\sigma^\circ_{\text{HV}}})
    
</center>  -->


 * PRVI (Polarimetric Radar Vegetation Index) : This functionality computes the polarimetric Radar vegetation index for full polarimetric SAR data. The required input and the computed output are as follows:
    
    ````python
        input : input_T3/C3_folder, window_size
        output: PRVI_FP.bin
    ````
    The formlation of PRVI interms of degree of polarization and cross-pol backscatter intensity can be expressed as follows: 
    <center>

    ![grvi](https://latex.codecogs.com/svg.latex?\Large&space;\text{PRVI}_{fp}=(1-\text{DOP}_{fp})\sigma^\circ_{\text{XY}})
        
    </center> 

    where DOP<sub>fp</sub> 3D Barakt degree of polarization and can be expressed as shown below. Further details on the PRVI can be found in [[1]](#1)

 * DOP (Degree of Polarization):  This functionality computes the 3D Barakat degree of polarization for full polarimetric SAR data. The required input and the computed output are as follows:

    ````python
        input : input_T3/C3_folder, window_size
        output: DOP_FP.bin
    ````
    <center>

    ![dopfp](https://latex.codecogs.com/svg.latex?\Large&space;\text{DOP}_{fp}=\sqrt{1-\frac{27\times\text{det([T3])}}{\text{(Trace[T3])}^3}})
        
    </center> 
    Further details on the Barakat Degree of polarization can be found in [[10]](#10)





**Compact-pol functions**
-------------------------
Compact-pol functionalities require the SAR data in the form of 2x2 covariance matrix (C2). A typical file structures of C2 matrix is as follows:

<center>

| C2 matrix files |              |
|-----------------|--------------|
| C11.bin         | C11.hdr      |
| C12_real.bin    | C12_real.hdr |
| C12_imag.bin    | C12_imag.hdr |
| C22.bin         | C22.hdr      |

</center>
<br>

 * Model Free 3-Component decomposition for Full-pol data (MF3CF) 
    
    <center>

    ![mcp](https://latex.codecogs.com/svg.latex?\Large&space;m_{\text{CP}}=\sqrt{1-\frac{4|\mathbf{C2}|}{\big(\mathrm{Trace}(\mathbf{C2})\big)^2}};\qquad{}\tan\theta_{\text{CP}}=\frac{m_{\text{CP}}{S_0}\left(\text{OC}-\text{SC}\right)}{\text{OC}\times\text{SC}+m_{\text{CP}}^{2}{S_0}^{2}})
    
    </center> 


    <center>

    ![cppowers](https://latex.codecogs.com/svg.latex?\Large&space;\\\P_{d}^{\text{CP}}=\frac{m_{\text{FP}}{S_0}}{2}{\left(1-\sin2\theta_{\text{CP}}\right)};\\\P_{v}^{\text{CP}}={S_0}\left(1-m_{\text{CP}}\right);\\\P_{s}^{\text{CP}}=\frac{m_{\text{CP}}{S_0}}{2}\left(1+\sin2\theta_{\text{CP}}\right))
    
    </center>
    <center>

    ![sparm](https://latex.codecogs.com/svg.latex?\Large&space;\\\S_0=\text{C11+C22};\qquad{}S_1=\text{C11-C22};\\\S_2=\text{C12+C21};\qquad{}S_3=\pm\text{j(C12-C21)})
    
    ![cparm](https://latex.codecogs.com/svg.latex?\Large&space;\text{SC}=\frac{S_0-S_3}{2};\qquad{}\text{OC}=\frac{S_0+S_3}{2};)
    
    </center>


````python
    input : input_C2_folder, window_size, tau
    output: Ps_CP.bin, Pd_CP.bin, Pv_CP.bin, Theta_CP.bin
````

 * Compact-pol Radar Vegetation Index (CpRVI) 
    <center>
    
    ![grvi](https://latex.codecogs.com/svg.latex?\Large&space;\text{CpRVI}=\left(1-\dfrac{3}{2}\text{GD}_{\text{ID}}\right)\Big(\frac{p}{q}\Big)^{2(\frac{3}{2}\text{GD}_{\text{ID}})})

    ![pqcp](https://latex.codecogs.com/svg.latex?\Large&space;p=\text{min\\{SC,OC\\}},q=\text{max\\{SC,OC\\}})
    
    ![cparm](https://latex.codecogs.com/svg.latex?\Large&space;\text{SC}=\frac{S_0-S_3}{2};\qquad{}\text{OC}=\frac{S_0+S_3}{2};)
    

    ![sparm](https://latex.codecogs.com/svg.latex?\Large&space;\\\S_0=\text{C11+C22};\qquad{}S_1=\text{C11-C22};\\\S_2=\text{C12+C21};\qquad{}S_3=\pm\text{j(C12-C21)})


    </center> 

where, GD<sub>ID</sub> is the geodesic distance between Kennaugh matrices (K) of the observed and the ideal depolarizer, p,q are minimum and maximum values of SC and OC which are functions of stocks parameters (S<sub>0</sub>, S<sub>1</sub>, S<sub>2</sub>, and S<sub>3</sub>). A detailed explanation of CpRVI is available in [[6]](#6).


````python
    input : input_C2_folder, window_size
    output: CpRVI.bin
````



 * improved S-&Omega; decomposition (iS-&Omega;) 
    
    The stokes paramters can be written in terms of the covariance matrx (C2) elements as follows:
    
    <center>
    
    ![sparm](https://latex.codecogs.com/svg.latex?\Large&space;\\\S_0=\text{C11+C22};\qquad{}S_1=\text{C11-C22};\\\S_2=\text{C12+C21};\qquad{}S_3=\pm\text{j(C12-C21)})    

    </center> 

    Then, the parameters Same-sense Circular (SC) and Opposite-sense Circular (OC) can be expressed as follows:

    <center>

    ![cparm](https://latex.codecogs.com/svg.latex?\Large&space;\text{SC}=\frac{S_0-S_3}{2};\qquad{}\text{OC}=\frac{S_0+S_3}{2};)
    </center>    
    <!-- <center>

    ![cparm](https://latex.codecogs.com/svg.latex?\Large&space;\Vec{\mathbf{S}}=\begin{bmatrix}S_{0}\\\S_{1}\\\S_{2}\\\\S_{3}\end{bmatrix}=\begin{bmatrix}C_{11}+C_{22}\\\C_{11}-C_{22}\\\C_{12}+C_{21}\\\pm\left(C_{12}-C_{21}\right)\end{bmatrix})

    </center> -->

    <center>

    ![cparm](https://latex.codecogs.com/svg.latex?\Large&space;\text{SC/OC}<1;\qquad{}\qquad{}\qquad{}\text{SC/OC}>1\\\P_s=\Omega\left(S_{r0}-\left(1-\Omega\right)\text{SC}\right);\qquad{}P_s=\Omega\left(1-\Omega\right)\text{OC}\\\P_d=\Omega\left(1-\Omega\right)\text{SC};\qquad{}P_d=\Omega\left(S_{r0}-\left(1-\Omega\right)\text{OC}\right))

    ![cparm](https://latex.codecogs.com/svg.latex?\Large&space;P_v=S_{r0}\left(1-\Omega\right))

    
    </center>    

````python
    input : input_C2_folder, window_size, tau, psi, chi
    output: Ps_iSOmega.bin, Pd_iSOmega.bin,Pv_iSOmega.bin
````

 * Degree of Polarization (DOP<sub>cp</sub>) 
    
    The conventional degree of polastokes paramters can be written in terms of the covariance matrx (C2) elements as follows:
    
    <center>
    
    ![sparm](https://latex.codecogs.com/svg.latex?\Large&space;\\\S_0=\text{C11+C22};\qquad{}S_1=\text{C11-C22};\\\S_2=\text{C12+C21};\qquad{}S_3=\pm\text{j(C12-C21)})    

    </center>  

    then,

    <center>

    ![dopcp](https://latex.codecogs.com/svg.latex?\Large&space;\text{DOP}_{cp}=\frac{\sqrt{S^2_1+S^2_2+S^2_3}}{S_0})
        
    </center> 

    [[11]](#11)

````python
    input : input_c2_folder, window_size, tau
    output: DOP_CP.bin
````





**Dual-pol**
------------
Dual-pol functionalities require the SAR data in the form of 2x2 covariance matrix (C2). A typical file structures of C2 matrix is as follows:

<center>

| C2 matrix files |              |
|-----------------|--------------|
| C11.bin         | C11.hdr      |
| C12_real.bin    | C12_real.hdr |
| C12_imag.bin    | C12_imag.hdr |
| C22.bin         | C22.hdr      |

</center>
<br>

 * Dual-pol Radar Vegetation Index (DpRVI) 

    <center>

    ![grvi](https://latex.codecogs.com/svg.latex?\Large&space;\text{DpRVI}=1-\Big(\frac{\lambda_1}{\lambda_1+\lambda_2}\Big)\sqrt{1-\frac{4\times\text{det([C2])}}{\text{(Trace[C2])}^2}})
    
    </center> 

where, C2 is co-variance matrix,  and  &lambda;<sub>1</sub> and &lambda;<sub>2</sub> are the eigen values of C2 matrix in descending order.

````python
    input : input_C2_folder, window_size
    output: DpRVI.bin
````

 * Radar Vegetation Index (RVI<sub>dp</sub>) 
<center>

![grvi](https://latex.codecogs.com/svg.latex?\Large&space;\text{RVI}_{dp}=\frac{4\sigma^\circ_{\text{XY}}}{\sigma^\circ_{\text{XX}}+\sigma^\circ_{\text{XY}}})
    
</center> 


````python
    input : input_c2_folder, window_size
    output: RVI_dp.bin
````

 * Polarimetric Radar Vegetation Index (PRVI<sub>dp</sub>) 
<center>

![grvi](https://latex.codecogs.com/svg.latex?\Large&space;\text{PRVI}_{dp}=(1-\text{DOP}_{dp})\sigma^\circ_{\text{XY}})
    
</center> 


````python
    input : input_c2_folder, window_size
    output: PRVI_dp.bin
````

 * Degree of Polarization (DOP<sub>dp</sub>) 
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

<a id="11">[11]</a> 
Lee, J.S. and Pottier, E., 2009. Polarimetric radar imaging: from basics to applications. CRC press.