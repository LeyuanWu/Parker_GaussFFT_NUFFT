# Parker_GaussFFT_NUFFT
Data and code for "Modified Parker's method for gravitational forward and inverse modeling using general polyhedral models"
______________________________________
## Abstract
We propose a new modified Parker's method for efficient gravitational forward modeling and inversion using general polyhedral models. We have made several important modifications to the classical method, including: 1) The new method is now applicable to a general polyhedron represented by triangulated surface or tetrahedral mesh, and with arbitrarily variable 3D density distribution. 2) An optimal Fourier-domain sampling strategy is used to improve the numerical accuracy of the new algorithm significantly. 3) A simple and effective automatic layering technique is introduced to accelerate the convergence rate of Parker's method. The method is demonstrated using both synthetic and real polyhedral models, including a sphere model approximated by a polyhedron, two asteroids, a digital elevation model in the Himalaya region, and the Yucca Flat basin model in Nevada. The numerical results show that, compared with analytical solutions of polyhedral models in the space domain, the modified Parker's method can improve the computational efficiency by several orders of magnitude while obtaining almost the same simulation results. The difference is well below existing instrumentation error level. By embedding the new forward algorithm into an iterative process, it can be used for fast inversion of density interfaces. Our new method is suitable for the efficient modeling and inversion of gravitational potential, gravitational vector, and gravitational gradient tensor caused by polyhedral models with a large number of faces, representing geological abnormal bodies, asteroids, and single or multilayer density interface models with triangulated surfaces.
_______________________________________
## Part 0: SUMMARY
This file contains data and MATLAB codes to reproduce the numerical results in the manuscript entitled 
"Modified Parker's method for gravitational forward and inverse modeling using general polyhedral models",
which has been submitted to Journal of Geophysical Research-Solid Earth

- Contact: Leyuan Wu <br>
- College of Science, Zhejiang University of Technology <br>
- Hangzhou, 310023, China <br>
- leyuanwu@zjut.edu.cn & leyuanwu@zju.edu.cn <br>

Last edited 2021-09-24
_______________________________________
## Part 1: INSTALLATION

The code is based on Non-Uniform Fast Fourier Transform (NUFFT) algorithms. It can be implemented using either 
the Courant Mathematics and Computing Laboratory (CMCL) NUFFT software package, which is freely available at   <br>
https://cims.nyu.edu/cmcl/nufft/nufft.html  <br>
or the Flatiron Institute Nonuniform Fast Fourier Transform (FINUFFT) software package, which is freely available at   <br>
https://finufft.readthedocs.io/en/latest/   <br>
Download either of the two NUFFT packages (or both of them for comparison) and load it into your MATLAB path.

### NOTE:
In our experience, the CMCL NUFFT software package can be used directly, while the FINUFFT software package
may have to be compiled locally (on your own computer) to make it work. 
Please visit https://finufft.readthedocs.io/en/latest/  for more details on the installation of the FINUFFT software package.
The numerical results provided in this manuscript is based on the FINUFFT software package. If the CMCL NUFFT 
software package is used instead, the numerical accuracy is about the same, while the time costs would be several times
longer than using the FINUFFT software package, depending on the specific model sizes.
_______________________________________
## Part 2: Source code file list

*List of source code files to implement the Fourier forward algorithm and neccessary numerical comparisons.*

    - phi_xyz_PolyhedronSI_Parker.m  ---- M-file to calculate the complete gravity caused by a polyhedron with constant density
    - phi_xyz_PolyhedronVI_Parker.m  ---- M-file to calculate the complete gravity caused by a polyhedron with variable density
    - phi_xyz_DEMSI_Parker.m  ---- M-file to calculate the complete gravity caused by a DEM with constant density
    - phi_xyz_sphere.m  ---- M-file to calculate the external complete gravity caused by a uniform sphere
    - sifft2.m  ---- M-file to calculate the 2D inverse shift FFT of a 2D spectrum F
    - sifft_X.m  ---- M-file to calculate the 1D inverse shift FFT along each row of a matrix F
    - sifft_Y.m  ---- M-file to calculate the 1D inverse shift FFT along each column of a matrix F
    - lgwt.m  ---- M-file to produces the Legendre-Gauss weights and nodes for computing definite integrals
    - DEM2TriMesh.m  ---- M-file to convert a DEM into a triangular mesh polyhedral model
    - GEOTIFF_READ.m  ---- M-file to read geotiff using imread and assign map info from infinfo.
    - triquad.m  ---- M-file to calculate the N^2 Gaussian nodes and weights for a bunch of 3D triangles simultaneously
    - tetraquad.m ---- M-file to calculate the N^3 Gaussian nodes and weights for a bunch of tetrahedrons simultaneously
    - sta_min_max_mean_rms.m ---- M-file to calculate statistics of the numerical results compared to the reference results
    - myplot.m ---- M-file to draw plots with different colors, markers and line styles

_______________________________________
## Part 3: Example code file list

*List of M-files and related .mat data files to carry out the numerical examples in the manuscript.*

- Example 1:
    - Sphere_Icosahedron_Const_VgT.m  ---- M-file to calculate the "constant sphere" model
|***Numerical results of example 1 are summarized in Figure 4, Figure 5 and Table 2 (Upper half)***

- Example 2:
    - Sphere_Icosahedron_Variable_RadialLinear_VgT.m ---- M-file to calculate the "variable sphere" model
|***Numerical results of example 2 are summarized in Table 2 (Lower half)***

- Example 3:
    - EROS_Const_Werner_GEEC_Fourier_VgT.m  ---- M-file to calculate the "constant EROS" model
        - EROS.mat   ---- mat-file containing the vertices and faces of the EROS model
        - EROS_Const2670_11272_22540_z0m6d3_dx0d1_dy0d1_VgT_Werner.mat  ---- mat-file containing the VecWerSch results
        - EROS_Const2670_11272_22540_z0m6d3_dx0d1_dy0d1_gT_Saraswati.mat ---- mat-file containing the GEEC forward results
***Numerical results of example 3 are summarized in Figure 7, Figure 8 and Table 3 (Upper half)***

- Example 4:
    - Bennu_Variable_RenZY_Fourier_gT.m ---- M-file to calculate the " variable Bennu" model
        - Bennu.mat   ---- mat-file containing the vertices and faces of the Bennu model
        - Bennus49152_pq1.2_295272.node   ---- node list of the tetrahedral mesh representation of the Bennu model
        - Bennus49152_pq1.2_1432834.ele   ---- element list of the tetrahedral mesh representation of the Bennu model
        - Bennu_Variable_25350_49152_z0m0d28_dx0d01_dy0d01_gT_RenZY.mat    ---- mat-file containing the GraPly results
***Numerical results of example 4 are summarized in Table 3 (Lower half)***

- Example 5:
    - /SRTMv4d1_Werner_PrismSum_Fourier_VgT.m ---- M-file to calculate the " constant DEM" model
        - srtm_54_07.tif   ---- .tif SRTM v4.1 DEM data
        - srtm_54_07_lon_p86d0_p88d0_lat_p27d0_p29d0_ne_Werner.mat    ---- mat-file containing the VecWerSch results
        - srtm_54_07_lon_p86d0_p88d0_lat_p27d0_p29d0_ne_PrismSum.mat  ---- mat-file containing the Prism-Summation results
***Numerical results of example 5 are summarized in Figure 11, Figure 12 and Table 4***
_______________________________________
## Part 4: Other related data and code:

- Different polyhedral-shaped models of the asteroid 101955 Bennu are freely available at the Web site  <br>
https://www.asteroidmission.org/updated-bennu-shape-model-3d-files

- Multiple polyhedral models of the asteroid 433 EROS can be downloaded at   <br>
https://sbn.psi.edu/pds/resource/nearbrowse.html

- Digital terrain data SRTM v4.1 are downloaded from http://dwtkns.com/srtm/

- The CMCL NUFFT software package is freely available at the Web site https://cims.nyu.edu/cmcl/nufft/nufft.html

- The FINUFFT software package is freely available at the Web site https://finufft.readthedocs.io/en/latest/

- The GraPly software package is freely available at the Web site https://github.com/zhong-yy/GraPly

- The GEEC software package is freely available at the Web site https://github.com/anitasaraswati/GEEC
______________________________________
## References: 
- [1]   Leslie Greengard and June-Yub Lee, 2004, Accelerating the Nonuniform Fast Fourier Transform. SIAM Review, 46:3, 443-454
- [2]   Lee, J.Y. and Greengard, L., 2005. The type 3 nonuniform FFT and its applications. Journal of Computational Physics, 206(1), pp.1-5.
- [3]   Alexander H. Barnett, Jeremy Magland, and Ludvig af Klinteberg, 2019. A Parallel Nonuniform Fast Fourier Transform Library Based on an "Exponential of Semicircle" Kernel. SIAM Journal on Scientific Computing 41:5, C479-C504
- [4]   Barnett, A.H., 2021. Aliasing error of the exp⁡(b1-z2) kernel in the nonuniform fast Fourier transform. Applied and Computational Harmonic Analysis, 51, pp.1-16.
- [5]   Wu, L., and G. Tian, 2014, High-precision Fourier forward modeling of potential fields: Geophysics, 79, no. 5, G59-G68
- [6]   Zhengyong Ren, Yiyuan Zhong, Chaojian Chen, Jingtian Tang, Thomas Kalscheuer, Hansruedi Maurer, and Yang Li. 2018. Gravity Gradient Tensor of Arbitrary 3D Polyhedral Bodies with up to Third-Order Polynomial Horizontal and Vertical Mass Contrasts. Surveys in Geophysics, 39(5):901-935, (Open access) DOI: https://doi.org/10.1007/s10712-018-9467-1
- [7]   Saraswati, A. T., Cattin, R., Mazzotti, S., & Cadio, C. 2019. New analytical solution and associated software for computing full-tensor gravitational field due to irregularly shaped bodies. Journal of Geodesy, 1-17, doi:10.1007/s00190-019-01309-y
- [8]   Greg von Winckel (2021). Gaussian Quadrature for Triangles (https://www.mathworks.com/matlabcentral/fileexchange/9230-gaussian-quadrature-for-triangles), MATLAB Central File Exchange. Retrieved June 5, 2021.
- [9]   Greg von Winckel (2021). Gauss Quadrature for Tetrahedra (https://www.mathworks.com/matlabcentral/fileexchange/9389-gauss-quadrature-for-tetrahedra), MATLAB Central File Exchange. Retrieved June 5, 2021.
- [10] Greg von Winckel (2021). Legendre-Gauss Quadrature Weights and Nodes (https://www.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes), MATLAB Central File Exchange. Retrieved September 23, 2021.
- [11] Yushin Ahn (2021). geotiff reader (https://www.mathworks.com/matlabcentral/fileexchange/29425-geotiff-reader), MATLAB Central File Exchange. Retrieved June 5, 2021.
