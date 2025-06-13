The MATLAB functions included in this repository are for the purpose of reconstructing the planar pressure field from two-component planar velocity data. The code was developed at the University of Southampton.

Instantaneous pressure field may be reconstructed using a Poisson solver using the functions "EU_FDM" and "TH_FDM". The former is used for when the time-derivatives are known and the latter for when Taylor's hypothesis must be invoked to estimate the time derviates.

The mean pressure may be reconstructed using a Poisson solver using the function "RANS_FDM".

A secondary option is to instead using the "OMNIPOL" function, which using a single-iteration omnidirectional integration with an arbitrary number of paths.

NOTE: Boundary conditions need to be set within the functions themselves. See "RANS_FDM" for an example of choosing different boundary conditions within the code.

If using this code in a published medium, please cite the following articles that used successive iterations of the pressure reconstruction code made available here:

1. van der Kindere JW, Laskari A, Ganapathisubramani B, de Kat R. Pressure from 2D snapshot PIV. Exp Fluids. 2019;60(2):32. doi: 10.1007/s00348-019-2678-5.
2. Ferreira, M. A., & Ganapathisubramani, B. (2020). PIV-based pressure estimation in the canopy of urban-like roughness. Experiments in Fluids, 61(3), 70.
3. Carter, Douglas W., and Bharathram Ganapathisubramani. "Low-Order Modeling and Sensor-Based Prediction of Stalled Airfoils at Moderate Reynolds Number." AIAA Journal (2023): 1-13.

If you are interested in results obtained from using this code for examining physical problems, then please see these articles: 

1. Ferreira, M. A., and B. Ganapathisubramani. "Scale interactions in velocity and pressure within a turbulent boundary layer developing over a staggered-cube array." Journal of Fluid Mechanics 910 (2021): A48.
2. Carter, D. W., & Ganapathisubramani, B. (2023). Data-driven determination of low-frequency dipole noise mechanisms in stalled airfoils. Experiments in Fluids, 64(2), 41.



