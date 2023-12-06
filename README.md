# Werner_et.al._submitted_2023

# RingSegmentation and smoothing.ipynb
Python code for tracking cytokinetic ring closure, generating output polygone with coordinate data and and overlay movie. Followed by coordinate smoothing and calculation of ingression speed of individual ring segments. 
Inpput is an image stack of an end on rig that closes over time
Outputs: - images sequence of .tiff files of input file overlayed with fit polygone (blue) perfect fit circle (green) and position of centre of perfect fit circle
        - .jug dile containing raw coordinates of fit polygone
        - Excel file of smoothed coordinates XYs.xlsx
        - Excel file of smoothed speed data speeedDiffCS.xlsx 
# Wavelet_analysis.m
Matlab code speedDiffCS.xlxs out put from RingSegmentation and smoothing.ipynb to perform wavelet synchosqueeze transform on speed data and fit the output to an ANH model. Followed by calculating an plotting period distribution of ANH fit output adjusted for angular frequency
Input: speedDiffCS.xlxs from RingSegmentation and smoothing.ipynb
output: - HDF5 file containing zero padded smoothed speed data (speedAdj), instantenious periods (IpsAll) and corresponding amplitudes (AmpsAll)
        - weighted histogram of period distribution where amplitude is adjusted by the angular frequency

# pulsed input model shorter.ipynb
mathematical model describing the cortex as an active actomyosin gel under the control of a reaction-diffusion model of RhoA activity, both of which were calibrated by in vivo measurements from C. elegans zygotes (Staddon, M.F., Munro, E.M., and Banerjee, S. (2022). Pulsatile contractions and pattern formation in excitable actomyosin cortex. PLoS Comput Biol 18, e1009981.) Model was adapted to include a pulsatile stimulus. Variabel inputs are the pulsatile stimulus (ds), the strain (sigma) and the duration (duration).

# closing_ring_model.ipynb
adaptation of the previous with periodic boundaries so that the domain was allowed to close. This new model consists of a series of active viscoelastic elements, arranged in a ring, that are pulled inward due to actomyosin tension. variable inputs are stimulus (Sr) strain (sigma0), ring radius (radius), time before ring closure can occour (t_radial) and friction in the radial direction relative to tangential (gamma_radial)
