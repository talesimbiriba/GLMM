# GLMM
This toolbox contains several scripts and functions for MATLAB, to unmix hyperspectral data using the Generalized Linear Mixing Model (GLMM).

This code was built over the ELMM code provided by Lucas Drumetz and his collegues. 

The use and redistribution of this code is allowed for noncomercial purposes as long as the copyrights presented here are replicated.

Author: Tales Imbiriba.

Date: April 2018.

Ref.:

[1] Imbiriba, T., Borsoi, R. A., Bermudez, J. C. M. (2018). Generalized linear mixing model accounting for endmember variability. 2018 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP).


The contents include:

- GLMM_ADMM.m: Function performing the unmixing with the GLMM
- GLMM_RealData.m: example of use of the function on a real hyperspectral dataset
- real_data_1.mat: real dataset used (crop of the DFC 2013 data)
- endmembers_houston.mat: reference endmember matrix used in the demo
- FCLSU.m : function performing the standard fully constrained least squared unmixing.
- CLSU.m : function performing the standard partially constrained least squared unmixing.
- SCLSU.m: function performing a scaled version of CLSU, which follows a particular case of the GLMM
- soft.m: soft thresholding, proximal operator of the L1 norm
- vector_soft_col.m: vector soft thresholding, proximal operator of the L2 norm
- pca_viz.m: function projecting data and endmembers on the space spanned by the first three principal components of the data and displaying a scatterplot
- rescale.m: function rescaling hyperspectral data between 0 and 1
- RMSEAndSTDForMatrix.m: dummy function to compute RMSE for matrix or vectors.
