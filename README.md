# roxsi_largescale_dataproc

This repository contains code for processing the data in the large-scale arrays of ROXSI. Other libraries are necessary for specific tasks, but the bulk of the data processing should be self-sufficient in this repository. 

After downloading the repository, the data parent directory must be included in the function data_dirpath.m so that the rest of the code knows where the data are.

In quality control (QC) analyses, colormaps in figures may be set with the cmocean package, which can be [downloaded here](https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps), or see the [GitHub page here](https://github.com/matplotlib/cmocean).