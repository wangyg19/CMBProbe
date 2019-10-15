# CMBProbe
This is the Python package for AC discrepancy maps in the paper

>Ian H. Sloan, Quoc T. Le Gia, Yu Guang Wang, Robert S. Womersley and Jan Hamann. [A new Probe of non-Gaussianity for Maps of the CMB](). arXiv preprint arXiv:, 2019.

Author: Dr. Yu Guang Wang (yuguang.wang@unsw.edu.au).

## Abstract
We introduce a newmathematical tool (a direction-dependent probe) to analyse the randomness
of purported isotropic Gaussian random fields on the sphere. We apply the probe to assess
the full-sky cosmic microwave background (CMB) temperature maps produced by the Planck
collaboration (PR2 2015 and PR3 2018), with special attention to the inpainted maps. To
study the randomness of the fields represented by each map we use the autocorrelation of the
sequence of probe coefficients (which are just the full-sky Fourier coefficients $a_{\ell 0}$ if the z axis
is taken in the probe direction). If the field is isotropic and Gaussian then the probe coefficients
for a given direction should be realisations of uncorrelated scalar Gaussian variables. We find
that for most of the maps there are many directions for which this is not the case.We make a first
attempt at justifying the features of the temperature maps that are contributing to the apparent
lack of randomness. In the case of Commander 2015 we mimic an aspect of the observed
behavior with a model field that is Gaussian but not isotropic. In contrast, the non-inpainted
2018 SEVEM map (which has visible equatorial pollution) is modelled by an isotropic Gaussian
random field plus a non-random needlet-like structure located near the galactic center.

## Citation 
If you use our codes and datasets, please cite:
```
@article{CMBProbe,
  title={A new Probe of non-Gaussianity for Maps of the CMB},
  author={Sloan, Ian H. and Le Gia, Quoc T. and and Wang, Yu Guang and Womersley, Robert S. and Hamann, Jan},
  journal={arXiv preprint arXiv:},
  year={2019}
}
```
## Environment Configuration and Package Installation
The code has been tested in Python 3.6 environment. This package relies on healpy package
* [healpy](https://healpy.readthedocs.io/en/latest/): Zonca, A., Singer, L., Lenz, D., Reinecke, M., Rosset, C., Hivon, E., and Gorski, K. (2019). "[Healpy: equal area pixelization and spherical harmonics transforms for data on the sphere in python. J. Open Source Softw., 4:1298.](https://joss.theoj.org/papers/10.21105/joss.01298)".

## Functions and Folders
* **utils**: This folder contains some basic tools/resources/auxiliary functions used for implementing our main functions, including
   1. SD: A folder that contains six examples of [symmetric spherical design points](https://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/ss.html), which are used in our paper. 
   2. tangent_field: A folder that contains several functions for generating three vector fields and their visualization used in our paper. These functions come from E. J. Fuselier and G. B. Wright as they use in "[*Stability and error estimates for vector field interpolation and decomposition on the sphere with RBFs. SIAM Journal on Numerical Analysis, 47(5):3213-39*](https://epubs.siam.org/doi/abs/10.1137/080730901)".
   3. m_map: A [mapping package](https://www.eoas.ubc.ca/~rich/map.html#ack) for Matlab. We have used some functions of this package for visualization of tangent fields. 
   4. QpS2.m: A function used for computing the weights and quadrature nodes (for a given degree and a specific type of quadrature rule) in both Cartesian and spherical coordinates. 

* **Setup.m**: The script for downloading NFFT package (compatible with the user's operating system and computing environment) and unzipping and installing the package in the current folder of FaVeST. If the installation is successful, users can test FaVeST in the demos and examples. The **Setup.m** and **Demo.m** have been tested on **Ubuntu 16.04.6, macOS High Sierra and Mojave, Windows7，8，10**. Please make sure that your operating system meets the requirement of NFFT package.</span>

* **Demo.m**: It tests **FaVeST_fwd.m** and **FaVeST_adj.m** on a tangent field. It is used to test whether users have successfully configured NFFT packages by **Setup.m**. 

* **FaVeST_fwd.m**: Main function for implementing forward FFTs computing Fourier coefficients associated with a quadrature rule: T - tangent field samples; L - degree for vector spherical harmonic; X,w - quadrature rule used for evaluating FFT. See Algorithm 1 in our paper.

* **FaVeST_adj.m**: Main function for implementing adjoint FFTs for vector spherical harmonic expansion with given inputs: alm -  Fourier coefficients for divergent-free part; blm - Fourier coefficients of curl-free part; X - evaluation points on the sphere.


* **Fig2a,2b,2c.m**, **Fig3a,3b,3c.m**, **Table1.m**, **Table2_Fig4.m**: These routines are used to reproduce the numerical results of the corresponding figures and tables of our paper.

## Demo
After running **Setup.m** successfully (meaning that NFFT package folder appears in **FaVeST** folder), users can obtain the following visualization results that is the same as Fig.2(a) in our paper. Then, good news~all is in order now, uers can try different settings by running **Demo.m**, or test our simulation programs by running other m-scripts (e.g. **Fig3a.m**). 

<img src="https://github.com/mingli-ai/FaVeST/blob/master/images/vf_1_gl.png" width="250"><img src="https://github.com/mingli-ai/FaVeST/blob/master/images/vf_1_rec_gl.png" width="250"><img src="https://github.com/mingli-ai/FaVeST/blob/master/images/vf_1_err_gl.png" width="250">


## Acknowledgement
Some of the results in this paper have been derived using the [HEALPix package](https://healpix.sourceforge.io/) [Gorski et al. (2005)](https://arxiv.org/abs/astro-ph/0409513). The authors acknowledge support from the Australian Research Council under Discovery Project DP180100506.

## Notes
The package **FaVeST** may be used for any research purposes under the following conditions:
* The user must acknowledge the use of **FaVeST** in publications resulting from the use of the functions/tools.
* The user may not redistribute **FaVeST** without separate permission.
* The user may not use this information for any commercial or revenue-bearing purposes without first obtaining permission from us.
