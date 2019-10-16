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

## Environment Configuration and Package Installation
The code has been tested in Python 3.6 environment. The CMBProbe relies on healpy package
* [healpy](https://healpy.readthedocs.io/en/latest/): Zonca, A., Singer, L., Lenz, D., Reinecke, M., Rosset, C., Hivon, E., and Gorski, K. (2019). "[Healpy: equal area pixelization and spherical harmonics transforms for data on the sphere in python. J. Open Source Softw., 4:1298.](https://joss.theoj.org/papers/10.21105/joss.01298)".
* To enhance the computational efficiency, we recommend users to run the codes in HPC facility.

## Functions and Folders
* **utils**: This folder contains some basic tools/resources/auxiliary functions used for implementing our main functions, including
   1. SD: A folder that contains six examples of [symmetric spherical design points](https://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/ss.html), which are used in our paper. 
   2. tangent_field: A folder that contains several functions for generating three vector fields and their visualization used in our paper. These functions come from E. J. Fuselier and G. B. Wright as they use in "[*Stability and error estimates for vector field interpolation and decomposition on the sphere with RBFs. SIAM Journal on Numerical Analysis, 47(5):3213-39*](https://epubs.siam.org/doi/abs/10.1137/080730901)".
   3. m_map: A [mapping package](https://www.eoas.ubc.ca/~rich/map.html#ack) for Matlab. We have used some functions of this package for visualization of tangent fields. 
   4. QpS2.m: A function used for computing the weights and quadrature nodes (for a given degree and a specific type of quadrature rule) in both Cartesian and spherical coordinates. 

* **Setup.m**: The script for downloading NFFT package (compatible with the user's operating system and computing environment) and unzipping and installing the package in the current folder of FaVeST. If the installation is successful, users can test FaVeST in the demos and examples. The **Setup.m** and **Demo.m** have been tested on **Ubuntu 16.04.6, macOS High Sierra and Mojave, Windows7，8，10**. Please make sure that your operating system meets the requirement of NFFT package.</span>

* **Demo.m**: The picture below shows the AC discrepancy map (right) with resolution Nside = 256 **ACD.py** for the Commander 2015 CMB temparature map (left) from [Planck Legacy Archive](https://pla.esac.esa.int/#maps).


## Demo
When the dependent packages are installed successfully, users can obtain the following Mollweide projection view of the AC discrepancy maps with resolution Nside = 1024 using **ACD.py** for the Commander 2015 and SEVEM 2018 CMB temparature maps from [Planck Legacy Archive](https://pla.esac.esa.int/#maps).

<img src="https://github.com/wangyg19/CMBProbe/blob/master/ACD_Commander2015_Nside1024_notitle.png" width="450" class="center"><img src="https://github.com/wangyg19/CMBProbe/blob/master/ACD_SEVEM2018_Nside1024_notitle.png" width="450" class="center">


## Acknowledgements
Some of the results in this paper have been derived using the [HEALPix package](https://healpix.sourceforge.io/) [Gorski et al. (2005)](https://arxiv.org/abs/astro-ph/0409513). The authors acknowledge support from the Australian Research Council under Discovery Project DP180100506.

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
## Notes
The package **CMBProbe** may be used for any research purposes under the following conditions:
* The user must acknowledge the use of **CMBProbe** in publications resulting from the use of the functions/tools.
* The user may not redistribute **CMBProbe** without separate permission.
* The user may not use this information for any commercial or revenue-bearing purposes without first obtaining permission from us.
