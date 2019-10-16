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
Install healpy using
>pip install healpy
* To enhance the computational efficiency, we recommend users to run the codes in HPC facility.

## Program routines
* 

* **ACDcommander15.py, ACDsevem18.py** compute the AC discrepancy maps with resolution Nside = 256 for the Commander 2015 and SEVEM 2018 CMB temparature maps from [Planck Legacy Archive](https://pla.esac.esa.int/#maps). They use the Fourier coefficients estimated by **coeff_cmb.py**


## Demo
When the dependent packages are installed successfully, users can obtain the following Mollweide projection view of the AC discrepancy maps with resolution Nside = 1024 using **ACD.py** for the Commander 2015 (left) and SEVEM 2018 CMB (right) temparature maps from [Planck Legacy Archive](https://pla.esac.esa.int/#maps).

  <img src="https://github.com/wangyg19/CMBProbe/blob/master/ACD_Commander2015_Nside1024_notitle.png" alt="acd_commander15" width="420"><img src="https://github.com/wangyg19/CMBProbe/blob/master/ACD_SEVEM2018_Nside1024_notitle.png" alt="acd_sevem18" width="420">


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
