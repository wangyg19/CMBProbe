# CMBProbe
This is the Python package for AC discrepancy maps in the paper

>Ian H. Sloan, Quoc T. Le Gia, Yu Guang Wang, Robert S. Womersley and Jan Hamann. [A New Probe of non-Gaussianity for CMB Maps](). arXiv preprint arXiv:, 2019.

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

## Dependent Package
The codes have been tested in Python 3.6 environment. The CMBProbe relies on healpy package
* [healpy](https://healpy.readthedocs.io/en/latest/): Zonca, A., Singer, L., Lenz, D., Reinecke, M., Rosset, C., Hivon, E., and Gorski, K. (2019). "[healpy: equal area pixelization and spherical harmonics transforms for data on the sphere in python. J. Open Source Softw., 4:1298.](https://joss.theoj.org/papers/10.21105/joss.01298)".
Install healpy using
```
pip install healpy
```
* To enhance the computational efficiency, we recommend users to run the codes in HPC facility.

## Program Routines
* **coeff_cmb.py** estimates the Fourier coefficients from CMB map (which is in .fits format). 

* **autocorr.py** computes autocorrelation for a sequence.

* **ACDcommander15.py, ACDsevem18.py** compute the AC discrepancy maps with resolution Nside = 256 for the Commander 2015 and SEVEM 2018 CMB temparature maps from [Planck Legacy Archive](https://pla.esac.esa.int/#maps). They use the Fourier coefficients estimated by **coeff_cmb.py**.

* **ACD_plt.py** plots the AC discrepancy map and saves it to .png format.

* **colormap.py** configures the colormap for CMB and AC discrepancy maps using platter **cmbmap.mat**.

* **cmbmap.mat** is the colormap (platter) for CMB maps and AC discrepancy maps.

* **Some Parameters** in the programs: 
* map_type   -  one of four CMB map methods: commander, nilc, sevem, smica. 
* vs         -  year version of CMB map, either '2015' or '2018'.
* plt_q      -  flag to choose inpainted or non-inpainted CMB maps, either 'I_inpainted' or 'I_noninpainted'. 
* L          -  max degree of Fourier coefficients used in computing AC discrepancy. 
* Nside_acd  -  Nside for the AC discrepancy map; lag_acmap - max lag for AC discrepancy.


## Demo
When the dependent packages are installed successfully, users can obtain the following Mollweide projection view of the AC discrepancy maps with resolution Nside = 256 using the following routines for the Commander 2015 (left) and SEVEM 2018 (right) CMB temparature maps from [Planck Legacy Archive](https://pla.esac.esa.int/#maps).

<img src="https://github.com/wangyg19/CMBProbe/blob/master/ACD_Commander2015_Nside1024_notitle.png" alt="acd_commander15" width="430"><img src="https://github.com/wangyg19/CMBProbe/blob/master/ACD_SEVEM2018_Nside1024_notitle.png" alt="acd_sevem18" width="430">

For example, run
```
python coeff_cmb.py
```
to generate Fourier coefficients from CMB map "COM_CMB_IQU-commander_2048_R3.01_full.fits", which can be downloaded from [Planck Legacy Archive](https://pla.esac.esa.int/#maps).

Run
```
python ACDcommander15.py
```
to compute AC discrepancies at HEALPix points at resolutioin Nside=256 (with about 786,432 points).

Then run
```
python ACDplt.py
```
to generate the AC discrepancy map for Commander 2015.

All the generated data are in .mat format and stored in the same folder where the program scripts are in. For initial test, try Nside for AC discrepancy map equal to 8 and degree L = 20.


## Acknowledgements
Some of the results in this paper have been derived using the [HEALPix package](https://healpix.sourceforge.io/) [Gorski et al. (2005)](https://arxiv.org/abs/astro-ph/0409513). The authors acknowledge support from the Australian Research Council under Discovery Project DP180100506.

## Citation 
If you use our codes and datasets, please cite:
```
@article{CMBProbe,
  title={A new Probe of non-Gaussianity for CMB Maps},
  author={Sloan, Ian H. and Le Gia, Quoc T. and and Wang, Yu Guang and Womersley, Robert S. and Hamann, Jan},
  journal={arXiv preprint arXiv:},
  year={2019}
}
```
## Notes
If you find any problem in programs, please contact Yu Guang Wang (yuguang.wang@unsw.edu.au).

The package **CMBProbe** may be used for any research purposes under the following conditions:
* The user must acknowledge the use of **CMBProbe** in publications resulting from the use of the functions/tools.
* The user may not redistribute **CMBProbe** without separate permission.
* The user may not use this information for any commercial or revenue-bearing purposes without first obtaining permission from us.
