# SNEmcee

**Dependencies:** <br>
Python libraries: <br>
corner <br>
emcee <br>
pandas <br>
numpy <br>
scipy <br>
matplotlib <br>

The required Python libraries can be installed using:
```bash
pip install corner emcee "pandas<2" "numpy<2" scipy matplotlib
```

**Input supernova data:** <br>
Supernova data files (observations) should be stored in the directory ‘SN_data’, with the following names: <br>
- `<SN_name>_veloc` - Fe II 5169 velocity measurements (km/s).
- `<SN_name>_mag` - Multi-band photometric measurements.
- `<SN_name>_lum` - Bolometric luminosity measurements (erg/s).

See the example files for SN2017eaw in `SN_data` for the expected file formats.  <br>

In addition, the file `distances.csv` should be stored in `SN_data` and contain the following columns: <br>
- `SN_name` - Supernova name.
- `distance` - Distance in Mpc.
- `distance_err` - Distance uncertainty in Mpc.


**SNEC model grid:** <br>
SNEmcee uses a pre-computed grid of 2,469 hydrodynamical explosion models calculated with SNEC, using KEPLER red supergiant progenitors with and without circumstellar material (CSM). <br>
The grid should be downloaded from: <br>
https://drive.google.com/file/d/1-7ki0ADkq2_Y8ffyzWMT42_VuvXQK69p/view?usp=sharing
<br>
And then added to the main directory and extracted using
`tar -xf SNEC_models.tar.gz` <br>
The result should be the directory `SNEmcee/SNEC_models` which contains 2,469 sub-directories with names like 'M9.0_Ni0.001_E0.1_Mix2.0_R1000_K10' <br>
Each model subdirectory contains the predicted bolometric light curve, multi-band magnitudes, Fe II 5169 velocities, and effective temperature evolution. <br>

**Fitting parameters:** <br>
SNEmcee fits eight parameters: six parameters defining the SNEC model grid and two additional fitting parameters (S and t). <br>

| Parameter | Description | Grid values | Units |
| --- | --- | --- | --- |
| Mzams | Progenitor ZAMS mass | 9, 10, 11, 13, 15, 17 | Solar masses |
| Ni | Nickel-56 mass | 0.001, 0.02, 0.07, 0.12 | Solar masses |
| E | Explosion energy | 0.1, 0.3, 0.5, 0.9, 1.3, 1.7 | 10^51 erg |
| R | CSM thickness | 0, 500, 1000, 2000 | Solar radii |
| K | CSM density parameter | 0, 10, 30, 60 | 10^17 g/cm |
| Mix | Nickel-56 mixing | 2, 8 | Solar masses |
| S | Luminosity scaling factor | - | Dimensionless |
| t | Explosion time shift | - | Days |

The parameter R in the code represents the CSM thickness (the outer CSM radius minus the progenitor radius). <br>
The scaling factor S accounts for uncertainties in the distance estimate, and t accounts for uncertainties in the assumed explosion time. <br>
SNEmcee interpolates between the discrete grid points to continuously explore the parameter space. <br>
For fits without CSM, R and K are fixed to zero. <br>

**Running SNEmcee:** <br>
To run SNEmcee on your SN of choice, use the following command: <br>

`python3 run_mcmc_snec_for_all_SNe.py -S <SN name> -t <fitting_type> -s <number of steps> -b <burn-in steps> -w <number of walkers> -o <output directory> -c <csm> -n <normalization> -l <luminosity threshold> -p <nonuniform priors>` <br>

- `-S` SN name [required. for example: SN2017eaw]
- `-t` fitting_type = lum, veloc, mag, lum-veloc, mag-veloc, lum-mag, lum-veloc-mag [default: lum]
- `-s` number of steps = <int> [default: 500]
- `-b` burn-in steps = <int> [default: 300]
- `-w` number of walkers = <int> [default: 30]
- `-o` output directory name [default: `output_<current time>`]
- `-c` csm = with, without, twostep, twostep-carryover [default: with]
- `-n` normalization = True, False [default: False]
- `-l` luminosity_threshold = True, False [default: False]
- `-p` nonuniform_priors: Use `-p True` to enable a Gaussian prior on S. Omit this argument to disable the optional Gaussian prior [default: None].

When normalization is enabled, each log-likelihood contribution is divided by the number of data points in the corresponding dataset. <br>
The temperature threshold (10^3.5 K) is always applied to velocity and magnitude fits, and is optionally applied to bolometric luminosity fits using the `-l` argument. <br>

**Example:** <br>
To fit the bolometric light curve and Fe II 5169 velocities of SN2017eaw with CSM, using 1,000 steps, 400 burn-in steps and 50 walkers, run: <br>

`python3 run_mcmc_snec_for_all_SNe.py -S SN2017eaw -t lum-veloc -s 1000 -b 400 -w 50 -o example -c with`

**CSM fitting options:** <br>
SNEmcee supports four CSM fitting modes, specified using the `-c` argument: <br>

- `with`: All eight parameters are fitted simultaneously, including the CSM parameters R and K.
- `without`: Only six parameters are fitted, with R and K fixed to zero.
- `twostep`: The fitting is performed in two stages. In the first stage, the models are fitted without CSM. In the second stage, CSM is included, using the final walker positions from the first stage as initial guesses for the six non-CSM parameters. The CSM parameters are initialized randomly within their allowed ranges.
- `twostep-carryover`: Similar to `twostep`, but the posterior distributions of the six non-CSM parameters obtained after burn-in in the first stage are approximated by degree-12 polynomial functions and used as priors in the second stage.

For both two-stage fitting modes, the final sampler contains only the results of the second stage. <br>
The specified number of MCMC steps is applied to each stage separately. <br>

**SNEmcee outputs:** <br>
Output files will be stored in `mcmc_results/<output directory>/<number of steps>step/<SN_name>_<fitting_type>_csm-<csm>_normalized<normalization>_TreshLum<luminosity_threshold>/`. <br>

Output files will include: <br>

flat_sampler.csv: <br>
MCMC sampler output table - a 2D table, where columns are the parameters, and rows are walkers and steps (including burn-in steps). <br>
The columns are ordered as follows: Mzams, Ni, E, R, K, Mix, S, t, without column headers. <br>
For fits without CSM, R and K are set to zero. <br>
For two-stage fits, only the samples from the second stage are saved. <br>

final_results.csv: <br>
Summary statistics for each parameter calculated from all walkers at the final MCMC step: mean and lower and upper uncertainties derived from the 16th and 84th percentiles. <br>

run_parameters.csv: <br>
Record of the input arguments used in that run. <br>

log_likelihoods.csv: <br>
Log-posterior values from the MCMC sampler, including burn-in steps and excluding non-finite values. For two-stage fits, only the values from the second stage are saved. <br>

For `twostep-carryover` fits, figures showing the first-stage parameter distributions and their polynomial approximations are also saved. <br>

**Running script for comparison plots:** <br>
After running run_mcmc_snec_for_all_SNe.py, you can use plot_snec_fits_wrapper.py to plot paper-ready figures - for example, to compare the effect of different likelihood functions (different fitting schemes), or the inclusion of CSM, on the posterior distributions: <br>

`python3 plot_snec_fits_wrapper.py -S <SN name> -s <number of steps> -f <figure_type> -o <output directory> -t <fitting type> -c <csm> -n <normalization> -l <luminosity threshold>` <br>


-S SN name [required. for example: SN2017eaw]<br>
-s number of steps = `<int>` [default: 500]<br>
-f figure type = lum, veloc, mag (or any combination of those, separated by -), corner, walkers, or one of the comparison plots listed below [default: lum] <br>
-o output directory name [default: `output_<current time>`] <br>
-t fitting_type = lum, veloc, mag, lum-veloc, mag-veloc, lum-mag, lum-veloc-mag [default: lum] <br>
-c csm = with, without, twostep, twostep-carryover [default: without]<br>
-n normalization = True, False [default: False]<br>
-l luminosity_threshold = True, False [default: False]<br>

**Figure types:** <br>
- `lum`: Bolometric light curve fit.
- `veloc`: Fe II 5169 velocity fit.
- `mag`: Multi-band magnitude fits.
- `corner`: Corner plot of the posterior distributions.
- `walkers`: Chain plots for each parameter.

**Comparison plots:** <br>
- `csm_comparison`
- `lum-mag-veloc_comparison`
- `lum-mag-veloc-Tthresh_comparison`
- `lum-veloc-normalized_comparison`
- `lum-veloc-twostep_comparison`
- `martinez_comparison`

The fitting settings should match those used in the original MCMC run. <br>
The required MCMC runs must be completed before generating comparison plots. <br>
The `martinez_comparison` option also requires additional model files from Martinez et al. (2020). <br>
The `martinez_comparison` option expects MCMC results generated with normalization and the luminosity threshold disabled. <br>

**Example:** <br>
`python3 plot_snec_fits_wrapper.py -S SN2017eaw -s 1000 -f lum-veloc -o example -t lum-veloc -c with` <br>


**Plotting outputs:** <br>
Output files will be stored in `figures/<output directory>/<number of steps>step/`. <br>
The Martinez comparison plots are stored in `figures/<output directory>/`. <br>

Output files will include summary figures, depending on the selected figure type: <br>
Fits of models to the input observation data for that SN. <br>
Comparison of model fits for each run-type compared, including the input observation data for that SN. <br>
Corner plot of posterior distributions for each parameter. <br>
Comparison corner plot with the posterior distributions for each run type overlaid on one another. <br>
Chain plots for each parameter, with the burn-in region indicated. <br>

**Citation:** <br>
If you use SNEmcee in your research, please cite our accompanying paper: <br>

Mohsen-Tanev, S., Shitrit, N., Bracha, S., & Arcavi, I. <br>
SNEmcee: MCMC Fitting of Numerical Progenitor and Explosion Models to Type IIP Supernova Light Curves and Velocities. <br>
