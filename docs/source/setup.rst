Setting up your runs
=====

Currently, there are two main analyses that Astroclimes can carry out: the first is to fit telluric model spectra to a sample of spectroscopic observations to find the best-fit values for the abundances of certain molecular species; the second is to use these best-fit telluric models to remove the telluric lines from this same sample of spectroscopic observations. There are two separate setup scripts that run each of the analyses, ``setup.py`` and ``setup_telrem.py``, respectively.

Run to find best-fit molecular abundance values
-----
``setup.py`` will run ``main.py``, which will loop over the list of spectroscopic observations and run an MCMC to find the best-fit abundance values for the selected molecules. To learn about the step-by-step of how `main.py` works, check [INCLUDE REFERENCE HERE]. Before running ``setup.py``, there are a number of variables inside it that you need to specify, which are listed and explained below:

``n_CPUs``: number of CPUs to use (choose more than 1 if you wish to use multiprocessing in the MCMC).

``home_directory``: path to parent directory where files and subdirectories will be saved.

``MCMC_directory``: path to subdirectory where results and plots from the MCMC will be saved.

``filename_mcmc_results``: file name (with path) where the MCMC results will be written.

``atm_profs_directory``: path to subdirectory where atmospheric profiles are stored.

``filename_em_line_spec``: file name (with path) of the model emission line spectra.

``obs_directory``: path to subdirectory where observations are stored.

``list_science_spectra``: list containing the file names of the spectroscopic observations to be analysed.

``instrument``: string containing the instrument name. Current instruments supported are: CARMENES (NIRPS and ESPRESSO in prep.). For a guide on how to tailor Astroclimes to your specific instrument, check [INCLUDE REFERENCE HERE].

``R_instrument``: resolution of the selected instrument.

``n_orders``: number of spectral orders the instrument data is divided into.

``filename_site_values``: file name (with path) where the information about the observations, the observing site and the weather conditions will be written. As a default, it will be saved in ``MCMC_directory``.

``include_stelmod``: boolean variable specifying whether or not to include a stellar model in the analysis.

``filename_stelmod_lam``: name of the file (with path) containing the wavelength distribution of the stellar model.

``filename_stelmod_spec``: name of the file (with path) containing the flux of the stellar model.

``molecs``: list containing names of the molecules to be included in the modelling via line-by-line absorption.

``molecs_for_cia``: list containing the reaction pairs to be included in the modellign via collision-induced absorption (CIA).

``stelpars``: Python object to store the stellar parameters. Here, the parameters needed are the star's projected rotational velocity 
:math:`vsini` and a linear limb-darkening coefficient :math:`\varepsilon`. The effective temperature :math:`T_\text{eff}`, surface gravity :math:`\log{g}` and metallicity FeH are not used in any of the calculations but are defined here for ease of access because they are needed to determine the best stellar model and limb darkening coefficient to use.

``orbpars``: Python object to store the orbital parameters. Here, the parameters needed are the systemic velocity :math:`V_\text{sys}`, the stellar RV semi-amplitude :math:`K_s`, the planet's orbital period `P` and mid-transit time `T_0`.

``scale_profs``: boolean variable specifying whether or not to scale the atmospheric profiles to match the values measured by the observatory's weather station (as reported in the observation's FITS header).

``vel_step``: velocity step to convert the telluric and stellar model to a constant :math:`\Delta`:math:`\lambda`/:math:`\lambda`.

``DMF_O2``: atmospheric dry air mole fraction of O2.

``free_molecs``: list of molecules whose abundances will be free parameters for the MCMC (these have to be included in ``molecs``).

``spec_orders``: spectral orders to be included in the analysis.

Once all of the parameters are specified, you can run the script by doing:

>>> python setup.py

Run to remove telluric lines
-----
In its current version, ``setup_telrem.py`` will do more than just remove telluric and stellar lines from the selected spectra. It will first run an automatic check for bad observations based on their measured SNR and airmass (whose limits can be customly specified) and then carry out several injection and retrieval tests for Astroclimes and Principal Component Analysis (PCA; as described in `Giacobbe et al. (2021) <https://www.nature.com/articles/s41586-021-03381-x>`_), where the level of erosion of the planetary signal is determined with an MCMC. IMPORTANT: there are two different MCMCs referenced here. The first MCMC, which is run by ``setup.py``, will be referred to as the "molecular abundance MCMC", whereas the second MCMC, run by ``setup_telrem.py``, will be referred to as the "telluric removal MCMC". Just as before, there are a number of variables inside ``setup_telrem.py`` that you need to specify, which are listed and explained below. Also note that ``setup.py`` must be run beforehnd, as its data products are used by ``setup_telrem.py``.

``n_CPUs``: number of CPUs to use (choose more than 1 if you wish to use multiprocessing in the telluric removal MCMC).

``home_directory``: path to parent directory where files and subdirectories will be saved.

``MCMC_directory``: path to subdirectory where results and plots from the molecular abundance MCMC were saved.

``filename_mcmc_results``: file name (with path) where the molecular abundance MCMC results were be written.

``atm_profs_directory``: path to subdirectory where atmospheric profiles are stored.

``obs_directory``: path to subdirectory where observations are stored.

``list_science_spectra``: list containing the file names of the spectroscopic observations to be analysed.

``instrument``: string containing the instrument name. Current instruments supported are: CARMENES (NIRPS and ESPRESSO in prep.). For a guide on how to tailor Astroclimes to your specific instrument, check [INCLUDE REFERENCE HERE].

``R_instrument``: resolution of the selected instrument.

``n_orders``: number of spectral orders the instrument data is divided into.

``n_pixels``: number of pixels each spectral order of the instrument has.

``main_telrem_directory``: path to subdirectory where results from this telluric removal analysis will be sored.

``filename_site_values``: file name (with path) where the information about the observations, the observing site and the weather conditions will be written. As a default, it will be saved in ``main_telrem_directory``. This is exactly the same file as before, so if one so chooses, there is no need to create another copy, just specify the path to the previously created one.

``molecs``: list containing names of the molecules to be included in the modelling via line-by-line absorption.

``molecs_for_cia``: list containing the reaction pairs to be included in the modellign via collision-induced absorption (CIA).

``bad_spec_orders_idxs``: indexes identifying the bad spectral orders to be excluded from the analysis (due to densely populated saturated telluric water lines or any other reason).

``SNR_lim``: measured SNR limit below which observations will be discarded.

``airmass_lim``: meaasured airmass limit above which observations will be discarded.

``add_bad_obs``: any additional observations that the user wants to remove that may not have been included in the automatic SNR and airmass screening.

``create_SNR_night_log``: boolean variable specifying whether or not to create a text file with the measured SNR for each observation of that night.

``make_airmass_SNR_plot``: boolean variable specifying whether or not to plot the measured SNR and airmass distribution as a function of time.

``scale_factors``: list of scale factors which will be used to multiply the model planetary signal.

``nc_PCA``: number of PCA components.

``filename_telrem_SNR_AC``: name of the file where the calculated Astroclimes SNRs will be written. These SNRs are not the measured SNRs from the observational spectra, but rather the SNRs of the detected planetary signal for each injected planetary signal strength. As a default, it will be saved in ``main_telrem_directory``.

``filename_telrem_SNR_PCA``: name of the file where the calculated PCA SNRs will be written. These SNRs are not the measured SNRs from the observational spectra, but rather the SNRs of the detected planetary signal for each injected planetary signal strength. As a default, it will be saved in ``main_telrem_directory``.

``filename_telrem_MCMC_results_AC``: name of the file where the telluric removal MCMC results for Astroclimes will be written. As a default, it will be saved in ``main_telrem_directory``.

``filename_telrem_MCMC_results_PCA``: name of the file where the telluric removal MCMC results for PCA will be written. As a default, it will be saved in ``main_telrem_directory``.

``filename_planet_lam``: name of the file (with path) containing the wavelength distribution of the planet model.

``filename_planet_spec``: name of the file (with path) containing the flux of the planet model.

``stelpars``: Python object to store the stellar parameters. Here, the parameters needed are the star's effective temperature :math:`T_\text{eff}` and radius `R`.

``orbpars``: Python object to store the orbital parameters. Here, the parameters needed are the systemic velocity :math:`V_\text{sys}`, the planet RV semi-amplitude :math:`K_p`, the planet radius :math:`R_p`, the planet's orbital period `P` and mid-transit time `T_0`.

``plot_telrem_steps_AC``: boolean variable specifying whether or not to plot the step-by-step of the telluric and stellar line removal with Astroclimes. As a default, this plot (and the two plots below) will be made for spectral order 12 (in the CARMENES spectral order distribution, this is an order which contains several typically non-saturated water lines).

``plot_telrem_steps_PCA``: boolean variable specifying whether or not to plot the step-by-step of the telluric and stellar line removal with PCA.

``plot_indiv_PCA_comps``: boolean variable specifying whether or not to plot the contributions of each PCA component. As a default, this will plot the first 6 PCA components, regardless of how many are included in the analysis. 

``run_MCMC_AC``: boolean variable specifying whether or not to run the telluric removal MCMC for Astroclimes.

``run_MCMC_PCA``: boolean variable specifying whether or not to run the telluric removal MCMC for PCA.

Once all of the parameters are specified, you can run the script by doing:

>>> python setup_telrem.py

SLURM runs
-----
If you have the ability of running code in a server that supports a batch system like SLURM, you can speed up the analysis process by running the SLURM-tailored versions of ``setup.py`` and ``setup_telrem.py``. They are basically the same, with just additional/alternative command line arguments, and have to be run through the bash scripts ``slurm_job_setup.sh`` and ``slurm_job_setup_telrem.sh``, respectively, where you can specify things such as how may CPUs you want per task and how much memory to allocate per CPU. You may run the former as shown below:

>>> sbatch --array=0-99%10 slurm_job_setup.sh

The ``--array`` option will run ``slurm_job_setup.sh`` 100 times, with task IDs ranging from 0 to 99, in batches of 10 (of course, you may change that to suit your case). The task ID is used to index which file from your list of observations will be analysed. That is one of the command line arguments given to ``setup_slurm.py``, the other being the number of CPUs, which must be specified in the bash script.

As for ``slurm_job_setup_telrem.sh``, you may run it as below:

>>> sbatch --array=0-15%16 slurm_job_setup_telrem.sh S{run_MCMC_AC} S{run_MCMC_PCA} S{nc_PCA}

In this case, the task IDs are used to index the scale factor among a pre-determined list of scale factors. This is the first command line argument that ``setup_telrem_slurm.py`` takes, the second one again being the number of CPUs, which is specified in the bash file. The remaining command line arguments that the Python script takes are the ones given to the bash script above. ``S{run_MCMC_AC}`` and ``S{run_MCMC_PCA}`` can be either 0 (False) or 1 (True) and ``S{nc_PCA}`` can be any integer. 
