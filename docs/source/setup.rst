Setting up your runs
=====

Currently, there are two main analyses that Astroclimes can carry out: the first is to fit telluric model spectra to a sample of spectroscopic observations to find the best-fit values for the abundances of certain molecular species; the second is to use these best-fit telluric models to remove the telluric lines from this same sample of spectroscopic observations. There are two separate setup scripts that run each of the analyses, ``setup.py`` and ``setup_telrem.py``, respectively.

Run to find best-fit molecular abundance values
-----
Before running ``setup.py``, there are a number of variables inside it that you need to specify, which are listed and explained below:


``home_directory``: path to parent directory where files and subdirectories will be saved

``MCMC_directory``: path to subdirectory where results and plots from the MCMC will be saved

``filename_mcmc_results``: file name (with path) where the MCMC results will be written

``atm_profs_directory``: path to subdirectory where atmospheric profiles are stored

``filename_em_line_spec``: file name (with path) of the model emission line spectra

``obs_directory``: path to subdirectory where observations are stored

``list_science_spectra``: list containing the file names of the spectroscopic observations to be analysed

``instrument``: string containing the instrument name. Current instruments supported are: CARMENES (NIRPS and ESPRESSO in prep.). For a guide on how to tailor Astroclimes to your specific instrument, check [INCLUDE REFERENCE HERE]

``R_instrument``: resolution of the selected instrument

``n_orders``: number of spectral orders the instrument data is divided into

``filename_site_values``: file name where the information about the observations, the site and the weather conditions will be written

``include_stelmod``: boolean variable specifying whether or not to include a stellar model in the analysis

``filename_stelmod_lam``: name of the file (with path) containing the wavelength distribution of the stellar model

``filename_stelmod_spec``: name of the file (with path) containing the flux of the stellar model

``molecs``: list containing names of the molecules to be included in the modelling via line-by-line absorption

``molecs_for_cia``: list containing the reaction pairs to be included in the modellign via collision-induced absorption (CIA)

``stelpars``: Python object to store the stellar parameters. Here, the parameters needed are the star's projected rotational velocity 
:math:`vsini` and a linear limb-darkening coefficient :math:`\varepsilon`. The effective temperature :math:`T_{eff}`, surface gravity :math:`\log{g}` and metallicity `FeH` are not used in any of the calculations but are defined here for ease of access because they are needed to determine the best stellar model and limb darkening coefficient to use

``orbpars``: Python object to store the stellar parameters. Here, the parameters needed are the systemic velocity :math:`V_{sys}`, the stellar RV semi-amplitude :math:`K_s`, the planet's orbital period `P` and mid-transit time `T0`.

``scale_profs``: boolean variable specifying whether or not to scale the atmospheric profiles to match the values measured by the observatory's weather station (as reported in the observation's FITS header)

``vel_step``: velocity step to convert the telluric and stellar model to a constant :math:`\Delta`:math:`\lambda`/:math:`\lambda`

``DMF_O2``: atmospheric dry air mole fraction of O2

``free_molecs``: list of molecules whose abundances will be free parameters for the MCMC (these have to be included in ``molecs``)

``spec_orders``: spectral orders to be included in the analysis

Once all of the parameters are specified, you can run the file by doing:

>>> python setup.py

This will run the script `main.py` on the selected spectroscopic observations to find the best-fit values for the abundances of the selected molecular species. To learn about the step-by-step of how `main.py` works, check [INCLUDE REFERENCE HERE].
