# =====================================================================================
# Basic packages
# =====================================================================================
import glob
import numpy as np
import os
from subprocess import call

# =====================================================================================
# Scripts
# =====================================================================================
import main
import funcs
import objects

'''
	Setup script for finding the best-fit molecular abundance values with Astroclimes
'''

## Specify number of CPUs to use
n_CPUs = 6

## Define your home directory
#home_directory = '/home/marceloaron/MarceloAron/PhD/thesis_work/Astroclimes/'
home_directory = '/media/marceloaron/New Volume/PhD/thesis_work/'
#home_directory = '/storage/astro2/phrgmq/Astroclimes/'

## Define the directory where the MCMC results will be stored
MCMC_directory = home_directory+'MCMC_results/tauboo_rerun_test/'

## Checking if MCMC directory already exists, and if not, creating it
if not os.path.isdir(MCMC_directory):
	call(['mkdir', MCMC_directory])

## Define the file name (with path) where the MCMC results will be stored
filename_mcmc_results = MCMC_directory+'results_mcmc.txt'

## Checking if the file with the MCMC results already exists, and if not, creating it with the header
if not os.path.isfile(filename_mcmc_results):
	txt = open(filename_mcmc_results, 'a')
	txt.write("## Time \t Ground_CO2 (ppm) \t Ground_CO2_error (ppm) \t Ground_CH4 (ppm) \t Ground_CH4_error (ppm) \t Ground_H2O (ppm) \t Ground_H2O_error (ppm) \t Ground_O2 (ppm) \t Ground_O2_error (ppm)\t Ground_N2 (ppm) \t Ground_N2_error (ppm)\t X_CO2 (ppm) \t X_CO2_error (ppm) \t X_CH4 (ppm) \t X_CH4_error (ppm)\t X_H2O (ppm) \t X_H2O_error (ppm)\t X_O2 (ppm) \t X_O2_error (ppm)\t X_N2 (ppm) \t X_N2_error (ppm)\t logp \t logp_error\n")
	txt.close()

## Define path to directory where relevant atmospheric profile files are stored
atm_profs_directory = home_directory+'atmosphere_profiles/GGG2020/fp/al/'

## Define path for emission line spectra file
filename_em_line_spec = home_directory+'auxiliary_files/skytable.fits'

## Define the directory path to where observations are stored
#obs_directory = home_directory+'data/CARMENES/'
obs_directory = '/home/marceloaron/MarceloAron/PhD/thesis_work/Astroclimes/data/CARMENES/'

## Creating list of the observation file names to be analysed
## Here we are specifying that we only want the NIR files from fibre A, which is a syntax specific to CARMENES, 
## other instruments might have different file name conventions
list_science_spectra = glob.glob(obs_directory+'*nir_A.fits')
list_science_spectra.sort()
list_science_spectra = list_science_spectra[:1]

## Pick your instrument (must supply its resolution)
instrument = 'CARMENES'
R_instrument = 80400	# CARMENES NIR resolution
n_orders = 28 			# Number of CARMENES NIR spectral orders

## Define path to site values file (if it doesn't exist, it will be created inside the specified directory)
filename_site_values = MCMC_directory+'site_values.txt'
if not os.path.isfile(filename_site_values):
	funcs.create_site_values_file(list_science_spectra, dirname=MCMC_directory)

## Set to False if you do not wish to include a model stellar spectra in the analysis
include_stelmod = True

## Define path to files containing the model stellar spectra wavelength and flux, respectively (the models used as default here have them in separate files)
## Here, we are using PHOENIX synthetic spectra from the Göttingen Spectral Library, available at https://phoenix.astro.physik.uni-goettingen.de/
## Download the template that best suits the star in question and also update the stellar parameters below
filename_stelmod_lam = home_directory+'data/PHOENIX_SPECTRA/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'
filename_stelmod_spec = home_directory+'data/PHOENIX_SPECTRA/lte06500-4.50+0.5.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'

## List of desired molecules to include in the modelling
molecs = ['CO2', 'CH4', 'H2O', 'O2', 'N2'] 		# Molecules included via line-by-line absorption
molecs_for_cia = ['O2-air'] 					# Molecules included via collision-induced absorption (CIA)

## Define stellar parameters through the StellarParams() object
stelpars = objects.StellarParams()
stelpars.Teff = 6465 	# Effective temperature, in K, from Soubiran et al. (2022)
stelpars.logg = 4.32 	# Surface gravity, in cm/s^2, from Soubiran et al. (2022)
stelpars.FeH = 0.28 	# Metallicity, in dex, from Soubiran et al. (2022)
stelpars.vsini = 16.8 	# Projected rotational velocity, in km/s, from Luck et al. (2017)
stelpars.eps = 0.328 	# Linear LD coefficient in the J band for a star with log(g) = 4.5 cm/s^2 and Teff = 6500 K, from Claret et al. (1995)

## Defining the orbital parameters through the OrbitParams() object
orbpars = objects.OrbitParams()
orbpars.Vsys = -16.4 		# Systemic velocity, in km/s, from Brogi et al. (2012)
orbpars.Ks = 0.4664 		# Star RV semi-amplitude, in km/s, from Brogi et al. (2012)
orbpars.P = 3.312433		# Orbital period, in days, from Brogi et al. (2012)
orbpars.T0 = 2455652.108	# Planet mid-transit time, in BJD, from Brogi et al

## Set to True if you wish to scale the pressure, temperature and humidity profiles to match the site measurements
scale_profs = False

## Choose a velocity step (in km/s) to convert the wavelength distribution from the stellar model and the cross-sections grid to a constant Delta_lambda/lambda
## The resulting wavelengt distribution will have a resolution of R = c/vel_step
vel_step = 1.0

## Define value for the O2 dry air mole fraction
DMF_O2 = 0.2095		# From Laughner et al. (2023a)

## Choose the molecules whose abundances will be a free parameter in the MCMC
free_molecs = ['CO2', 'CH4', 'H2O', 'O2']

## Choose which spectral orders to include in the analysis
spec_orders_idxs = np.arange(n_orders)
bad_spec_orders_idxs = [8, 9, 10, 17, 18, 19, 20, 21, 22] 		## Orders that are too densely populated with saturated telluric water lines
spec_orders = np.delete(spec_orders_idxs, bad_spec_orders_idxs)

## Running main.py
main.run_main(atm_profs_directory=atm_profs_directory,
			  MCMC_directory=MCMC_directory,
			  filename_mcmc_results=filename_mcmc_results, 
			  filename_em_line_spec=filename_em_line_spec, 
			  list_science_spectra=list_science_spectra, 
			  include_stelmod=include_stelmod, 
			  filename_stelmod_lam=filename_stelmod_lam, 
			  filename_stelmod_spec=filename_stelmod_spec, 
			  molecs=molecs, 
			  molecs_for_cia=molecs_for_cia, 
			  free_molecs=free_molecs, 
			  spec_orders=spec_orders, 
			  instrument=instrument, 
			  R_instrument=R_instrument, 
			  stelpars=stelpars, 
			  orbpars=orbpars, 
			  scale_profs=scale_profs, 
			  vel_step=vel_step, 
			  DMF_O2=DMF_O2,
			  n_CPUs=n_CPUs)