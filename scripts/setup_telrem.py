# =====================================================================================
# Basic packages
# =====================================================================================
import glob
import os
from subprocess import call

# =====================================================================================
# Scripts
# =====================================================================================
import objects
import pca_functions
import funcs
import telrem

'''
	Setup script for removing telluric and stellar lines and doing injection and recovery tests with Astroclimes and PCA
'''

for nc_PCA in [1,2,3,4,5,6,10]:
	## Some constants
	R_Sun = 6.957*1e8	# Solar radius, in m (from NASA Sun Fact Sheet)
	R_Jup = 6.9911*1e7	# Jupiter radius, in m (from NASA Jupiter Fact Sheet)

	## Specify number of CPUs to use
	n_CPUs = 1

	## Define your home directory
	#home_directory = '/Users/marceloaronkeniger/PhD/thesis_work/Astroclimes/'
	#home_directory = '/home/marceloaron/MarceloAron/PhD/thesis_work/Astroclimes/'
	home_directory = '/media/marceloaron/New Volume/PhD/thesis_work/'
	#home_directory = '/storage/astro2/phrgmq/Astroclimes/'

	## Directory path for the "main analysis" MCMC (the one to determine the best-fit molecular abundances)
	MCMC_directory = home_directory+'MCMC_results/tauboo/'

	## Define the file name where the MCMC results are stored (by default, this should be inside the Astroclimes/MCMC_results/ directory)
	filename_mcmc_results = MCMC_directory+'results_mcmc.txt'
	
	## Directory path for the atmospheric profiles
	atm_profs_directory = home_directory+'atmosphere_profiles/GGG2020/fp/al/'

	## Define the directory path to where observations are stored
	#obs_directory = home_directory+'data/CARMENES/'
	obs_directory = '/home/marceloaron/MarceloAron/PhD/thesis_work/Astroclimes/data/CARMENES/'

	## Creating list of the observation file names to be analysed
	## Here we are specifying that we only want the NIR files from fibre A, which is a syntax specific to CARMENES, 
	## other instruments might have different file name conventions
	list_science_spectra = glob.glob(obs_directory+'*nir_A.fits')
	list_science_spectra.sort()

	## Pick your instrument and fill out the relevant information
	instrument = 'CARMENES'
	R_instrument = 80400	# CARMENES NIR resolution
	n_orders = 28 			# Number of CARMENES NIR spectral orders
	n_pixels = 4080			# Number of pixels (number of wavelength points) per order

	## Define directory where the telluric removal results will be stored
	main_telrem_directory = home_directory+'telluric_removal/TauBoo/2018_03_26/'
	if not os.path.isdir(main_telrem_directory):
		call(['mkdir', main_telrem_directory])

	## Define path to site values file (if it doesn't exist, it will be created inside the specified directory)
	filename_site_values = main_telrem_directory+'site_values.txt'
	if not os.path.isfile(filename_site_values):
		funcs.create_site_values_file(list_science_spectra, dirname=main_telrem_directory)

	## List of desired molecules to include in the modelling
	molecs = ['CO2', 'CH4', 'H2O', 'O2', 'N2'] 		# Molecules included via line-by-line absorption
	molecs_for_cia = ['O2-air'] 					# Molecules included via collision-induced absorption (CIA)

	## Choose bad spectral orders to be removed from the analysis
	bad_spec_orders_idxs = [8, 9, 10, 17, 18, 19, 20, 21, 22] 		## Orders that are too densely populated with saturated telluric water lines

	## We will now run a function that determines "bad" observations based on their measured SNR (from the FITS header) and airmass
	SNR_lim = 100
	airmass_lim = 1.7
	## Manually add any extra bad observations that may not be included in the automatic SNR and airmass flagging
	add_bad_obs = [13,32,82] 		# These are the indexes of the spectra from the original sample, before removing any observations
	create_SNR_night_log = False 	# Optional: creating a file with the exposure times and measured spectra SNR
	make_airmass_SNR_plot = True 	# Optional: plotting the airmass and SNR over time
	bad_obs_idxs = telrem.filter_bad_obs(instrument=instrument, 
										 list_science_spectra=list_science_spectra, 
										 n_orders=n_orders, 
										 bad_spec_orders_idxs=bad_spec_orders_idxs, 
										 filename_site_values=filename_site_values, 
										 SNR_lim=SNR_lim, 
										 airmass_lim=airmass_lim, 
										 add_bad_obs=add_bad_obs, 
										 main_telrem_directory=main_telrem_directory, 
										 create_SNR_night_log=create_SNR_night_log,
										 make_airmass_SNR_plot=make_airmass_SNR_plot)

	## Range of planet signal strength multipliers for the injected planetary signal
	scale_factors = [1,2,3,3.25,3.50,3.75,4,4.25,4.50,4.75,5,6,7,8,9]

	## Number of PCA components
	#nc_PCA = 5

	## Filename to store the SNRs of each run (this is not the same SNR as the one used for the bad data clip above!)
	filename_telrem_SNR_AC = 'SNR_logs_ac.txt'
	filename_telrem_SNR_PCA = f'SNR_logs_pca_nc_{nc_PCA}.txt'

	## Checking if the SNR log files exist and if not, creating them and adding the header
	if not os.path.isfile(main_telrem_directory+filename_telrem_SNR_AC):
		txt = open(main_telrem_directory+filename_telrem_SNR_AC, 'a')

		txt.write('## Orders removed: ')
		for i in range(len(bad_spec_orders_idxs)):
			txt.write(f'{bad_spec_orders_idxs[i]+1}, ')

		txt.write('\n')
		txt.write('## Spectra removed: ')
		for i in range(len(bad_obs_idxs)):
			txt.write(f'{bad_obs_idxs[i]+1}, ')

		txt.write('\n')
		txt.write("## Injected planet strength \t SNR_inj_pos \t vRest_inj_pos (km/s) \t kpVec_inj_pos (km/s) \t SNR_max \t vRest_max (km/s) \t kpVec_max (km/s)\n")
		txt.close()

	if not os.path.isfile(main_telrem_directory+filename_telrem_SNR_PCA):
		txt = open(main_telrem_directory+filename_telrem_SNR_PCA, 'a')

		txt.write('## Orders removed: ')
		for i in range(len(bad_spec_orders_idxs)):

			txt.write(f'{bad_spec_orders_idxs[i]+1}, ')
		txt.write('\n')
		txt.write('## Spectra removed: ')
		for i in range(len(bad_obs_idxs)):

			txt.write(f'{bad_obs_idxs[i]+1}, ')
		txt.write('\n')
		txt.write("## Injected planet strength \t nc_PCA \t SNR_inj_pos \t vRest_inj_pos (km/s) \t kpVec_inj_pos (km/s) \t SNR_max \t vRest_max (km/s) \t kpVec_max (km/s)\n")
		txt.close()

	## Filename to store the MCMC results of each run (this is not the same MCMC as the results_mcmc.txt file from above!)
	filename_telrem_MCMC_results_AC = 'results_telrem_mcmc_ac.txt'
	filename_telrem_MCMC_results_PCA = f'results_telrem_mcmc_pca_nc_{nc_PCA}.txt'

	## Check if the MCMC results files exist and if not, creating them and adding the header
	if not os.path.isfile(main_telrem_directory+filename_telrem_MCMC_results_AC):
		txt = open(main_telrem_directory+filename_telrem_MCMC_results_AC, 'a')
		txt.write("## Injected planet strength \t sf \t sf_error \t log_sf \t log_sf_error \t Kp (km/s) \t Kp_error (km/s) \t Vsys (km/s) \t Vsys_error (km/s) \t logp \t logp_error\n")
		txt.close()

	if not os.path.isfile(main_telrem_directory+filename_telrem_MCMC_results_PCA):
		txt = open(main_telrem_directory+filename_telrem_MCMC_results_PCA, 'a')
		txt.write("## Injected planet strength \t nc_PCA \t sf \t sf_error \t log_sf \t log_sf_error \t Kp (km/s) \t Kp_error (km/s) \t Vsys (km/s) \t Vsys_error (km/s) \t logp \t logp_error\n")
		txt.close()

	## Checking if the data cubes containing the wavelength, observational spectra and model telluric spectra already exist, and if not, creating them

	## This creates both the wavelength and observational spectra cube, if not already existent (this only needs to be done once)
	if not os.path.isfile(main_telrem_directory+'all_nights_lam.npy') or not os.path.isfile(main_telrem_directory+'all_nights_spec.npy'):
		all_nights_lam, all_nights_spec = pca_functions.get_obs_cube(n_orders, n_pixels, list_science_spectra, generate=True, save_dirname=main_telrem_directory)

	## This creates the model spectra cube, if not already existent (this only needs to be done once)
	if not os.path.isfile(main_telrem_directory+'all_nights_spec_mod.npy'):
		all_nights_spec_mod = pca_functions.get_mod_cube(n_orders, n_pixels, list_science_spectra, MCMC_directory, atm_profs_directory, filename_mcmc_results, molecs, molecs_for_cia, generate=True, save_dirname=main_telrem_directory)

	## Specify filenames where the planet model wavelength and flux are kept (models here come from GENESIS code from Gandhi & Madhusudhan 2017)
	filename_planet_lam = home_directory+'telluric_removal/planet_models/lam.txt'
	filename_planet_flux = home_directory+'telluric_removal/planet_models/h2o_-3.0_ch4_-20.0_hcn_-20.0.txt'

	## Define stellar parameters through the StellarParams() object
	stelpars = objects.StellarParams()
	stelpars.Teff = 6465 	# Effective temperature, in K, from Soubiran et al. (2022)
	stelpars.R = 1.42*R_Sun	# Stellar radius in m, from Borsa et al. (2015), also used in Webb et al. (2022), consistent with Valenti & Fisher (2005), who report R_star = 1.419 R_Sun

	## Defining the orbital parameters through the OrbitParams() object
	orbpars = objects.OrbitParams()
	orbpars.Vsys = -16.4 		# Systemic velocity, in km/s, from Brogi et al. (2012)
	orbpars.Kp = 110 			# Planet RV semi-amplitude in km/s, from Brogi et al. (2012), also consistent with Lockwood et al. (2014; 111 km/s)
	orbpars.P = 3.312433		# Orbital period, in days, from Brogi et al. (2012)
	orbpars.T0 = 2455652.108	# Planet mid-transit time, in BJD, from Brogi et al
	orbpars.Rp = 1.2*R_Jup	# Planet radius in m, from Webb et al. (2022)

	## Decide if you want to make certain optional plots
	plot_telrem_steps_AC = True
	plot_telrem_steps_PCA = True
	plot_indiv_PCA_comps = True

	## Decide if you want to run the MCMC for AC or PCA (this MCMC is to determine the erosion of the injected signal)
	run_MCMC_AC = False
	run_MCMC_PCA = False

	telrem.run_telrem(home_directory=home_directory, 
					  MCMC_directory=MCMC_directory,
					  atm_profs_directory=atm_profs_directory, 
				  	  filename_mcmc_results=filename_mcmc_results, 
				  	  main_telrem_directory=main_telrem_directory, 
				  	  list_science_spectra=list_science_spectra, 
				  	  filename_site_values=filename_site_values, 
				  	  molecs=molecs, 
				  	  molecs_for_cia=molecs_for_cia, 
				  	  R_instrument=R_instrument,
				  	  n_orders=n_orders,
				  	  n_pixels=n_pixels, 
				  	  bad_spec_orders_idxs=bad_spec_orders_idxs, 
				  	  bad_obs_idxs=bad_obs_idxs, 
				  	  add_bad_obs=add_bad_obs, 
				  	  scale_factors=scale_factors, 
				  	  nc_PCA=nc_PCA, 
				  	  filename_telrem_SNR_AC=filename_telrem_SNR_AC, 
				  	  filename_telrem_SNR_PCA=filename_telrem_SNR_PCA, 
				  	  filename_telrem_MCMC_results_AC=filename_telrem_MCMC_results_AC, 
				  	  filename_telrem_MCMC_results_PCA=filename_telrem_MCMC_results_PCA,  
				  	  filename_planet_lam=filename_planet_lam, 
				  	  filename_planet_flux=filename_planet_flux, 
				  	  stelpars=stelpars, 
				  	  orbpars=orbpars, 
				  	  plot_telrem_steps_AC=plot_telrem_steps_AC, 
				  	  plot_telrem_steps_PCA=plot_telrem_steps_PCA, 
				  	  plot_indiv_PCA_comps=plot_indiv_PCA_comps, 
				  	  run_MCMC_AC=run_MCMC_AC, 
				  	  run_MCMC_PCA=run_MCMC_PCA,
				  	  n_CPUs=n_CPUs)