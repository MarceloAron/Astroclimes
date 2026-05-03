# =====================================================================================
# Basic packages
# =====================================================================================
import numpy as np
import glob
import os
from subprocess import call

# =====================================================================================
# Scripts
# =====================================================================================
import objects
import funcs
import funcs_telrem
import retrv_exo_atm

'''
	Setup script for removing telluric and stellar lines Astroclimes and PCA and verifying the detection of an exoplanetary atmospheric signal
	You may choose to use data from multiple nights 
	Detection is assessed by comparison with pre-computed model or by carrying out a full retrieval analysis by sampling the parameter space (not yet implemented)
'''

## Some constants
R_Sun = 6.957*1e8	# Solar radius, in m (from NASA Sun Fact Sheet)
R_Jup = 6.9911*1e7	# Jupiter radius, in m (from NASA Jupiter Fact Sheet)

## Define your home directory
home_directory = '/home/Astroclimes/'

## Directory path for the atmospheric profiles
## If GGG2020 atmospheric profiles are being used, this path must point to the directory where
## the vertical/, vmrs-vertical/ and maps-vertical/ directories are
atm_profs_directory = home_directory+'atmosphere_profiles/GGG2020/fp/al/'

## Define the directory path to where observations are stored
obs_directory = home_directory+'data/CARMENES/nir/tauboo/'

## Defining the nights, format is YYYY_MM_DD, and night refers to the observing night 
## (e.g. if the observations were taken at 20180327T01h00m00s, the corresponding night would be 2018_03_26)
nights = ['2018_03_26', '2018_05_11', '2019_03_12', '2019_03_15', '2019_04_11']

## Compiling the necessary files for each observing night to be combined in the analysis

## Directory path for the "main analysis" MCMC (the one to determine the best-fit molecular abundances)
MCMC_directories = [home_directory+f'MCMC_results/tauboo/{nights[0]}/',
					home_directory+f'MCMC_results/tauboo/{nights[1]}/',
					home_directory+f'MCMC_results/tauboo/{nights[2]}/',
					home_directory+f'MCMC_results/tauboo/{nights[3]}/',
					home_directory+f'MCMC_results/tauboo/{nights[4]}/']

## Define the file name where the MCMC results are stored (by default, this should be inside the Astroclimes/MCMC_results/ directory)
filenames_mcmc_results = [	MCMC_directories[0]+'results_mcmc.txt',
							MCMC_directories[1]+'results_mcmc.txt',
							MCMC_directories[2]+'results_mcmc.txt',
							MCMC_directories[3]+'results_mcmc.txt',
							MCMC_directories[4]+'results_mcmc.txt']

## Creating list of the observation file names to be analysed
## Here we are specifying that we only want the NIR files from fibre A, which is a syntax specific to CARMENES, 
## other instruments might have different file name conventions
list_science_spectra_night_1 = glob.glob(obs_directory+'*201803*nir_A.fits')
list_science_spectra_night_1.sort()
list_science_spectra_night_2 = glob.glob(obs_directory+'*201805*nir_A.fits')
list_science_spectra_night_2.sort()
list_science_spectra_night_3 = glob.glob(obs_directory+'*20190313*nir_A.fits')
list_science_spectra_night_3.sort()
list_science_spectra_night_4 = glob.glob(obs_directory+'*20190316*nir_A.fits')
list_science_spectra_night_4.sort()
list_science_spectra_night_5 = glob.glob(obs_directory+'*20190411*nir_A.fits')+glob.glob(obs_directory+'*20190412*nir_A.fits')
list_science_spectra_night_5.sort()
lists_science_spectra = [list_science_spectra_night_1,
						 list_science_spectra_night_2,
						 list_science_spectra_night_3,
						 list_science_spectra_night_4,
						 list_science_spectra_night_5]

## Pick your instruments and fill out the relevant information
instruments = ['CARMENES','CARMENES','CARMENES','CARMENES','CARMENES']
R_instruments = [80400,80400,80400,80400,80400]
n_orders = [28,28,28,28,28] 	
n_pixels = [4080,4080,4080,4080,4080]

## Define stellar parameters through the StellarParams() object. Here we have parameters for Tau Bootis
stelpars = objects.StellarParams()
stelpars.Teff = 6465 	# Effective temperature, in K, from Soubiran et al. (2022)
stelpars.R = 1.42*R_Sun	# Stellar radius in m, from Borsa et al. (2015), also used in Webb et al. (2022), consistent with Valenti & Fisher (2005), who report R_star = 1.419 R_Sun

## Defining the orbital parameters through the OrbitParams() object. Here we have parameters for Tau Bootis b
orbpars = objects.OrbitParams()
orbpars.Vsys = -16.4 		# Systemic velocity, in km/s, from Brogi et al. (2012)
orbpars.Kp = 110 			# Planet RV semi-amplitude in km/s, from Brogi et al. (2012), also consistent with Lockwood et al. (2014; 111 km/s) and Panwar et al. (2024; 110.9)
orbpars.P = 3.312433		# Orbital period, in days, from Brogi et al. (2012). There is also P = 3.312454 days from Justesen & Albrecht (2019)
orbpars.T0 = 2455652.108	# Planet mid-transit time, in HJD, from Brogi et al. (2012). There is also T0 = 2456402.3797 JD from Justesen & Albrecht (2019)
orbpars.Rp = 1.2*R_Jup		# Planet radius in m, from Webb et al. (2022)

## If I choose to use the T0 from Justesen & Albrecht (2019), I need to convert their reported T0 from MJD to BJD
## Plus, whenever calculating the phases, I need to add a 0.5 shift because their T0 is given in the rest frame of the star, no the planet!
#from astropy.coordinates import SkyCoord
#from astropy import units as u
#from astropy.time import Time
#orbpars.P = 3.312454 		# Orbital period, in days, from Justesen & Albrecht (2019)
#T0_JA19_MJD = 56401.8797	# Planet mid-transit time, in MJD, from Justesen & Albrecht (2019)
#T0_JA19_JD = T0_JA19_MJD + 2400000.5 	# Converting to JD
#RA = 206.8156
#DEC = 17.4569
#coords = SkyCoord(RA, DEC, unit=(u.deg,u.deg))
#ltt_bary_T0_JA19 = Time(T0_JA19_JD, format='jd', location=(0,0)).light_travel_time(coords)
#T0_JA19_BJD = T0_JA19_JD + ltt_bary_T0_JA19.jd # Converting from JD to BJD
#orbpars.T0 = T0_JA19_BJD
#
#orbpars.Vsys = -11.5 		# Systemic velocity, in km/s, from Webb et al. (2022)
#orbpars.Kp = 106 			# Planet RV semi-amplitude in km/s, from Webb et al. (2022)

## Define directory where the telluric removal results will be stored
main_retrv_exo_directory = home_directory+'telluric_removal/TauBoo/combined_nights/'
if not os.path.isdir(main_retrv_exo_directory):
	call(['mkdir', main_retrv_exo_directory])

## Define path to site values file (if it doesn't exist, it will be created inside the specified directory)
filenames_site_values = [main_retrv_exo_directory+f'site_values_{nights[0]}.txt',
						 main_retrv_exo_directory+f'site_values_{nights[1]}.txt',
						 main_retrv_exo_directory+f'site_values_{nights[2]}.txt',
						 main_retrv_exo_directory+f'site_values_{nights[3]}.txt',
						 main_retrv_exo_directory+f'site_values_{nights[4]}.txt']

for i in range(len(filenames_site_values)):
	if not os.path.isfile(filenames_site_values[i]):
		funcs.create_site_values_file(lists_science_spectra[i], dirname=main_retrv_exo_directory)

## List of desired molecules to include in the modelling
molecs = ['CO2', 'CH4', 'H2O', 'O2', 'N2'] 		# Molecules included via line-by-line absorption
molecs_for_cia = ['O2-air'] 					# Molecules included via collision-induced absorption (CIA)

## Choose bad spectral orders to be removed from the analysis (e.g. those that are too densely populated with saturated telluric water lines)
bad_spec_orders_idxs = [[8, 9, 10, 17, 18, 19, 20, 21, 22],
						[8, 9, 10, 17, 18, 19, 20, 21, 22],
						[8, 9, 10, 17, 18, 19, 20, 21, 22],
						[8, 9, 10, 17, 18, 19, 20, 21, 22],
						[8, 9, 10, 17, 18, 19, 20, 21, 22]]

## We will now run a function that determines "bad" observations based on their measured SNR (from the FITS header) and airmass
SNR_lims = [100,100,100,100,100]
airmass_lims = [1.7,1.7,1.7,1.7,1.7]
## Manually add any extra bad observations that may not be included in the automatic SNR and airmass flagging
add_bad_obs = [[13,32,82],[],[],[],[]] 	# These are the indexes of the spectra from the original sample, before removing any observations
create_SNR_night_log = False 			# Optional: creating a file with the exposure times and measured spectra SNR
make_obs_conds_plot = True 				# Optional: plotting the observing conditions (currently plots airmass, humidity, PWV and instrument SNR)
bad_obs_idxs = funcs_telrem.filter_bad_obs(	nights=nights,
											instruments=instruments, 
											lists_science_spectra=lists_science_spectra, 
											n_orders=n_orders, 
											bad_spec_orders_idxs=bad_spec_orders_idxs, 
											filenames_site_values=filenames_site_values, 
											filenames_mcmc_results=filenames_mcmc_results, 
											SNR_lims=SNR_lims, 
											airmass_lims=airmass_lims, 
											add_bad_obs=add_bad_obs, 
											main_directory=main_retrv_exo_directory,
											P=orbpars.P,
											T0=orbpars.T0, 
											create_SNR_night_log=create_SNR_night_log,
											make_obs_conds_plot=make_obs_conds_plot)

## Number of PCA components
nc_PCAs = [3,14,5,5,5]

## Choose the "deep line threshold", which defines the points below which transmission will be masked out of the analysis when calculating the CCF
deep_line_threshold = 0.2

## Filename (with path) to store the SNRs of each run (this is not the same SNR as the one used for the bad data clip above!)
filename_retrv_exo_SNR_AC = main_retrv_exo_directory+'SNR_logs_ac.txt'
filename_retrv_exo_SNR_PCA = main_retrv_exo_directory+'SNR_logs_pca.txt'

## Checking if the SNR log files exist and if not, creating them and adding the header
if not os.path.isfile(filename_retrv_exo_SNR_AC):
	txt = open(filename_retrv_exo_SNR_AC, 'a')

	txt.write('## Orders removed: ')
	for i in range(len(nights)):
		txt.write(str(list(np.array(bad_spec_orders_idxs[i])+1))+'\t')

	txt.write('\n')
	txt.write('## Spectra removed: ')
	for i in range(len(nights)):
		txt.write(str(list(bad_obs_idxs[i]+1))+'\t')

	txt.write('\n')
	txt.write("## SNR_inj_pos \t vRest_inj_pos (km/s) \t kpVec_inj_pos (km/s) \t SNR_max \t vRest_max (km/s) \t kpVec_max (km/s)\n")
	txt.close()

if not os.path.isfile(filename_retrv_exo_SNR_PCA):
	txt = open(filename_retrv_exo_SNR_PCA, 'a')

	txt.write('## Orders removed: ')
	for i in range(len(nights)):
		txt.write(str(list(np.array(bad_spec_orders_idxs[i])+1))+'\t')

	txt.write('\n')
	txt.write('## Spectra removed: ')
	for i in range(len(nights)):
		txt.write(str(list(bad_obs_idxs[i]+1))+'\t')

	txt.write('\n')
	txt.write("## nc_PCA \t SNR_inj_pos \t vRest_inj_pos (km/s) \t kpVec_inj_pos (km/s) \t SNR_max \t vRest_max (km/s) \t kpVec_max (km/s)\n")
	txt.close()

## Checking if the data cubes containing the wavelength, observational spectra and model telluric spectra already exist, and if not, creating them

## This creates both the wavelength and observational spectra cube, if not already existent (this only needs to be done once per night)
for i,night in enumerate(nights):
	if not os.path.isfile(main_retrv_exo_directory+f'all_lam_{night}.npy') or not os.path.isfile(main_retrv_exo_directory+f'all_spec_{night}.npy'):
		all_lam, all_spec = funcs_telrem.get_obs_cube(n_orders[i], n_pixels[i], lists_science_spectra[i], generate=True, save_dirname=main_retrv_exo_directory, save_lam_filename=f'all_lam_{night}.npy', save_spec_filename=f'all_spec_{night}.npy')

	## This creates the model spectra cube, if not already existent (this step takes a while, but only needs to be done once)
	if not os.path.isfile(main_retrv_exo_directory+f'all_spec_mod_{night}.npy'):
		all_spec_mod = funcs_telrem.get_mod_cube(n_orders[i], n_pixels[i], lists_science_spectra[i], MCMC_directories[i], atm_profs_directory, filenames_mcmc_results[i], molecs, molecs_for_cia, generate=True, save_dirname=main_retrv_exo_directory, save_spec_mod_filename=f'all_spec_mod_{night}.npy')

## Decide if you want to make certain optional plots
make_gen_cube_plots = False
plot_telrem_steps_AC = False
plot_telrem_steps_PCA = False
plot_indiv_PCA_comps = False

## Checking if the plots directory already exists and if not, creating it
dirname_for_plots = main_retrv_exo_directory+'plots/'
if not os.path.isdir(dirname_for_plots):
	call(['mkdir', dirname_for_plots])

## Checking if the directory where the NumPy objects will be saved already exists and if not, creating it
numpy_objs_directory = main_retrv_exo_directory+'analysis_files/'
if not os.path.isdir(numpy_objs_directory):
	call(['mkdir', numpy_objs_directory])

all_nights_cubes = retrv_exo_atm.telrem_each_night( nights=nights,
													main_retrv_exo_directory=main_retrv_exo_directory, 
													filenames_site_values=filenames_site_values, 
													n_orders=n_orders,
													bad_spec_orders_idxs=bad_spec_orders_idxs, 
													bad_obs_idxs=bad_obs_idxs, 
													nc_PCAs=nc_PCAs, 
													filename_retrv_exo_SNR_AC=filename_retrv_exo_SNR_AC, 
													filename_retrv_exo_SNR_PCA=filename_retrv_exo_SNR_PCA, 
													orbpars=orbpars, 
													make_gen_cube_plots=make_gen_cube_plots,
													plot_telrem_steps_AC=plot_telrem_steps_AC, 
													plot_telrem_steps_PCA=plot_telrem_steps_PCA, 
													plot_indiv_PCA_comps=plot_indiv_PCA_comps,
													deep_line_threshold=deep_line_threshold,
													dirname_for_plots=dirname_for_plots,
													numpy_objs_directory=numpy_objs_directory)

## If you've already run the function retrv_exo_atm.telrem_each_night(), no need to run it again, just unpack the values
night_1_cube = np.load(numpy_objs_directory+f'cubes_{nights[0]}.npy')
night_2_cube = np.load(numpy_objs_directory+f'cubes_{nights[1]}.npy')
night_3_cube = np.load(numpy_objs_directory+f'cubes_{nights[2]}.npy')
night_4_cube = np.load(numpy_objs_directory+f'cubes_{nights[3]}.npy')
night_5_cube = np.load(numpy_objs_directory+f'cubes_{nights[4]}.npy')

## Now we have to create the variable "all_nights_cubes" to match the format given by retrv_exo_atm.telrem_each_night() and as taken by retrv_exo_atm.get_planet_sig_all_nights()
all_nights_cubes = []
for ii in range(len(night_1_cube)):
	if ii == 6:
		all_nights_cubes.append([night_1_cube[ii].astype(bool), 
								 night_2_cube[ii].astype(bool),
								 night_3_cube[ii].astype(bool),
								 night_4_cube[ii].astype(bool),
								 night_5_cube[ii].astype(bool)])
	else:
		all_nights_cubes.append([night_1_cube[ii], 
								 night_2_cube[ii],
								 night_3_cube[ii],
								 night_4_cube[ii],
								 night_5_cube[ii]])

## Specify filenames where the planet model wavelength and flux are kept (models here come from GENESIS code from Gandhi & Madhusudhan 2017)
filename_planet_lam = home_directory+'telluric_removal/planet_models/lam.txt' 							# Wavelength in microns
filename_planet_flux = home_directory+'telluric_removal/planet_models/h2o_-3.0_ch4_-20.0_hcn_-20.0.txt' # Flux in W m^-2 m^-1 (W/m)

## Defining the parameters for the straight lines that mark the "signal" and "noise" region in the SNR plots.
## y = a*x + b
SNR_lines_params = [[1.8,20], [-2.2,20], [3.0,20], [-3.5,20], [2.0,20]]

all_nights_results = retrv_exo_atm.get_planet_sig_all_nights(	nights=nights, 
																all_nights_cubes=all_nights_cubes, 
																main_retrv_exo_directory=main_retrv_exo_directory, 
																filenames_site_values=filenames_site_values, 
																R_instruments=R_instruments, 
																n_orders=n_orders, 
																bad_spec_orders_idxs=bad_spec_orders_idxs, 
																bad_obs_idxs=bad_obs_idxs, 
																nc_PCAs=nc_PCAs, 
																filename_retrv_exo_SNR_AC=filename_retrv_exo_SNR_AC, 
																filename_retrv_exo_SNR_PCA=filename_retrv_exo_SNR_PCA, 
																filename_planet_lam=filename_planet_lam, 
																filename_planet_flux=filename_planet_flux, 
																stelpars=stelpars, 
																orbpars=orbpars, 
																SNR_lines_params=SNR_lines_params,
																dirname_for_plots=dirname_for_plots,
																numpy_objs_directory=numpy_objs_directory)

## If you've already run the function retrv_exo_atm.get_planet_sig_all_nights(), no need to run it again, just unpack the values
all_nights_phases = []
all_nights_RVs = []
all_nights_tr = []
all_nights_ll = []
all_nights_planet_trail = []
all_nights_vmap = []
all_nights_SNR = []

all_nights_tr_PCA = []
all_nights_ll_PCA = []
all_nights_planet_trail_PCA = []
all_nights_vmap_PCA = []
all_nights_SNR_PCA = []

for night in nights:
	print(night)
	phases = np.load(numpy_objs_directory+f'phases_{night}.npy')
	RVs = np.load(numpy_objs_directory+f'RVs_{night}.npy')
	tr = np.load(numpy_objs_directory+f'tr_{night}.npy')
	ll = np.load(numpy_objs_directory+f'll_{night}.npy')
	planet_trail_mean_rem = np.load(numpy_objs_directory+f'planet_trail_{night}.npy')
	vmap = np.load(numpy_objs_directory+f'vmap_{night}.npy')
	SNR = np.load(numpy_objs_directory+f'SNR_{night}.npy')
	tr_PCA = np.load(numpy_objs_directory+f'tr_PCA_{night}.npy')
	ll_PCA = np.load(numpy_objs_directory+f'll_PCA_{night}.npy')
	planet_trail_PCA_mean_rem = np.load(numpy_objs_directory+f'planet_trail_PCA_{night}.npy')
	vmap_PCA = np.load(numpy_objs_directory+f'vmap_PCA_{night}.npy')
	SNR_PCA = np.load(numpy_objs_directory+f'SNR_PCA_{night}.npy')

	all_nights_phases.append(phases)
	all_nights_RVs.append(RVs)
	all_nights_tr.append(tr)
	all_nights_ll.append(ll)
	all_nights_planet_trail.append(planet_trail_mean_rem)
	all_nights_vmap.append(vmap)
	all_nights_SNR.append(SNR)

	all_nights_tr_PCA.append(tr_PCA)
	all_nights_ll_PCA.append(ll_PCA)
	all_nights_planet_trail_PCA.append(planet_trail_PCA_mean_rem)
	all_nights_vmap_PCA.append(vmap_PCA)
	all_nights_SNR_PCA.append(SNR_PCA)

## Now we have to create the variable "all_nights_results" to match the format given by retrv_exo_atm.get_planet_sig_all_nights() and as taken by retrv_exo_atm.combine_sig_multi_nights()
all_nights_results = [	all_nights_phases,
						all_nights_RVs,
						all_nights_tr,
						all_nights_ll,
						all_nights_planet_trail,
						all_nights_vmap,
						all_nights_SNR,
						all_nights_tr_PCA,
						all_nights_ll_PCA,
						all_nights_planet_trail_PCA,
						all_nights_vmap_PCA,
						all_nights_SNR_PCA]

dirname_for_plots = main_retrv_exo_directory+'plots/'
filename_suffix = '_all_nights'
retrv_exo_atm.combine_sig_multi_nights(all_nights_results, dirname_for_plots, filename_suffix)
