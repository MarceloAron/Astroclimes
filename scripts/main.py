# =====================================================================================
# Basic packages
# =====================================================================================
import numpy as np
import os
from subprocess import call
import emcee
import copy
from multiprocessing import Pool

# =====================================================================================
# Scripts
# =====================================================================================
import funcs
import fit
import plots

## Color-blind friendly colors: blue, orange, green, pink, brown, purple, gray, red and yellow
CBcols= ['#377EB8', '#FF7F00', '#4DAF4A',
		 '#F781BF', '#A65628', '#984EA3',
		 '#999999', '#E41A1C', '#DEDE00']

## Some necessary constants
kB = 1.380649*1e-23	# Boltzmann constant, in m^2 kg s^-2 K^-1
R = 8.31446 		# Gas constant, in J mol^-1 K^-1
N_A = 6.022*1e23	# Avogadro constant, in molecules/mole
c = 299792458		# Speed of light, in m/s
DU = 2.69*1e20		# Dobson unit, in molecules/m^2

def run_main(atm_profs_directory, MCMC_directory, filename_mcmc_results, filename_em_line_spec, list_science_spectra, filename_site_values, include_stelmod, filename_stelmod_lam, filename_stelmod_spec, molecs, molecs_for_cia, free_molecs, spec_orders, instrument, R_instrument, stelpars, orbpars, scale_profs, vel_step, DMF_O2, n_CPUs=1):

	R_regrid = (c*1e-3)/vel_step	# Resolution of our regridded model 
	if include_stelmod:
		## Unpacking and preparing the stellar model spectra to be used in the analysis
		lam_stelmod, spec_stelmod = funcs.get_phoenix_spectra(filename_stelmod_lam, filename_stelmod_spec, vel_step)
		lam_stelmod_broad, spec_stelmod_broad = funcs.broaden_phoenix_spectra(lam_stelmod, spec_stelmod, vel_step, stelpars.vsini, stelpars.eps)	# Broadening to include rotational broadening
		lam_stelmod_broad, spec_stelmod_broad = funcs.inst_broaden_gauss(lam_stelmod_broad, spec_stelmod_broad, R_regrid, R_instrument)				# Broadening to include instrumental broadening

	## Looping over all spectra
	for filename_spectra in list_science_spectra:
		## Creating the folder to store the individual results if not already existent
		dirname = MCMC_directory+filename_spectra.split('/')[-1].split('.')[0]+'/'
		if not os.path.isdir(dirname):
			call(['mkdir', dirname])

		## Getting all of the observational data (NaNs and infs are already removed)
		hdr_obs, all_lam_obs, all_spec_obs = funcs.get_spectral_data(filename_spectra, instrument=instrument)

		## Unpacking data from FITS file header, tailored per instrument
		date_obs, BJD, hgt_site, P_site, T_site, relhum_site, airmass, V_BERV = funcs.get_header_data(hdr_obs, instrument)
		
		## Unpacking information from the GGG2020 atmospheric profile
		hgt_atm, P_atm, T_atm, Xs_atm = funcs.get_ggg2020_atm_profiles(date_obs, atm_profs_directory, molecs, hgt_site)

		## If choosen, scaling the pressure, temperature and humidity profiles to match the site measurements
		if scale_profs:
			T_atm = T_site*(T_atm/T_atm[-1])
			P_atm = P_site*(P_atm/P_atm[-1])
			Xs_atm['H2O'] = relhum_site*(Xs_atm['H2O']/Xs_atm['H2O'][-1])
		
		## Plotting the initial atmospheric profile (pressure, temperature, CO2 abundance and humidity as a function of height)
		plots.plot_atm_prof([[hgt_atm, 1e-2*10**P_atm, T_atm, Xs_atm['CO2'], Xs_atm['H2O']]],
							[''],
							savefig=True, dirname=dirname, filename='atm_profile')
		
		## Plotting the initial molecular abundance profiles
		plots.plot_molecs_atm_prof([[hgt_atm, Xs_atm['CO2'], Xs_atm['CH4'], Xs_atm['H2O'], Xs_atm['O2'], Xs_atm['N2']]],
							[''],
							savefig=True, dirname=dirname, filename='molecs_atm_profile')

		## Converting the relative abundances in ppm to number densities using equation (9) from Laughner et al. (2023a)
		Cn_atm = {}
		Sn_atm = {}
		nideal = (10**P_atm)/(kB*T_atm)
		for key in Xs_atm:
			Cn_atm[key] = 1e-6*Xs_atm[key]*nideal/(1+1e-6*Xs_atm['H2O'])
			
			## Integrating the number densities over the atmosphere column to get the total column amount of each molecule (in molecules/m^2, converted to Dobson units)
			Sn_atm[key] = np.trapz(Cn_atm[key][::-1], x=hgt_atm[::-1])/DU

		## Calculating the column averaged volume mixing ratios (the expression can be found in Morino et al. 2011, Laughner et al. 2023a and OCO-2 Retrieval Algorithm)
		Xs_atm_CA_init = [1e6*DMF_O2*(Sn_atm[key]/Sn_atm['O2']) for key in Xs_atm]	## Initial column-averaged dry air mole fractions
		Xs_atm_ground_init = [Xs_atm[key][-1] for key in Xs_atm]					## Initial ground-level dry air mole fractions

		## Keeping only the desired orders, lam_obs/spec_obs are lists that contain the selected orders
		## We store them on a list to allow for different orders to have different numbers of data points
		## (in case any are removed for being NaNs, infs or any emission lines)
		lam_obs = [all_lam_obs[x] for x in spec_orders]
		spec_obs = [all_spec_obs[x] for x in spec_orders]

		lam_ranges = [[np.min(lam_obs[i]), np.max(lam_obs[i])] for i in range(len(lam_obs))]
		
		## Merging consecutive orders - this is done simply to limit the spectral range we work with in get_cross_secs_dic() in the interest of saving computational time
		merge_lam_ranges = []
		low_end = lam_ranges[0][0]
		for i in range(len(lam_ranges)):
			try:
				high_end = lam_ranges[i][1]
				next_low_end = lam_ranges[i+1][0]
				next_high_end = lam_ranges[i+1][1]
				if 1.01*high_end > 0.99*next_low_end:
					#print(f'Merging orders {spec_orders[i]+1} and {spec_orders[i+1]+1}')
					high_end = next_high_end
				else:
					merge_lam_ranges.append([low_end, high_end])
					low_end = lam_ranges[i+1][0]
			except IndexError:
				merge_lam_ranges.append([low_end, high_end])
				#merge_lam_ranges.append([lam_ranges[i][0], lam_ranges[i][1]])
		## Adding the last order if it hasn't been merged with the previous ones
		if merge_lam_ranges[-1][1] < lam_ranges[-1][1]:
			merge_lam_ranges.append([lam_ranges[-1][0], lam_ranges[-1][1]])
		
		## Getting the emission sky model, already interpolated to the observational spectra wavelength grid
		spec_skymodel = funcs.get_skymodel_spectra(filename_em_line_spec, lam_obs)
		
		## Using the emission line spectra to mask out the points associated with emission lines
		## Anything above the defined threshold will be considered an emission line (thresholds are determined by visual inspection for each order)
		em_lines_lims = [1000, 1000, 1000, 5000, 7000, 20000, 8000, 20000, 30000, 25000, 130000, 500, 1000, 5000, 100000, 100000, 4000, 40000, 1000, 1000, 15000, 70000, 100000, 15000, 40000, 4000, 2000, 3000]
		em_lines_lims = [em_lines_lims[x] for x in spec_orders]
		for i in range(len(lam_obs)):
			em_lines_mask = spec_skymodel[i] > em_lines_lims[i]
			plots.plot_em_line_mask(lam_obs[i], spec_obs[i], spec_skymodel[i], em_lines_mask, savefig=True, dirname=dirname, filenames=[f'em_lines_mask_ord_{spec_orders[i]+1}', f'spec_em_lines_mask_ord_{spec_orders[i]+1}'])
			lam_obs[i] = lam_obs[i][~em_lines_mask]
			spec_obs[i] = spec_obs[i][~em_lines_mask]

		## Saving the non-normalised spectra in a .txt file in case it is needed for troubleshooting and debugging
		for i in range(len(spec_orders)):
			at = (lam_obs[i]*1e6, spec_obs[i])
			np.savetxt(dirname+f'not_norm_spec_ord_{spec_orders[i]+1}.txt', np.vstack(at).T)

		## Getting the information from the cross-section files
		molecs_cross_secs = funcs.get_cross_secs_dic(molecs, merge_lam_ranges, vel_step)

		## Interpolating the cross sections
		int_cross_secs = funcs.interpolate_cross_secs(P_atm, T_atm, hgt_atm, molecs, molecs_cross_secs, airmass) 

		## Getting the information from the CIA files
		lam = molecs_cross_secs['CO2']['lam']
		molecs_cias = funcs.get_cia_dict(lam, molecs_for_cia)

		## Interpolating the CIA
		int_cias = funcs.interpolate_cia(P_atm, T_atm, hgt_atm, molecs_for_cia, molecs_cias, airmass)

		## Calculating Rayleigh and aerosol scattering
		tau_rayleigh = funcs.rayleigh_scatter(lam*1e6, (10**P_site)*1e-2, hgt_site*1e-3)
		tau_aerosol = funcs.aerosol_scatter(lam*1e6, airmass)
		
		## Calculating the star's radial velocity on the observer's reference frame and shifting the model stellar spectra to said frame
		phases = ((BJD - orbpars.T0)%orbpars.P)/orbpars.P 									# Orbital phases (for Tau Boo b!)
		RV_star = orbpars.Vsys - V_BERV + orbpars.Ks*np.sin(2*np.pi*(phases+0.5)) 			# Phases are calculated for planet, so we add 0.5 to get star orbital phase
		lam_stelmod_shift = lam_stelmod_broad*(1+RV_star/(c*1e-3))
		
		## Interpolating the model stellar spectra to the observational wavelength distribution
		if include_stelmod:
			spec_stelmod_int = [np.interp(lam_obs[i], lam_stelmod_shift, spec_stelmod_broad) for i in range(len(lam_obs))]
			spec_stelmod_norm = [spec_stelmod_int[i]/np.max(spec_stelmod_int[i]) for i in range(len(spec_stelmod_int))]
		else:
			spec_stelmod_norm = [np.ones(len(spec_obs[i])) for i in range(len(lam_obs))] 	## This will just define the model stellar spectra as an array of ones

		## Window sizes for the median filters (one for each order included)
		window_sizes_1 = [0.10]*len(spec_orders)
		window_sizes_2 = [0.01]*len(spec_orders)

		## Getting the normalisation mask that will be used to normalise the observational and model data
		sel_CO2_ground = 8*1e21		# Arbitrary CO2 ground level for dummy model spectra used in the normalisation process, in molecules/m^3
		sel_H2O_ground = 1.5*1e23 	# Arbitrary H2O ground level for dummy model spectra used in the normalisation process, in molecules/m^3
		norm_mask = funcs.get_normalisation_mask(lam, lam_obs, int_cross_secs, int_cias, tau_rayleigh, tau_aerosol, Cn_atm, molecs, molecs_for_cia, R_regrid, R_instrument, sel_CO2_ground, sel_H2O_ground, window_sizes_1, spec_stelmod_norm, spec_orders, plot_bool=True, dirname=dirname)
		spec_obs_norm = funcs.normalise_spectra(lam_obs, spec_obs, norm_mask, window_sizes_2, spec_orders, plot_bool=True, dirname=dirname)

		## Computing the initial model spectra to save it and later compare it to the best-fit model spectra
		mod_to_obs_spec_init = funcs.spec_handle(lam, lam_obs, int_cross_secs, int_cias, tau_rayleigh, tau_aerosol, Cn_atm, molecs, molecs_for_cia, R_regrid, R_instrument)
		spec_mod_init_with_stellar_lines = [mod_to_obs_spec_init[i]*spec_stelmod_norm[i] for i in range(len(mod_to_obs_spec_init))]
		spec_mod_norm_init = funcs.normalise_spectra(lam_obs, spec_mod_init_with_stellar_lines, norm_mask, window_sizes_2, spec_orders)
		
		## Determining cutoff mask, which will remove points below a certain transmissivity from the calculation of the log likelihood
		## This is done because we think the deep lines might be skewing the results
		cutoff_mask = []
		for spec_norm in spec_obs_norm:
			cutoff_mask.append(spec_norm > 0.2)

		## Setting up the MCMC	
		profiles = []
		mcmc_params = []
		for fm in free_molecs:
			profiles.append(Cn_atm[fm]/Cn_atm[fm][-1])
			mcmc_params.append(Cn_atm[fm][-1])
		mcmc_params = np.array(mcmc_params)

		if 'O2' in free_molecs:
			free_molecs_minus_O2 = copy.deepcopy(free_molecs)
			free_molecs_minus_O2.remove('O2')
		else:
			free_molecs_minus_O2 = free_molecs

		mcmc_par_names = ['Cn_'+fm for fm in free_molecs]
		pars_for_plots = ['g_'+fm for fm in free_molecs]
		pars_for_plots.append('logp')
		labels = ['g_'+fm+' (ppm)' for fm in free_molecs]
		labels.append('logp')
		pars_for_plots_2 = ['X_'+fm+'_CA' for fm in free_molecs_minus_O2]
		labels_2 = ['X_'+fm+' (ppm)' for fm in free_molecs_minus_O2]
		
		pos = mcmc_params + 0.1*mcmc_params*np.random.randn(10, len(mcmc_params))
		nwalkers, ndim = pos.shape
		max_steps = 5000
		check_each_n_steps = 100
		index = 0
		autocorr = np.zeros(max_steps)
		old_tau = -np.inf
		taufile = open(dirname+'autocorr.txt', 'a')

		backfile = dirname+'backupmcmc.h5'
		backend = emcee.backends.HDFBackend(backfile)
		backend.reset(nwalkers, ndim)

		## Saving the necessary parameters into a temporary .npy file so they can be unpacked at the start of the MCMC 
		## instead of giving them as parameters for the log_probability function (this saves up computational time)
		MCMC_glob_pars = np.array([lam_obs, spec_obs_norm, R_regrid, R_instrument, lam, molecs, molecs_for_cia, Cn_atm, int_cross_secs, int_cias, tau_rayleigh, tau_aerosol, norm_mask, window_sizes_2, cutoff_mask, spec_stelmod_norm, free_molecs, profiles], dtype=object)
		temp_npy_file = f'temp_MCMC_{list_science_spectra.index(filename_spectra)}.npy' 
		np.save(temp_npy_file, MCMC_glob_pars, allow_pickle=True)

		fit.unpack_global_vars(temp_npy_file)

		## Running the MCMC
		with Pool(processes=n_CPUs) as pool:
			#sampler = emcee.EnsembleSampler(nwalkers, ndim, fit.log_probability, 
			#								args=(lam_obs, spec_obs_norm, R_regrid, R_instrument, lam, molecs, molecs_for_cia, Cn_atm, int_cross_secs, int_cias, tau_rayleigh, tau_aerosol, norm_mask, window_sizes_2, cutoff_mask, spec_stelmod_norm, free_molecs, profiles),
			#								backend=backend)
			sampler = emcee.EnsembleSampler(nwalkers, ndim, fit.log_probability, 
											backend=backend, pool=pool)
			for sample in sampler.sample(pos, iterations=max_steps, progress=True):
				if sampler.iteration % check_each_n_steps:
					continue
				
				tau = sampler.get_autocorr_time(tol=0)
				autocorr[index] = np.mean(tau)
				taufile.write('{} \t {:.4f} \n'.format(sampler.iteration,autocorr[index]))
				index += 1

				converged = (np.all(tau*100 < sampler.iteration)) & (np.all(np.abs(old_tau - tau)/tau < 0.01))
				if converged:
					print('MCMC converged after {} steps'.format(sampler.iteration))
					break
				old_tau = tau
			taufile.close()

		call(['rm', temp_npy_file])
			
		#sampler = funcs.get_mcmc_result(dirname)
		samples = sampler.get_chain()
		logp = sampler.get_log_prob()
		if abs(100*np.std(logp[-1,:])/np.median(logp[-1,:])) > 10:
			bad_walkers = np.where(logp[-1,:] < (np.median(logp[-1,:]) - np.std(logp[-1,:])))[0]
		else:
			bad_walkers = []
		samples = np.delete(samples, bad_walkers, axis=1)
		logp = np.delete(logp, bad_walkers, axis=1)
		nsteps = np.shape(samples)[0]
		discard, thin = np.max([np.where(abs(100*np.std(logp, axis=1)/np.median(logp, axis=1)) < 0.1)[0][0], int(nsteps/2)]), 1
		flat_samples = sampler.get_chain(discard=discard, thin=thin, flat=True)
		logp = sampler.get_log_prob()
		logpflat = sampler.get_log_prob(discard=discard, thin=thin, flat=True)
		samples, flat_samples, logp, logpflat = funcs.delete_bad_walkers(samples, flat_samples, logp, logpflat, bad_walkers, nsteps, nwalkers, discard, thin)

		## Unpacking the results from the MCMC
		mcmc_median = np.median(flat_samples, axis=0)
		mcmc_std = np.std(flat_samples, axis=0)
		
		Cn_atm_f = copy.deepcopy(Cn_atm)
		for fm in free_molecs:
			Cn_atm_f[fm] = mcmc_median[free_molecs.index(fm)]*profiles[free_molecs.index(fm)]

		## Saving the molecular profiles (it becomes easier to generate the final models without having to rerun a bunch of stuff)
		np.save(dirname+'Cn_atm_f.npy', Cn_atm_f, allow_pickle=True)

		## From the new number density profiles (Cn_atm_f, which have the same shape as the old ones, just scaled),
		## I calculate the final DMFs using a rearranged form of equation (9) from Laughner et al. (2023a)
		Xs_atm_f = dict.fromkeys(molecs)
		Sn_atm_f = dict.fromkeys(molecs)
		Xs_atm_f['H2O'] = 1e6/(nideal/Cn_atm_f['H2O'] - 1)	## First I separately calculate the posterior H2O DMF, in ppm (hence the 1e6 factor!)
		for key in Xs_atm:
			Xs_atm_f[key] = 1e6*Cn_atm_f[key]*(1+1e-6*Xs_atm_f['H2O'])/nideal
			Sn_atm_f[key] = np.trapz(Cn_atm_f[key][::-1], x=hgt_atm[::-1])/DU
		
		## Plotting the final atmospheric profile (pressure, temperature, CO2 abundance and humidity as a function of height)
		plots.plot_atm_prof([[hgt_atm, 1e-2*10**P_atm, T_atm, Xs_atm_f['CO2'], Xs_atm_f['H2O']]],
							[''],
							savefig=True, dirname=dirname, filename='atm_profile_f')
		
		## Plotting the final molecular abundance profiles
		plots.plot_molecs_atm_prof([[hgt_atm, Xs_atm_f['CO2'], Xs_atm_f['CH4'], Xs_atm_f['H2O'], Xs_atm_f['O2'], Xs_atm_f['N2']]],
							[''],
							savefig=True, dirname=dirname, filename='molecs_atm_profile_f')

		## Now I will define some quantities necessary to make the chain and corner plots
		Xs = 1e6*samples*(1+1e-6*Xs_atm_f['H2O'][-1])/nideal[-1]
		Cn = np.einsum('ijk,kl->ijkl',samples, profiles)
		Sn = np.trapz(Cn[:,:,:, ::-1], x=hgt_atm[::-1])/DU

		Xs_flat = 1e6*flat_samples*(1+1e-6*Xs_atm_f['H2O'][-1])/nideal[-1]
		Cn_flat = np.einsum('ik,kl->ikl',flat_samples, profiles)
		Sn_flat = np.trapz(Cn_flat[:,:, ::-1], x=hgt_atm[::-1])/DU
		
		smps = np.concatenate((Xs, logp[:, :, np.newaxis]), axis=2)
		flat_smps = np.concatenate((Xs_flat, logpflat[:, np.newaxis]), axis=1)
		
		try:
			O2_index = mcmc_par_names.index('Cn_O2')
			smps_2 = 1e6*DMF_O2*np.delete(Sn, O2_index, axis=2)/Sn[:,:,O2_index,np.newaxis]
			flat_smps_2 = 1e6*DMF_O2*np.delete(Sn_flat, O2_index, axis=1)/Sn_flat[:,O2_index,np.newaxis]
		except ValueError:
			smps_2 = 1e6*DMF_O2*Sn/Sn_atm_f['O2']
			flat_smps_2 = 1e6*DMF_O2*Sn_flat/Sn_atm_f['O2']

		init_vals = [Xs_atm_ground_init[free_molecs.index(fm)] for fm in free_molecs]
		init_vals.append(logpflat[0])
		init_vals_2 = [Xs_atm_CA_init[free_molecs.index(fm)] for fm in free_molecs_minus_O2]
		plots.make_chain_plot(smps, len(pars_for_plots), labels, savefig=True, dirname=dirname)
		plots.make_chain_plot(smps_2, len(pars_for_plots_2), labels_2, savefig=True, dirname=dirname, filename='chain_2')
		plots.make_corner_plot(len(pars_for_plots), flat_smps, logpflat, labels, init_vals, savefig=True, dirname=dirname)
		plots.make_corner_plot(len(pars_for_plots_2), flat_smps_2, logpflat, labels_2, init_vals_2, savefig=True, dirname=dirname, filename='corner_2')
		plots.make_autocorr_plot(check_each_n_steps, dirname=dirname)

		## Defining the variables to print in the results file
		## These are the ground level values and their uncertainties
		final_g_vals = np.array([Xs_atm_f[key][-1] for key in Xs_atm_f])
		final_g_uncs = np.zeros(len(molecs))
		for fm in free_molecs:
			final_g_uncs[molecs.index(fm)] = 1e6*mcmc_std[free_molecs.index(fm)]*(1+1e-6*Xs_atm_f['H2O'][-1])/nideal[-1]

		CA_vals = np.median(flat_smps_2, axis=0)
		CA_errs = np.std(flat_smps_2, axis=0)
		final_X_vals = np.array([1e6*DMF_O2*(Sn_atm_f[key]/Sn_atm_f['O2']) for key in Xs_atm_f])
		final_X_uncs = np.zeros(len(molecs))
		for fm in free_molecs_minus_O2:
			final_X_vals[molecs.index(fm)] = CA_vals[free_molecs_minus_O2.index(fm)]
			final_X_uncs[molecs.index(fm)] = CA_errs[free_molecs_minus_O2.index(fm)]

		logp_med, u_logp_med = np.median(logpflat, axis=0), np.std(logpflat, axis=0)

		text =	(f"{date_obs.value[:-4]:19s} \t " +
				f"{final_g_vals[0]:10.4f} \t {final_g_uncs[0]:6.4f} \t " +
				f"{final_g_vals[1]:10.4f} \t {final_g_uncs[1]:6.4f} \t " +
				f"{final_g_vals[2]:10.4f} \t {final_g_uncs[2]:6.4f} \t " +
				f"{final_g_vals[3]:10.4f} \t {final_g_uncs[3]:6.4f} \t " +
				f"{final_g_vals[4]:10.4f} \t {final_g_uncs[4]:6.4f} \t " +
				f"{final_X_vals[0]:10.4f} \t {final_X_uncs[0]:6.4f} \t " +
				f"{final_X_vals[1]:10.4f} \t {final_X_uncs[1]:6.4f} \t " +
				f"{final_X_vals[2]:10.4f} \t {final_X_uncs[2]:6.4f} \t " +
				f"{final_X_vals[3]:10.4f} \t {final_X_uncs[3]:6.4f} \t " +
				f"{final_X_vals[4]:10.4f} \t {final_X_uncs[4]:6.4f} \t " +
				f"{logp_med:10.4f} \t {u_logp_med:10.4f}\n")
		txt = open(MCMC_directory+filename_mcmc_results, 'a')
		txt.write(text)
		txt.close()

		## Computing and plotting the model with the results from the MCMC
		xlims = [[0.9675,0.9700], [0.9775,0.9800], [1.0060,1.0100], [1.0120,1.0130], [1.0325,1.0340], [1.0575,1.0600], [1.0675,1.0700], [1.0825,1.0850], [1.115,1.118], [1.125,1.128], [1.160,1.163], [1.1794,1.1799], [1.203,1.206], [1.210,1.213], [1.255,1.258], [1.279, 1.281], [1.305,1.311], [1.320,1.325], [1.345,1.350], [1.395,1.400], [1.425,1.432], [1.445,1.450], [1.485,1.490], [1.515,1.520], [1.5713,1.5718], [1.610,1.612], [1.655,1.660], [1.705,1.710]]
		xlim = [xlims[x] for x in spec_orders]
		mod_to_obs_spec = funcs.spec_handle(lam, lam_obs, int_cross_secs, int_cias, tau_rayleigh, tau_aerosol, Cn_atm_f, molecs, molecs_for_cia, R_regrid, R_instrument)
		spec_mod_with_stellar_lines = [mod_to_obs_spec[i]*spec_stelmod_norm[i] for i in range(len(mod_to_obs_spec))]
		spec_mod_norm = funcs.normalise_spectra(lam_obs, spec_mod_with_stellar_lines, norm_mask, window_sizes_2, spec_orders)
		for i in range(len(spec_orders)):
			zoom = True
			plots.plot_tel_lines(lam_obs[i]*1e6, [spec_obs_norm[i], spec_mod_norm_init[i], spec_mod_norm[i]], labels=['Obs. spectra', 'Initial model', 'Best-fit model'], zoom=zoom, xlim=xlim[i], cols=['gray', CBcols[0], CBcols[1]], savefig=True, dirname=dirname, filename=f'model_spectra_ord_{spec_orders[i]+1}')
			zoom = False
			plots.plot_tel_lines(lam_obs[i]*1e6, [spec_obs_norm[i], spec_mod_norm_init[i], spec_mod_norm[i]], labels=['Obs. spectra', 'Initial model', 'Best-fit model'], zoom=zoom, xlim=xlim[i], cols=['gray', CBcols[0], CBcols[1]], savefig=True, dirname=dirname, filename=f'model_spectra_zoom_ord_{spec_orders[i]+1}')
			at = (lam_obs[i]*1e6, spec_obs_norm[i], spec_mod_norm[i], spec_mod_norm_init[i])
			np.savetxt(dirname+f'spec_ord_{spec_orders[i]+1}.txt', np.vstack(at).T)

		print(f"Progress: {list_science_spectra.index(filename_spectra)+1}/{len(list_science_spectra)}")
