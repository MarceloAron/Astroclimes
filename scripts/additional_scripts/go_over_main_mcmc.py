'''
Script to go over the results from the molecular abundance MCMC for whatever reason
'''
# =====================================================================================
# Basic packages
# =====================================================================================
import numpy as np
import os
from subprocess import call
import emcee
import copy
from multiprocessing import Pool
import sys

# =====================================================================================
# Scripts
# =====================================================================================
sys.path.append('/home/marceloaron/MarceloAron/Astroclimes/scripts/')
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

def run_main_lite(atm_profs_directory, MCMC_directory, filename_mcmc_results, filename_em_line_spec, list_science_spectra, include_stelmod, filename_stelmod_lam, filename_stelmod_spec, molecs, molecs_for_cia, free_molecs, spec_orders, instrument, R_instrument, stelpars, orbpars, scale_profs, vel_step, DMF_O2, n_CPUs=1):

	R_regrid = (c*1e-3)/vel_step	# Resolution of our regridded model 

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
		
		pos = mcmc_params + 0.1*mcmc_params*np.random.randn(10, len(mcmc_params))
		nwalkers, ndim = pos.shape
		max_steps = 5000

		sampler = funcs.get_mcmc_result(dirname)
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

		## From the new number density profiles (Cn_atm_f, which have the same shape as the initial ones, just scaled),
		## I calculate the final DMFs using a rearranged form of equation (9) from Laughner et al. (2023a) - these can have different shapes than the initial ones
		Xs_atm_f = dict.fromkeys(molecs)
		Sn_atm_f = dict.fromkeys(molecs)
		Xs_atm_f['H2O'] = 1e6/(nideal/Cn_atm_f['H2O'] - 1)	## First I separately calculate the posterior H2O DMF, in ppm (hence the 1e6 factor!)
		for key in Xs_atm:
			Xs_atm_f[key] = 1e6*Cn_atm_f[key]*(1+1e-6*Xs_atm_f['H2O'])/nideal
			Sn_atm_f[key] = np.trapz(Cn_atm_f[key][::-1], x=hgt_atm[::-1])/DU
		
		## Now I will define some quantities necessary to make the chain and corner plots
		Xs = 1e6*samples*(1+1e-6*Xs_atm_f['H2O'][-1])/nideal[-1]
		Cn = np.einsum('ijk,kl->ijkl',samples, profiles)
		Sn = np.trapz(Cn[:,:,:, ::-1], x=hgt_atm[::-1])/DU

		Xs_flat = 1e6*flat_samples*(1+1e-6*Xs_atm_f['H2O'][-1])/nideal[-1]
		Cn_flat = np.einsum('ik,kl->ikl',flat_samples, profiles)
		Sn_flat = np.trapz(Cn_flat[:,:, ::-1], x=hgt_atm[::-1])/DU

		H2O_index = mcmc_par_names.index('Cn_H2O')
		H2O_prof_f = Xs_atm_f['H2O']/Xs_atm_f['H2O'][-1]
	
		for_PWV = np.einsum('ij,l->ijl',Xs[:,:,H2O_index], H2O_prof_f)
		PWVs = np.zeros(shape=np.shape(Xs[:,:,H2O_index]))
		for i in range(np.shape(PWVs)[0]):
			for j in range(np.shape(PWVs)[1]):
				PWVs[i,j] = funcs.get_PWV(hgt_atm*airmass, P_atm, T_atm, for_PWV[i,j,:])

		for_PWV_flat = np.einsum('i,l->il',Xs_flat[:,H2O_index], H2O_prof_f)
		PWVs_flat = np.zeros(len(Xs_flat[:,H2O_index]))
		for i in range(len(PWVs_flat)):
			PWVs_flat[i] = funcs.get_PWV(hgt_atm*airmass, P_atm, T_atm, for_PWV_flat[i,:])
		
		smps = np.concatenate((Xs, logp[:, :, np.newaxis]), axis=2)
		flat_smps = np.concatenate((Xs_flat, logpflat[:, np.newaxis]), axis=1)
		
		try:
			O2_index = mcmc_par_names.index('Cn_O2')
			smps_2 = 1e6*DMF_O2*np.delete(Sn, O2_index, axis=2)/Sn[:,:,O2_index,np.newaxis]
			flat_smps_2 = 1e6*DMF_O2*np.delete(Sn_flat, O2_index, axis=1)/Sn_flat[:,O2_index,np.newaxis]
		except ValueError:
			smps_2 = 1e6*DMF_O2*Sn/Sn_atm_f['O2']
			flat_smps_2 = 1e6*DMF_O2*Sn_flat/Sn_atm_f['O2']

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

		final_PWV_val = np.median(PWVs_flat, axis=0)
		final_PWV_unc = np.std(PWVs_flat, axis=0)
		#final_PWV_val = funcs.get_PWV(hgt_atm, P_atm, T_atm, Xs_atm_f['H2O'])
		#final_PWV_unc = funcs.get_PWV(hgt_atm, P_atm, T_atm, final_g_uncs[H2O_index]*np.ones(len(hgt_atm)))

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
				f"{final_PWV_val:10.4f} \t {final_PWV_unc:6.4f} \t " +
				f"{logp_med:10.4f} \t {u_logp_med:10.4f}\n")
		txt = open(filename_mcmc_results, 'a')
		txt.write(text)
		txt.close()

		print(f"Progress: {list_science_spectra.index(filename_spectra)+1}/{len(list_science_spectra)}")
