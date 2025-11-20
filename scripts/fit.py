# =====================================================================================
# Basic packages
# =====================================================================================
import numpy as np
from scipy import interpolate
import sys
import copy
from astropy.time import Time
import time

# =====================================================================================
# Scripts
# =====================================================================================
import funcs

def unpack_global_vars(temp_npy_file):
	#lam_obs, spec_obs, R_regrid, R_obs, lam_model, molecs, molecs_cia, Cn_atm, int_cross_secs, int_cia, tau_rayleigh, tau_aerosol, norm_mask, norm_window_sizes, cutoff_mask, stellar_spectra, free_molecs, free_molec_profiles

	global lam_obs
	global spec_obs
	global R_regrid
	global R_obs
	global lam_model
	global molecs
	global molecs_cia
	global Cn_atm
	global int_cross_secs
	global int_cia
	global tau_rayleigh
	global tau_aerosol
	global norm_mask
	global norm_window_sizes
	global cutoff_mask
	global stellar_spectra
	global free_molecs
	global free_molec_profiles

	glob_pars = np.load(temp_npy_file, allow_pickle=True)

	lam_obs = glob_pars[0] 
	spec_obs = glob_pars[1] 
	R_regrid = glob_pars[2] 
	R_obs = glob_pars[3] 
	lam_model = glob_pars[4] 
	molecs = glob_pars[5] 
	molecs_cia = glob_pars[6] 
	Cn_atm = glob_pars[7] 
	int_cross_secs = glob_pars[8] 
	int_cia = glob_pars[9] 
	tau_rayleigh = glob_pars[10] 
	tau_aerosol = glob_pars[11] 
	norm_mask = glob_pars[12] 
	norm_window_sizes = glob_pars[13] 
	cutoff_mask = glob_pars[14] 
	stellar_spectra = glob_pars[15] 
	free_molecs = glob_pars[16]
	free_molec_profiles = glob_pars[17]

def log_likelihood(params):

	## Updating the molecular abundance values
	Cn_atm_tmp = copy.deepcopy(Cn_atm)	## Have to make a copy to avoid the Cn_atm from being changed throughout the MCMC
	for i, ff in enumerate(free_molecs):
		Cn_atm_tmp[ff] = params[i]*free_molec_profiles[i]

	## Generating the model spectra and normalising it
	spec_mod = funcs.spec_handle(lam_model, lam_obs, int_cross_secs, int_cia, tau_rayleigh, tau_aerosol, Cn_atm_tmp, molecs, molecs_cia, R_regrid, R_obs)
	spec_mod_with_stellar_lines = [spec_mod[i]*stellar_spectra[i] for i in range(len(spec_mod))]
	spec_mod_norm = funcs.normalise_spectra(lam_obs, spec_mod_with_stellar_lines, norm_mask, window_sizes=norm_window_sizes)

	## Rule to not use too deep lines in the calculation of the cross correlation and log likelihood,
	## as we believe they might be skewing the results
	spec_obs_cut = funcs.cut_deep_lines(spec_obs, cutoff_mask)
	spec_mod_norm_cut = funcs.cut_deep_lines(spec_mod_norm, cutoff_mask)
	
	## Calculating the cross correlation and the log likelihood
	cc = funcs.cross_correlate(spec_obs_cut, spec_mod_norm_cut)
	log_like = funcs.get_loglike(spec_obs_cut, spec_mod_norm_cut, cc)

	for i in range(np.shape(cc)[0]):
		## Check against NaNs and infinity values
		if cc[i] > 1 or np.isfinite(cc[i]) == False or np.isnan(cc[i]) == True:
			return -np.inf

	return np.sum(log_like)

def log_priors(params):

	## Since our parameters are abundances, it makes no physical sense for them to be negative,
	## so that is the only prior we impose
	for par in params:
		if par < 0.: return -np.inf

	return 0

def log_probability(params):
	lp = log_priors(params)
	
	if not np.isfinite(lp):
		return -np.inf
	if np.isnan(lp):
		print("Log prior is NaN")
		return -np.inf
	
	ll = log_likelihood(params)

	if not np.isfinite(ll):
		return -np.inf
	if np.isnan(ll):
		print("Log likelihood is NaN")
		return -np.inf
	return lp + ll