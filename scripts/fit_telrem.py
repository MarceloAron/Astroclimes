# =====================================================================================
# Basic packages
# =====================================================================================
import numpy as np

# =====================================================================================
# Scripts
# =====================================================================================
import pca_functions

def unpack_global_vars(temp_npy_lam, temp_npy_res):
	## Unpacking the wavelength distribution and residuals to be compared to the model globally
	## to reduce the overheads of parallelism
	global lam
	global residuals
	lam = np.load(temp_npy_lam)
	residuals = np.load(temp_npy_res)

def get_logL(obs,mod,rule):
	no, nf, nx = obs.shape
	logLtemp = 0
	for io in range(no):
		for j in range(nf):
			fV = obs[io,j][rule[io,j]]
			gV = mod[io,j][rule[io,j]]
			nx = len(fV)
			idV = np.ones(nx)
			fV -= np.dot(fV, idV) / nx
			gV -= np.dot(gV, idV) / nx
			sf2 = np.dot(fV,fV) / nx
			sg2 = np.dot(gV,gV) / nx
			R = np.dot(fV,gV) / nx
			logLtemp -= 0.5 * nx * np.log(sf2+sg2-2*R)
	return logLtemp

def get_logL_save_dims(obs,mod,rule):
	no, nf, nx = obs.shape
	logLtemp = np.zeros(shape=(no,nf))
	for io in range(no):
		for j in range(nf):
			fV = obs[io,j][rule[io,j]]
			gV = mod[io,j][rule[io,j]]
			nx = len(fV)
			idV = np.ones(nx)
			fV -= np.dot(fV, idV) / nx
			gV -= np.dot(gV, idV) / nx
			sf2 = np.dot(fV,fV) / nx
			sg2 = np.dot(gV,gV) / nx
			R = np.dot(fV,gV) / nx
			logLtemp[io,j] = -0.5 * nx * np.log(sf2+sg2-2*R)
	return logLtemp

def log_likelihood(params, planet_mod_lam, planet_mod, pf, V_BERV, ip, rule_cube):
	Vsys, Kp, logsf = params

	## planet_mod has +1 added to set the planet level relative to the star, so we need to subtract it 
	## before applying the scaling factor and then add it back before calcualting the model sequence
	## We define a new variable because otherwise planet_mod won't be the same in each iteration, which 
	## shouldn't be the case
	fMod = planet_mod - 1
	fMod *= 10**logsf
	fMod += 1
	modseq = pca_functions.create_model_sequence(planet_mod_lam, fMod, lam, pf, Kp, Vsys, V_BERV, ip)

	log_like = get_logL(modseq, residuals, rule_cube)

	return log_like

def log_priors(params):
	Vsys, Kp, logsf = params
	
	if Vsys < -50 or Vsys > 50 or Kp < 10 or Kp > 300 or logsf < -2 or logsf > 2:
		return -np.inf

	return 0

def log_probability(params, planet_mod_lam, planet_mod, pf, V_BERV, ip, rule_cube):
	lp = log_priors(params)
	
	if not np.isfinite(lp):
		return -np.inf
	if np.isnan(lp):
		print("Log prior is NaN")
		return -np.inf
	
	ll = log_likelihood(params, planet_mod_lam, planet_mod, pf, V_BERV, ip, rule_cube)

	if not np.isfinite(ll):
		return -np.inf
	if np.isnan(ll):
		print("Log likelihood is NaN")
		return -np.inf
	return lp + ll