## Basic packages
import numpy as np
from scipy import interpolate
import copy
from scipy.interpolate import splrep, splev
from astropy.time import Time
import matplotlib.pyplot as plt
from scipy import stats

## Scripts
import funcs

## Some constants
kB = 1.380649*1e-23	# Boltzmann constant, in m^2 kg s^-2 K^-1
R = 8.31446 		# Gas constant, in J mol^-1 K^-1
N_A = 6.022*1e23	# Avogadro constant, in molecules/mole
c = 299792458		# Speed of light, in m/s
#DU = 2.69*1e16		# Dobson unit, in molecules/cm^2
DU = 2.69*1e20		# Dobson unit, in molecules/m^2
h = 6.626e-34        # Planck's constant

## Color-blind friendly colors: blue, orange, green, pink, brown, purple, gray, red and yellow
CBcols= ['#377EB8', '#FF7F00', '#4DAF4A',
		 '#F781BF', '#A65628', '#984EA3',
		 '#999999', '#E41A1C', '#DEDE00']

def filter_bad_obs(nights, instruments, lists_science_spectra, n_orders, bad_spec_orders_idxs, filenames_site_values, filenames_mcmc_results, SNR_lims, airmass_lims, add_bad_obs, main_directory, P, T0, create_SNR_night_log=False, make_obs_conds_plot=False, plots_cols=CBcols):
	## Before carrying on with the telluric line removal, we check for bad observations
	SNRs_all_nights = []
	airmasses_all_nights = []
	rr_bad_obs_all_nights = []
	bad_obs_idxs_all_nights = []
	for ii in range(len(nights)):
		## First, loop over the observational sample to get the measured instrumental SNR from each spectrum
		if create_SNR_night_log:
			txt = open(main_directory+f'night_logs_{nights[ii]}.txt', 'w')
			txt.write('Date obs \t Target \t Exp time \t SNR (mean all orders excluding bad ones)\n')
		SNRs = np.zeros(shape=(len(lists_science_spectra[ii]),n_orders[ii]))
		for filename_spectra in lists_science_spectra[ii]:
			## Getting all of the spectral data data (NaNs and infs are already removed)
			hdr_obs, all_lam_obs, all_spec_obs = funcs.get_spectral_data(filename_spectra, instrument=instruments[ii])
			## These header keys are specific for CARMENES, if you are dealing with another instrument, these might be different
			for jj in range(n_orders[ii]):
				SNRs[lists_science_spectra[ii].index(filename_spectra),jj] = 0.5*(hdr_obs[f'HIERARCH CARACAL FOX SNR {jj}']+hdr_obs[f'HIERARCH CARACAL FOX SNR {jj+n_orders[ii]}'])
			SNRs_nobad_ords = np.delete(SNRs, bad_spec_orders_idxs[ii], axis=1)
			SNRs_nobad_ords_mean = np.mean(SNRs_nobad_ords, axis=1)
			if create_SNR_night_log:
				tgt_name = hdr_obs['OBJECT']
				date_obs = Time(hdr_obs['DATE-OBS'], format='isot', scale='utc')
				exptime = hdr_obs['EXPTIME']
				txt.write(f'{date_obs.value:>20} \t {tgt_name:15} \t {exptime:.2f}s \t {SNRs_nobad_ords_mean[lists_science_spectra[ii].index(filename_spectra)]:<6.2f}\n')
		if create_SNR_night_log:
			txt.close()

		SNRs_all_nights.append(SNRs_nobad_ords_mean)
		## Second, get the airmasses from the site_values.txt file
		airmasses = np.loadtxt(filenames_site_values[ii], usecols=(9), unpack=True)

		## Determine bad observations from their SNR and airmasses, whose limits are given as arguments to the function
		rr_bad_obs = (SNRs_nobad_ords_mean < SNR_lims[ii]) | (airmasses > airmass_lims[ii])

		## Manually adding any other observations deemed bad that were not included in the automatic flagging above
		rr_bad_obs[add_bad_obs[ii]] = True

		bad_obs_idxs = np.arange(len(lists_science_spectra[ii]))[rr_bad_obs]
		print('Night:',nights[ii])
		print('Orders removed:',np.array(bad_spec_orders_idxs[ii])+1)
		print('Spectra removed:',bad_obs_idxs+1)

		airmasses_all_nights.append(airmasses)
		rr_bad_obs_all_nights.append(rr_bad_obs)
		bad_obs_idxs_all_nights.append(bad_obs_idxs)

	if make_obs_conds_plot:
		## Plotting the airmass, relative humidity, PWV and instrument SNR over time, highlighting observations flagged as bad
		fig = plt.figure(figsize=(8,6)) 
		ax1 = plt.subplot2grid((4,1), (0,0), rowspan=1, colspan=1)
		ax2 = plt.subplot2grid((4,1), (1,0), rowspan=1, colspan=1)
		ax3 = plt.subplot2grid((4,1), (2,0), rowspan=1, colspan=1)
		ax4 = plt.subplot2grid((4,1), (3,0), rowspan=1, colspan=1)

		for ii in range(len(nights)):
			BJD, rel_hum = np.loadtxt(filenames_site_values[ii], usecols=(4,7), unpack=True)
			phases = ((BJD - T0)%P)/P
			PWV, u_PWV = np.loadtxt(filenames_mcmc_results[ii], usecols=(21,22), unpack=True)

			ax1.plot(phases, airmasses_all_nights[ii], 'o', mec=plots_cols[ii], mfc='None', label=f'Night {ii+1}')
			ax1.plot(phases[rr_bad_obs_all_nights[ii]], airmasses_all_nights[ii][rr_bad_obs_all_nights[ii]], 'x', c='black')

			ax2.plot(phases, rel_hum, 'o', mec=plots_cols[ii], mfc='None')
			ax2.plot(phases[rr_bad_obs_all_nights[ii]], rel_hum[rr_bad_obs_all_nights[ii]], 'x', c='black')

			ax3.errorbar(phases, PWV, yerr=u_PWV, fmt='o', mec=plots_cols[ii], mfc='None', zorder=4)
			ax3.plot(phases[rr_bad_obs_all_nights[ii]], PWV[rr_bad_obs_all_nights[ii]], 'x', c='black', zorder=5)

			ax4.plot(phases, SNRs_all_nights[ii], 'o', mec=plots_cols[ii], mfc='None')
			ax4.plot(phases[rr_bad_obs_all_nights[ii]], SNRs_all_nights[ii][rr_bad_obs_all_nights[ii]], 'x', c='black')

		ax1.legend(loc='best')
		ax1.axhline(airmass_lims[0], ls='dashed', c='gray')
		ax4.axhline(SNR_lims[0], ls='dashed', c='gray')

		ax1.set_xticks([])
		ax1.set_ylabel('Airmass', fontsize=14)
		ax1.invert_yaxis()

		ax2.set_xticks([])
		ax2.set_ylim(30,102)
		ax2.set_ylabel('Rel. hum. (\%)', fontsize=14)

		ax3.set_xticks([])
		#ax3.set_ylim(1,9)
		ax3.set_ylabel('PWV (mm)', fontsize=14)

		ax4.tick_params(axis='x', labelsize=14)
		ax4.set_ylabel('Inst. SNR', fontsize=14)
		ax4.set_xlabel('Orbital phase', fontsize=14)

		fig.align_ylabels()
		plt.tight_layout()
		plt.subplots_adjust(hspace=0)
		plt.savefig(main_directory+'snr_airmass_hum_pwv_dist.pdf', bbox_inches='tight')
		plt.savefig(main_directory+'snr_airmass_hum_pwv_dist.jpeg', dpi=300, bbox_inches='tight')
		#plt.show()
		plt.close()

	if len(bad_obs_idxs_all_nights) == 1:
		return bad_obs_idxs_all_nights[0]
	else:
		return bad_obs_idxs_all_nights


''' FULL PRINCIPAL COMPONENT ANALYSIS as in Giacobbe et al. (2021) '''
def standardise_data(mat):
	''' 
	Standardise one order of data ('spec') before PCA. In IDL this would 
	happen inside the function PCOMP.pro with the option /STADARDIZE. 
	'''
	nf, nx = mat.shape
	fStd = mat.copy()
	for i in range(nx):
		fStd[:,i] -= np.mean(fStd[:,i])
		# This is the biased stdev (normalised by nx rather than nx-1)
		# It needs changing to match CORRELATE.pro
		fStd[:,i] /= np.std(fStd[:,i])
	
	## Alternative
	## fStd = (fIn - np.mean(fIn, axis=0))/np.std(fIn, axis=0)
	return fStd

def get_eigenv_via_PCA(mat,nc,reconstruct=False):
	''' 
	Run the SVD on the matrix 'mat' (nf, nx) and return the first 'nc' eigen-
	vectors in addition to a column of 1s. If 'reconstruct=True', it also returns
	the reconstructed matrix with the other (nf-nc) components. 
	'''
	nf, nx = mat.shape
	xMat = np.ones((nf,nc+1))
	u, s, vh = np.linalg.svd(mat, full_matrices=False)
	xMat[:,1:] = u[:,0:nc]
	if reconstruct:
		ss = s.copy()
		ss[0:nc] = 0.0
		ww = np.diag(ss)
		res = u.dot(ww.dot(vh))
		return xMat, res
	else:
		return xMat

def linear_regression(X,Y):
	''' 
	Calculate the multi-variate linear regression fit between the matrix of 
	eigenvectors X (nf, ncomponents) and the spectra Y (nf, nx) 
	''' 
	XT = X.T
	term1 = np.linalg.inv(np.dot(XT,X))
	term2 = np.dot(term1,XT)
	beta = np.dot(term2,Y)
	return np.dot(X,beta)

def clean_cube(mat):
	rr = (np.isfinite(mat)) & (mat > 0.05)
	
	rra = rr[0]
	for nobs in range(1,np.shape(rr)[0]):
		rra *= rr[nobs,:]

	## Alternative
	#rra = np.prod(rr, axis=0, dtype=bool)

	clean_mat = mat[:,rra]

	return clean_mat, rra

def run_pca(all_lam, all_spec, nc, generate=True, save_dirname='../', save_filename='all_pca_fit'):
	if generate:
		## Running the PCA
		nor, nobs, nwav = np.shape(all_lam)
		all_pca_fit = np.empty((0,nobs,nwav))
		for n in range(nor):
			## Cleaning up the data of NaN's and inf's for a certain order, as well as points with very low transmission
			fIn, rra = clean_cube(all_spec[n,:])
			fStd = standardise_data(fIn)
			xMat = get_eigenv_via_PCA(fStd, nc)
			pca_fit = linear_regression(xMat, fIn)
			pca_full = np.ones(shape=(nobs,nwav))
			pca_full[:,rra] = pca_fit
			all_pca_fit = np.concatenate((all_pca_fit, pca_full[np.newaxis,:,:]), axis=0)
		np.save(save_dirname+save_filename+f'_nc_{nc}.npy', all_pca_fit)

		return all_pca_fit
	else:
		all_pca_fit = np.load(save_dirname+save_filename+f'_nc_{nc}.npy')

		return all_pca_fit

def get_obs_cube(nord, nwav, filenames_spectra, instrument='CARMENES', generate=True, save_dirname='../', save_lam_filename='all_lam', save_spec_filename='all_spec'):
	if generate:
		all_lam = np.empty((nord,0,nwav))
		all_spec = np.empty((nord,0,nwav))

		for filename_spectra in filenames_spectra:
			## Getting all of the observational spectral data (NaNs and infs not removed)
			hdr_obs, all_lam_obs, all_spec_obs = funcs.get_spectral_data(filename_spectra, instrument=instrument, clean_NaNs=False)
			
			## Concatenating this data to create a "cube" with all of the available observations, dimension (nOrders, nObs, nWavelength)
			all_lam = np.concatenate((all_lam, all_lam_obs[:,np.newaxis]), axis=1)
			all_spec = np.concatenate((all_spec, all_spec_obs[:,np.newaxis]), axis=1)
			
			print(f"Progress: {filenames_spectra.index(filename_spectra)+1}/{len(filenames_spectra)}")

		np.save(save_dirname+save_lam_filename+'.npy', all_lam)
		np.save(save_dirname+save_spec_filename+'.npy', all_spec)


		return all_lam, all_spec
	else:
		all_lam = np.load(save_dirname+save_lam_filename+'.npy')
		all_spec = np.load(save_dirname+save_spec_filename+'.npy')

		return all_lam, all_spec

def get_mod_cube(nord, nwav, filenames_spectra, MCMC_directory, atm_profs_directory, filename_mcmc_results, molecs, molecs_for_cia, instrument='CARMENES', R_instrument=80400, generate=True, save_dirname='../', save_spec_mod_filename='all_spec_mod'):
	if generate:
		ground_CO2, u_ground_CO2, ground_CH4, u_ground_CH4, ground_H2O, u_ground_H2O, ground_O2, u_ground_O2, X_CO2, u_X_CO2, X_CH4, u_X_CH4, X_H2O, u_X_H2O = np.loadtxt(filename_mcmc_results, usecols=(1,2,3,4,5,6,7,8,11,12,13,14,15,16), unpack=True)

		hdr_obs, all_lam_obs, all_spec_obs = funcs.get_spectral_data(filenames_spectra[0], instrument=instrument, clean_NaNs=False)
		
		# Choosing the spectral orders to include in the analysis
		spec_orders = [x for x in range(len(all_lam_obs))]

		## Keeping only the desired orders
		lam_obs = []
		lam_obs = [all_lam_obs[x] for x in spec_orders]

		lam_ranges = [[np.min(lam_obs[0]), np.max(lam_obs[-1])]]

		## Getting the information from the cross-section files
		vel_step = 1.0 					# This is the velocity step (in km/s) to convert the wavelength distribution
		R_regrid = (c*1e-3)/vel_step	# Resolution of our regridded model 
		molecs_cross_secs = funcs.get_cross_secs_dic(molecs, lam_ranges, vel_step)

		all_spec_mod = np.empty((nord,0,nwav))

		for filename_spectra in filenames_spectra:
			MCMC_directory_indiv = MCMC_directory+filename_spectra.split('/')[-1].split('.')[0]+'/'
			
			## Getting all of the observational data (NaNs and infs are already removed)
			hdr_obs, all_lam_obs, all_spec_obs = funcs.get_spectral_data(filename_spectra, instrument=instrument, clean_NaNs=False)
			
			## Unpacking data from FITS file header, tailored per instrument
			date_obs, BJD, hgt_site, P_site, T_site, relhum_site, airmass, V_BERV = funcs.get_header_data(hdr_obs, instrument)
			
			## Getting the atmospheric profile
			hgt_atm, P_atm, T_atm, Xs_atm = funcs.get_ggg2020_atm_profiles(date_obs, atm_profs_directory, molecs, hgt_site)

			## Interpolating the cross sections
			int_cross_secs = funcs.interpolate_cross_secs(P_atm, T_atm, hgt_atm, molecs, molecs_cross_secs, airmass) 

			## Getting the information from the CIA files
			lam = molecs_cross_secs['CO2']['lam']
			molecs_cias = funcs.get_cia_dict(lam, molecs_for_cia)

			## Interpolating the CIA
			int_cias = funcs.interpolate_cia(P_atm, T_atm, hgt_atm, molecs_for_cia, molecs_cias, airmass)

			## Calculating Rayleigh and aerosol scattering
			c_rayleigh = funcs.rayleigh_scatter(lam*1e6, (10**P_site)*1e-2, hgt_site*1e-3)
			c_aerosol = funcs.aerosol_scatter(lam*1e6, airmass)

			## Unpacking the final molecular profiles from that speficid run
			Cn_atm = np.load(MCMC_directory_indiv+'Cn_atm_f.npy', allow_pickle=True)
			Cn_atm = Cn_atm.item() 	## This is necessary due to the way the variable is saved with np.save()

			## Computing the model with the results from the MCMC
			mod_to_obs_spec = funcs.spec_handle(lam, lam_obs, int_cross_secs, int_cias, c_rayleigh, c_aerosol, Cn_atm, molecs, molecs_for_cia, R_regrid, R_instrument)
			all_spec_mod = np.concatenate((all_spec_mod, mod_to_obs_spec[:,np.newaxis]), axis=1)

			print(f"Progress: {filenames_spectra.index(filename_spectra)+1}/{len(filenames_spectra)}")

		all_spec_mod = np.array(all_spec_mod, dtype=float)
		np.save(save_dirname+save_spec_mod_filename+'.npy', all_spec_mod, allow_pickle=True)

		return all_spec_mod
	else:
		all_spec_mod = np.load(save_dirname+save_spec_mod_filename+'.npy', allow_pickle=True)

		return all_spec_mod

def mask_cube(cube, mask_val, rule_cube=None, create_rule_cube=True, deep_line_threshold=0.2):
	'''
	This function either creates a "rule cube" for an array of dimensions (nOrds, nObs, nWave), 
	where the rule is True if points are finite and above a certain threshold (calculated from the median), 
	and then assigns a certain value ("mask_val") to the points that are False. 
	Alternatively, it simply masks the False values, if a rule cube is already provided.
	'''
	cc = copy.deepcopy(cube)
	nor, nobs, nwav = np.shape(cube)
	if create_rule_cube:
		rule_cube = np.ones(shape=(nor,nobs,nwav), dtype=bool)
		for i in range(nor):
			rra = np.ones(shape=(nwav), dtype=bool)
			for j in range(nobs):
				rr = (np.isfinite(cc[i,j,:])) & (cc[i,j,:] > deep_line_threshold*np.nanmedian(cc[i,j,:]))
				rra *= rr
			rule_cube[i,:,:] = rra
	
	for i in range(nor):
		for j in range(nobs):
			cc[i,j,:][~rule_cube[i,j,:]] = mask_val

	return cc, rule_cube

def blackbody(T,wl):
	''' 
	Computing black-body radiation via Planck's function.
	INPUTS:
	- T: temperature in Kelvin
	- wl: wavelengths in meters

	OUTPUTS:
	- blackbody flux (W m-2 m-1), that is power per unit area and unit wavelength 
	'''
	t1 = 2.0*h*c**2/wl**5
	t2 = h*c/(wl*kB*T)
	return t1 / (np.exp(t2) - 1.0)

def get_gaussian_ip(res, wl):
	''' 
	Gaussian instrument profile calculation given target resolution 'res' and
	wavelength vector 'wl'. Wavelengths MUST BE already spaced at constant steps in
	velocity (or in delta(lambda)/lambda), that is constant resolving power R. 
	The code will not check for the spacing, so be aware! 
	'''
	lamPix = 0.5*(wl[1:]+wl[0:-1])
	dlamPix = (wl[1:]-wl[0:-1])
	dvPix = np.mean(dlamPix / lamPix)    # Dropping a factor 'c' (in common with dvFWHM below)
	dvFWHM = 1.0 / res
	fwhmPix = dvFWHM / dvPix
	sigmaPix = 0.5 * fwhmPix / np.sqrt(2.0 * np.log(2.0))
	hker = int(4.0 * fwhmPix)
	nker = hker * 2 + 1
	xker = np.arange(nker) - hker
	yker = np.exp(-0.5*(xker/sigmaPix)**2)
	yker /= yker.sum()

	return yker, hker

def create_model_sequence(wMod, fMod, wData, ph, kp, vSys, vBary, ip):
	''' 
	It creates a model spectral sequence based on a set of radial velocities. 

	INPUTS:

	- wMod, fMod: model wavelength and fluxes (in the right units), MUST BE GRIDDED AT CONSTANT R
	- wData: data wavelength, needs to be in (nOrder,nPixel) format
	- ph: orbital phases obtained as the fractional part of (BJD-T0)/P
	- kp: planet orbital RV semi-amplitude (km/s)
	- vSys: systemic velocity (km/s)
	- vBary: barycentric velocity in the telluric frame, so -VBERV if extracted from CARMENES headers 
	- ip: the gaussian instrumental profile given the wavelengths of the model, if in doubt use get_gaussian_kernel()

	OUTPUT: 

	- specSeq, i.e. the entire spectral sequence for the night with format (nOrders, nFrames, nPixels)
	
	'''
	nf = len(ph)
	no, nx = wData.shape
	specSeq = np.zeros((no,nf,nx))
	# Measure length of IP
	hlen = int((len(ip)-1)/2)
	fModConv = np.convolve(fMod,ip,mode='same')
	cs = splrep(wMod[hlen:-hlen],fModConv[hlen:-hlen],s=0)
	# Compute planet RVs vs time
	rvPl = vSys + vBary + kp*np.sin(2.0 * np.pi * ph)
	for j in range(nf):
		wShift = wData * (1 - rvPl[j]/(c*1e-3))
		fShift = splev(wShift, cs, der=0)
		specSeq[:,j,:] = fShift.copy()
	return specSeq

def create_planet_signal(all_lam, filename_planet_lam, filename_planet_flux, filename_site_values, Res_obs, T_star, R_star, R_planet, scale, T0, P, Kp, Vsys, generate=True, save_dirname='../', save_base_filename='planet_mod'):

	if generate:
		## Unpacking the wavelength and flux of the model planetary signal
		ll = np.loadtxt(filename_planet_lam, skiprows=1) 		# In microns
		ll = ll*1e-6											# Converting from microns to m
		F_planet = np.loadtxt(filename_planet_flux, skiprows=1) # In W m^-2 m^-1

		## Regridding the wavelength so it has constant Delta_lambda/lambda (resolution here needs to be high)
		vel_step = (c*1e-3)/500000
		ll_regrid = funcs.regrid_wavelength(ll, vel_step)
		int_F_planet = interpolate.interp1d(ll, F_planet, bounds_error=False, fill_value=0)
		F_planet_regrid = int_F_planet(ll_regrid)

		## Calculating the Gaussian kernel that will be used to include instrumental broadening
		ip, ttt = get_gaussian_ip(Res_obs, ll_regrid)

		## Calculating the stellar flux with the blackbody radiation equation
		B_star = blackbody(T_star, ll_regrid)
		F_star = np.pi*B_star

		## Calculating the total flux (planet signal is calculated relative to the star, so star is 1, anything extra is planet)
		fMod = 1 + scale*(F_planet_regrid/F_star)*(R_planet/R_star)**2

		## Barycentric Earth radial velocity is needed to shift the planetary signal to the right reference frame
		BJD, V_BERV = np.loadtxt(filename_site_values, usecols=(4,15), unpack=True)
		V_BERV = -V_BERV ## V_BERV is the radial velocity of the Earth with respect to the barycentre, but we want the opposite, hence the signal change

		phases = ((BJD - T0)%P)/P

		modseq = create_model_sequence(ll_regrid, fMod, all_lam[:,0,:], phases, Kp, Vsys, V_BERV, ip)
	
		np.save(save_dirname+'all_'+save_base_filename+'.npy', modseq)
		np.save(save_dirname+''+save_base_filename+'_lam.npy', ll_regrid)
		np.save(save_dirname+''+save_base_filename+'.npy', fMod)

		return modseq, ll_regrid, fMod

	else:
		modseq = np.load(save_dirname+'all_'+save_base_filename+'.npy')
		ll_regrid = np.load(save_dirname+''+save_base_filename+'_lam.npy')
		fMod = np.load(save_dirname+''+save_base_filename+'.npy')

		return modseq, ll_regrid, fMod

def get_fixed_ccf_cube(spec, wl, rv, cs, rule_cube):
	'''
	Computing Pearson's correlation coefficient between a data cube and a model
	spectrum shifted by a range of radial velocities. It produces a cross-correlation
	cube with dimensions (nOrders,nSpectra,nRV).
	Additionally, it computes the log-likelihood value for each point of the CC
	cube via the Brogi & Line (2019) formula.
	INPUTS:
		spec:		the data cube with dimensions (nOrders, nSpectra, nPixels)
		wl  :		the wavelength solution of 'spec' (nOrders, nPixels)
		rv  :		the radial velocity lag vector in km/s
		cs  :		the output of scipy.interpolate.splrep() on the model spectrum to be
					correlated with the data.
		rule_cube:	rule cube to mask out bad points (NaNs and deep lines)
	OUTPUTS:
		tr:	the 'trace' cube (dimensions (nOrders, nSpectra, nRV)) containing the
			Pearson's correlation coefficient per order, spectrum, and RV
		ll:	the log-likelihood cube corresponding to 'tr' 
	'''
	no, nf, nx = spec.shape
	ncc = len(rv)
	tr = np.zeros((no,nf,ncc))
	ll = np.zeros((no,nf,ncc))
	for ic in range(ncc):
		beta = -rv[ic]/(c*1e-3)
		wShift = wl*np.sqrt((1+beta)/(1-beta))
		fShift = splev(wShift,cs,der=0,ext=2)
		for io in range(no):
			for j in range(nf):
				gVec = fShift[io,rule_cube[io,j,:]].copy()
				nx = len(gVec)
				N = float(nx)
				idMat = np.ones(nx)
				gVec -= (gVec @ idMat)/N
				sg2 = (gVec @ gVec)/N
				fVec = spec[io,j,rule_cube[io,j,:]].copy()
				fVec -= (fVec @ idMat)/N
				sf2 = (fVec @ fVec)/N
				R = (fVec @ gVec)/N
				tr[io,j,ic] = R/np.sqrt(sf2 * sg2)
				ll[io,j,ic] = -0.5*N*np.log(sf2+sg2-2*R)
	return tr, ll

def get_velocity_map(ccf,vIn,ph,vObs,vRest,kpVec):
	''' 
	Getting the Vrest-Kp matrix from the CCF cube summed over all the valid 
	orders (not observations?). Note that it can also be used to sum up logL functions in the same way, 
	although the robustness of interpolating logL values has not been tested yet. 
	INPUTS:
	- vIn   :	RV lag vector used to compute the CCFs = 'rv' in get_fixed_ccf_cube()
	- ph    :	vector of orbital phases 
	- vObs  :	combined barycentric (w.r.t. the Earth, not the opposite) and systemic 
				velocities. 
	- vRest :	output rest-frame velocities -- MUST BE NARROWER than vIn to avoid 
				extrapolation of CCF values.
	- kpVec : 	vector of Kp values to be tested, in km/s
	RETURNS: 
	    diag:	the (Vrest-Kp) matrix with y-axis corresponding to 'kpVec' and x-axis 
				to vRest. 
	'''
	nf, ncc = ccf.shape
	nkp = len(kpVec)
	nvs = len(vRest)
	diag = np.zeros((nkp,nvs))
	trail = np.zeros((nf,nvs))
	for ik in range(nkp):
		vPl = vObs + kpVec[ik]*np.sin(2.0 * np.pi * ph)
		for j in range(nf):
			xin = vIn - vPl[j]
			int_fit = interpolate.interp1d(xin,ccf[j,:])
			trail[j,:] = int_fit(vRest)
		diag[ik,] = np.sum(trail,axis=0)
	return diag

def mask_data_post_pca(cube, rule_cube, sigVal=5.0, dynamicSigma=False):
	''' 
	Masking the noisy columns in the data after telluric removal
	'''
	no, nf, nx = cube.shape
	masked = cube.copy()
	rule_cube_masked = rule_cube.copy()
	for io in range(no):
		varVec = np.var(cube[io,],axis=0)
		ivalid = varVec > 0
		nvalid = ivalid.sum()
		if dynamicSigma:
			sigVal = stats.norm.isf(0.5/nvalid)
		medVar = np.median(varVec[ivalid])        
		masked[io,:,varVec > sigVal*medVar] = 0.0
		rule_cube_masked[io,:,varVec > sigVal*medVar] = False
	return masked, rule_cube_masked

def model_reprocessing(lam, spec, planet_model):
	a