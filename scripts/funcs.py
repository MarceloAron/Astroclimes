# =====================================================================================
# Basic packages
# =====================================================================================
import numpy as np
import astropy.io.fits as pyfits
import glob
from astropy.time import Time
import pandas as pd
from scipy import interpolate
import h5py
import emcee
import copy
#import time

# =====================================================================================
# Scripts
# =====================================================================================
import plots
import gen_funcs

## Some constants
kB = 1.380649*1e-23	# Boltzmann constant, in m^2 kg s^-2 K^-1
R = 8.31446 		# Gas constant, in J mol^-1 K^-1
N_A = 6.022*1e23	# Avogadro constant, in molecules/mole
c = 299792458		# Speed of light, in m/s
#DU = 2.69*1e16		# Dobson unit, in molecules/cm^2
DU = 2.69*1e20		# Dobson unit, in molecules/m^2

#home_directory = '/Users/marceloaronkeniger/PhD/thesis_work/Astroclimes/'
#home_directory = '/home/marceloaron/MarceloAron/PhD/thesis_work/Astroclimes/'
#home_directory = '/storage/astro2/phrgmq/Astroclimes/'
home_directory = '/media/marceloaron/New Volume/PhD/thesis_work/'

plots_directory = home_directory+'plots/'

def get_spectral_data(filename, instrument='CARMENES', clean_NaNs=True):
	'''
	Unpack spectral data from a given instrument. By default, this function cleans the spectra of NaNs and infs
	'''
	if instrument == 'CARMENES':
		f = pyfits.open(filename)
		hdr = f[0].header
		data1 = f[1].data 	# This is the flux
		data4 = f[4].data 	# This is the wavelength, in Angstrom

		if not clean_NaNs:
			return hdr, data4*1e-10, data1
		else:
			## Getting rid of NaNs and infs and converting the wavelength from Angstrom to m 
			all_lam_carmenes = []
			all_spec_carmenes = []
			for n in range(data4.shape[0]):
				rr = np.isfinite(data4[n,:]) & np.isfinite(data1[n,:])
				all_spec_carmenes.append(data1[n,:][rr])
				all_lam_carmenes.append(data4[n,:][rr]*1e-10)

			return hdr, all_lam_carmenes, all_spec_carmenes
	
	elif instrument == 'NIRPS':
		f = pyfits.open(filename)
		hdr = f[0].header
		data = f[1].data
		lam = data['WAVE'][0]*1e-10 	# This is the wavelength, converted from  Angstrom to m
		try:
			flux = data['FLUX'][0]
		except KeyError:
			flux = data['FLUX_EL'][0]

		## Getting rid of NaNs and infs
		if not clean_NaNs:
			return hdr, lam, data
		else:
			rr = np.isfinite(lam) & np.isfinite(flux)
			lam_clean = lam[rr]
			flux_clean = flux[rr]

			return hdr, lam_clean, flux_clean

	else:
		raise ValueError('Instrument name not included in list. Please pick another one. Current available options are: CARMENES, NIRPS')

def get_header_data(hdr, instrument):
	ESO_instruments = ['ESPRESSO', 'NIRPS']
	if instrument == 'CARMENES':
		date_obs = Time(hdr['DATE-OBS'], format='isot', scale='utc') 		## UTC at start of exposure
		BJD = hdr['HIERARCH CARACAL BJD'] + 2400000							# Barycentric Julian Day at middle of exposure
		P_site = hdr['HIERARCH CAHA GEN AMBI PRESSURE'] 					# Pressure in hPa
		P_site = np.log10(P_site*100) 										# Converting the pressure from hPa to log(Pa)
		T_site = hdr['HIERARCH CAHA GEN AMBI TEMPERATURE']+273.15 			# Temperature in ºC, converted to K
		hgt_site = hdr['HIERARCH CAHA TEL GEOELEV']							# Height above sea level, in m
		relhum_site = hdr['HIERARCH CAHA GEN AMBI RHUM']					# Relative humidity, in %
		relhum_site_ppmv = relhum2ppmv(T_site, P_site, relhum_site)			# Converting relative humidity to ppmv
		airmass = hdr['AIRMASS'] 											# Airmass at the start of exposure
		try:
			V_BERV = hdr['HIERARCH CARACAL BERV']					 		# Barycentric correction (in km/s), from BarCor
		except KeyError:
			print('V_BERV not available, setting dummy value!')
			V_BERV = -999.99
	elif instrument in ESO_instruments:
		date_obs = Time(hdr['DATE-OBS'], format='isot', scale='utc') 									## UTC at start of exposure
		BJD = hdr['HIERARCH ESO QC BJD'] 																# Barycentric Julian Day (TDB) at middle of exposure
		P_site = 0.5*(hdr['HIERARCH ESO TEL AMBI PRES START']+hdr['HIERARCH ESO TEL AMBI PRES END']) 	# Pressure in hPa (mean between start and end of exposure)
		P_site = np.log10(P_site*100) 																	# Converting the pressure from hPa to log(Pa)
		T_site = hdr['HIERARCH ESO TEL AMBI TEMP']+273.15 												# Temperature in ºC, converted to K
		hgt_site = hdr['HIERARCH ESO TEL GEOELEV']														# Height above sea level, in m
		relhum_site = hdr['HIERARCH ESO TEL AMBI RHUM']													# Relative humidity, in %
		relhum_site_ppmv = relhum2ppmv(T_site, P_site, relhum_site)										# Converting relative humidity to ppmv
		airmass = 0.5*(hdr['HIERARCH ESO TEL AIRM END']+hdr['HIERARCH ESO TEL AIRM START'])				# Mean of airmass at start and end of exposure
		try:
			V_BERV = hdr['HIERARCH ESO QC BERV'] 														# Barycentric correction, in km/s
		except KeyError:
			print('V_BERV not available, setting dummy value!')
			V_BERV = -999.99
	else:
		raise ValueError('Instrument not included in list. Please pick another one. Current available options are: CARMENES, NIRPS, ESPRESSO')

	return date_obs, BJD, hgt_site, P_site, T_site, relhum_site_ppmv, airmass, V_BERV

def create_site_values_file(files, instrument='CARMENES', dirname='../', filename='site_values.txt'):
	ESO_instruments = ['ESPRESSO', 'NIRPS']
	txt = open(dirname+filename, 'w')
	txt.write('## Target name \t RA \t DEC \t Obs. start date \t BJD at mid-exp \t P(hPa) \t T (K) \t Humidity (%) \t Humidity (ppm) \t Airmass \t Azimuth (deg) \t Elevation (deg) \t Obs. latitude (deg) \t Obs. longitude (deg) \t Obs. altitude (m) \t Barycentric velocity (km/s) \n')
	for i,fname in enumerate(files):
		f = pyfits.open(fname)
		hdr = f[0].header
		if instrument == 'CARMENES':
			tgt_name = hdr['OBJECT'].replace(" ", "") 							# Target name, erasing spaces
			if len(tgt_name) == 0: tgt_name = 'None' 							# Setting target name to "None" if not present
			RA = hdr['RA'] 														# Right ascension, in degrees (J2000)
			DEC = hdr['DEC'] 													# Declination, in degrees (J2000)
			date_obs = Time(hdr['DATE-OBS'], format='isot', scale='utc') 		# UTC at start of exposure
			BJD = hdr['HIERARCH CARACAL BJD'] + 2400000							# Barycentric Julian Day at (presumably) middle of exposure, from BarCor
			temperature = hdr['HIERARCH CAHA GEN AMBI TEMPERATURE'] + 273.15 	# Temperature in ºC, converted to K
			pressure = hdr['HIERARCH CAHA GEN AMBI PRESSURE']					# Pressure in hPa
			rhum = hdr['HIERARCH CAHA GEN AMBI RHUM']							# Relative humidity, in %
			rhum_ppmv = relhum2ppmv(temperature, pressure, rhum)				# Converting relative humidity to ppmv
			airmass = hdr['AIRMASS'] 											# Airmass at the start of exposure
			AZ = hdr['HIERARCH CAHA TEL POS AZ_START'] 							# Telescope azimuth at exposure start, in degrees
			ELEV = hdr['HIERARCH CAHA TEL POS EL_START'] 						# Telescope elevation at exposure start, in degrees
			LAT = hdr['HIERARCH CAHA TEL GEOLAT'] 								# Geographical latitude, in degrees
			LON = hdr['HIERARCH CAHA TEL GEOLON']								# Geographical longitude, in degrees
			GEOELEV = hdr['HIERARCH CAHA TEL GEOELEV'] 							# Height above sea level, in m
			try:
				V_BERV = hdr['HIERARCH CARACAL BERV']					 		# Barycentric correction (in km/s), from BarCor
			except KeyError:
				print('V_BERV not available, setting dummy value!')
				V_BERV = -999.99
		elif instrument in ESO_instruments:
			tgt_name = hdr['OBJECT'].replace(" ", "") 														# Target name, erasing spaces
			if len(tgt_name) == 0: tgt_name = 'None' 														# Setting target name to "None" if not present
			RA = hdr['RA'] 																					# Right ascension, in degrees (J2000)
			DEC = hdr['DEC'] 																				# Declination, in degrees (J2000)
			date_obs = Time(hdr['DATE-OBS'], format='isot', scale='utc') 									# UTC at start of exposure
			BJD = hdr['HIERARCH ESO QC BJD']																# Barycentric Julian Day (TDB) at middle of exposure
			temperature = hdr['HIERARCH ESO TEL AMBI TEMP']+273.15 											# Temperature in ºC, converted to K
			pressure = 0.5*(hdr['HIERARCH ESO TEL AMBI PRES START']+hdr['HIERARCH ESO TEL AMBI PRES END']) 	# Pressure in hPa
			rhum = hdr['HIERARCH ESO TEL AMBI RHUM']														# Relative humidity, in %
			rhum_ppmv = relhum2ppmv(temperature, pressure, rhum)											# Converting relative humidity to ppmv
			airmass = 0.5*(hdr['HIERARCH ESO AIRM END']+hdr['HIERARCH ESO AIRM START'])						# Mean of airmass at start and end of exposure
			AZ = hdr['HIERARCH ESO TEL AZ'] 																# Telescope azimuth at exposure start, in degrees (S=0, W=90)
			ELEV = hdr['HIERARCH ESO TEL ALT'] 																# Telescope elevation at exposure start, in degrees
			LAT = hdr['HIERARCH ESO TEL GEOLAT'] 															# Geographical latitude, in degrees (+=North)
			LON = hdr['HIERARCH ESO TEL GEOLON']															# Geographical longitude, in degrees (+=East)
			GEOELEV = hdr['HIERARCH ESO TEL GEOELEV'] 														# Height above sea level, in m
			try:
				V_BERV = hdr['HIERARCH ESO QC BERV']					 									# Barycentric correction (in km/s)
			except KeyError:
				print('V_BERV not available, setting dummy value!')
				V_BERV = -999.99
		else:
			raise ValueError('Instrument name not included in list. Please pick another one. Current available options are: CARMENES')

		txt.write(f'{tgt_name:15s} \t {RA:15.4f} \t {DEC:15.4f} \t {date_obs.value:19s} \t {BJD:15.7f} \t {pressure:15.4f} \t {temperature:15.4f} \t {rhum:15.4f} \t {rhum_ppmv:15.4f} \t {airmass:15.4f} \t {AZ:15.4f} \t {ELEV:15.4f} \t {LAT:15.4f} \t {LON:15.4f} \t {GEOELEV:15.4f} \t {V_BERV:15.4f}\n')
		print(f"Progres:{i+1}/{len(files)}")
	txt.close()

def get_skymodel_spectra(filename, lam_obs, spec_type='emission'):
	if type(lam_obs) == list:
		spec_skymodel = []
		for i in range(len(lam_obs)):
			if filename[-4:] == '.txt':
				lam, spec = np.loadtxt(filename, unpack=True)
				lam = lam*1e-9 # Converting from nm to um
				spec_skymodel.append(np.interp(lam_obs[i], lam, spec))
			elif filename[-5:] == '.fits':
				f = pyfits.open(filename)
				data = f[1].data
				lam = data['lam']
				lam = lam*1e-9 # Converting from nm to um
				if spec_type == 'emission':
					spec = data['flux']
				elif spec_type == 'absorption':
					spec = data['trans']
				spec_skymodel.append(np.interp(lam_obs[i], lam, spec))
	else:
		f = pyfits.open(filename)
		data = f[1].data
		lam = data['lam']
		lam = lam*1e-9 # Converting from nm to um
		if spec_type == 'emission':
			spec = data['flux']
		elif spec_type == 'absorption':
			spec = data['trans']
		spec_skymodel = np.interp(lam_obs, lam, spec)

	return spec_skymodel

def get_phoenix_spectra(filename_lam, filename_spec, vel_step):
	lam = pyfits.open(filename_lam)[0].data*1e-10 				# Converting wavelength from Angstrom to m
	spec = pyfits.open(filename_spec)[0].data*1e-7/(1e-2) 		# Flux in erg s^-1 cm^-1 converted to W/m

	## Regridding the wavelength distribution to have a constant Delta_lambda/lambda
	lam_regrid = regrid_wavelength(lam, vel_step)

	## Interpolating the flux to this new wavelength grid
	interp = interpolate.interp1d(lam, spec)
	spec_regrid = interp(lam_regrid)

	return lam_regrid, spec_regrid

def rot_profile(Dv, vsini, eps):
	'''
	Calculating the rotation profile given by Equation (18.14) in Gray (2005)
	eps = linear limb darkening coefficient when the continuum intensity is given by Ic = I0c*(1 - eps + eps*cos(theta)),
	where theta is the angle between the line of sight and the incoming flux 
	theta = 0 deg corresponds to the disk center and theta = 90 deg corresponds to the limb
	Claret et al (1995) LD coefficients for R, I, J, H, K bands for log(g) = 4.5 cm/s^2 and Teff = 6500 K: 0.516, 0.440, 0.328, 0.263, 0.231
	'''
	c1 = 2*(1-eps)/(np.pi*vsini*(1-eps/3))
	c2 = (np.pi*eps)/(2*np.pi*vsini*(1-eps/3))
	G = c1*(1 - (Dv/vsini)**2)**(1/2) + c2*(1 - (Dv/vsini)**2)
	G[np.where(np.isnan(G))] = 0

	return G

def broaden_phoenix_spectra(lam, spec, vel_step=1, vsini=0, eps=0, plot_fig=False, dirname=plots_directory, filename='phoenix_broad'):
	Dv = np.arange(-int(vsini)-5,int(vsini)+5+vel_step, step=vel_step)
	Dv_fine = np.empty(0)
	for i in range(len(Dv)):
		nsteps = int((Dv[i]+vel_step/2 - (Dv[i]-vel_step/2))/(vel_step/10) + 1)
		fine = np.linspace(Dv[i]-vel_step/2, Dv[i]+vel_step/2, num=nsteps)
		Dv_fine = np.concatenate((Dv_fine, fine))
	G_fine = rot_profile(Dv_fine, vsini, eps)
	G = np.mean(np.reshape(G_fine, (len(Dv),int(len(Dv_fine)/len(Dv)))), axis=1)
	num = int(len(G)/2)
	lam_broad = lam[num:-num]
	spec_broad = np.convolve(spec, G, mode='same')[num:-num]

	return lam_broad, spec_broad

def separate_spectra(lam, spec, lam_ranges):
	lams = []
	specs = []
	for lam_range in lam_ranges:
		rr = (lam >= lam_range[0]) & (lam <= lam_range[1])
		l, s = lam[rr], spec[rr]
		lams.append(l)
		specs.append(s)

	return lams, specs

def relhum2ppmv(T, pressure, relhum):
	'''
	Function to convert relative humidity in % to ppmv using the approach from Murphy & Koop (2005)
	Credit: David Armstrong
	T must be in K
	Pressure must be in mbar or hPa
	Relative humidity must be in %
	'''

	logT = np.log(T)

	## Convert GDAS relative humidity to ppmv
	## Vapour pressure over water, equation (10) from Murphy & Koop (2005), valid for 123K < T < 332K
	log_ew = 54.842763 - 6763.22 / T - 4.210 * logT + 0.000367 * T + np.tanh(0.0415 * (T - 218.8)) * (53.878 - 1331.22 / T - 9.44523 * logT + 0.014025 * T)
	## Vapour pressure over hexagonal ixe, equation (7) from Murphy & Koop (2005), valid for T > 110K
	log_ei = 9.550426 - 5723.265 / T + 3.53068 * logT - 0.00728332 * T

	## Vapour pressure
	log_e = np.min(np.array([log_ew,log_ei]),axis=0)
	p_sat = np.exp(log_e) / 100.

	p_h2o = relhum / 100. * p_sat

	ppmv = p_h2o / pressure * 1e6

	return ppmv

def get_CAMS_data(filename, keys):
	f = h5py.File(filename)

	T0 = Time('1970-01-01T00:00:00.000', format='isot', scale='utc').jd
	tts = f['valid_time'][:] 				## Time in seconds since 1970-01-01 00:00:00 UTC
	tt = tts/(60*60*24) + T0 				## Converting time to days and then adding it to T0 to get the date of observation, in JD
	date_obs = Time(tt, format='jd').isot 	## Converting time to isot format

	lat = f['latitude'][:] 		## Latitude, in degrees
	lon = f['longitude'][:] 	## Longitude, in degrees

	pars = []
	for key in keys:
		pars.append(f[key][:])
	#xCO2 = f['tcco2'][:] 		# Column averaged CO2 dry air mole fraction, in ppm
	#xCH4 = f['tcch4'][:] 		# Column averaged CH4 dry air mole fraction, in ppm

	return f, date_obs, lat, lon, pars#xCO2, xCH4

def get_ggg2020_atm_profiles(date_obs, dirname, molecs, hgt_bot=0.):
	'''
	Function to unpack the data from the atmospheric profiles generated by the GGG2020 algorithm
	Data is stored in .map, .mod and .vmr files
	'''
	
	## First I identify which set of data is closer in time to the observation
	ff = glob.glob(dirname+'maps-vertical/*.map')
	ff.sort()
	profs_times = []
	for i,f in enumerate(ff):
		t = f.split('/')[-1].split('_')[-1][:-4]
		year, month, day, hour = t[:4], t[4:6], t[6:8], t[8:10]
		profs_times.append(year+'-'+month+'-'+day+'T'+hour+':00:00.000')

	profs_times_jd = Time(profs_times, format='isot').jd
	date_obs_jd = Time(date_obs, format='isot').jd

	idx = np.argmin(abs(profs_times_jd-date_obs_jd))
	ds = ff[idx].split('/')[-1].split('_')[-1][:-4]

	## Then I select the appropriate files
	map_file = glob.glob(dirname+f'maps-vertical/*{ds}*.map')[0]
	#mod_file = glob.glob(dirname+f'vertical/*{ds}*.mod')[0] 		## Not needed
	vmr_file = glob.glob(dirname+f'vmrs-vertical/*{ds}*.vmr')[0]

	## Unpacking height, temperature and pressure from the .map file
	## Depending on how the .vmr and .map files are computed, they may or may not have an extra line stating the CO source,
	## so you may have to change the number of skipped rows when unpacking their values
	hgt_atm, T_atm, P_atm = np.loadtxt(map_file, skiprows=12, usecols=(0,1,2), unpack=True, delimiter=',')
	hgt_atm = hgt_atm*1e3 		# Converting from km to m
	P_atm = np.log10(P_atm*1e2)	# Converting from hPa to log(Pa)

	## Unpacking the molecular abundances as a function of height from the .vmr file and putting them in a dictionary
	## These are given in volume mixing ratios, which I believe here are equivalent to DRY air mole fractions
	Xs_atm = {}
	vmr_header = np.loadtxt(vmr_file, skiprows=8, max_rows=1, dtype=str, unpack=True).tolist()
	for molec in molecs:
		Xs_atm[molec] = 1e6*np.loadtxt(vmr_file, skiprows=9, usecols=(vmr_header.index(molec)), unpack=True)

	## If the bottom height is above ground level (that is, if the data comes from a certain altitude),
	## Then we interpolate the profile to this level and set it as the bottom of the profile, removing any
	## points below
	if hgt_bot > 0.:
		P_bot = np.interp(hgt_bot, hgt_atm, P_atm)
		T_bot = np.interp(hgt_bot, hgt_atm, T_atm)
		Xs_bot = {}
		for X in Xs_atm:
			Xs_bot[X] = np.interp(hgt_bot, hgt_atm, Xs_atm[X])

		rr = hgt_atm > hgt_bot
		hgt_atm = np.append(hgt_bot, hgt_atm[rr])
		P_atm = np.append(P_bot, P_atm[rr])
		T_atm = np.append(T_bot, T_atm[rr])
		for X in Xs_atm:
			Xs_atm[X] = np.append(Xs_bot[X], Xs_atm[X][rr])

	## Inverting the variables so that the first entry in the array is at the top of the atmosphere
	hgt_atm = hgt_atm[::-1]
	P_atm = P_atm[::-1]
	T_atm = T_atm[::-1]
	for X in Xs_atm:
		Xs_atm[X] = Xs_atm[X][::-1]

	return hgt_atm, P_atm, T_atm, Xs_atm

def interp_abundances(Xs_mipas, hgt_mipas, hgt_gdas):
	Xs_gdas = {}
	for X in Xs_mipas:
		Xs_gdas[X] = np.interp(hgt_gdas, hgt_mipas[::-1], Xs_mipas[X][::-1])
		## Previous implementation using scipy.interpolate.interp1d
		#xinterp = interpolate.interp1d(hgt_mipas, Xs_mipas[X], fill_value='extrapolate')
		#Xs_gdas[X] = interp(hgt_gdas)

	return Xs_gdas

def regrid_wavelength(lam, vel_step):
	## Converting the wavelength distribution to have a constant Delta_lambda/lambda
	## (which I think is the same as having a constant velocity step)
	## vel_step has to be given in km/s
	r = vel_step/(c*1e-3)
	ndata = int(np.log(lam.max()/lam.min())/np.log(1+r))
	w = lam.min()*(1+r)**np.arange(ndata)

	return w

def get_cross_secs_dic(molecs, lam_ranges, vel_step):
	molecs_cross_secs = {}

	for molec in molecs:
		f = h5py.File(home_directory+'auxiliary_files/opacities/lowT/'+molec.lower()+'.hdf5')
		cross_secs = f[molec.lower()]['cross_sec'][:]	# Cross section, in log(m^2/molecule)
		P = f[molec.lower()]['P'][:]					# Pressure, in log(Pa)
		T = f[molec.lower()]['T'][:]					# Temperature, in K
		lam = f[molec.lower()]['lam'][:]				# Wavelength, in m

		## Reducing our sample to include only the wavelengths associated with the observational spectra 
		## (with some extra on the edges to avoid issues with the convolution later)
		new_lam = np.empty(0)
		new_cross_secs = np.empty((0,11,11))
		regrid_lam = np.empty(0)
		for lam_range in lam_ranges:
			#rr = (lam >= 0.99*lam_range[0]) & (lam <= 1.01*lam_range[1])
			rr = (lam >= 0.999*lam_range[0]) & (lam <= 1.001*lam_range[1])	# Had to lower the range because consecutive orders were overlapping
			new_lam = np.append(new_lam, lam[rr])
			new_cross_secs = np.append(new_cross_secs, cross_secs[rr], axis=0)
			w = regrid_wavelength(lam[rr], vel_step)
			regrid_lam = np.append(regrid_lam, w)
		
		int_cross_secs = np.zeros(shape=(len(regrid_lam), len(P), len(T)))
		for i in range(len(P)):
			for j in range(len(T)):
				int_cross_sec = interpolate.interpn((new_lam, P, T), new_cross_secs, (regrid_lam, P[i], T[j]), bounds_error=True, fill_value=None)
				int_cross_secs[:,i,j] = int_cross_sec
		
		molecs_cross_secs[molec] = { 'cross_sec': int_cross_secs,
									 'lam': regrid_lam,
									 'P': P,
									 'T': T}
		
	return molecs_cross_secs

def interpolate_cross_secs(P_atm, T_atm, hgt_atm, molecs, molecs_cross_secs, airmass):
	tau_each_X = dict.fromkeys(molecs, {})
	cos_theta = 1/airmass
	Pint = 0.5*P_atm[1:] + 0.5*P_atm[:-1]
	Tint = 0.5*(T_atm[1:] + T_atm[:-1])
	r = abs(hgt_atm[1:]-hgt_atm[:-1])/cos_theta
	for molec in molecs:
		P = molecs_cross_secs[molec]['P']
		T = molecs_cross_secs[molec]['T']
		lam = molecs_cross_secs[molec]['lam']
		cross_sec = molecs_cross_secs[molec]['cross_sec']

		tau_each_layer = np.zeros(shape=(len(P_atm)-1,len(lam)))
		for i in range(len(P_atm)-1):
			int_cross_sec = interpolate.interpn((lam,P,T), cross_sec, (lam,Pint[i],Tint[i]), bounds_error=True)
			tau_each_layer[i,:] = r[i]*10**(int_cross_sec)
		tau_each_X[molec] = np.array(tau_each_layer)

	return tau_each_X

def get_spectra(molecs, Xs_atm, tau_each_X):
	tau = 0
	for molec in molecs:
		X_atm = Xs_atm[molec]
		X = 0.5*(X_atm[1:] + X_atm[:-1])
		tau += np.einsum('i,ij->j',X,tau_each_X[molec])
	
	return tau

def get_cia_dict(lam, molecs):
	molecs_cias = {}
	for molec in molecs:
		fs = glob.glob(home_directory+f'auxiliary_files/opacities/cia/{molec}/*')
		fs.sort()
		molecs_cias[molec] = {}
		molecs_cias[molec]['Ts'] = {}
		for f in fs:
			try:
				wvn, k = np.loadtxt(f, usecols=(0,1), unpack=True)
				ll = 1e-2*(1/wvn)
				int_cias = interpolate.interp1d(ll, k, bounds_error=False, fill_value=0)
				nk = int_cias(lam)
				molecs_cias[molec]['noT'] = {	'lam': lam,
												'cia': nk}
			except IsADirectoryError:
				nf = glob.glob(f+'/*')
				nf.sort()
				Ts = np.zeros(len(nf))
				nks = np.zeros(shape=(len(lam), len(Ts)))
				for i,n in enumerate(nf):
					wvn, k = np.loadtxt(n, usecols=(0,1), skiprows=1, unpack=True)
					ll = 1e-2*(1/wvn)
					int_cias = interpolate.interp1d(ll, k, bounds_error=False, fill_value=0)
					nk = int_cias(lam)
					nks[:,i] = nk
					Ts[i] = np.loadtxt(n, usecols=(4), max_rows=1, unpack=True)
				lamrange = f.split('/')[-1]
				molecs_cias[molec]['Ts'][lamrange] = {	'lam': lam,
														'T': Ts,
														'cia': nks}

	return molecs_cias

def interpolate_cia(P_atm, T_atm, hgt_atm, molecs, molecs_cias, airmass):
	tau_each_X = dict.fromkeys(molecs, {})
	cos_theta = 1/airmass
	Tint = 0.5*(T_atm[1:] + T_atm[:-1])
	r = abs(hgt_atm[1:]-hgt_atm[:-1])/cos_theta
	for molec in molecs:
		lam = molecs_cias[molec]['noT']['lam']
		for key in molecs_cias[molec]:
			if key == 'noT':
				lam = molecs_cias[molec]['noT']['lam']
				cia = molecs_cias[molec]['noT']['cia']

				tau_each_X_noT = cia[np.newaxis, :]*(r[:, np.newaxis]*1e2)*(1e-6)**2
			elif key == 'Ts':
				tau_each_X_T = np.zeros(shape=(len(P_atm)-1,len(lam)))
				for skey in molecs_cias[molec][key]:
					lam = molecs_cias[molec]['Ts'][skey]['lam']
					T = molecs_cias[molec]['Ts'][skey]['T']
					cia = molecs_cias[molec]['Ts'][skey]['cia']

					tau_each_layer = np.zeros(shape=(len(P_atm)-1,len(lam)))
					for i in range(len(P_atm)-1):
						try:
							int_cia = interpolate.interpn((lam,T), cia, (lam,Tint[i]), bounds_error=True)
						except ValueError:
							T_idx = np.argmin(abs(T - Tint[i]))
							int_cia = cia[:,T_idx]
						tau_each_layer[i,:] = int_cia*(r[i]*1e2)*(1e-6)**2
					tau_each_X_T += np.array(tau_each_layer)
		tau_each_X[molec] = tau_each_X_noT + tau_each_X_T

	return tau_each_X

def get_spectra_cia(molecs, Xs_atm, tau_each_X):
	tau = 0
	for molec in molecs:
		if molec.split('-')[1] == 'air':
			X_atm_1 = Xs_atm[molec.split('-')[0]]
			X_1 = 0.5*(X_atm_1[1:] + X_atm_1[:-1])
			X_atm_2 = Xs_atm['O2'] + Xs_atm['N2']# + Xs_atm['add']
			X_2 = 0.5*(X_atm_2[1:] + X_atm_2[:-1])
		else:
			X_atm_1 = Xs_atm[molec.split('-')[0]]
			X_1 = 0.5*(X_atm_1[1:] + X_atm_1[:-1])
			X_atm_2 = Xs_atm[molec.split('-')[1]]
			X_2 = 0.5*(X_atm_2[1:] + X_atm_2[:-1])
		tau += np.einsum('i,i,ij->j',X_1, X_2, tau_each_X[molec])
	
	return tau

def gauss(mu, sigma, fwhm):
	#x = np.arange(-5*round(fwhm), 5*round(fwhm)+1)
	x = np.arange(-int(4*fwhm), int(4*fwhm)+1)	# To make it the same as Matteo's implementation
	gauss_kernel = np.exp(-(x-mu)**2/(2*sigma**2))

	return gauss_kernel/np.sum(gauss_kernel)

def inst_broaden_gauss(lam, spec, R, R_instrument):
	'''
	Function to broaden a spectra to include instrumental broadeing based on a Gaussian instrumental profile
	Since this involves a convolution, the edges of the wavelength distribution and resulting spectra are trimmed
	'''
	g_fwhm = R/R_instrument
	sigma = g_fwhm/(2*np.sqrt(2*np.log(2)))
	gg = gauss(0, sigma, g_fwhm)
	gnum = int(len(gg)/2)
	broad_spec = np.convolve(spec, gg, mode='same')

	return lam[gnum:-gnum], broad_spec[gnum:-gnum]

def spec_handle(lam, lam_obs, int_cross_secs, int_cias, tau_rayleigh, tau_aerosol, Cn_atm, molecs, molecs_for_cia, R_regrid, R_obs):
	'''
	Function to normalise the spectra (includes instrumental broadening and interpolation to observed wavelength distribution)
	'''

	## Generating the spectra
	tau_lbl = get_spectra(molecs, Cn_atm, int_cross_secs)

	## Generating the CIA spectra
	tau_cia = get_spectra_cia(molecs_for_cia, Cn_atm, int_cias)

	## Plugging all of the effects together
	full_spec = np.exp(-(tau_lbl + tau_cia + tau_rayleigh + tau_aerosol*np.log(10)))

	## Broadening the model spectra to include instrumental brodening based on the resolution of the spectrograph used
	lam, broad_spec = inst_broaden_gauss(lam, full_spec, R_regrid, R_obs)

	## Interpolating the model spectra so it is in the same wavelength distribution as the observational spectra
	mod_to_obs_spl = interpolate.splrep(lam, broad_spec, k=1)
	mod_to_obs_spec = [interpolate.splev(lam_obs[i], mod_to_obs_spl) for i in range(len(lam_obs))]

	return np.array(mod_to_obs_spec, dtype=object)

def get_normalisation_mask(lam, lam_obs, int_cross_secs, int_cias, tau_rayleigh, tau_aerosol, Cn_atm, molecs, molecs_for_cia, R_regrid, R_obs, select_CO2_ground, select_H2O_ground, window_sizes, stellar_spectra, orders=[], plot_bool=False, filename='norm_mask', dirname=plots_directory):

	## First I'll change the H2O and CO2 abundances to arbitrary values so that the lines are not too deep
	## and the continuum baseline is not lowered by the vast number of lines
	Cn_CO2_profile = Cn_atm['CO2']/Cn_atm['CO2'][-1]
	Cn_H2O_profile = Cn_atm['H2O']/Cn_atm['H2O'][-1]
	Cn_atm_dummy = copy.deepcopy(Cn_atm)
	Cn_atm_dummy['CO2'] = Cn_CO2_profile*select_CO2_ground
	Cn_atm_dummy['H2O'] = Cn_H2O_profile*select_H2O_ground

	## Now I'll generate the model spectra (this function already broadens and interpolates the spectra
	## to the observed wavelength distribution)
	mod_to_obs_spec = spec_handle(lam, lam_obs, int_cross_secs, int_cias, tau_rayleigh, tau_aerosol, Cn_atm_dummy, molecs, molecs_for_cia, R_regrid, R_obs)

	## Combining the model telluric spectra with the model stellar spectra
	mod_to_obs_spec = [mod_to_obs_spec[i]*stellar_spectra[i] for i in range(len(mod_to_obs_spec))]

	## Normalising the model spectra
	norm_mask = []
	for i in range(len(mod_to_obs_spec)):
		win = int(2*len(mod_to_obs_spec[i])*window_sizes[i] + 1)
		filt = pd.Series(mod_to_obs_spec[i]).rolling(window=win, min_periods=1, center=True).median()	# Median filter with window size defined by argument
		spec_clip = mod_to_obs_spec[i]/filt - 1 														# Dividing spectra by median filter and subtracting one so it is centred on 0
		MAD = np.median(np.abs(spec_clip))
		idx_out = (spec_clip > -1*MAD)# & (spec_clip < 5*MAD)
		norm_mask.append(idx_out)
		if plot_bool:
			plots.plot_norm_mask(lam_obs[i]*1e6, mod_to_obs_spec[i], filt, spec_clip, idx_out, savefig=True, filename=filename+f'_ord_{orders[i]+1}', dirname=dirname)
		
	return norm_mask

def normalise_spectra(lam, spec, norm_mask, window_sizes, orders=[], plot_bool=False, filename='norm', dirname=plots_directory):
	spec_norm = []
	for i in range(len(spec)):
		win = int(2*len(spec[i])*window_sizes[i] + 1)
		filt = pd.Series(spec[i][norm_mask[i]]).rolling(window=win, min_periods=1, center=True).median()	# Second median filter, now with the lines removed
		filt_interp = np.interp(lam[i], lam[i][norm_mask[i]], filt)
		spec_norm.append(spec[i]/filt_interp)
		if plot_bool:
			plots.plot_normalisation(lam[i]*1e6, spec[i], norm_mask[i], filt_interp, savefig=True, filename=filename+f'_ord_{orders[i]+1}', dirname=dirname)

	return np.array(spec_norm, dtype=object)

def cross_correlate(spec_obs, spec_model):
	cc = []
	for i in range(len(spec_obs)):
		smod = spec_model[i]
		sobs = spec_obs[i]
		id_matrix = np.ones(len(smod))
		sobs_norm = np.dot(sobs, id_matrix)/len(sobs)
		smod_norm = np.dot(smod, id_matrix)/len(smod)
		sobs_msub = sobs - sobs_norm
		smod_msub = smod - smod_norm
		cc.append(np.dot(sobs_msub, smod_msub)/np.sqrt(np.dot(smod_msub, smod_msub)*np.dot(sobs_msub, sobs_msub)))

	return cc

def get_loglike(spec_obs, spec_model, cc=None, method='brogiline'):
	log_like = []
	for i in range(len(spec_obs)):
		smod = spec_model[i] 
		sobs = spec_obs[i]
		if method == 'zucker':
			## Calculating the log likelihood from the cross-correlation following Zucker (2003)
			log_like.append(-0.5*len(sobs)*np.log(1-cc[i]**2))
		elif method == 'brogiline':
			## Calculating the log likelihood from the cross-correlation following Brogi & Line (2019)
			id_matrix = np.ones(len(sobs))
			sobs_msub = sobs - np.dot(sobs, id_matrix)/len(sobs)
			sf2 = np.dot(sobs_msub, sobs_msub)/len(sobs_msub)
			smod_msub = smod - np.dot(smod, id_matrix)/len(smod)
			sg2 = np.dot(smod_msub, smod_msub)/len(smod_msub)
			R = np.dot(sobs_msub, smod_msub)/len(sobs_msub)#cc[i]*np.sqrt(sf2*sg2)
			log_like.append(-0.5*len(sobs)*np.log(sf2 - 2*R + sg2))

	return log_like

def get_PWV(hgt_atm, P_atm, T_atm, Xs_H2O):
	R = 8.31446 		# Gas constant, in J mol^-1 K^-1
	mmol_H2O = 0.0182	# Mole mass of water, in kg
	rhoH2O = 1e3 		# Density of water, in kg/m^3
	
	## Calculating the precipitable water vapour (PWV) using eq. (9) from Smette et al. (2015)
	## Xs must be in ppmv, P in Pa, T in K and hgt in km
	PWV = (mmol_H2O/(rhoH2O*R))*np.trapz(Xs_H2O[::-1]*10**P_atm[::-1]/T_atm[::-1], x=hgt_atm[::-1]*1e-3)

	return PWV

def wet_2_dry(Xs_atm):
	Xs_dry = copy.deepcopy(Xs_atm)
	for key in Xs_dry:
		Xs_dry[key] = 1e6*Xs_atm[key]/(1e6 - Xs_atm['H2O'])
	Xs_dry['H2O'] = 0

	return Xs_dry

def rayleigh_scatter(lam, P, H):
	## lam has to be in micrometer, P in hPa and H in km
	cons = (P/1013.25)*(0.00864 + 6.5*10**(-6)*H)
	var = lam**(-(3.916+0.074*lam + 0.050/lam))
	tau = cons*var
	
	return tau

def aerosol_scatter(lam, airmass):
	## lam has to be in micrometer
	k0 = 0.013
	alpha = -1.38
	k = k0*lam**(alpha)
	tau = 0.4*k*airmass

	return tau

def delete_bad_walkers(samples, flat_samples, logp, logpflat, bad_walkers, max_steps, nwalkers, discard, thin):
	samples = np.delete(samples, bad_walkers, axis=1)
	idxs1 = np.arange(0,int((max_steps-discard)*nwalkers/thin),1)
	idxs2 = np.array([], dtype=int)
	for bad in bad_walkers:
		idxs2 = np.append(idxs2,np.arange(bad,int((max_steps-discard)*nwalkers/thin),nwalkers))
	flat_samples = flat_samples[np.delete(idxs1,idxs2),:]
	logp = np.delete(logp, bad_walkers, axis=1)
	logpflat = logpflat[np.delete(idxs1,idxs2)]

	return samples, flat_samples, logp, logpflat

def get_mcmc_result(dirname):
	## Global variables
	global sampler

	sampler = emcee.backends.HDFBackend(dirname+"backupmcmc.h5", read_only=True)
	#mcmcobjfile = open(dirname+'mcmcobject','rb') 
	#sampler = pickle.load(mcmcobjfile)
	
	return sampler

def cut_deep_lines(specs, cutoff_mask):
	specs_cut = []
	for i in range(len(specs)):
		specs_cut.append(specs[i][cutoff_mask[i]])

	return specs_cut

def calc_earth_dist(coords_1, coords_2):
	''' 
	Function to calculate the distance between two points on Earth based on their latitude and longitude
	Uses the Haversine formula
	coors - (lat,lon)
	'''
	r = 6378 	## Equatorial radius of the Earth, in km
	lat_1 = coords_1[0]*np.pi/180 	## Getting latitude of point 1 and converting to radians
	lon_1 = coords_1[1]*np.pi/180 	## Getting longitude of point 1 and converting to radians
	lat_2 = coords_2[0]*np.pi/180 	## Getting latitude of point 2 and converting to radians
	lon_2 = coords_2[1]*np.pi/180 	## Getting longitude of point 2 and converting to radians

	t1 = np.sin((lat_2-lat_1)/2)**2
	t2 = np.cos(lat_1)*np.cos(lat_2)*np.sin((lon_2-lon_1)/2)**2
	
	d = 2*r*np.arcsin((np.sqrt(t1+t2))) 	## Distance, in km

	return d

def get_spec_from_txt(dirnames, order=None):
	data = []
	for i, dirname in enumerate(dirnames):
		print(f'Progress: {i+1}/{len(dirnames)}')
		all_lam = np.empty(0)
		all_spec_obs = np.empty(0)
		all_spec_mod = np.empty(0)
		all_spec_mod_init = np.empty(0)
		if order:
			spec_files = glob.glob(dirname+f'/spec_ord_{order}.txt')
		else:
			spec_files = glob.glob(dirname+'/spec_ord*.txt')
		spec_files.sort()
		for file in spec_files:
			try: 
				lam, spec_obs, spec_mod = np.loadtxt(file, unpack=True)
				all_lam = np.concatenate((all_lam, lam))
				all_spec_obs = np.concatenate((all_spec_obs, spec_obs))
				all_spec_mod = np.concatenate((all_spec_mod, spec_mod))
			except ValueError:
				lam, spec_obs, spec_mod, spec_mod_init = np.loadtxt(file, unpack=True)
				all_lam = np.concatenate((all_lam, lam))
				all_spec_obs = np.concatenate((all_spec_obs, spec_obs))
				all_spec_mod = np.concatenate((all_spec_mod, spec_mod))
				all_spec_mod_init = np.concatenate((all_spec_mod_init, spec_mod_init))
		data.append([all_lam, all_spec_obs, all_spec_mod, all_spec_mod_init])

	return data

def get_not_norm_spec_from_txt(dirnames, order=None):
	data = []
	for i, dirname in enumerate(dirnames):
		print(f'Progress: {i+1}/{len(dirnames)}')
		all_lam = np.empty(0)
		all_spec_obs = np.empty(0)
		if order:
			spec_files = glob.glob(dirname+f'/not_norm_spec_ord_{order}.txt')
		else:
			spec_files = glob.glob(dirname+'/not_norm_spec_ord*.txt')
		spec_files.sort()
		for file in spec_files:
			lam, spec_obs = np.loadtxt(file, unpack=True)
			all_lam = np.concatenate((all_lam, lam))
			all_spec_obs = np.concatenate((all_spec_obs, spec_obs))
		data.append([all_lam, all_spec_obs])

	return data

def time2phase(time,per,T0):
	phase = ((time-T0)%per)/per
	#for ii in range(len(phase)):
	#	if phase[ii] > 0.5: phase[ii] = phase[ii] - 1
	return phase