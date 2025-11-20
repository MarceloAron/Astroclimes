class StellarParams(object):
	'''
	Stellar parameters:
		Teff	: float		- effective temperature, in K
		logg	: float		- surface gravity, in cm/s2
		FeH		: float		- metallicity [Fe/H], in dex
		vsini	: float		- projected rotational velocity, in km/s
		LD		: string	- limb darkening law used: 
								'uni'	- uniform, no LD;
								'lin'	- linear, default;
								'quad'	- quadratic;
								'nl'	- non-linear;
		eps		: float		- linear limb darkening coefficient; 
		R 		: float 	- stellar radius, in m

	Set the parameters by calling
		stelparams = StellarParams()
		stelparams.Teff = 5000.0

	Default is a sun-like star in terms of Teff, logg, and [Fe/H]. 
	'''
	def __init__(self):
		self.Teff = 5750
		self.logg = 4.5
		self.FeH = 0.0
		self.vsini = 10.0
		self.LD = 'lin'
		self.eps = 0.1
		self.R = 6.957*1e8	# Solar radius, in m (from NASA Sun Fact Sheet)

class OrbitParams(object):
	'''
	Orbital parameters:
		Vsys	: float		- systemic velocity, in km/s - positive means redshift, negative means blueshift
		Ks		: float		- stellar RV semi-amplitude, in km/s
		Kp		: float		- planet RV semi-amplitude, in km/s
		P		: float		- orbital period, in days
		T0		: float		- planet mid-transit time, in BJD
		Rp		: float		- planet radius, in m

	Set the parameters by calling
		orbparams = OrbitParams()
		orbparams.Vsys = 10

	Default uses dummy values
	'''
	def __init__(self):
		self.Vsys = 10.0
		self.Ks = 0.5
		self.Kp = 100.
		self.P = 3.
		self.T0 = 2450000.
		self.Rp = 6.9911*1e7	# Jupiter radius, in m (from NASA Jupiter Fact Sheet)

class FitParams(object):
	## Currently not used by the algorithm, may be implemented in future versions
	'''
	Fit parameters:
		g_CO2: float - ground level abundance of CO2, in molecules/m^3.
		g_CH4: float - ground level abundance of CH4, in molecules/m^3.
		g_H2O: float - ground level abundance of H2O, in molecules/m^3.
		g_O2: float - ground level abundance of O2, in molecules/m^3.

	Set the parameters by calling
		params = FitParams()
		params.gCO2 = 5e21
	'''
	def __init__(self):
		self.gCO2 = 1e21
		self.gCH4 = 1e19
		self.gH2O = 1e23
		self.gO2 = 1e24
		self.gN2 = 1e25