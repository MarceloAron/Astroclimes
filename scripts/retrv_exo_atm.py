# =====================================================================================
# Basic packages
# =====================================================================================
import numpy as np
from scipy import interpolate
import os
from subprocess import call
import emcee
from multiprocessing import Pool

# =====================================================================================
# Scripts
# =====================================================================================
import funcs
import funcs_telrem
import plots
#import fit_retrv_exo_atm

## Some constants
R_Sun = 6.957*1e8	# Solar radius, in m (from NASA Sun Fact Sheet)
R_Jup = 6.9911*1e7	# Jupiter radius, in m (from NASA Jupiter Fact Sheet)

## Color-blind friendly colors: blue, orange, green, pink, brown, purple, gray, red and yellow
CBcols= ['#377EB8', '#FF7F00', '#4DAF4A',
		 '#F781BF', '#A65628', '#984EA3',
		 '#999999', '#E41A1C', '#DEDE00']

def telrem_each_night(nights, main_retrv_exo_directory, filenames_site_values, n_orders, bad_spec_orders_idxs, bad_obs_idxs, nc_PCAs, filename_retrv_exo_SNR_AC, filename_retrv_exo_SNR_PCA, orbpars, make_gen_cube_plots, plot_telrem_steps_AC, plot_telrem_steps_PCA, plot_indiv_PCA_comps, deep_line_threshold, dirname_for_plots, numpy_objs_directory):
	'''
	This function carries out the telluric and stellar line removal for each of the selected nights, using both Astroclimes and PCA
	Optionally, you can also plot the observational and model data cubes, as well as the residuals, along with the step-by-step
	removal process from Astroclimes and PCA and the contribution of each individual PCA component
	'''
	## Defining variables that will store information for all nights
	all_nights_lam = []
	all_nights_spec = []
	all_nights_spec_mod = []
	all_nights_pca_fit = []
	all_nights_res_ac = []
	all_nights_res_pca = []
	all_nights_rule_cube = []

	## Looping over all nights
	for i, night in enumerate(nights):
		
		## Unpacking the observational wavelength and spectra cube (no need to give any sensible arguments when generate=False, besides save_dirname and save_lam_filename and save_spec_filename if file does not have default name)
		all_lam, all_spec = funcs_telrem.get_obs_cube(None, None, None, generate=False, save_dirname=main_retrv_exo_directory, save_lam_filename=f'all_lam_{night}', save_spec_filename=f'all_spec_{night}')

		## Unpacking the model spectra cube (no need to give any sensible arguments when generate=False, besideds save_dirname and save_spec_mod_filename if file does not have default name)
		all_spec_mod = funcs_telrem.get_mod_cube(None, None, None, None, None, None, None, None, generate=False, save_dirname=main_retrv_exo_directory, save_spec_mod_filename=f'all_spec_mod_{night}')

		## Deleting unwanted orders (due to too many saturated telluric lines)
		orders = np.arange(1,n_orders[i]+1)
		good_orders = np.delete(orders, bad_spec_orders_idxs[i])
		all_lam = np.delete(all_lam, bad_spec_orders_idxs[i], axis=0)
		all_spec = np.delete(all_spec, bad_spec_orders_idxs[i], axis=0)
		all_spec_mod = np.delete(all_spec_mod, bad_spec_orders_idxs[i], axis=0)

		## Deleting unwanted observations (due to bad SNR, bad weather, etc.)
		all_lam = np.delete(all_lam, bad_obs_idxs[i], axis=1)
		all_spec = np.delete(all_spec, bad_obs_idxs[i], axis=1)
		all_spec_mod = np.delete(all_spec_mod, bad_obs_idxs[i], axis=1)

		## Calculating the "median divided" observational AND model flux cubes before using our telluric models to remove the telluric lines
		all_spec_mod_corr = all_spec_mod.copy()
		all_spec_corr = all_spec.copy()
		for jj in range(np.shape(all_lam)[0]):
			all_spec_mod_corr[jj,:,:] = all_spec_mod[jj,:,:]/np.nanmedian(all_spec_mod[jj,:,:], axis=1)[:,np.newaxis]
			all_spec_corr[jj,:,:] = all_spec[jj,:,:]/np.nanmedian(all_spec[jj,:,:], axis=1)[:,np.newaxis]

		## Masking deep lines
		all_spec_corr_masked, rule_cube = funcs_telrem.mask_cube(all_spec_corr, mask_val=0.0, deep_line_threshold=deep_line_threshold)
		all_spec_corr_masked_for_plots, rule_cube_for_plots = funcs_telrem.mask_cube(all_spec_corr, mask_val=np.nan, deep_line_threshold=deep_line_threshold)

		## Removing the telluric and stellar lines from the observational spectra using Astroclimes
		## Dividing the observational spectra by the telluric model spectra should remove the telluric lines
		res_ac = all_spec_corr_masked/all_spec_mod_corr

		## Now we divide the residuals from above by their mean over all nights (i.e. a master stellar template), which should remove the stellar lines
		## Note that this step only works if one has several observations taken on the same night, and even still some of the planetary signal may be removed as well
		res_ac_mean = np.nanmean(res_ac, axis=1)[:,np.newaxis,:] # Mean of the telluric line removed spectra over all observations, to remove stellar lines
		#res_ac_median = np.nanmedian(res_ac, axis=1)[:,np.newaxis,:]	# Median of the telluric line removed spectra over all observations, to remove stellar lines
		res_ac_mean_rem = res_ac/res_ac_mean
		
		## Getting the PCA "cube" for both the case with and without the injected planetary signal
		nc_PCA = nc_PCAs[i]
		if not os.path.isfile(main_retrv_exo_directory+f'analysis_files/all_pca_fit_{night}_nc_{nc_PCA}'):
			all_pca_fit = funcs_telrem.run_pca(all_lam, all_spec_corr, nc_PCA, generate=True, save_dirname=main_retrv_exo_directory, save_filename=f'analysis_files/all_pca_fit_{night}')
		else:
			all_pca_fit = funcs_telrem.run_pca(0, 0, nc_PCA, generate=False, save_dirname=main_retrv_exo_directory, save_filename=f'analysis_files/all_pca_fit_{night}')

		## Removing telluric and stellar lines using PCA
		res_pca = all_spec_corr/all_pca_fit

		T0 = orbpars.T0
		P = orbpars.P	
		
		BJD = np.loadtxt(filenames_site_values[i], usecols=(4), unpack=True)
		BJD = np.delete(BJD, bad_obs_idxs[i], axis=0)
		phases = ((BJD - T0)%P)/P

		if make_gen_cube_plots:

			## Choosing orders to plot
			plot_orders = good_orders

			## Plotting the observational spectra as it comes, with its median removed, and with its median removed and with deep lines masked, respectively
			plots.plot_cubes(all_lam, all_spec, list(good_orders), plot_orders, 'Observational cube', y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'obs_flux_cube_{night}')
			plots.plot_cubes(all_lam, all_spec_corr, list(good_orders), plot_orders, 'Observational cube (median removed)', vmin=0, vmax=1.1, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'obs_flux_cube_medrem_{night}')
			plots.plot_cubes(all_lam, all_spec_corr_masked_for_plots, list(good_orders), plot_orders, 'Observational cube (median removed, deep lines masked)', vmin=0, vmax=1.1, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'obs_flux_cube_medrem_deepmask_{night}')

			## Plotting the Astroclimes model cube and residuals, the mean of the residuals, and the residuals with the mean removed, respectively
			plots.plot_cubes(all_lam, all_spec_mod_corr, list(good_orders), plot_orders, 'Model cube (Astroclimes)', vmin=0, vmax=1.1, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'mod_flux_cube_{night}')
			plots.plot_cubes(all_lam, res_ac, list(good_orders), plot_orders, 'Residuals (Astroclimes)', vmin=0.90, vmax=1.05, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'res_ac_flux_cube_{night}')
			plots.plot_cubes(all_lam, np.concatenate(([res_ac_mean]*np.shape(all_lam)[1]), axis=1), list(good_orders), plot_orders, 'Residuals (Astroclimes, mean)', vmin=0.8, vmax=1.05, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'res_ac_flux_cube_mean_{night}')
			plots.plot_cubes(all_lam, res_ac_mean_rem, list(good_orders), plot_orders, 'Residuals (Astroclimes, mean removed)', vmin=0.9, vmax=1.05, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'res_ac_flux_cube_mean_rem_{night}')

			## Plotting the PCA model and PCA residuals without and with the injected planetary signal, respectively
			plots.plot_cubes(all_lam, all_pca_fit, list(good_orders), plot_orders, 'PCA cube', y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'pca_flux_cube_{night}_nc_{nc_PCA}')
			plots.plot_cubes(all_lam, res_pca, list(good_orders), plot_orders, 'Residuals (PCA)', vmin=0.90, vmax=1.05, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'pca_res_flux_cube_{night}_nc_{nc_PCA}')

		chosen_order_idx = list(good_orders).index(12)
		lam = all_lam[chosen_order_idx,0,:]
		if plot_telrem_steps_AC:
			## Plotting the steps taken on the telluric removal process
			step_1, rule_cube = funcs_telrem.mask_cube(all_spec, np.nan, rule_cube=rule_cube, create_rule_cube=False)
			step_1 = step_1[chosen_order_idx,:,:]
			step_2, rule_cube = funcs_telrem.mask_cube(all_spec_corr, np.nan, rule_cube=rule_cube, create_rule_cube=False)
			step_2 = step_2[chosen_order_idx,:,:]
			step_3, rule_cube = funcs_telrem.mask_cube(all_spec_mod_corr, np.nan, rule_cube=rule_cube, create_rule_cube=False)
			step_3 = step_3[chosen_order_idx,:,:]
			step_4, rule_cube = funcs_telrem.mask_cube(res_ac, np.nan, rule_cube=rule_cube, create_rule_cube=False)
			step_4 = step_4[chosen_order_idx,:,:]
			#step_4p5 = np.concatenate(([res_ac_mean]*np.shape(all_lam)[1]), axis=1)[chosen_order_idx,:,:]
			#step_4p5, rule_cube = funcs_telrem.mask_cube(np.concatenate(([res_ac_mean]*np.shape(all_lam)[1]), axis=1), np.nan, rule_cube=rule_cube, create_rule_cube=False)
			#step_4p5 = step_4p5[chosen_order_idx,:,:]
			step_5, rule_cube = funcs_telrem.mask_cube(res_ac_mean_rem, np.nan, rule_cube=rule_cube, create_rule_cube=False)
			step_5 = step_5[chosen_order_idx,:,:]
			step_6_PCA, rule_cube = funcs_telrem.mask_cube(res_pca, np.nan, rule_cube=rule_cube, create_rule_cube=False)
			step_6_PCA = step_6_PCA[chosen_order_idx,:,:]
			cube_steps = [step_1, step_2, step_3, step_4, step_5, step_6_PCA]
			vmins = [None, None, None, np.nanmedian(step_4.flatten())-3*np.nanstd(step_4.flatten()), np.nanmedian(step_5.flatten())-3*np.nanstd(step_5.flatten()), np.nanmedian(step_5.flatten())-3*np.nanstd(step_5.flatten())]
			vmaxs = [None, None, None, np.nanmedian(step_4.flatten())+3*np.nanstd(step_4.flatten()), np.nanmedian(step_5.flatten())+3*np.nanstd(step_5.flatten()), np.nanmedian(step_5.flatten())+3*np.nanstd(step_5.flatten())]
			plots.plot_cube_steps(lam, cube_steps, vmins=vmins, vmaxs=vmaxs, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'removal_steps_ac_ord_{good_orders[chosen_order_idx]}_with_PCA_{night}')

		if plot_telrem_steps_PCA:
			## Now for PCA
			step_1, rule_cube = funcs_telrem.mask_cube(all_spec, np.nan, rule_cube=rule_cube, create_rule_cube=False)
			step_1 = step_1[chosen_order_idx,:,:]
			step_2, rule_cube = funcs_telrem.mask_cube(all_spec_corr, np.nan, rule_cube=rule_cube, create_rule_cube=False)
			step_2 = step_2[chosen_order_idx,:,:]
			step_3_PCA, rule_cube = funcs_telrem.mask_cube(all_pca_fit, np.nan, rule_cube=rule_cube, create_rule_cube=False)
			step_3_PCA = step_3_PCA[chosen_order_idx,:,:]
			step_4_PCA, rule_cube = funcs_telrem.mask_cube(res_pca, np.nan, rule_cube=rule_cube, create_rule_cube=False)
			step_4_PCA = step_4_PCA[chosen_order_idx,:,:]
			cube_steps = [step_1, step_2, step_3_PCA, step_4_PCA]
			vmins = [None, None, None, np.nanmedian(step_4_PCA.flatten())-3*np.nanstd(step_4_PCA.flatten())]
			vmaxs = [None, None, None, np.nanmedian(step_4_PCA.flatten())+3*np.nanstd(step_4_PCA.flatten())]
			plots.plot_cube_steps(lam, cube_steps, vmins=vmins, vmaxs=vmaxs, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'removal_steps_pca_nc_{nc_PCA}_ord_{good_orders[chosen_order_idx]}_{night}')

		if plot_indiv_PCA_comps:
			## Plotting the difference between each PCA component
			all_pca_fit_1 = funcs_telrem.run_pca(all_lam, all_spec_corr, 1, generate=True, save_dirname=main_retrv_exo_directory, save_filename='all_pca_fit')
			all_pca_fit_2 = funcs_telrem.run_pca(all_lam, all_spec_corr, 2, generate=True, save_dirname=main_retrv_exo_directory, save_filename='all_pca_fit')
			all_pca_fit_3 = funcs_telrem.run_pca(all_lam, all_spec_corr, 3, generate=True, save_dirname=main_retrv_exo_directory, save_filename='all_pca_fit')
			all_pca_fit_4 = funcs_telrem.run_pca(all_lam, all_spec_corr, 4, generate=True, save_dirname=main_retrv_exo_directory, save_filename='all_pca_fit')
			all_pca_fit_5 = funcs_telrem.run_pca(all_lam, all_spec_corr, 5, generate=True, save_dirname=main_retrv_exo_directory, save_filename='all_pca_fit')
			all_pca_fit_6 = funcs_telrem.run_pca(all_lam, all_spec_corr, 6, generate=True, save_dirname=main_retrv_exo_directory, save_filename='all_pca_fit')
			#all_pca_fit_7 = funcs_telrem.run_pca(all_lam, all_spec_corr, 7, generate=True, save_dirname=main_retrv_exo_directory, save_filename='all_pca_fit')
			#all_pca_fit_8 = funcs_telrem.run_pca(all_lam, all_spec_corr, 8, generate=True, save_dirname=main_retrv_exo_directory, save_filename='all_pca_fit')
			all_pca_fit_9 = funcs_telrem.run_pca(all_lam, all_spec_corr, 9, generate=True, save_dirname=main_retrv_exo_directory, save_filename='all_pca_fit')
			all_pca_fit_10 = funcs_telrem.run_pca(all_lam, all_spec_corr, 10, generate=True, save_dirname=main_retrv_exo_directory, save_filename='all_pca_fit')
			comp_1 = all_pca_fit_1[chosen_order_idx,:,:] 
			comp_2 = all_pca_fit_2[chosen_order_idx,:,:] - all_pca_fit_1[chosen_order_idx,:,:]
			comp_3 = all_pca_fit_3[chosen_order_idx,:,:] - all_pca_fit_2[chosen_order_idx,:,:]
			comp_4 = all_pca_fit_4[chosen_order_idx,:,:] - all_pca_fit_3[chosen_order_idx,:,:]
			comp_5 = all_pca_fit_5[chosen_order_idx,:,:] - all_pca_fit_4[chosen_order_idx,:,:]
			comp_6 = all_pca_fit_6[chosen_order_idx,:,:] - all_pca_fit_5[chosen_order_idx,:,:]
			#comp_7 = all_pca_fit_7[chosen_order_idx,:,:] - all_pca_fit_6[chosen_order_idx,:,:]
			#comp_8 = all_pca_fit_8[chosen_order_idx,:,:] - all_pca_fit_7[chosen_order_idx,:,:]
			#comp_9 = all_pca_fit_9[chosen_order_idx,:,:] - all_pca_fit_8[chosen_order_idx,:,:]
			comp_10 = all_pca_fit_10[chosen_order_idx,:,:] - all_pca_fit_9[chosen_order_idx,:,:]
			pca_components = [all_spec_corr[chosen_order_idx,:,:], comp_1, comp_2, comp_3, comp_4, comp_5, comp_6, comp_10]
			vmins = [None, None, np.nanmedian(comp_2.flatten())-3*np.nanstd(comp_2.flatten()), np.nanmedian(comp_3.flatten())-3*np.nanstd(comp_3.flatten()), np.nanmedian(comp_4.flatten())-3*np.nanstd(comp_4.flatten()), np.nanmedian(comp_5.flatten())-3*np.nanstd(comp_5.flatten()), np.nanmedian(comp_6.flatten())-3*np.nanstd(comp_6.flatten()), np.nanmedian(comp_10.flatten())-3*np.nanstd(comp_10.flatten())]
			vmaxs = [None, None, np.nanmedian(comp_2.flatten())+3*np.nanstd(comp_2.flatten()), np.nanmedian(comp_3.flatten())+3*np.nanstd(comp_3.flatten()), np.nanmedian(comp_4.flatten())+3*np.nanstd(comp_4.flatten()), np.nanmedian(comp_5.flatten())+3*np.nanstd(comp_5.flatten()), np.nanmedian(comp_6.flatten())+3*np.nanstd(comp_6.flatten()), np.nanmedian(comp_10.flatten())+3*np.nanstd(comp_10.flatten())]
			plots.plot_cube_steps(lam, pca_components, vmins=vmins, vmaxs=vmaxs, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'pca_components_ord_{good_orders[chosen_order_idx]}_{night}')

		all_nights_lam.append(all_lam)
		all_nights_spec.append(all_spec)
		all_nights_spec_mod.append(all_spec_mod)
		all_nights_pca_fit.append(all_pca_fit)
		all_nights_res_ac.append(res_ac_mean_rem)
		all_nights_res_pca.append(res_pca)
		all_nights_rule_cube.append(rule_cube)

		nightly_cube = [all_lam,
						all_spec,
						all_spec_mod,
						all_pca_fit,
						res_ac_mean_rem,
						res_pca,
						rule_cube]
    	
		np.save(numpy_objs_directory+f'cubes_{night}.npy', nightly_cube)

	return all_nights_lam, all_nights_spec, all_nights_spec_mod, all_nights_pca_fit, all_nights_res_ac, all_nights_res_pca, all_nights_rule_cube

def get_planet_sig_all_nights(nights, all_nights_cubes, main_retrv_exo_directory, filenames_site_values, R_instruments, n_orders, bad_spec_orders_idxs, bad_obs_idxs, nc_PCAs, filename_retrv_exo_SNR_AC, filename_retrv_exo_SNR_PCA, filename_planet_lam, filename_planet_flux, stelpars, orbpars, SNR_lines_params, dirname_for_plots, numpy_objs_directory):
	'''
	This function calculates the detection strength of a planetary signal for each of the nights individually and saved the results so they can be combined later
	'''

	all_nights_lam = all_nights_cubes[0]
	#all_nights_spec = all_nights_cubes[1]
	#all_nights_spec_mod = all_nights_cubes[2]
	#all_nights_pca_fit = all_nights_cubes[3]
	all_nights_res_ac = all_nights_cubes[4]
	all_nights_res_pca = all_nights_cubes[5]
	all_nights_rule_cube = all_nights_cubes[6]

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

	for i, night in enumerate(nights):

		lam = all_nights_lam[i]
		res_ac = all_nights_res_ac[i]
		res_pca = all_nights_res_pca[i]
		rule_cube = all_nights_rule_cube[i]

		## Creating a model planetary signal
		T_star = stelpars.Teff
		R_star = stelpars.R
		R_planet = orbpars.Rp
		T0 = orbpars.T0
		P = orbpars.P	
		Kp = orbpars.Kp			
		Vsys = orbpars.Vsys
		## Planet cube already comes without the bad orders, so there's no need to delete them
		planet_cube, planet_mod_lam, planet_model = funcs_telrem.create_planet_signal(lam, filename_planet_lam, filename_planet_flux, filenames_site_values[i], 
																						R_instruments[i], T_star, R_star, R_planet, 1, T0, P, Kp, Vsys, generate=True, save_dirname=main_retrv_exo_directory)
		
		## Deleting unwanted observations (due to bad SNR, bad weather, etc.)
		planet_cube = np.delete(planet_cube, bad_obs_idxs[i], axis=1)

		## Producing the planet's trail and velocity map
		cs = interpolate.splrep(planet_mod_lam, planet_model)
		rv = np.arange(-300,301,1.5)
		BJD, V_BERV = np.loadtxt(filenames_site_values[i], usecols=(4,15), unpack=True)
		BJD = np.delete(BJD, bad_obs_idxs[i], axis=0)
		V_BERV = np.delete(V_BERV, bad_obs_idxs[i], axis=0)
		V_BERV = -V_BERV ## V_BERV is the radial velocity of the Earth with respect to the Sun, but we want the opposite
		phases = ((BJD - T0)%P)/P
		#phases = phases - 0.5
		vRest = np.arange(-75,76,1.5)
		kpVec = np.arange(0,251)
		vObs = Vsys + V_BERV
		RVs = vObs + Kp*np.sin(2*np.pi*phases)

		## Calculating the CCF and log likelihood between the residuals and the model planetary signal
		## For the Astroclimes removal
		#res_ac, rule_cube_ac = funcs_telrem.mask_data_post_pca(res_ac, rule_cube)
		tr, ll = funcs_telrem.get_fixed_ccf_cube(res_ac, lam[:,0,:], rv, cs, rule_cube)
		planet_trail = np.nansum(tr, axis=0)
		planet_trail_mean_rem = planet_trail - np.mean(planet_trail, axis=1)[:,np.newaxis]
		vmap = funcs_telrem.get_velocity_map(planet_trail,rv,phases,vObs,vRest,kpVec)

		## For the PCA removal
		#res_pca, rule_cube_pca = funcs_telrem.mask_data_post_pca(res_pca, rule_cube)
		tr_PCA, ll_PCA = funcs_telrem.get_fixed_ccf_cube(res_pca, lam[:,0,:], rv, cs, rule_cube)
		planet_trail_PCA = np.nansum(tr_PCA, axis=0)
		planet_trail_PCA_mean_rem = planet_trail_PCA - np.mean(planet_trail_PCA, axis=1)[:,np.newaxis]
		vmap_PCA = funcs_telrem.get_velocity_map(planet_trail_PCA,rv,phases,vObs,vRest,kpVec)

		## Choosing orders to plot
		orders = np.arange(1,n_orders[0]+1)
		good_orders = np.delete(orders, bad_spec_orders_idxs[i])
		plot_orders = good_orders

		## Plotting the model planet cube
		plots.plot_cubes(lam, planet_cube, list(good_orders), plot_orders, 'Model planetary signal', vmin=0, vmax=1.1, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'mod_planet_cube_{night}')
		
		## Plotting the planet trail (as it as mean removed) and velocity map for Astroclimes and PCA (reminder that we only calculated these for the case with an injected signal), respectively
		#plots.plot_squares(planet_trail, extent=[np.min(rv),np.max(rv),np.min(phases),np.max(phases)], xlabel='RV (km/s)', ylabel='Orbital phase', cmap='gray', savefig=True, dirname=dirname_for_plots, filename=f'planet_trail_{night}')
		plots.plot_squares(planet_trail_mean_rem, extent=[np.min(rv),np.max(rv),np.min(phases),np.max(phases)], xlabel='RV (km/s)', ylabel='Orbital phase', vmin=np.min([np.min(planet_trail_mean_rem),np.min(planet_trail_PCA_mean_rem)]), vmax=np.max([np.max(planet_trail_mean_rem),np.max(planet_trail_PCA_mean_rem)]), cmap='gray', savefig=True, dirname=dirname_for_plots, filename=f'planet_trail_{night}')
		plots.plot_squares(vmap, extent=[np.min(vRest),np.max(vRest),0,len(vmap)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', plot_cmap=True, savefig=True, dirname=dirname_for_plots, filename=f'vmap_{night}')

		#plots.plot_squares(planet_trail_PCA, extent=[np.min(rv),np.max(rv),np.min(phases),np.max(phases)], xlabel='RV (km/s)', ylabel='Orbital phase', cmap='gray', savefig=True, dirname=dirname_for_plots, filename=f'planet_trail_PCA_nc_{nc_PCAs[i]}_{night}')
		plots.plot_squares(planet_trail_PCA_mean_rem, extent=[np.min(rv),np.max(rv),np.min(phases),np.max(phases)], xlabel='RV (km/s)', ylabel='Orbital phase', vmin=np.min([np.min(planet_trail_mean_rem),np.min(planet_trail_PCA_mean_rem)]), vmax=np.max([np.max(planet_trail_mean_rem),np.max(planet_trail_PCA_mean_rem)]), cmap='gray', savefig=True, dirname=dirname_for_plots, filename=f'planet_trail_PCA_{night}_nc_{nc_PCAs[i]}')
		plots.plot_squares(vmap_PCA, extent=[np.min(vRest),np.max(vRest),0,len(vmap)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', plot_cmap=True, savefig=True, dirname=dirname_for_plots, filename=f'vmap_PCA_nc_{nc_PCAs[i]}_{night}')
		
		## Calculating the SNR for the velocity map

		## Defining two lines with roughly the same inclination (determined by visual inspection) 
		## as the planetary signal on the velocity map, one slightly above (l1) and the other slightly below (l2)
		l1 = SNR_lines_params[i][0]*vRest + Kp + SNR_lines_params[i][1]
		l2 = SNR_lines_params[i][0]*vRest + Kp - SNR_lines_params[i][1]
		
		## Creating a boolean mask where True represents the values outside the region
		## delimited by the two lines defined above, along with two arrays containing only the values
		## on the "noise region" and on the "signal region" (mostly for plotting/verification purposes)
		vmap_mask = np.zeros_like(vmap, dtype=bool)
		vmap_noise_region = np.zeros_like(vmap)
		vmap_signal_region = np.zeros_like(vmap)
		vmap_PCA_noise_region = np.zeros_like(vmap_PCA)
		vmap_PCA_signal_region = np.zeros_like(vmap_PCA)
		for jj in range(len(vRest)):
			rr = (kpVec > l1[jj]) | (kpVec < l2[jj])
			vmap_mask[rr,jj] = True
			vmap_noise_region[rr,jj] = vmap[rr,jj]
			vmap_signal_region[~rr,jj] = vmap[~rr,jj]
			vmap_PCA_noise_region[rr,jj] = vmap_PCA[rr,jj]
			vmap_PCA_signal_region[~rr,jj] = vmap_PCA[~rr,jj]		
		
		## The noise level is calculated as the mean of the points outside the region delimited by the two lines
		noise = np.mean(vmap[vmap_mask])
		noise_std = np.std(vmap[vmap_mask])
		noise_PCA = np.mean(vmap_PCA[vmap_mask])
		noise_PCA_std = np.std(vmap_PCA[vmap_mask])
		
		## The SNR is simply the vmap values minus this noise divided by the standard deviation of the noise, and we calculate 
		## the maximum inside the region delimited by the two lines. We know the planet signal is inside this region, so we want 
		## to quantify it. For cases when the injected signal is small, the global maximum can be elsewhere, which is why we 
		## explicitly calculate the maximum inside the region
		SNR = (vmap-noise)/noise_std
		SNR_max = np.max(SNR[~vmap_mask])
		SNR_PCA = (vmap_PCA-noise_PCA)/noise_PCA_std
		SNR_PCA_max = np.max(SNR_PCA[~vmap_mask])
		
		## Coordinates in the vRest X kpVec space where the maximum occurs
		vRest_max_idx = np.where(SNR == SNR_max)[1][0]
		kpVec_max_idx = np.where(SNR == SNR_max)[0][0]
		vRest_max_idx_PCA = np.where(SNR_PCA == SNR_PCA_max)[1][0]
		kpVec_max_idx_PCA = np.where(SNR_PCA == SNR_PCA_max)[0][0]
		
		## Coordinates and value of the exact location of the injected signal
		kpVec_inj_sig_pos_idx = np.where(kpVec == Kp)[0][0]
		vRest_inj_sig_pos_idx = np.where(vRest == 0)[0][0]
		SNR_inj_sig_pos = SNR[kpVec_inj_sig_pos_idx,vRest_inj_sig_pos_idx]
		SNR_PCA_inj_sig_pos = SNR_PCA[kpVec_inj_sig_pos_idx,vRest_inj_sig_pos_idx]
		
		plots.plot_SNR_square(	SNR, extent=[np.min(vRest),np.max(vRest),np.min(kpVec),np.max(kpVec)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', 
								vRest=vRest, kpVec=kpVec, l1=l1, l2=l2, SNR_max=SNR_max, vRest_max_idx=vRest_max_idx, kpVec_max_idx=kpVec_max_idx, 
								SNR_inj_sig_pos=SNR_inj_sig_pos, vRest_inj_sig_pos_idx=vRest_inj_sig_pos_idx, kpVec_inj_sig_pos_idx=kpVec_inj_sig_pos_idx,
								vmin=np.min([np.min(SNR),np.min(SNR_PCA)]), vmax=np.max([np.max(SNR),np.max(SNR_PCA)]),
								savefig=True, dirname=dirname_for_plots, filename=f'SNR_{night}')
	
		plots.plot_SNR_square(	SNR_PCA, extent=[np.min(vRest),np.max(vRest),np.min(kpVec),np.max(kpVec)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', 
								vRest=vRest, kpVec=kpVec, l1=l1, l2=l2, SNR_max=SNR_PCA_max, vRest_max_idx=vRest_max_idx_PCA, kpVec_max_idx=kpVec_max_idx_PCA, 
								SNR_inj_sig_pos=SNR_PCA_inj_sig_pos, vRest_inj_sig_pos_idx=vRest_inj_sig_pos_idx, kpVec_inj_sig_pos_idx=kpVec_inj_sig_pos_idx, 
								vmin=np.min([np.min(SNR),np.min(SNR_PCA)]), vmax=np.max([np.max(SNR),np.max(SNR_PCA)]),
								savefig=True, dirname=dirname_for_plots, filename=f'SNR_PCA_{night}_nc_{nc_PCAs[i]}')
		
		text = (f"{SNR_inj_sig_pos:10.4f} \t {vRest[vRest_inj_sig_pos_idx]:10.4f} \t {kpVec[kpVec_inj_sig_pos_idx]:10.4f} \t " +
				f"{SNR_max:10.4f} \t {vRest[vRest_max_idx]:10.4f} \t {kpVec[kpVec_max_idx]:10.4f}\n")
		txt = open(filename_retrv_exo_SNR_AC, 'a')
		txt.write(text)
		txt.close()

		text =	(f"{nc_PCAs[i]:<5.0f} \t {SNR_PCA_inj_sig_pos:10.4f} \t {vRest[vRest_inj_sig_pos_idx]:10.4f} \t {kpVec[kpVec_inj_sig_pos_idx]:10.4f} \t " +
				 f"{SNR_PCA_max:10.4f} \t {vRest[vRest_max_idx_PCA]:10.4f} \t {kpVec[kpVec_max_idx_PCA]:10.4f}\n")
		txt = open(filename_retrv_exo_SNR_PCA, 'a')
		txt.write(text)
		txt.close()

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

		#for i in range(len(nights)):
		np.save(numpy_objs_directory+f'phases_{night}.npy', phases)
		np.save(numpy_objs_directory+f'RVs_{night}.npy', RVs)
		np.save(numpy_objs_directory+f'tr_{night}.npy', tr)
		np.save(numpy_objs_directory+f'll_{night}.npy', ll)
		np.save(numpy_objs_directory+f'planet_trail_{night}.npy', planet_trail_mean_rem)
		np.save(numpy_objs_directory+f'vmap_{night}.npy', vmap)
		np.save(numpy_objs_directory+f'SNR_{night}.npy', SNR)
		np.save(numpy_objs_directory+f'tr_PCA_{night}.npy', tr_PCA)
		np.save(numpy_objs_directory+f'll_PCA_{night}.npy', ll_PCA)
		np.save(numpy_objs_directory+f'planet_trail_PCA_{night}.npy', planet_trail_PCA_mean_rem)
		np.save(numpy_objs_directory+f'vmap_PCA_{night}.npy', vmap_PCA)
		np.save(numpy_objs_directory+f'SNR_PCA_{night}.npy', SNR_PCA)

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

	return all_nights_results

def combine_sig_multi_nights(all_nights_results, dirname_for_plots, filename_suffix):

	all_nights_phases = all_nights_results[0]
	all_nights_RVs = all_nights_results[1]
	#all_nights_tr = all_nights_results[2]
	#all_nights_ll = all_nights_results[3]
	all_nights_planet_trail = all_nights_results[4]
	all_nights_vmap = all_nights_results[5]
	#all_nights_SNR = all_nights_results[6]

	#all_nights_tr_PCA = all_nights_results[7]
	#all_nights_ll_PCA = all_nights_results[8]
	all_nights_planet_trail_PCA = all_nights_results[9]
	all_nights_vmap_PCA = all_nights_results[10]
	#all_nights_SNR_PCA = all_nights_results[11]

	rv = np.arange(-300,301,1.5)
	vRest = np.arange(-75,76,1.5)
	kpVec = np.arange(0,251)

	## Combining the signal from all nights by summing their velocity maps
	comb_vmap = np.sum(all_nights_vmap, axis=0)
	comb_vmap_PCA = np.sum(all_nights_vmap_PCA, axis=0)
	
	plots.plot_squares(comb_vmap, extent=[np.min(vRest),np.max(vRest),0,len(comb_vmap)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', plot_cmap=True, savefig=True, dirname=dirname_for_plots, filename='combined_vmap'+filename_suffix)
	plots.plot_squares(comb_vmap_PCA, extent=[np.min(vRest),np.max(vRest),0,len(comb_vmap_PCA)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', plot_cmap=True, savefig=True, dirname=dirname_for_plots, filename='combined_vmap_PCA'+filename_suffix)
	plots.plot_squares(comb_vmap_PCA, extent=[np.min(vRest),np.max(vRest),0,len(comb_vmap_PCA)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', vmin=-4, vmax=10, plot_cmap=True, savefig=True, dirname=dirname_for_plots, filename='combined_vmap_PCA_same_scale'+filename_suffix)
	
	#comb_phases = np.concatenate((all_nights_phases[0],all_nights_phases[1],all_nights_phases[2],all_nights_phases[3],all_nights_phases[4],all_nights_phases[5]))
	#comb_planet_trails = np.concatenate((all_nights_planet_trail[0],all_nights_planet_trail[1],all_nights_planet_trail[2],all_nights_planet_trail[3],all_nights_planet_trail[4],all_nights_planet_trail[5]))
	#comb_planet_trails_PCA = np.concatenate((all_nights_planet_trail_PCA[0],all_nights_planet_trail_PCA[1],all_nights_planet_trail_PCA[2],all_nights_planet_trail_PCA[3],all_nights_planet_trail_PCA[4],all_nights_planet_trail_PCA[5]))
	#plots.plot_squares(comb_planet_trails, extent=[np.min(rv),np.max(rv),np.min(comb_phases),np.max(comb_phases)], xlabel='RV (km/s)', ylabel='Orbital phase', vmin=np.min([np.min(comb_planet_trails),np.min(comb_planet_trails_PCA)]), vmax=np.max([np.max(comb_planet_trails),np.max(comb_planet_trails_PCA)]), savefig=True, dirname=dirname_for_plots, filename='combined_planet_trail')
	#plots.plot_planet_trail_multi_nights(comb_planet_trails, extent=[np.min(rv),np.max(rv),0,len(comb_phases)], xlabel='RV (km/s)', ylabel='Orbital phase', comb_phases=comb_phases, indiv_phases=all_nights_phases, RVs=all_nights_RVs, savefig=True, dirname=dirname_for_plots, filename='combined_planet_trail')
	#plots.plot_squares(comb_planet_trails_PCA, extent=[np.min(rv),np.max(rv),np.min(comb_phases),np.max(comb_phases)], xlabel='RV (km/s)', ylabel='Orbital phase', vmin=np.min([np.min(comb_planet_trails),np.min(comb_planet_trails_PCA)]), vmax=np.max([np.max(comb_planet_trails),np.max(comb_planet_trails_PCA)]), savefig=True, dirname=dirname_for_plots, filename='combined_planet_trail_PCA')
	#plots.plot_planet_trail_multi_nights(comb_planet_trails_PCA, extent=[np.min(rv),np.max(rv),0,len(comb_phases)], xlabel='RV (km/s)', ylabel='Orbital phase', comb_phases=comb_phases, indiv_phases=all_nights_phases, RVs=all_nights_RVs, savefig=True, dirname=dirname_for_plots, filename='combined_planet_trail_PCA')
	#plots.plot_planet_trail_multi_nights(all_nights_planet_trail_PCA[3], extent=[np.min(rv),np.max(rv),0,len(all_nights_phases[3])], xlabel='RV (km/s)', ylabel='Orbital phase', comb_phases=all_nights_phases[3], indiv_phases=all_nights_phases[3:4], RVs=all_nights_RVs[3:4], savefig=False)

	## For (properly?) plotting the planet trail of all nights when there is an overlap in phases
	#comb_phases_interp = np.linspace(np.min(comb_phases), np.max(comb_phases), 1000)
	#comb_planet_trails_interp = np.ones(shape=(len(comb_phases_interp),len(rv)))
	#for kk in range(np.shape(comb_planet_trails_interp)[1]):
	#	comb_planet_trails_interp[:,kk] = interpolate.interpn((np.sort(comb_phases),rv), comb_planet_trails[np.argsort(comb_phases),:], (comb_phases_interp, rv[kk]))
	#plots.plot_planet_trail_multi_nights(comb_planet_trails_interp, extent=[np.min(rv),np.max(rv),0,len(comb_phases_interp)], xlabel='RV (km/s)', ylabel='Orbital phase', comb_phases=comb_phases_interp, indiv_phases=all_nights_phases, RVs=all_nights_RVs, savefig=True, dirname=dirname_for_plots, filename='combined_planet_trail')

	comb_vmap_mask = np.zeros_like(comb_vmap, dtype=bool)
	comb_vmap_noise_region = np.zeros_like(comb_vmap)
	comb_vmap_PCA_noise_region = np.zeros_like(comb_vmap_PCA)
	kpVec_lims = [80,140]
	vRest_lims = [-30,30]
	for ii in range(len(vRest)):
		rr_kpVec = (kpVec < kpVec_lims[0]) | (kpVec >kpVec_lims[1])
		comb_vmap_mask[rr_kpVec,ii] = True
		comb_vmap_noise_region[rr_kpVec,ii] = comb_vmap[rr_kpVec,ii]
		comb_vmap_PCA_noise_region[rr_kpVec,ii] = comb_vmap_PCA[rr_kpVec,ii]
	
	for jj in range(len(kpVec)):
		rr_vRest = (vRest < vRest_lims[0]) | (vRest > vRest_lims[1])
		comb_vmap_mask[jj,rr_vRest] = True
		comb_vmap_noise_region[jj,rr_vRest] = comb_vmap[jj,rr_vRest]
		comb_vmap_PCA_noise_region[jj,rr_vRest] = comb_vmap_PCA[jj,rr_vRest]
	
	## The noise level is calculated as the mean of the points outside the delimited region
	noise = np.mean(comb_vmap[comb_vmap_mask])
	noise_std = np.std(comb_vmap[comb_vmap_mask])
	noise_PCA = np.mean(comb_vmap_PCA[comb_vmap_mask])
	noise_PCA_std = np.std(comb_vmap_PCA[comb_vmap_mask])
	
	## The SNR is simply the comb_vmap values minus this noise divided by the standard deviation of the noise
	comb_SNR = (comb_vmap-noise)/noise_std
	comb_SNR_max = np.max(comb_SNR[~comb_vmap_mask])
	comb_SNR_PCA = (comb_vmap_PCA-noise_PCA)/noise_PCA_std
	comb_SNR_PCA_max = np.max(comb_SNR_PCA[~comb_vmap_mask])

	plots.plot_multi_SNR_squares([comb_SNR, comb_SNR_PCA], extent=[np.min(vRest),np.max(vRest),np.min(kpVec),np.max(kpVec)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', vRest=vRest, kpVec=kpVec, kpVec_lims=kpVec_lims, vRest_lims=vRest_lims, SNR_max=[comb_SNR_max,comb_SNR_PCA_max], vmin=np.min([np.min(comb_SNR),np.min(comb_SNR_PCA)]), vmax=np.max([np.max(comb_SNR),np.max(comb_SNR_PCA)]), savefig=True, dirname=dirname_for_plots, filename='comb_SNRs_same_plot'+filename_suffix)
	plots.plot_multi_SNR_squares([comb_SNR, comb_SNR_PCA], extent=[np.min(vRest),np.max(vRest),np.min(kpVec),np.max(kpVec)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', vRest=vRest, kpVec=kpVec, kpVec_lims=kpVec_lims, vRest_lims=vRest_lims, SNR_max=[comb_SNR_max,comb_SNR_PCA_max], vmin=-6, vmax=6, savefig=True, dirname=dirname_for_plots, filename='comb_SNRs_same_plot'+filename_suffix+'_same_scale')
	plots.plot_multi_vmap_squares([comb_vmap, comb_vmap_PCA], extent=[np.min(vRest),np.max(vRest),np.min(kpVec),np.max(kpVec)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', vRest=vRest, kpVec=kpVec, vmin=np.min([np.min(comb_vmap),np.min(comb_vmap_PCA)]), vmax=np.max([np.max(comb_vmap),np.max(comb_vmap_PCA)]), savefig=True, dirname=dirname_for_plots, filename='comb_vmaps_same_plot'+filename_suffix)
