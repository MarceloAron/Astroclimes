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
import fit_inj_rec

## Some constants
R_Sun = 6.957*1e8	# Solar radius, in m (from NASA Sun Fact Sheet)
R_Jup = 6.9911*1e7	# Jupiter radius, in m (from NASA Jupiter Fact Sheet)

## Color-blind friendly colors: blue, orange, green, pink, brown, purple, gray, red and yellow
CBcols= ['#377EB8', '#FF7F00', '#4DAF4A',
		 '#F781BF', '#A65628', '#984EA3',
		 '#999999', '#E41A1C', '#DEDE00']

def run_inj_rec(main_inj_rec_directory, filename_site_values, R_instrument, n_orders, bad_spec_orders_idxs, bad_obs_idxs, scale_factors, nc_PCA, filename_inj_rec_SNR_AC, filename_inj_rec_SNR_PCA, filename_inj_rec_MCMC_results_AC, filename_inj_rec_MCMC_results_PCA, filename_planet_lam, filename_planet_flux, stelpars, orbpars, plot_telrem_steps_AC, plot_telrem_steps_PCA, plot_indiv_PCA_comps, run_MCMC_AC, run_MCMC_PCA, deep_line_threshold, n_CPUs=1):

	## Loop over all scale factors
	for n, scale_factor in enumerate(scale_factors):
		## Checking if the subdirectory associated with this signal strength multiplier already exists and if not, creating it
		dirname = main_inj_rec_directory+'sf_'+str(scale_factor).replace('.','_')+'x/'
		if not os.path.isdir(dirname):
			call(['mkdir', dirname])

		## Unpacking the observational wavelength and spectra cube (no need to give any sensible arguments when generate=False, besides save_dirname and save_lam_filename and save_spec_filename if file does not have default name)
		all_lam, all_spec = funcs_telrem.get_obs_cube(None, None, None, generate=False, save_dirname=main_inj_rec_directory)

		## Unpacking the model spectra cube (no need to give any sensible arguments when generate=False, besides save_dirname and save_spec_mod_filename if file does not have default name)
		all_spec_mod = funcs_telrem.get_mod_cube(None, None, None, None, None, None, None, None, generate=False, save_dirname=main_inj_rec_directory)
		#all_spec_mod = np.load(main_inj_rec_directory+'all_spec_mod_molecfit.npy', allow_pickle=True)
		
		## Creating a model planetary signal
		T_star = stelpars.Teff
		R_star = stelpars.R
		R_planet = orbpars.Rp
		T0 = orbpars.T0
		P = orbpars.P	
		Kp = orbpars.Kp			
		Vsys = orbpars.Vsys
		planet_cube, planet_mod_lam, planet_model = funcs_telrem.create_planet_signal(	all_lam, filename_planet_lam, filename_planet_flux, filename_site_values, 
																						R_instrument, T_star, R_star, R_planet, scale_factor, T0, P, Kp, Vsys, generate=True, save_dirname=dirname)

		## Injecting the model planetary signal into the observational cube
		all_spec_wps = all_spec*planet_cube

		## In the analysis that follows, we will be working both with the observational spectra with and without an injected planetary signal, 
		## and all variables associated with the former have the "wps" subscript.

		## Deleting unwanted orders (due to too many saturated telluric lines)
		orders = np.arange(1,n_orders+1)
		good_orders = np.delete(orders, bad_spec_orders_idxs)
		all_lam = np.delete(all_lam, bad_spec_orders_idxs, axis=0)
		all_spec = np.delete(all_spec, bad_spec_orders_idxs, axis=0)
		all_spec_wps = np.delete(all_spec_wps, bad_spec_orders_idxs, axis=0)
		all_spec_mod = np.delete(all_spec_mod, bad_spec_orders_idxs, axis=0)
		planet_cube = np.delete(planet_cube, bad_spec_orders_idxs, axis=0)

		## Deleting unwanted observations (due to bad SNR, bad weather, etc.)
		all_lam = np.delete(all_lam, bad_obs_idxs, axis=1)
		all_spec = np.delete(all_spec, bad_obs_idxs, axis=1)
		all_spec_wps = np.delete(all_spec_wps, bad_obs_idxs, axis=1)
		all_spec_mod = np.delete(all_spec_mod, bad_obs_idxs, axis=1)
		planet_cube = np.delete(planet_cube, bad_obs_idxs, axis=1)

		## Calculating the "median divided" observational AND model flux cubes before using our telluric models to remove the telluric lines
		all_spec_mod_corr = all_spec_mod.copy()
		all_spec_corr = all_spec.copy()
		all_spec_wps_corr = all_spec_wps.copy()
		for i in range(np.shape(all_lam)[0]):
			all_spec_mod_corr[i,:,:] = all_spec_mod[i,:,:]/np.nanmedian(all_spec_mod[i,:,:], axis=1)[:,np.newaxis]
			all_spec_corr[i,:,:] = all_spec[i,:,:]/np.nanmedian(all_spec[i,:,:], axis=1)[:,np.newaxis]
			all_spec_wps_corr[i,:,:] = all_spec_wps[i,:,:]/np.nanmedian(all_spec_wps[i,:,:], axis=1)[:,np.newaxis]

		## Masking deep lines
		all_spec_corr_masked, rule_cube = funcs_telrem.mask_cube(all_spec_corr, mask_val=0.0, deep_line_threshold=deep_line_threshold)
		all_spec_wps_corr_masked, rule_cube_wps = funcs_telrem.mask_cube(all_spec_wps_corr, mask_val=0.0, deep_line_threshold=deep_line_threshold)
		
		rule_cube_pca = rule_cube.copy() 			# Explicitly setting these in case I want to use the post-detrending mask later, after which the rule cube won't be the same for Astroclimes and PCA
		rule_cube_wps_pca = rule_cube_wps.copy() 	# Explicitly setting these in case I want to use the post-detrending mask later, after which the rule cube won't be the same for Astroclimes and PCA

		## Removing the telluric and stellar lines from the observational spectra using Astroclimes
		## Dividing the observational spectra by the telluric model spectra should remove the telluric lines
		res_ac = all_spec_corr_masked/all_spec_mod_corr
		res_ac_wps = all_spec_wps_corr_masked/all_spec_mod_corr	

		## Now we divide the residuals from above by their mean over all nights, which should remove the stellar lines
		## Note that this step only works if one has several observations taken on the same night, and even still some of
		## the planetary signal may be removed as well
		res_ac_mean = np.nanmean(res_ac, axis=1)[:,np.newaxis,:]
		res_ac_mean_rem = res_ac/res_ac_mean

		res_ac_wps_mean = np.nanmean(res_ac_wps, axis=1)[:,np.newaxis,:] 	# Mean of the telluric line removed spectra over all observations, to remove stellar lines
		res_ac_wps_mean_rem = res_ac_wps/res_ac_wps_mean
		
		## Getting the PCA "cube" for both the case with and without the injected planetary signal
		all_pca_fit = funcs_telrem.run_pca(all_lam, all_spec_corr, nc_PCA, generate=True, save_dirname=dirname)
		all_pca_fit_wps = funcs_telrem.run_pca(all_lam, all_spec_wps_corr, nc_PCA, generate=True, save_dirname=dirname, save_filename='all_pca_fit_wps')
		## Uncomment the two lines below if you do not wish to save the PCA models
		#call(['rm', dirname+f'all_pca_fit_nc_{nc}.npy'])
		#call(['rm', dirname+f'all_pca_fit_wps_nc_{nc}.npy'])

		## Removing telluric and stellar lines using PCA
		res_pca = all_spec_corr/all_pca_fit
		res_pca_wps = all_spec_wps_corr/all_pca_fit_wps

		## The difference between the residuals with and without the injected planetary signal should yield precisely the planetary signal
		diff_res_ac = res_ac_wps_mean_rem-res_ac_mean_rem
		diff_res_pca = res_pca_wps-res_pca

		## Producing the planet's trail and velocity map (only for the case with the injected planetary signal)
		cs = interpolate.splrep(planet_mod_lam, planet_model)
		rv = np.arange(-300,301,1.5)
		BJD, V_BERV = np.loadtxt(filename_site_values, usecols=(4,15), unpack=True)
		V_BERV = np.delete(V_BERV, bad_obs_idxs, axis=0)
		V_BERV = -V_BERV ## V_BERV is the radial velocity of the Earth with respect to the Sun, but we want the opposite
		BJD = np.delete(BJD, bad_obs_idxs, axis=0)
		phases = ((BJD - T0)%P)/P
		vRest = np.arange(-75,76,1.5)
		kpVec = np.arange(0,251)
		vObs = Vsys + V_BERV

		## For the Astroclimes removal
		#res_ac_wps_mean_rem, rule_cube_wps = funcs_telrem.mask_data_post_pca(res_ac_wps_mean_rem, rule_cube_wps) ## Applying post-detrending mask and updating the rule cube
		tr, ll = funcs_telrem.get_fixed_ccf_cube(res_ac_wps_mean_rem, all_lam[:,0,:], rv, cs, rule_cube_wps)
		planet_trail = np.nansum(tr, axis=0)
		planet_trail_mean_rem = planet_trail - np.mean(planet_trail, axis=1)[:,np.newaxis]
		vmap = funcs_telrem.get_velocity_map(planet_trail,rv,phases,vObs,vRest,kpVec)

		## For the Astroclimes removal (without injected planet signal)
		#tr, ll = funcs_telrem.get_fixed_ccf_cube(res_ac_mean_rem, all_lam[:,0,:], rv, cs, rule_cube)
		#planet_trail_nps = np.nansum(tr, axis=0)
		#planet_trail_nps_mean_rem = planet_trail_nps - np.mean(planet_trail_nps, axis=1)[:,np.newaxis]
		#vmap_nps = funcs_telrem.get_velocity_map(planet_trail_nps,rv,phases,vObs,vRest,kpVec)

		## For the PCA removal
		#res_pca_wps, rule_cube_wps_pca = funcs_telrem.mask_data_post_pca(res_pca_wps, rule_cube_wps) ## Applying post-detrending mask and updating the rule cube
		tr_PCA, ll_PCA = funcs_telrem.get_fixed_ccf_cube(res_pca_wps, all_lam[:,0,:], rv, cs, rule_cube_wps_pca)
		planet_trail_PCA = np.nansum(tr_PCA, axis=0)
		planet_trail_PCA_mean_rem = planet_trail_PCA - np.mean(planet_trail_PCA, axis=1)[:,np.newaxis]
		vmap_PCA = funcs_telrem.get_velocity_map(planet_trail_PCA,rv,phases,vObs,vRest,kpVec)

		## For the PCA removal (without injected planet signal)
		#res_pca, rule_cube_pca = funcs_telrem.mask_data_post_pca(res_pca, rule_cube) ## Applying post-detrending mask and updating the rule cube
		tr_PCA_nps, ll_PCA_nps = funcs_telrem.get_fixed_ccf_cube(res_pca, all_lam[:,0,:], rv, cs, rule_cube_pca)
		planet_trail_PCA_nps = np.nansum(tr_PCA_nps, axis=0)
		planet_trail_PCA_nps_mean_rem = planet_trail_PCA_nps - np.mean(planet_trail_PCA_nps, axis=1)[:,np.newaxis]
		vmap_PCA_nps = funcs_telrem.get_velocity_map(planet_trail_PCA_nps,rv,phases,vObs,vRest,kpVec)
		
		''' Uncomment if you'd like to include injection and recovery tests with Molecfit (you'll need your own Molecfit models)
		## Including Molecfit
		all_spec_mod_molecfit = np.load(main_inj_rec_directory+'all_spec_mod_molecfit.npy', allow_pickle=True)
		all_spec_mod_molecfit = np.delete(all_spec_mod_molecfit, bad_spec_orders_idxs, axis=0)
		all_spec_mod_molecfit = np.delete(all_spec_mod_molecfit, bad_obs_idxs, axis=1)
		all_spec_mod_molecfit_corr = all_spec_mod_molecfit.copy()
		for i in range(np.shape(all_lam)[0]):
			all_spec_mod_molecfit_corr[i,:,:] = all_spec_mod_molecfit[i,:,:]/np.nanmedian(all_spec_mod_molecfit[i,:,:], axis=1)[:,np.newaxis]
		
		## Removing the telluric and stellar lines from the observational spectra using Molecfit
		## Dividing the observational spectra by the telluric model spectra should remove the telluric lines
		res_molecfit = all_spec_corr_masked/all_spec_mod_molecfit_corr
		res_molecfit_wps = all_spec_wps_corr_masked/all_spec_mod_molecfit_corr	
		
		## Now we divide the residuals from above by their mean over all nights, which should remove the stellar lines
		## Note that this step only works if one has several observations taken on the same night, and even still some of
		## the planetary signal may be removed as well
		res_molecfit_mean = np.nanmean(res_molecfit, axis=1)[:,np.newaxis,:]
		res_molecfit_mean_rem = res_molecfit/res_molecfit_mean
		
		res_molecfit_wps_mean = np.nanmean(res_molecfit_wps, axis=1)[:,np.newaxis,:] 	# Mean of the telluric line removed spectra over all observations, to remove stellar lines
		res_molecfit_wps_mean_rem = res_molecfit_wps/res_molecfit_wps_mean
		diff_res_molecfit = res_molecfit_wps_mean_rem-res_molecfit_mean_rem
		
		## For the Molecfit removal
		#res_ac_wps_mean_rem, rule_cube_wps = funcs_telrem.mask_data_post_pca(res_ac_wps_mean_rem, rule_cube_wps) ## Applying post-detrending mask and updating the rule cube
		tr_molecfit, ll_molecfit = funcs_telrem.get_fixed_ccf_cube(res_molecfit_wps_mean_rem, all_lam[:,0,:], rv, cs, rule_cube_wps)
		planet_trail_molecfit = np.nansum(tr_molecfit, axis=0)
		planet_trail_molecfit_mean_rem = planet_trail_molecfit - np.mean(planet_trail_molecfit, axis=1)[:,np.newaxis]
		vmap_molecfit = funcs_telrem.get_velocity_map(planet_trail_molecfit,rv,phases,vObs,vRest,kpVec)
		'''

		## Choosing orders to plot
		plot_orders = good_orders

		## Checking if the plots directory already exists and if not, creating it
		dirname_for_plots = dirname+'plots/'
		if not os.path.isdir(dirname_for_plots):
			call(['mkdir', dirname_for_plots])
		
		## Plotting the observational spectra as it comes, with its median removed, and with its median removed and with deep lines masked, respectively
		plots.plot_cubes(all_lam, all_spec, list(good_orders), plot_orders, 'Observational cube', y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='obs_flux_cube')
		plots.plot_cubes(all_lam, all_spec_corr, list(good_orders), plot_orders, 'Observational cube (median removed)', vmin=0, vmax=1.1, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='obs_flux_cube_medrem')
		plots.plot_cubes(all_lam, all_spec_corr_masked, list(good_orders), plot_orders, 'Observational cube (median removed, deep lines masked)', vmin=0, vmax=1.1, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='obs_flux_cube_medrem_deepmask')

		## Same thing as above but for the case with an injected planetary signal (the difference won't be visible)
		plots.plot_cubes(all_lam, all_spec_wps, list(good_orders), plot_orders, 'Observational cube', y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='wps_obs_flux_cube')
		plots.plot_cubes(all_lam, all_spec_wps_corr, list(good_orders), plot_orders, 'Observational cube (median removed)', vmin=0, vmax=1.1, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='wps_obs_flux_cube_medrem')
		plots.plot_cubes(all_lam, all_spec_wps_corr_masked, list(good_orders), plot_orders, 'Observational cube (median removed, deep lines masked)', vmin=0, vmax=1.1, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='wps_obs_flux_cube_medrem_deepmask')

		## Plotting the model planet cube and the Astroclimes model cube
		plots.plot_cubes(all_lam, planet_cube, list(good_orders), plot_orders, 'Model planetary signal', vmin=0, vmax=1.1, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='mod_planet_cube')
		plots.plot_cubes(all_lam, all_spec_mod_corr, list(good_orders), plot_orders, 'Model cube (Astroclimes)', vmin=0, vmax=1.1, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='mod_flux_cube')
		#plots.plot_cubes(all_lam, all_spec_mod_molecfit_corr, list(good_orders), plot_orders, 'Model cube (Molecfit)', vmin=0, vmax=1.1, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='mod_flux_cube_molecfit')

		## Plotting the Astroclimes residuals, the mean of the residuals, and the residuals with the mean removed, respectively
		plots.plot_cubes(all_lam, res_ac, list(good_orders), plot_orders, 'Residuals (Astroclimes)', vmin=0.90, vmax=1.05, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='res_ac_flux_cube')
		plots.plot_cubes(all_lam, res_ac_mean, list(good_orders), plot_orders, 'Residuals (Astroclimes, mean)', vmin=0.8, vmax=1.05, savefig=True, dirname=dirname_for_plots, filename='res_ac_flux_cube_mean')
		plots.plot_cubes(all_lam, res_ac_mean_rem, list(good_orders), plot_orders, 'Residuals (Astroclimes, mean removed)', vmin=0.9, vmax=1.05, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='res_ac_flux_cube_mean_rem')

		## Same thing as above but for the case with an injected planetary signal (the difference probably won't be visible again here, as the planet signal is buried in noise)
		plots.plot_cubes(all_lam, res_ac_wps, list(good_orders), plot_orders, 'Residuals (Astroclimes)', vmin=0.90, vmax=1.05, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='wps_res_ac_flux_cube')
		plots.plot_cubes(all_lam, res_ac_wps_mean, list(good_orders), plot_orders, 'Residuals (Astroclimes, mean)', vmin=0.8, vmax=1.05, savefig=True, dirname=dirname_for_plots, filename='wps_res_ac_flux_cube_mean')
		plots.plot_cubes(all_lam, res_ac_wps_mean_rem, list(good_orders), plot_orders, 'Residuals (Astroclimes, mean removed)', vmin=0.9, vmax=1.05, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='wps_res_ac_flux_cube_mean_rem')

		## Plotting the PCA model and PCA residuals without and with the injected planetary signal, respectively
		plots.plot_cubes(all_lam, all_pca_fit, list(good_orders), plot_orders, 'PCA cube', y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'pca_flux_cube_nc_{nc_PCA}')
		plots.plot_cubes(all_lam, res_pca, list(good_orders), plot_orders, 'Residuals (PCA)', vmin=0.90, vmax=1.05, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'pca_res_flux_cube_nc_{nc_PCA}')
		plots.plot_cubes(all_lam, res_pca_wps, list(good_orders), plot_orders, 'Residuals (PCA)', vmin=0.90, vmax=1.05, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'wps_pca_res_flux_cube_nc_{nc_PCA}')

		## Plotting the planet trail (as it as mean removed) and velocity map for Astroclimes and PCA (reminder that we only calculated these for the case with an injected signal), respectively
		#plots.plot_squares(planet_trail, extent=[np.min(rv),np.max(rv),0,len(phases)], xlabel='RV (km/s)', ylabel='Orbital phase', y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='planet_trail_sf_'+str(scale_factor).replace('.','_')+'x')
		plots.plot_squares(planet_trail_mean_rem, extent=[np.min(rv),np.max(rv),0,len(phases)], xlabel='RV (km/s)', ylabel='Orbital phase', y_phase=phases, vmin=np.min([np.min(planet_trail_mean_rem),np.min(planet_trail_PCA_mean_rem)]), vmax=np.max([np.max(planet_trail_mean_rem),np.max(planet_trail_PCA_mean_rem)]), savefig=True, dirname=dirname_for_plots, filename='planet_trail_mean_rem_sf_'+str(scale_factor).replace('.','_')+'x')
		plots.plot_squares(vmap, extent=[np.min(vRest),np.max(vRest),0,len(vmap)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', savefig=True, dirname=dirname_for_plots, filename='vmap_sf_'+str(scale_factor).replace('.','_')+'x')
		#plots.plot_squares(vmap_nps, extent=[np.min(vRest),np.max(vRest),0,len(vmap_nps)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', savefig=True, dirname=dirname_for_plots, filename='vmap_nps_sf_'+str(scale_factor).replace('.','_')+'x')

		#plots.plot_squares(planet_trail_PCA, extent=[np.min(rv),np.max(rv),0,len(phases)], xlabel='RV (km/s)', ylabel='Orbital phase', y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'planet_trail_PCA_nc_{nc_PCA}_sf_'+str(scale_factor).replace('.','_')+'x')
		plots.plot_squares(planet_trail_PCA_mean_rem, extent=[np.min(rv),np.max(rv),0,len(phases)], xlabel='RV (km/s)', ylabel='Orbital phase', y_phase=phases, vmin=np.min([np.min(planet_trail_mean_rem),np.min(planet_trail_PCA_mean_rem)]), vmax=np.max([np.max(planet_trail_mean_rem),np.max(planet_trail_PCA_mean_rem)]), savefig=True, dirname=dirname_for_plots, filename=f'planet_trail_PCA_nc_{nc_PCA}_mean_rem_sf_'+str(scale_factor).replace('.','_')+'x')
		plots.plot_squares(vmap_PCA, extent=[np.min(vRest),np.max(vRest),0,len(vmap_PCA)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', savefig=True, dirname=dirname_for_plots, filename=f'vmap_PCA_nc_{nc_PCA}_sf_'+str(scale_factor).replace('.','_')+'x')
		
		plots.plot_multi_squares([planet_trail_mean_rem, planet_trail_PCA_mean_rem], extent=[np.min(rv),np.max(rv),0,len(phases)], xlabel='RV (km/s)', ylabel='Orbital phase', y_phase=phases, vmin=np.min([np.min(planet_trail_mean_rem),np.min(planet_trail_PCA_mean_rem)]), vmax=np.max([np.max(planet_trail_mean_rem),np.max(planet_trail_PCA_mean_rem)]), savefig=True, dirname=dirname_for_plots, filename='all_planet_trails_mean_rem_sf_'+str(scale_factor).replace('.','_')+'x')

		## Plotting the difference between the residuals with and without the injected planetary signal, for Astroclimes and PCA, respectively
		plots.plot_cubes(all_lam, diff_res_ac, list(good_orders), plot_orders, 'Difference residuals (Astroclimes)', vmin=0, vmax=1.1, savefig=True, dirname=dirname_for_plots, filename='diff_res_ac_flux_cube')
		plots.plot_cubes(all_lam, diff_res_pca, list(good_orders), plot_orders, 'Difference residuals (PCA)', vmin=0, vmax=1.1, savefig=True, dirname=dirname_for_plots, filename=f'diff_res_pca_flux_cube_nc_{nc_PCA}')

		## Plotting the injected signal minus the residuals difference for Astroclimes and PCA, respectively
		plots.plot_cubes(all_lam, planet_cube - diff_res_ac, list(good_orders), plot_orders, 'Injected planetary signal - residuals difference (Astroclimes)', vmin=0, vmax=1.1, savefig=True, dirname=dirname_for_plots, filename='planet_signal_minus_diff_res_ac_flux_cube')
		plots.plot_cubes(all_lam, planet_cube - diff_res_pca, list(good_orders), plot_orders, 'Injected planetary signal - residuals difference (PCA)', vmin=0, vmax=1.1, savefig=True, dirname=dirname_for_plots, filename=f'planet_signal_minus_diff_res_pca_flux_cube_nc_{nc_PCA}')
		
		'''
		## Plotting involving Molecfit
		plots.plot_cubes(all_lam, res_molecfit, list(good_orders), plot_orders, 'Residuals (Molecfit)', vmin=0.90, vmax=1.05, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='res_molecfit_flux_cube')
		plots.plot_cubes(all_lam, res_molecfit_mean, list(good_orders), plot_orders, 'Residuals (Molecfit, mean)', vmin=0.8, vmax=1.05, savefig=True, dirname=dirname_for_plots, filename='res_molecfit_flux_cube_mean')
		plots.plot_cubes(all_lam, res_molecfit_mean_rem, list(good_orders), plot_orders, 'Residuals (Molecfit, mean removed)', vmin=0.9, vmax=1.05, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='res_molecfit_flux_cube_mean_rem')
		plots.plot_cubes(all_lam, res_molecfit_wps, list(good_orders), plot_orders, 'Residuals (Molecfit)', vmin=0.90, vmax=1.05, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='wps_res_molecfit_flux_cube')
		plots.plot_cubes(all_lam, res_molecfit_wps_mean, list(good_orders), plot_orders, 'Residuals (Molecfit, mean)', vmin=0.8, vmax=1.05, savefig=True, dirname=dirname_for_plots, filename='wps_res_molecfit_flux_cube_mean')
		plots.plot_cubes(all_lam, res_molecfit_wps_mean_rem, list(good_orders), plot_orders, 'Residuals (Molecfit, mean removed)', vmin=0.9, vmax=1.05, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename='wps_res_molecfit_flux_cube_mean_rem')
		#plots.plot_squares(planet_trail_mean_rem, extent=[np.min(rv),np.max(rv),0,len(phases)], xlabel='RV (km/s)', ylabel='Orbital phase', y_phase=phases, vmin=np.min([np.min(planet_trail_mean_rem),np.min(planet_trail_PCA_mean_rem),np.min(planet_trail_molecfit_mean_rem)]), vmax=np.max([np.max(planet_trail_mean_rem),np.max(planet_trail_PCA_mean_rem),np.max(planet_trail_molecfit_mean_rem)]), savefig=True, dirname=dirname_for_plots, filename='planet_trail_mean_rem_sf_'+str(scale_factor).replace('.','_')+'x')
		#plots.plot_squares(planet_trail_molecfit, extent=[np.min(rv),np.max(rv),0,len(phases)], xlabel='RV (km/s)', ylabel='Orbital phase', y_phase=phases, vmin=np.min([np.min(planet_trail_mean_rem),np.min(planet_trail_PCA_mean_rem),np.min(planet_trail_molecfit_mean_rem)]), vmax=np.max([np.max(planet_trail_mean_rem),np.max(planet_trail_PCA_mean_rem),np.max(planet_trail_molecfit_mean_rem)]), savefig=True, dirname=dirname_for_plots, filename='planet_trail_molecfit_sf_'+str(scale_factor).replace('.','_')+'x')
		plots.plot_squares(planet_trail_molecfit_mean_rem, extent=[np.min(rv),np.max(rv),0,len(phases)], xlabel='RV (km/s)', ylabel='Orbital phase', y_phase=phases, vmin=np.min([np.min(planet_trail_mean_rem),np.min(planet_trail_PCA_mean_rem),np.min(planet_trail_molecfit_mean_rem)]), vmax=np.max([np.max(planet_trail_mean_rem),np.max(planet_trail_PCA_mean_rem),np.max(planet_trail_molecfit_mean_rem)]), savefig=True, dirname=dirname_for_plots, filename='planet_trail_molecfit_mean_rem_sf_'+str(scale_factor).replace('.','_')+'x')
		plots.plot_squares(vmap_molecfit, extent=[np.min(vRest),np.max(vRest),0,len(vmap_molecfit)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', savefig=True, dirname=dirname_for_plots, filename='vmap_molecfit_sf_'+str(scale_factor).replace('.','_')+'x')
		#plots.plot_squares(planet_trail_PCA_mean_rem, extent=[np.min(rv),np.max(rv),0,len(phases)], xlabel='RV (km/s)', ylabel='Orbital phase', y_phase=phases, vmin=np.min([np.min(planet_trail_mean_rem),np.min(planet_trail_PCA_mean_rem),np.min(planet_trail_molecfit_mean_rem)]), vmax=np.max([np.max(planet_trail_mean_rem),np.max(planet_trail_PCA_mean_rem),np.max(planet_trail_molecfit_mean_rem)]), savefig=True, dirname=dirname_for_plots, filename=f'planet_trail_PCA_nc_{nc_PCA}_mean_rem_sf_'+str(scale_factor).replace('.','_')+'x')
		plots.plot_squares(planet_trail_PCA_nps_mean_rem, extent=[np.min(rv),np.max(rv),0,len(phases)], xlabel='RV (km/s)', ylabel='Orbital phase', y_phase=phases, vmin=np.min([np.min(planet_trail_mean_rem),np.min(planet_trail_PCA_mean_rem)]), vmax=np.max([np.max(planet_trail_mean_rem),np.max(planet_trail_PCA_mean_rem)]), savefig=True, dirname=dirname_for_plots, filename=f'planet_trail_PCA_nps_nc_{nc_PCA}_mean_rem_sf_'+str(scale_factor).replace('.','_')+'x')
		#plots.plot_multi_squares([planet_trail_mean_rem, planet_trail_molecfit_mean_rem, planet_trail_PCA_mean_rem], extent=[np.min(rv),np.max(rv),0,len(phases)], xlabel='RV (km/s)', ylabel='Orbital phase', y_phase=phases, vmin=np.min([np.min(planet_trail_mean_rem),np.min(planet_trail_PCA_mean_rem)]), vmax=np.max([np.max(planet_trail_mean_rem),np.max(planet_trail_PCA_mean_rem)]), savefig=True, dirname=dirname_for_plots, filename='all_planet_trails_mean_rem_sf_'+str(scale_factor).replace('.','_')+'x')
		plots.plot_cubes(all_lam, diff_res_molecfit, list(good_orders), plot_orders, 'Difference residuals (Molecfit)', vmin=0, vmax=1.1, savefig=True, dirname=dirname_for_plots, filename='diff_res_molecfit_flux_cube')
		plots.plot_cubes(all_lam, planet_cube - diff_res_molecfit, list(good_orders), plot_orders, 'Injected planetary signal - residuals difference (Molecfit)', vmin=0, vmax=1.1, savefig=True, dirname=dirname_for_plots, filename='planet_signal_minus_diff_res_molecfit_flux_cube')
		'''

		chosen_order_idx = list(good_orders).index(12)
		lam = all_lam[chosen_order_idx,0,:]
		if plot_telrem_steps_AC:
			## Plotting the steps taken on the telluric removal process
			step_1, rule_cube_wps = funcs_telrem.mask_cube(all_spec_wps, np.nan, rule_cube=rule_cube_wps, create_rule_cube=False)
			step_1 = step_1[chosen_order_idx,:,:]
			step_2, rule_cube_wps = funcs_telrem.mask_cube(all_spec_wps_corr, np.nan, rule_cube=rule_cube_wps, create_rule_cube=False)
			step_2 = step_2[chosen_order_idx,:,:]
			step_3, rule_cube_wps = funcs_telrem.mask_cube(all_spec_mod_corr, np.nan, rule_cube=rule_cube_wps, create_rule_cube=False)
			step_3 = step_3[chosen_order_idx,:,:]
			step_4, rule_cube_wps = funcs_telrem.mask_cube(res_ac_wps, np.nan, rule_cube=rule_cube_wps, create_rule_cube=False)
			step_4 = step_4[chosen_order_idx,:,:]
			#step_4p5 = np.concatenate(([res_ac_wps_mean]*99), axis=1)[chosen_order_idx,:,:]
			#step_4p5, rule_cube_wps = funcs_telrem.mask_cube(np.concatenate(([res_ac_wps_mean]*99), axis=1), np.nan, rule_cube=rule_cube_wps, create_rule_cube=False)
			#step_4p5 = step_4p5[chosen_order_idx,:,:]
			step_5, rule_cube_wps = funcs_telrem.mask_cube(res_ac_wps_mean_rem, np.nan, rule_cube=rule_cube_wps, create_rule_cube=False)
			step_5 = step_5[chosen_order_idx,:,:]
			step_6_PCA, rule_cube_wps = funcs_telrem.mask_cube(res_pca_wps, np.nan, rule_cube=rule_cube_wps, create_rule_cube=False)
			step_6_PCA = step_6_PCA[chosen_order_idx,:,:]
			cube_steps = [step_1, step_2, step_3, step_4, step_5, step_6_PCA]
			vmins = [None, None, None, np.nanmedian(step_4.flatten())-3*np.nanstd(step_4.flatten()), np.nanmedian(step_5.flatten())-3*np.nanstd(step_5.flatten()), np.nanmedian(step_5.flatten())-3*np.nanstd(step_5.flatten())]
			vmaxs = [None, None, None, np.nanmedian(step_4.flatten())+3*np.nanstd(step_4.flatten()), np.nanmedian(step_5.flatten())+3*np.nanstd(step_5.flatten()), np.nanmedian(step_5.flatten())+3*np.nanstd(step_5.flatten())]
			plots.plot_cube_steps(lam, cube_steps, vmins=vmins, vmaxs=vmaxs, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'removal_steps_ac_ord_{good_orders[chosen_order_idx]}_with_PCA')

			'''
			## Now with Molecfit residuals as well
			step_6_molecfit, rule_cube_wps = funcs_telrem.mask_cube(res_molecfit_wps_mean_rem, np.nan, rule_cube=rule_cube_wps, create_rule_cube=False)
			step_6_molecfit = step_6_molecfit[chosen_order_idx,:,:]
			cube_steps = [step_1, step_2, step_3, step_4, step_5, step_6_molecfit, step_6_PCA]
			vmins = [None, None, None, np.nanmedian(step_4.flatten())-3*np.nanstd(step_4.flatten()), np.nanmedian(step_5.flatten())-3*np.nanstd(step_5.flatten()), np.nanmedian(step_5.flatten())-3*np.nanstd(step_5.flatten()), np.nanmedian(step_5.flatten())-3*np.nanstd(step_5.flatten())]
			vmaxs = [None, None, None, np.nanmedian(step_4.flatten())+3*np.nanstd(step_4.flatten()), np.nanmedian(step_5.flatten())+3*np.nanstd(step_5.flatten()), np.nanmedian(step_5.flatten())+3*np.nanstd(step_5.flatten()), np.nanmedian(step_5.flatten())+3*np.nanstd(step_5.flatten())]
			plots.plot_cube_steps(lam, cube_steps, vmins=vmins, vmaxs=vmaxs, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'removal_steps_ac_ord_{good_orders[chosen_order_idx]}_with_PCA_and_molecfit')
			'''

		if plot_telrem_steps_PCA:
			## Now for PCA
			step_1, rule_cube_wps_pca = funcs_telrem.mask_cube(all_spec_wps, np.nan, rule_cube=rule_cube_wps_pca, create_rule_cube=False)
			step_1 = step_1[chosen_order_idx,:,:]
			step_2, rule_cube_wps_pca = funcs_telrem.mask_cube(all_spec_wps_corr, np.nan, rule_cube=rule_cube_wps_pca, create_rule_cube=False)
			step_2 = step_2[chosen_order_idx,:,:]
			step_3_PCA, rule_cube_wps_pca = funcs_telrem.mask_cube(all_pca_fit_wps, np.nan, rule_cube=rule_cube_wps_pca, create_rule_cube=False)
			step_3_PCA = step_3_PCA[chosen_order_idx,:,:]
			step_4_PCA, rule_cube_wps_pca = funcs_telrem.mask_cube(res_pca_wps, np.nan, rule_cube=rule_cube_wps_pca, create_rule_cube=False)
			step_4_PCA = step_4_PCA[chosen_order_idx,:,:]
			cube_steps = [step_1, step_2, step_3_PCA, step_4_PCA]
			vmins = [None, None, None, np.nanmedian(step_4_PCA.flatten())-3*np.nanstd(step_4_PCA.flatten())]
			vmaxs = [None, None, None, np.nanmedian(step_4_PCA.flatten())+3*np.nanstd(step_4_PCA.flatten())]
			#plots.plot_cube_steps(lam, cube_steps, vmins=vmins, vmaxs=vmaxs, savefig=True, dirname=dirname_for_plots, filename=f'removal_steps_pca_nc_{nc_PCA}_ord_{good_orders[chosen_order_idx]}')
			plots.plot_cube_steps(lam, cube_steps, vmins=vmins, vmaxs=vmaxs, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'removal_steps_pca_nc_{nc_PCA}_ord_{good_orders[chosen_order_idx]}')

		if plot_indiv_PCA_comps:
			## Plotting the difference between each PCA component
			all_pca_fit_wps_1 = funcs_telrem.run_pca(all_lam, all_spec_wps_corr, 1, generate=False, save_dirname=dirname, save_filename='all_pca_fit_wps')
			all_pca_fit_wps_2 = funcs_telrem.run_pca(all_lam, all_spec_wps_corr, 2, generate=False, save_dirname=dirname, save_filename='all_pca_fit_wps')
			all_pca_fit_wps_3 = funcs_telrem.run_pca(all_lam, all_spec_wps_corr, 3, generate=False, save_dirname=dirname, save_filename='all_pca_fit_wps')
			all_pca_fit_wps_4 = funcs_telrem.run_pca(all_lam, all_spec_wps_corr, 4, generate=False, save_dirname=dirname, save_filename='all_pca_fit_wps')
			all_pca_fit_wps_5 = funcs_telrem.run_pca(all_lam, all_spec_wps_corr, 5, generate=False, save_dirname=dirname, save_filename='all_pca_fit_wps')
			all_pca_fit_wps_6 = funcs_telrem.run_pca(all_lam, all_spec_wps_corr, 6, generate=False, save_dirname=dirname, save_filename='all_pca_fit_wps')
			#all_pca_fit_wps_7 = funcs_telrem.run_pca(all_lam, all_spec_wps_corr, 7, generate=False, save_dirname=dirname, save_filename='all_pca_fit_wps')
			#all_pca_fit_wps_8 = funcs_telrem.run_pca(all_lam, all_spec_wps_corr, 8, generate=False, save_dirname=dirname, save_filename='all_pca_fit_wps')
			all_pca_fit_wps_9 = funcs_telrem.run_pca(all_lam, all_spec_wps_corr, 9, generate=False, save_dirname=dirname, save_filename='all_pca_fit_wps')
			all_pca_fit_wps_10 = funcs_telrem.run_pca(all_lam, all_spec_wps_corr, 10, generate=False, save_dirname=dirname, save_filename='all_pca_fit_wps')
			comp_1 = all_pca_fit_wps_1[chosen_order_idx,:,:] 
			comp_2 = all_pca_fit_wps_2[chosen_order_idx,:,:] - all_pca_fit_wps_1[chosen_order_idx,:,:]
			comp_3 = all_pca_fit_wps_3[chosen_order_idx,:,:] - all_pca_fit_wps_2[chosen_order_idx,:,:]
			comp_4 = all_pca_fit_wps_4[chosen_order_idx,:,:] - all_pca_fit_wps_3[chosen_order_idx,:,:]
			comp_5 = all_pca_fit_wps_5[chosen_order_idx,:,:] - all_pca_fit_wps_4[chosen_order_idx,:,:]
			comp_6 = all_pca_fit_wps_6[chosen_order_idx,:,:] - all_pca_fit_wps_5[chosen_order_idx,:,:]
			#comp_7 = all_pca_fit_wps_7[chosen_order_idx,:,:] - all_pca_fit_wps_6[chosen_order_idx,:,:]
			#comp_8 = all_pca_fit_wps_8[chosen_order_idx,:,:] - all_pca_fit_wps_7[chosen_order_idx,:,:]
			#comp_9 = all_pca_fit_wps_9[chosen_order_idx,:,:] - all_pca_fit_wps_8[chosen_order_idx,:,:]
			comp_10 = all_pca_fit_wps_10[chosen_order_idx,:,:] - all_pca_fit_wps_9[chosen_order_idx,:,:]
			pca_components = [all_spec_wps_corr[chosen_order_idx,:,:], comp_1, comp_2, comp_3, comp_4, comp_5, comp_6, comp_10]
			vmins = [None, None, np.nanmedian(comp_2.flatten())-3*np.nanstd(comp_2.flatten()), np.nanmedian(comp_3.flatten())-3*np.nanstd(comp_3.flatten()), np.nanmedian(comp_4.flatten())-3*np.nanstd(comp_4.flatten()), np.nanmedian(comp_5.flatten())-3*np.nanstd(comp_5.flatten()), np.nanmedian(comp_6.flatten())-3*np.nanstd(comp_6.flatten()), np.nanmedian(comp_10.flatten())-3*np.nanstd(comp_10.flatten())]
			vmaxs = [None, None, np.nanmedian(comp_2.flatten())+3*np.nanstd(comp_2.flatten()), np.nanmedian(comp_3.flatten())+3*np.nanstd(comp_3.flatten()), np.nanmedian(comp_4.flatten())+3*np.nanstd(comp_4.flatten()), np.nanmedian(comp_5.flatten())+3*np.nanstd(comp_5.flatten()), np.nanmedian(comp_6.flatten())+3*np.nanstd(comp_6.flatten()), np.nanmedian(comp_10.flatten())+3*np.nanstd(comp_10.flatten())]
			#plots.plot_cube_steps(lam, pca_components, vmins=vmins, vmaxs=vmaxs, savefig=True, dirname=dirname_for_plots, filename=f'pca_components_ord_{good_orders[chosen_order_idx]}')
			plots.plot_cube_steps(lam, pca_components, vmins=vmins, vmaxs=vmaxs, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'pca_components_ord_{good_orders[chosen_order_idx]}')

		## Plotting the injected planetary model cube compared with the residuals difference for both Astroclimes and PCA
		step_1, rule_cube_wps = funcs_telrem.mask_cube(planet_cube, np.nan, rule_cube=rule_cube_wps, create_rule_cube=False)
		step_1 = step_1[chosen_order_idx,:,:] - np.nanmedian(step_1[chosen_order_idx,:,:])
		step_2, rule_cube_wps = funcs_telrem.mask_cube(diff_res_ac, np.nan, rule_cube=rule_cube_wps, create_rule_cube=False)
		step_2 = step_2[chosen_order_idx,:,:] - np.nanmedian(step_2[chosen_order_idx,:,:], axis=1)[:,np.newaxis]
		step_3, rule_cube_wps_pca = funcs_telrem.mask_cube(diff_res_pca, np.nan, rule_cube=rule_cube_wps_pca, create_rule_cube=False)
		step_3 = step_3[chosen_order_idx,:,:] - np.nanmedian(step_3[chosen_order_idx,:,:], axis=1)[:,np.newaxis]
		## Limiting the wavelength range
		rr = (lam*1e6 > 1.162) & (lam*1e6 < 1.165)
		cube_steps = [step_1[:,rr], step_2[:,rr], step_3[:,rr]]
		vmins = [-0.0003, -0.0003, -0.0003]
		vmaxs = [0.0002, 0.0002, 0.0002]
		#plots.plot_cube_steps(lam[rr], cube_steps, vmins=vmins, vmaxs=vmaxs, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'compare_planet_mod_res_diff_ac_pca_nc_{nc_PCA}_ord_{good_orders[chosen_order_idx]}')
		plots.plot_cube_steps_extra_panel(lam[rr], cube_steps, vmin=-3, vmax=2, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'compare_planet_mod_res_diff_ac_pca_nc_{nc_PCA}_ord_{good_orders[chosen_order_idx]}_extra_panel')

		'''
		## Now with Molecfit as well
		step_4, rule_cube_wps = funcs_telrem.mask_cube(diff_res_molecfit, np.nan, rule_cube=rule_cube_wps, create_rule_cube=False)
		step_4 = step_4[chosen_order_idx,:,:] - np.nanmedian(step_4[chosen_order_idx,:,:], axis=1)[:,np.newaxis]
		## Limiting the wavelength range
		rr = (lam*1e6 > 1.162) & (lam*1e6 < 1.165)
		cube_steps = [step_1[:,rr], step_2[:,rr], step_4[:,rr], step_3[:,rr]]
		vmins = [-0.0003, -0.0003, -0.0003, -0.0003]
		vmaxs = [0.0002, 0.0002, 0.0002, 0.0002]
		#plots.plot_cube_steps(lam[rr], cube_steps, vmins=vmins, vmaxs=vmaxs, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'compare_planet_mod_res_diff_ac_molecfit_pca_nc_{nc_PCA}_ord_{good_orders[chosen_order_idx]}')
		plots.plot_cube_steps_extra_panel_and_molecfit(lam[rr], cube_steps, vmin=-3, vmax=2, y_phase=phases, savefig=True, dirname=dirname_for_plots, filename=f'compare_planet_mod_res_diff_ac_molecfit_pca_nc_{nc_PCA}_ord_{good_orders[chosen_order_idx]}_extra_panel')
		'''

		## Calculating the SNR for the velocity map

		## Defining two lines with roughly the same inclination (determined by visual inspection) 
		## as the planetary signal on the velocity map, one slightly above (l1) and the other slightly below (l2)
		l1 = 1.8*vRest + Kp + 20 ## For the signal on March 26th 2018
		l2 = 1.8*vRest + Kp - 20 ## For the signal on March 26th 2018
		#l1 = -2.2*vRest + Kp + 20 ## For the signal on May 11th 2018
		#l2 = -2.2*vRest + Kp - 20 ## For the signal on May 11th 2018
		#l1 = 3.0*vRest + Kp + 20 ## For the signal on March 12th 2019
		#l2 = 3.0*vRest + Kp - 20 ## For the signal on March 12th 2019
		#l1 = -3.5*vRest + Kp + 20 ## For the signal on March 15th 2019
		#l2 = -3.5*vRest + Kp - 20 ## For the signal on March 15th 2019
		#l1 = 2.0*vRest + Kp + 20 ## For the signal on April 11th 2019
		#l2 = 2.0*vRest + Kp - 20 ## For the signal on April 11th 2019
		#l1 = 50.0*vRest + Kp + 400 ## For the signal on April 14th 2019
		#l2 = 50.0*vRest + Kp - 400 ## For the signal on April 14th 2019
		
		## Creating a boolean mask where True represents the values outside the region
		## delimited by the two lines defined above, along with two arrays containing only the values
		## on the "noise region" and on the "signal region" (mostly for plotting/verification purposes)
		vmap_mask = np.zeros_like(vmap, dtype=bool)
		vmap_noise_region = np.zeros_like(vmap)
		vmap_signal_region = np.zeros_like(vmap)
		vmap_PCA_noise_region = np.zeros_like(vmap_PCA)
		vmap_PCA_signal_region = np.zeros_like(vmap_PCA)
		vmap_PCA_nps_noise_region = np.zeros_like(vmap_PCA_nps)
		vmap_PCA_nps_signal_region = np.zeros_like(vmap_PCA_nps)
		for i in range(len(vRest)):
			rr = (kpVec > l1[i]) | (kpVec < l2[i])
			vmap_mask[rr,i] = True
			vmap_noise_region[rr,i] = vmap[rr,i]
			vmap_signal_region[~rr,i] = vmap[~rr,i]
			vmap_PCA_noise_region[rr,i] = vmap_PCA[rr,i]
			vmap_PCA_signal_region[~rr,i] = vmap_PCA[~rr,i]
			vmap_PCA_nps_noise_region[rr,i] = vmap_PCA_nps[rr,i]
			vmap_PCA_nps_signal_region[~rr,i] = vmap_PCA_nps[~rr,i]
		
		## The noise level is calculated as the mean of the points outside the region delimited by the two lines
		noise = np.mean(vmap[vmap_mask])
		noise_std = np.std(vmap[vmap_mask])
		noise_PCA = np.mean(vmap_PCA[vmap_mask])
		noise_PCA_std = np.std(vmap_PCA[vmap_mask])
		noise_PCA_nps = np.mean(vmap_PCA_nps[vmap_mask])
		noise_PCA_nps_std = np.std(vmap_PCA_nps[vmap_mask])
		
		## The SNR is simply the vmap values minus this noise divided by the standard deviation of the noise, and we calculate 
		## the maximum inside the region delimited by the two lines. We know the planet signal is inside this region, so we want 
		## to quantify it. For cases when the injected signal is small, the global maximum can be elsewhere, which is why we 
		## explicitly calculate the maximum inside the region
		SNR = (vmap-noise)/noise_std
		SNR_max = np.max(SNR[~vmap_mask])
		SNR_PCA = (vmap_PCA-noise_PCA)/noise_PCA_std
		SNR_PCA_max = np.max(SNR_PCA[~vmap_mask])

		SNR_PCA_nps = (vmap_PCA_nps-noise_PCA_nps)/noise_PCA_nps_std
		SNR_PCA_nps_max = np.max(SNR_PCA_nps[~vmap_mask])
		SNR_PCA_DeltaCCF = ((vmap_PCA - vmap_PCA_nps)-noise_PCA_nps)/noise_PCA_nps_std
		SNR_PCA_DeltaCCF_max = np.max(SNR_PCA_DeltaCCF[~vmap_mask])
		
		## Coordinates in the vRest X kpVec space where the maximum occurs (inside the "signal" region)
		vRest_max_idx = np.where(SNR == SNR_max)[1][0]
		kpVec_max_idx = np.where(SNR == SNR_max)[0][0]
		vRest_max_idx_PCA = np.where(SNR_PCA == SNR_PCA_max)[1][0]
		kpVec_max_idx_PCA = np.where(SNR_PCA == SNR_PCA_max)[0][0]
		
		vRest_max_idx_PCA_nps = np.where(SNR_PCA_nps == SNR_PCA_nps_max)[1][0]
		kpVec_max_idx_PCA_nps = np.where(SNR_PCA_nps == SNR_PCA_nps_max)[0][0]
		vRest_max_idx_PCA_DeltaCCF = np.where(SNR_PCA_DeltaCCF == SNR_PCA_DeltaCCF_max)[1][0]
		kpVec_max_idx_PCA_DeltaCCF = np.where(SNR_PCA_DeltaCCF == SNR_PCA_DeltaCCF_max)[0][0]
		
		## Coordinates and value of the exact location of the injected signal
		kpVec_inj_sig_pos_idx = np.where(kpVec == Kp)[0][0]
		vRest_inj_sig_pos_idx = np.where(vRest == 0)[0][0]
		SNR_inj_sig_pos = SNR[kpVec_inj_sig_pos_idx,vRest_inj_sig_pos_idx]
		SNR_PCA_inj_sig_pos = SNR_PCA[kpVec_inj_sig_pos_idx,vRest_inj_sig_pos_idx]
		
		SNR_PCA_nps_inj_sig_pos = SNR_PCA_nps[kpVec_inj_sig_pos_idx,vRest_inj_sig_pos_idx]
		SNR_PCA_DeltaCCF_inj_sig_pos = SNR_PCA_DeltaCCF[kpVec_inj_sig_pos_idx,vRest_inj_sig_pos_idx]
		
		''' With Molecfit
		vmap_molecfit_noise_region = np.zeros_like(vmap_molecfit)
		vmap_molecfit_signal_region = np.zeros_like(vmap_molecfit)
		for i in range(len(vRest)):
			vmap_molecfit_noise_region[rr,i] = vmap_molecfit[rr,i]
			vmap_molecfit_signal_region[~rr,i] = vmap_molecfit[~rr,i]
		noise_molecfit = np.mean(vmap_molecfit[vmap_mask])
		noise_molecfit_std = np.std(vmap_molecfit[vmap_mask])
		SNR_molecfit = (vmap_molecfit-noise_molecfit)/noise_molecfit_std
		SNR_molecfit_max = np.max(SNR_molecfit[~vmap_mask])
		vRest_max_idx_molecfit = np.where(SNR_molecfit == SNR_molecfit_max)[1][0]
		kpVec_max_idx_molecfit = np.where(SNR_molecfit == SNR_molecfit_max)[0][0]
		SNR_molecfit_inj_sig_pos = SNR_molecfit[kpVec_inj_sig_pos_idx,vRest_inj_sig_pos_idx]
		'''

		plots.plot_SNR_square(	SNR, extent=[np.min(vRest),np.max(vRest),np.min(kpVec),np.max(kpVec)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', 
								vRest=vRest, kpVec=kpVec, l1=l1, l2=l2, SNR_max=SNR_max, vRest_max_idx=vRest_max_idx, kpVec_max_idx=kpVec_max_idx, 
								SNR_inj_sig_pos=SNR_inj_sig_pos, vRest_inj_sig_pos_idx=vRest_inj_sig_pos_idx, kpVec_inj_sig_pos_idx=kpVec_inj_sig_pos_idx,
								vmin=np.min([np.min(SNR),np.min(SNR_PCA)]), vmax=np.max([np.max(SNR),np.max(SNR_PCA)]),
								#vmin=np.min([np.min(SNR),np.min(SNR_PCA),np.min(SNR_molecfit)]), vmax=np.max([np.max(SNR),np.max(SNR_PCA),np.max(SNR_molecfit)]),
								savefig=True, dirname=dirname_for_plots, filename='SNR_sf_'+str(scale_factor).replace('.','_')+'x')

		plots.plot_SNR_square(	SNR, extent=[np.min(vRest),np.max(vRest),np.min(kpVec),np.max(kpVec)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', 
								vRest=vRest, kpVec=kpVec, l1=l1, l2=l2, SNR_max=SNR_max, vRest_max_idx=vRest_max_idx, kpVec_max_idx=kpVec_max_idx, 
								SNR_inj_sig_pos=SNR_inj_sig_pos, vRest_inj_sig_pos_idx=vRest_inj_sig_pos_idx, kpVec_inj_sig_pos_idx=kpVec_inj_sig_pos_idx,
								vmin=-6, vmax=14,
								savefig=True, dirname=dirname_for_plots, filename='SNR_sf_'+str(scale_factor).replace('.','_')+'x_same_scale')
		
		plots.plot_SNR_square(	SNR_PCA, extent=[np.min(vRest),np.max(vRest),np.min(kpVec),np.max(kpVec)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', 
								vRest=vRest, kpVec=kpVec, l1=l1, l2=l2, SNR_max=SNR_PCA_max, vRest_max_idx=vRest_max_idx_PCA, kpVec_max_idx=kpVec_max_idx_PCA, 
								SNR_inj_sig_pos=SNR_PCA_inj_sig_pos, vRest_inj_sig_pos_idx=vRest_inj_sig_pos_idx, kpVec_inj_sig_pos_idx=kpVec_inj_sig_pos_idx, 
								vmin=np.min([np.min(SNR),np.min(SNR_PCA)]), vmax=np.max([np.max(SNR),np.max(SNR_PCA)]),
								#vmin=np.min([np.min(SNR),np.min(SNR_PCA),np.min(SNR_molecfit)]), vmax=np.max([np.max(SNR),np.max(SNR_PCA),np.max(SNR_molecfit)]),
								savefig=True, dirname=dirname_for_plots, filename=f'SNR_PCA_nc_{nc_PCA}_sf_'+str(scale_factor).replace('.','_')+'x')

		plots.plot_SNR_square(	SNR_PCA, extent=[np.min(vRest),np.max(vRest),np.min(kpVec),np.max(kpVec)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', 
								vRest=vRest, kpVec=kpVec, l1=l1, l2=l2, SNR_max=SNR_PCA_max, vRest_max_idx=vRest_max_idx_PCA, kpVec_max_idx=kpVec_max_idx_PCA, 
								SNR_inj_sig_pos=SNR_PCA_inj_sig_pos, vRest_inj_sig_pos_idx=vRest_inj_sig_pos_idx, kpVec_inj_sig_pos_idx=kpVec_inj_sig_pos_idx, 
								vmin=-6, vmax=14,
								savefig=True, dirname=dirname_for_plots, filename=f'SNR_PCA_nc_{nc_PCA}_sf_'+str(scale_factor).replace('.','_')+'x_same_scale')

		''' With Molecfit
		#plots.plot_SNR_square(	SNR_molecfit, extent=[np.min(vRest),np.max(vRest),np.min(kpVec),np.max(kpVec)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', 
		#						vRest=vRest, kpVec=kpVec, l1=l1, l2=l2, SNR_max=SNR_molecfit_max, vRest_max_idx=vRest_max_idx_molecfit, kpVec_max_idx=kpVec_max_idx_molecfit, 
		#						SNR_inj_sig_pos=SNR_molecfit_inj_sig_pos, vRest_inj_sig_pos_idx=vRest_inj_sig_pos_idx, kpVec_inj_sig_pos_idx=kpVec_inj_sig_pos_idx,
		#						vmin=np.min([np.min(SNR),np.min(SNR_PCA),np.min(SNR_molecfit)]), vmax=np.max([np.max(SNR),np.max(SNR_PCA),np.max(SNR_molecfit)]),
		#						savefig=True, dirname=dirname_for_plots, filename='SNR_molecfit_sf_'+str(scale_factor).replace('.','_')+'x')
		#
		#plots.plot_SNR_square(	SNR_molecfit, extent=[np.min(vRest),np.max(vRest),np.min(kpVec),np.max(kpVec)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', 
		#						vRest=vRest, kpVec=kpVec, l1=l1, l2=l2, SNR_max=SNR_molecfit_max, vRest_max_idx=vRest_max_idx_molecfit, kpVec_max_idx=kpVec_max_idx_molecfit, 
		#						SNR_inj_sig_pos=SNR_molecfit_inj_sig_pos, vRest_inj_sig_pos_idx=vRest_inj_sig_pos_idx, kpVec_inj_sig_pos_idx=kpVec_inj_sig_pos_idx,
		#						vmin=-6, vmax=14,
		#						savefig=True, dirname=dirname_for_plots, filename='SNR_molecfit_sf_'+str(scale_factor).replace('.','_')+'x_same_scale')
		'''

		plots.plot_SNR_square(	SNR_PCA_nps, extent=[np.min(vRest),np.max(vRest),np.min(kpVec),np.max(kpVec)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', 
								vRest=vRest, kpVec=kpVec, l1=l1, l2=l2, SNR_max=SNR_PCA_nps_max, vRest_max_idx=vRest_max_idx_PCA_nps, kpVec_max_idx=kpVec_max_idx_PCA_nps, 
								SNR_inj_sig_pos=SNR_PCA_nps_inj_sig_pos, vRest_inj_sig_pos_idx=vRest_inj_sig_pos_idx, kpVec_inj_sig_pos_idx=kpVec_inj_sig_pos_idx, 
								savefig=True, dirname=dirname_for_plots, filename=f'SNR_PCA_nps_nc_{nc_PCA}_sf_'+str(scale_factor).replace('.','_')+'x')

		plots.plot_SNR_square(	SNR_PCA_DeltaCCF, extent=[np.min(vRest),np.max(vRest),np.min(kpVec),np.max(kpVec)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', 
								vRest=vRest, kpVec=kpVec, l1=l1, l2=l2, SNR_max=SNR_PCA_DeltaCCF_max, vRest_max_idx=vRest_max_idx_PCA_DeltaCCF, kpVec_max_idx=kpVec_max_idx_PCA_DeltaCCF, 
								SNR_inj_sig_pos=SNR_PCA_DeltaCCF_inj_sig_pos, vRest_inj_sig_pos_idx=vRest_inj_sig_pos_idx, kpVec_inj_sig_pos_idx=kpVec_inj_sig_pos_idx, 
								savefig=True, dirname=dirname_for_plots, filename=f'SNR_PCA_DeltaCCF_nc_{nc_PCA}_sf_'+str(scale_factor).replace('.','_')+'x')
		
		#plots.plot_multi_SNR_squares([SNR, SNR_molecfit, SNR_PCA], extent=[np.min(vRest),np.max(vRest),np.min(kpVec),np.max(kpVec)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', vRest=vRest, kpVec=kpVec, l1=l1, l2=l2, SNR_max=[np.max(SNR), np.max(SNR_molecfit), np.max(SNR_PCA)], vmin=np.min([np.min(SNR),np.min(SNR_PCA),np.min(SNR_molecfit)]), vmax=np.max([np.max(SNR),np.max(SNR_PCA),np.max(SNR_molecfit)]), savefig=True, dirname=dirname_for_plots, filename='SNRs_sf_'+str(scale_factor).replace('.','_')+'x_same_plot')
		plots.plot_multi_SNR_squares([SNR, SNR_PCA], extent=[np.min(vRest),np.max(vRest),np.min(kpVec),np.max(kpVec)], xlabel='Planet rest frame velocities (km/s)', ylabel='Kp (km/s)', vRest=vRest, kpVec=kpVec, l1=l1, l2=l2, SNR_max=[np.max(SNR), np.max(SNR_PCA)], vmin=np.min([np.min(SNR),np.min(SNR_PCA)]), vmax=np.max([np.max(SNR),np.max(SNR_PCA)]), savefig=True, dirname=dirname_for_plots, filename='SNRs_sf_'+str(scale_factor).replace('.','_')+'x_same_plot')
		
		text = (f"{scale_factor:<5.2f} \t {SNR_inj_sig_pos:10.4f} \t {vRest[vRest_inj_sig_pos_idx]:10.4f} \t {kpVec[kpVec_inj_sig_pos_idx]:10.4f} \t " +
				f"{SNR_max:10.4f} \t {vRest[vRest_max_idx]:10.4f} \t {kpVec[kpVec_max_idx]:10.4f}\n")
		txt = open(filename_inj_rec_SNR_AC, 'a')
		txt.write(text)
		txt.close()

		''' With Molecfit
		#text = (f"{scale_factor:<5.2f} \t {SNR_molecfit_inj_sig_pos:10.4f} \t {vRest[vRest_inj_sig_pos_idx]:10.4f} \t {kpVec[kpVec_inj_sig_pos_idx]:10.4f} \t " +
		#		f"{SNR_molecfit_max:10.4f} \t {vRest[vRest_max_idx_molecfit]:10.4f} \t {kpVec[kpVec_max_idx_molecfit]:10.4f}\n")
		#txt = open(filename_inj_rec_SNR_AC.replace('_ac','_molecfit'), 'a')
		#txt.write(text)
		#txt.close()
		'''

		text =	(f"{scale_factor:<5.2f} \t {nc_PCA:<5.0f} \t {SNR_PCA_inj_sig_pos:10.4f} \t {vRest[vRest_inj_sig_pos_idx]:10.4f} \t {kpVec[kpVec_inj_sig_pos_idx]:10.4f} \t " +
				 f"{SNR_PCA_max:10.4f} \t {vRest[vRest_max_idx_PCA]:10.4f} \t {kpVec[kpVec_max_idx_PCA]:10.4f}\n")
		txt = open(filename_inj_rec_SNR_PCA, 'a')
		txt.write(text)
		txt.close()

		text =	(f"{scale_factor:<5.2f} \t {nc_PCA:<5.0f} \t {SNR_PCA_nps_inj_sig_pos:10.4f} \t {vRest[vRest_inj_sig_pos_idx]:10.4f} \t {kpVec[kpVec_inj_sig_pos_idx]:10.4f} \t " +
				 f"{SNR_PCA_nps_max:10.4f} \t {vRest[vRest_max_idx_PCA_nps]:10.4f} \t {kpVec[kpVec_max_idx_PCA_nps]:10.4f}\n")
		txt = open(filename_inj_rec_SNR_PCA.replace('pca_nc_','pca_nps_nc_'), 'a')
		txt.write(text)
		txt.close()

		text =	(f"{scale_factor:<5.2f} \t {nc_PCA:<5.0f} \t {SNR_PCA_DeltaCCF_inj_sig_pos:10.4f} \t {vRest[vRest_inj_sig_pos_idx]:10.4f} \t {kpVec[kpVec_inj_sig_pos_idx]:10.4f} \t " +
				 f"{SNR_PCA_DeltaCCF_max:10.4f} \t {vRest[vRest_max_idx_PCA_DeltaCCF]:10.4f} \t {kpVec[kpVec_max_idx_PCA_DeltaCCF]:10.4f}\n")
		txt = open(filename_inj_rec_SNR_PCA.replace('pca_nc_','pca_DeltaCCF_nc_'), 'a')
		txt.write(text)
		txt.close()
		
		## Running MCMCs to quantify how much Astroclimes and PCA change the injected planetary signal
		
		if run_MCMC_AC:
			## Setting up the MCMC for the Astroclimes residuals
			MCMC_telrem_directory_AC = dirname+'MCMC_telrem_AC/'
			if not os.path.isdir(MCMC_telrem_directory_AC):
				call(['mkdir', MCMC_telrem_directory_AC])

			logsf = 0.01
			mcmc_params = [Vsys, Kp, logsf]
			mcmc_params = np.array(mcmc_params)

			ip, ttt = funcs_telrem.get_gaussian_ip(R_instrument, planet_mod_lam)

			pars_for_plots = ['Vsys', 'Kp', 'log_sf']
			pars_for_plots.append('logp')
			labels = ['Vsys (km/s)', 'Kp (km/s)', 'log_sf']
			labels.append('logp')

			pos = mcmc_params + 0.1*mcmc_params*np.random.randn(10, len(mcmc_params))
			nwalkers, ndim = pos.shape
			max_steps = 5000
			check_each_n_steps = 100
			index = 0
			autocorr = np.zeros(max_steps)
			old_tau = -np.inf
			taufile = open(MCMC_telrem_directory_AC+'autocorr.txt', 'a')

			backfile = MCMC_telrem_directory_AC+'backupmcmc.h5'
			backend = emcee.backends.HDFBackend(backfile)
			backend.reset(nwalkers, ndim)

			## Saving the necessary parameters into a temporary .npy file so they can be unpacked at the start of the MCMC 
			## instead of giving them as parameters for the log_probability function (this saves up computational time)
			temp_npy_lam = f'temp_lam_MCMC_AC_sf_{scale_factor}x.npy'
			temp_npy_res = f'temp_res_MCMC_AC_sf_{scale_factor}x.npy'
			np.save(temp_npy_lam, all_lam[:,0,:])
			np.save(temp_npy_res, res_ac_wps_mean_rem)

			fit_inj_rec.unpack_global_vars(temp_npy_lam, temp_npy_res)

			## Running the MCMC
			with Pool(processes=n_CPUs) as pool:
				#sampler = emcee.EnsembleSampler(nwalkers, ndim, fit_inj_rec.log_probability, 
				#								args=(planet_mod_lam, planet_model, all_lam[:,0,:], res_ac_wps_mean_rem, phases, V_BERV, ip, rule_cube_wps),
				#								backend=backend, pool=pool)
				sampler = emcee.EnsembleSampler(nwalkers, ndim, fit_inj_rec.log_probability, 
								args=(planet_mod_lam, planet_model, phases, V_BERV, ip, rule_cube_wps),
								backend=backend, pool=pool)
				for sample in sampler.sample(pos, iterations=max_steps, progress=True):
					if sampler.iteration % check_each_n_steps:
						continue
					
					tau = sampler.get_autocorr_time(tol=0)
					autocorr[index] = np.mean(tau)
					taufile.write(f'{sampler.iteration} \t {autocorr[index]:.4f} \n')
					index += 1

					converged = (np.all(tau*100 < sampler.iteration)) & (np.all(np.abs(old_tau - tau)/tau < 0.01))
					if converged:
						print('MCMC converged after {} steps'.format(sampler.iteration))
						break
					old_tau = tau
			taufile.close()

			call(['rm', temp_npy_lam])
			call(['rm', temp_npy_res])
			
			#sampler = funcs.get_mcmc_result(MCMC_telrem_directory_AC)
			samples = sampler.get_chain()
			logp = sampler.get_log_prob()
			bad_walkers = np.where(logp[-1,:] < (np.median(logp[-1,:]) - 3*np.std(logp[-1,:])))[0]
			#bad_walkers = np.where(logp[-1,:] < (np.median(logp[-1,:]) - 1*np.std(logp[-1,:])))[0]
			nsteps = np.shape(samples)[0]
			discard, thin = int(nsteps/2), 1
			flat_samples = sampler.get_chain(discard=discard, thin=thin, flat=True)
			logp = sampler.get_log_prob()
			logpflat = sampler.get_log_prob(discard=discard, thin=thin, flat=True)
			samples, flat_samples, logp, logpflat = funcs.delete_bad_walkers(samples, flat_samples, logp, logpflat, bad_walkers, nsteps, nwalkers, discard, thin)

			## Unpacking the results from the MCMC
			mcmc_median = np.median(flat_samples, axis=0)
			mcmc_std = np.std(flat_samples, axis=0)

			smps = np.concatenate((samples, logp[:, :, np.newaxis]), axis=2)
			flat_smps = np.concatenate((flat_samples, logpflat[:, np.newaxis]), axis=1)

			init_vals = [Vsys, Kp, logsf]
			init_vals.append(logpflat[0])
			plots.make_chain_plot(smps, len(pars_for_plots), labels, savefig=True, dirname=MCMC_telrem_directory_AC)
			plots.make_corner_plot(len(pars_for_plots), flat_smps, logpflat, labels, init_vals, savefig=True, dirname=MCMC_telrem_directory_AC)
			plots.make_autocorr_plot(check_each_n_steps, dirname=MCMC_telrem_directory_AC)

			logp_med, u_logp_med = np.median(logpflat, axis=0), np.std(logpflat, axis=0)

			Vsys_f = mcmc_median[0]
			u_Vsys_f = mcmc_std[0]
			Kp_f = mcmc_median[1]
			u_Kp_f = mcmc_std[1]
			logsf_f = mcmc_median[2]
			u_logsf_f = mcmc_std[2]
			sf_f = 10**logsf_f
			u_sf_f = u_logsf_f*sf_f*np.log(10) 	## Propagration of uncertainty, from f(x) = 10^x

			text =	(f"{scale_factor:<10.2f} \t " +
					f"{sf_f:10.4f} \t {u_sf_f:6.4f} \t " +
					f"{logsf_f:10.4f} \t {u_logsf_f:6.4f} \t " +
					f"{Kp_f:10.4f} \t {u_Kp_f:6.4f} \t " +
					f"{Vsys_f:10.4f} \t {u_Vsys_f:6.4f} \t " +
					f"{logp_med:10.4f} \t {u_logp_med:10.4f}\n")
			txt = open(filename_inj_rec_MCMC_results_AC, 'a')
			txt.write(text)
			txt.close()
		
		''' With Molecfit
		run_MCMC_MOLECFIT = False
		if run_MCMC_MOLECFIT:
			## Setting up the MCMC for the Astroclimes residuals
			MCMC_telrem_directory_MOLECFIT = dirname+'MCMC_telrem_MOLECFIT/'
			if not os.path.isdir(MCMC_telrem_directory_MOLECFIT):
				call(['mkdir', MCMC_telrem_directory_MOLECFIT])

			logsf = 0.01
			mcmc_params = [Vsys, Kp, logsf]
			mcmc_params = np.array(mcmc_params)

			ip, ttt = funcs_telrem.get_gaussian_ip(R_instrument, planet_mod_lam)

			pars_for_plots = ['Vsys', 'Kp', 'log_sf']
			pars_for_plots.append('logp')
			labels = ['Vsys (km/s)', 'Kp (km/s)', 'log_sf']
			labels.append('logp')

			pos = mcmc_params + 0.1*mcmc_params*np.random.randn(10, len(mcmc_params))
			nwalkers, ndim = pos.shape
			max_steps = 5000
			check_each_n_steps = 100
			index = 0
			autocorr = np.zeros(max_steps)
			old_tau = -np.inf
			taufile = open(MCMC_telrem_directory_MOLECFIT+'autocorr.txt', 'a')

			backfile = MCMC_telrem_directory_MOLECFIT+'backupmcmc.h5'
			backend = emcee.backends.HDFBackend(backfile)
			backend.reset(nwalkers, ndim)

			## Saving the necessary parameters into a temporary .npy file so they can be unpacked at the start of the MCMC 
			## instead of giving them as parameters for the log_probability function (this saves up computational time)
			temp_npy_lam = f'temp_lam_MCMC_MOLECFIT_sf_{scale_factor}x.npy'
			temp_npy_res = f'temp_res_MCMC_MOLECFIT_sf_{scale_factor}x.npy'
			np.save(temp_npy_lam, all_lam[:,0,:])
			np.save(temp_npy_res, res_molecfit_wps_mean_rem)

			fit_inj_rec.unpack_global_vars(temp_npy_lam, temp_npy_res)

			## Running the MCMC
			with Pool(processes=n_CPUs) as pool:
				#sampler = emcee.EnsembleSampler(nwalkers, ndim, fit_inj_rec.log_probability, 
				#								args=(planet_mod_lam, planet_model, all_lam[:,0,:], res_molecfit_wps_mean_rem, phases, V_BERV, ip, rule_cube_wps),
				#								backend=backend, pool=pool)
				sampler = emcee.EnsembleSampler(nwalkers, ndim, fit_inj_rec.log_probability, 
								args=(planet_mod_lam, planet_model, phases, V_BERV, ip, rule_cube_wps),
								backend=backend, pool=pool)
				for sample in sampler.sample(pos, iterations=max_steps, progress=True):
					if sampler.iteration % check_each_n_steps:
						continue
					
					tau = sampler.get_autocorr_time(tol=0)
					autocorr[index] = np.mean(tau)
					taufile.write(f'{sampler.iteration} \t {autocorr[index]:.4f} \n')
					index += 1

					converged = (np.all(tau*100 < sampler.iteration)) & (np.all(np.abs(old_tau - tau)/tau < 0.01))
					if converged:
						print('MCMC converged after {} steps'.format(sampler.iteration))
						break
					old_tau = tau
			taufile.close()

			call(['rm', temp_npy_lam])
			call(['rm', temp_npy_res])
			
			#sampler = funcs.get_mcmc_result(MCMC_telrem_directory_MOLECFIT)
			samples = sampler.get_chain()
			logp = sampler.get_log_prob()
			bad_walkers = np.where(logp[-1,:] < (np.median(logp[-1,:]) - 3*np.std(logp[-1,:])))[0]
			#bad_walkers = np.where(logp[-1,:] < (np.median(logp[-1,:]) - 1*np.std(logp[-1,:])))[0]
			nsteps = np.shape(samples)[0]
			discard, thin = int(nsteps/2), 1
			flat_samples = sampler.get_chain(discard=discard, thin=thin, flat=True)
			logp = sampler.get_log_prob()
			logpflat = sampler.get_log_prob(discard=discard, thin=thin, flat=True)
			samples, flat_samples, logp, logpflat = funcs.delete_bad_walkers(samples, flat_samples, logp, logpflat, bad_walkers, nsteps, nwalkers, discard, thin)

			## Unpacking the results from the MCMC
			mcmc_median = np.median(flat_samples, axis=0)
			mcmc_std = np.std(flat_samples, axis=0)

			smps = np.concatenate((samples, logp[:, :, np.newaxis]), axis=2)
			flat_smps = np.concatenate((flat_samples, logpflat[:, np.newaxis]), axis=1)

			init_vals = [Vsys, Kp, logsf]
			init_vals.append(logpflat[0])
			plots.make_chain_plot(smps, len(pars_for_plots), labels, savefig=True, dirname=MCMC_telrem_directory_MOLECFIT)
			plots.make_corner_plot(len(pars_for_plots), flat_smps, logpflat, labels, init_vals, savefig=True, dirname=MCMC_telrem_directory_MOLECFIT)
			plots.make_autocorr_plot(check_each_n_steps, dirname=MCMC_telrem_directory_MOLECFIT)

			logp_med, u_logp_med = np.median(logpflat, axis=0), np.std(logpflat, axis=0)

			Vsys_f = mcmc_median[0]
			u_Vsys_f = mcmc_std[0]
			Kp_f = mcmc_median[1]
			u_Kp_f = mcmc_std[1]
			logsf_f = mcmc_median[2]
			u_logsf_f = mcmc_std[2]
			sf_f = 10**logsf_f
			u_sf_f = u_logsf_f*sf_f*np.log(10) 	## Propagration of uncertainty, from f(x) = 10^x

			text =	(f"{scale_factor:<10.2f} \t " +
					f"{sf_f:10.4f} \t {u_sf_f:6.4f} \t " +
					f"{logsf_f:10.4f} \t {u_logsf_f:6.4f} \t " +
					f"{Kp_f:10.4f} \t {u_Kp_f:6.4f} \t " +
					f"{Vsys_f:10.4f} \t {u_Vsys_f:6.4f} \t " +
					f"{logp_med:10.4f} \t {u_logp_med:10.4f}\n")
			txt = open(filename_inj_rec_MCMC_results_AC.replace('_ac','_molecfit'), 'a')
			txt.write(text)
			txt.close()
		'''
		
		if run_MCMC_PCA:
			## Now for PCA
			MCMC_telrem_directory_PCA = dirname+f'MCMC_telrem_PCA_nc_{nc_PCA}/'
			if not os.path.isdir(MCMC_telrem_directory_PCA):
				call(['mkdir', MCMC_telrem_directory_PCA])

			logsf = 0.01
			mcmc_params = [Vsys, Kp, logsf]
			mcmc_params = np.array(mcmc_params)

			ip, ttt = funcs_telrem.get_gaussian_ip(R_instrument, planet_mod_lam)

			pars_for_plots = ['Vsys', 'Kp', 'log_sf']
			pars_for_plots.append('logp')
			labels = ['Vsys (km/s)', 'Kp (km/s)', 'log_sf']
			labels.append('logp')

			pos = mcmc_params + 0.1*mcmc_params*np.random.randn(10, len(mcmc_params))
			nwalkers, ndim = pos.shape
			max_steps = 5000
			check_each_n_steps = 100
			index = 0
			autocorr = np.zeros(max_steps)
			old_tau = -np.inf
			taufile = open(MCMC_telrem_directory_PCA+'autocorr.txt', 'a')

			backfile = MCMC_telrem_directory_PCA+'backupmcmc.h5'
			backend = emcee.backends.HDFBackend(backfile)
			backend.reset(nwalkers, ndim)

			## Saving the necessary parameters into a temporary .npy file so they can be unpacked at the start of the MCMC 
			## instead of giving them as parameters for the log_probability function (this saves up computational time)
			temp_npy_lam = f'temp_lam_MCMC_PCA_nc_{nc_PCA}_sf_{scale_factor}x.npy'
			temp_npy_res = f'temp_res_MCMC_PCA_nc_{nc_PCA}_sf_{scale_factor}x.npy'
			np.save(temp_npy_lam, all_lam[:,0,:])
			np.save(temp_npy_res, res_pca_wps)

			fit_inj_rec.unpack_global_vars(temp_npy_lam, temp_npy_res)

			## Running the MCMC
			with Pool(processes=n_CPUs) as pool:
				#sampler = emcee.EnsembleSampler(nwalkers, ndim, fit_inj_rec.log_probability, 
				#								args=(planet_mod_lam, planet_model, all_lam[:,0,:], res_pca_wps, phases, V_BERV, ip, rule_cube_wps_pca),
				#								backend=backend)
				sampler = emcee.EnsembleSampler(nwalkers, ndim, fit_inj_rec.log_probability, 
												args=(planet_mod_lam, planet_model, phases, V_BERV, ip, rule_cube_wps_pca),
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

			call(['rm', temp_npy_lam])
			call(['rm', temp_npy_res])

			#sampler = funcs.get_mcmc_result(MCMC_telrem_directory_PCA)
			samples = sampler.get_chain()
			logp = sampler.get_log_prob()
			bad_walkers = np.where(logp[-1,:] < (np.median(logp[-1,:]) - 3*np.std(logp[-1,:])))[0]
			#bad_walkers = np.where(logp[-1,:] < (np.median(logp[-1,:]) - 1*np.std(logp[-1,:])))[0]
			nsteps = np.shape(samples)[0]
			discard, thin = int(nsteps/2), 1
			flat_samples = sampler.get_chain(discard=discard, thin=thin, flat=True)
			logp = sampler.get_log_prob()
			logpflat = sampler.get_log_prob(discard=discard, thin=thin, flat=True)
			samples, flat_samples, logp, logpflat = funcs.delete_bad_walkers(samples, flat_samples, logp, logpflat, bad_walkers, nsteps, nwalkers, discard, thin)

			## Unpacking the results from the MCMC
			mcmc_median = np.median(flat_samples, axis=0)
			mcmc_std = np.std(flat_samples, axis=0)

			smps = np.concatenate((samples, logp[:, :, np.newaxis]), axis=2)
			flat_smps = np.concatenate((flat_samples, logpflat[:, np.newaxis]), axis=1)

			init_vals = [Vsys, Kp, logsf]
			init_vals.append(logpflat[0])
			plots.make_chain_plot(smps, len(pars_for_plots), labels, savefig=True, dirname=MCMC_telrem_directory_PCA)
			plots.make_corner_plot(len(pars_for_plots), flat_smps, logpflat, labels, init_vals, savefig=True, dirname=MCMC_telrem_directory_PCA)
			plots.make_autocorr_plot(check_each_n_steps, dirname=MCMC_telrem_directory_PCA)
			
			logp_med, u_logp_med = np.median(logpflat, axis=0), np.std(logpflat, axis=0)

			Vsys_f_PCA = mcmc_median[0]
			u_Vsys_f_PCA = mcmc_std[0]
			Kp_f_PCA = mcmc_median[1]
			u_Kp_f_PCA = mcmc_std[1]
			logsf_f_PCA = mcmc_median[2]
			u_logsf_f_PCA = mcmc_std[2]
			sf_f_PCA = 10**logsf_f_PCA
			u_sf_f_PCA = u_logsf_f_PCA*sf_f_PCA*np.log(10) 	## Propagration of uncertainty, from f(x) = 10^x

			text =	(f"{scale_factor:<10.2f} \t {nc_PCA:<5.0f} \t" +
					f"{sf_f_PCA:10.4f} \t {u_sf_f_PCA:6.4f} \t " +
					f"{logsf_f_PCA:10.4f} \t {u_logsf_f_PCA:6.4f} \t " +
					f"{Kp_f_PCA:10.4f} \t {u_Kp_f_PCA:6.4f} \t " +
					f"{Vsys_f_PCA:10.4f} \t {u_Vsys_f_PCA:6.4f} \t " +
					f"{logp_med:10.4f} \t {u_logp_med:10.4f}\n")
			txt = open(filename_inj_rec_MCMC_results_PCA, 'a')
			txt.write(text)
			txt.close()
		
		print(f'Progress:{n+1}/{len(scale_factors)}')