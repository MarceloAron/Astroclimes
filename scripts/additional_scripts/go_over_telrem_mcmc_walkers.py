'''
Script to go over the results from the telluric removal MCMC to check and remove any walkers that may not have converged
'''
# =====================================================================================
# Basic packages
# =====================================================================================
import numpy as np
import sys

# =====================================================================================
# Scripts
# =====================================================================================
sys.path.append('/home/marceloaron/MarceloAron/Astroclimes/scripts/')
import funcs
import plots

home_directory = '/media/marceloaron/New Volume/PhD/thesis_work/'

scale_factors = [1,2,3,3.25,3.50,3.75,4,4.25,4.50,4.75,5,6,7,8,9,10]
nc_PCA = 2 	## Number of components on the PCA

main_telrem_directory = home_directory+'telluric_removal/TauBoo/2018_03_26/rerun_test/'

Kp = 110 	
Vsys = -16.4 

for n, scale_factor in enumerate(scale_factors):
	dirname = main_telrem_directory+'sf_'+str(scale_factor).replace('.','_')+'x/'

	dirname_MCMC_AC = dirname+'MCMC_telrem_AC/'

	logsf = 0.01
	mcmc_params = [Vsys, Kp, logsf]
	mcmc_params = np.array(mcmc_params)

	mcmc_par_names = ['Vsys', 'Kp', 'log_sf']
	pars_for_plots = ['Vsys', 'Kp', 'log_sf']
	pars_for_plots.append('logp')
	labels = ['Vsys (km/s)', 'Kp (km/s)', 'log_sf']
	labels.append('logp')

	pos = mcmc_params + 0.1*mcmc_params*np.random.randn(10, len(mcmc_params))
	nwalkers, ndim = pos.shape
	max_steps = 5000

	sampler = funcs.get_mcmc_result(dirname_MCMC_AC)
	samples = sampler.get_chain()
	logp = sampler.get_log_prob()
	nsteps = np.shape(samples)[0]
	discard, thin = int(nsteps/2), 1
	flat_samples = sampler.get_chain(discard=discard, thin=thin, flat=True)
	logp = sampler.get_log_prob()
	logpflat = sampler.get_log_prob(discard=discard, thin=thin, flat=True)
	#bad_walkers = np.where(logp[-1,:] < (np.median(logp[-1,:]) - 1*np.std(logp[-1,:])))[0]
	bad_walkers = np.where(logp[-1,:] < (np.median(logp[-1,:]) - 3*np.std(logp[-1,:])))[0]
	samples, flat_samples, logp, logpflat = funcs.delete_bad_walkers(samples, flat_samples, logp, logpflat, bad_walkers, nsteps, nwalkers, discard, thin)

	mcmc_median = np.median(flat_samples, axis=0)
	mcmc_std = np.std(flat_samples, axis=0)

	smps = np.concatenate((samples, logp[:, :, np.newaxis]), axis=2)
	flat_smps = np.concatenate((flat_samples, logpflat[:, np.newaxis]), axis=1)

	init_vals = [Vsys, Kp, logsf]
	init_vals.append(logpflat[0])
	plots.make_chain_plot(smps, len(pars_for_plots), labels, savefig=True, dirname=dirname_MCMC_AC)
	plots.make_corner_plot(len(pars_for_plots), flat_smps, logpflat, labels, init_vals, savefig=True, dirname=dirname_MCMC_AC)
	#plots.make_autocorr_plot(check_each_n_steps, dirname=dirname_MCMC_AC)

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
	txt = open(main_telrem_directory+'results_telrem_mcmc_ac.txt', 'a')
	txt.write(text)
	txt.close()

	dirname_MCMC_PCA = dirname+f'MCMC_telrem_PCA_nc_{nc_PCA}/'

	logsf = 0.01
	mcmc_params = [Vsys, Kp, logsf]
	mcmc_params = np.array(mcmc_params)

	mcmc_par_names = ['Vsys', 'Kp', 'log_sf']
	pars_for_plots = ['Vsys', 'Kp', 'log_sf']
	pars_for_plots.append('logp')
	labels = ['Vsys (km/s)', 'Kp (km/s)', 'log_sf']
	labels.append('logp')

	pos = mcmc_params + 0.1*mcmc_params*np.random.randn(10, len(mcmc_params))
	nwalkers, ndim = pos.shape
	max_steps = 5000

	sampler = funcs.get_mcmc_result(dirname_MCMC_PCA)
	samples = sampler.get_chain()
	logp = sampler.get_log_prob()
	nsteps = np.shape(samples)[0]
	discard, thin = int(nsteps/2), 1
	flat_samples = sampler.get_chain(discard=discard, thin=thin, flat=True)
	logp = sampler.get_log_prob()
	logpflat = sampler.get_log_prob(discard=discard, thin=thin, flat=True)
	#bad_walkers = np.where(logp[-1,:] < (np.median(logp[-1,:]) - 1*np.std(logp[-1,:])))[0]
	bad_walkers = np.where(logp[-1,:] < (np.median(logp[-1,:]) - 3*np.std(logp[-1,:])))[0]
	samples, flat_samples, logp, logpflat = funcs.delete_bad_walkers(samples, flat_samples, logp, logpflat, bad_walkers, nsteps, nwalkers, discard, thin)

	mcmc_median = np.median(flat_samples, axis=0)
	mcmc_std = np.std(flat_samples, axis=0)

	smps = np.concatenate((samples, logp[:, :, np.newaxis]), axis=2)
	flat_smps = np.concatenate((flat_samples, logpflat[:, np.newaxis]), axis=1)

	init_vals = [Vsys, Kp, logsf]
	init_vals.append(logpflat[0])
	plots.make_chain_plot(smps, len(pars_for_plots), labels, savefig=True, dirname=dirname_MCMC_PCA)
	plots.make_corner_plot(len(pars_for_plots), flat_smps, logpflat, labels, init_vals, savefig=True, dirname=dirname_MCMC_PCA)
	#plots.make_autocorr_plot(check_each_n_steps, dirname=dirname_MCMC_PCA)

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
	txt = open(main_telrem_directory+f'results_telrem_mcmc_pca_nc_{nc_PCA}.txt', 'a')
	txt.write(text)
	txt.close()
