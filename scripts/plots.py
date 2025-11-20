# =====================================================================================
# Basic packages
# =====================================================================================
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import corner

# =====================================================================================
# Scripts
# =====================================================================================
from gen_funcs import roundto1

## Color-blind friendly colors: blue, orange, green, pink, brown, purple, gray, red and yellow
CBcols= ['#377EB8', '#FF7F00', '#4DAF4A',
		 '#F781BF', '#A65628', '#984EA3',
		 '#999999', '#E41A1C', '#DEDE00']

## Added extra: blue, light blue, orange, light orange, green, light green, red, light red,
## purple, light purple, brown, light brown, pink, light pink, gray, light gray, beige-yellow,
## light beige-yellow, turquoise, light turquoise
CBcols2 = [	'#1F77B4', '#AEC7E8', '#FF7F0E', 
			'#FFBB78', '#2CA02C', '#98DF8A', 
			'#D62728', '#FF9896', '#9467BD', 
			'#C5B0D5', '#8C564B', '#C49C94', 
			'#E377C2', '#F7B6D2', '#7F7F7F', 
			'#C7C7C7', '#BCBD22', '#DBDB8D', 
			'#17BECF', '#9EDAE5']

rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)

## Define your home directory
home_directory = '/home/marceloaron/MarceloAron/PhD/thesis_work/Astroclimes/'

## Define default directory where plots will be saved
save_directory = home_directory+'plots/'

def plot_em_line_mask(lam, spec, spec_skymodel, em_lines_mask, savefig=False, dirname=save_directory, filenames=['em_lines_mask', 'spec_em_lines_mask']):
	lam_em = lam[em_lines_mask]
	plt.figure()
	plt.plot(lam*1e6, spec_skymodel)
	plt.plot(lam_em*1e6, spec_skymodel[em_lines_mask], 'o-')
	plt.xlabel('Wavelength ('+r'$\mu$'+'m)')
	plt.ylabel('Spectra')
	plt.tick_params(axis='x')
	plt.tick_params(axis='y')

	if savefig:		
		plt.savefig(dirname+filenames[0]+'.pdf')
		plt.savefig(dirname+filenames[0]+'.jpeg', dpi=300)

	#plt.show()
	plt.close()

	plt.figure()
	plt.plot(lam*1e6, spec)
	plt.plot(lam_em*1e6, spec[em_lines_mask], 'o-')
	plt.xlabel('Wavelength ('+r'$\mu$'+'m)')
	plt.ylabel('Spectra')
	plt.tick_params(axis='x')
	plt.tick_params(axis='y')

	if savefig:		
		plt.savefig(dirname+filenames[1]+'.pdf')
		plt.savefig(dirname+filenames[1]+'.jpeg', dpi=300)

	#plt.show()
	plt.close()

def plot_norm_mask(lam, spec, filt, spec_clip, idx_out, savefig=False, dirname=save_directory, filename='norm_mask'):
	fig, axes = plt.subplots(3, figsize=(18,14), sharex=True)

	axes[0].plot(lam, spec, c=CBcols[0])
	axes[0].plot(lam, filt, c='black')
	#axes[0].set_ylim(-0.1, 1.1)
	axins = axes[0].inset_axes([0.05, 0.1, 0.5, 0.5])
	axins.plot(lam, spec, c=CBcols[0])
	axins.plot(lam, filt, c='black')
	axins.set_xlim(np.mean(lam)+4e-3, np.mean(lam)+4e-3+2e-3)
	#axins.set_ylim(0.87, 1.02)
	axins.set_ylim(0.87*np.max(spec), 1.02*np.max(spec))
	axins.set_xticks(ticks=[])
	axins.set_yticks(ticks=[])
	mark_inset(axes[0], axins, loc1=2, loc2=1, fc='None', ec='k', zorder=10)
	axes[0].set_ylabel('Spectra', fontsize=26)
	axes[0].tick_params(axis='x', labelsize=26)
	axes[0].tick_params(axis='y', labelsize=26)

	MAD = np.median(np.abs(spec_clip))
	axes[1].axhline(1*MAD, ls='dashed', c='grey')
	axes[1].axhline(-1*MAD, ls='dashed', c='gray')
	axes[1].plot(lam, spec_clip, c=CBcols[0], label='Spectra')
	axes[1].plot(lam[idx_out], spec_clip[idx_out], c=CBcols[1], label='Spectra without lines')
	axins = axes[1].inset_axes([0.05, 0.1, 0.5, 0.5])
	axins.axhline(1*MAD, ls='dashed', c='grey')
	axins.axhline(-1*MAD, ls='dashed', c='gray')
	axins.plot(lam, spec_clip)
	axins.plot(lam[idx_out], spec_clip[idx_out], 'o-', ms=4.)
	axins.set_xlim(np.mean(lam)+4e-3, np.mean(lam)+4e-3+2e-3)
	axins.set_ylim(-0.1,0.05)
	axins.set_xticks(ticks=[])
	axins.set_yticks(ticks=[])
	mark_inset(axes[1], axins, loc1=2, loc2=1, fc='None', ec='k', zorder=10)
	axes[1].set_ylim(-1.0,0.1)
	axes[1].set_ylabel('Spectra', fontsize=26)
	axes[1].tick_params(axis='x', labelsize=26)
	axes[1].tick_params(axis='y', labelsize=26)

	axes[2].plot(lam, spec, c=CBcols[0], label='Spectra')
	axes[2].plot(lam[idx_out], spec[idx_out], c=CBcols[1], label='Spectra without lines')
	#axes[2].set_ylim(0., 1.1)
	axes[2].axhline(1.0, ls='dashed', c='grey')
	axins = axes[2].inset_axes([0.05, 0.1, 0.5, 0.5])
	axins.plot(lam, spec)
	axins.plot(lam[idx_out], spec[idx_out], 'o-', ms=4.)
	axins.axhline(1.0, ls='dashed', c='grey')
	axins.set_xlim(np.mean(lam)+4e-3, np.mean(lam)+4e-3+2e-3)
	#axins.set_ylim(0.87, 1.02)
	axins.set_ylim(0.87*np.max(spec), 1.02*np.max(spec))
	axins.set_xticks(ticks=[])
	axins.set_yticks(ticks=[])
	mark_inset(axes[2], axins, loc1=2, loc2=1, fc='None', ec='k', zorder=10)
	axes[2].set_xlabel('Wavelength ('+r'$\mu$'+'m)', fontsize=26)
	axes[2].set_ylabel('Spectra', fontsize=26)
	axes[2].tick_params(axis='x', labelsize=26)
	axes[2].tick_params(axis='y', labelsize=26)

	plt.subplots_adjust(hspace=0)

	if savefig:		
		plt.savefig(dirname+filename+'.pdf')
		plt.savefig(dirname+filename+'.jpeg', dpi=300)

	#plt.show()
	plt.close()

def plot_normalisation(lam, spec, norm_mask, filt, savefig=False, dirname=save_directory, filename='norm'):
	fig, axes = plt.subplots(2, figsize=(18,14), sharex=True)

	axes[0].plot(lam, spec, c=CBcols[0])
	axes[0].plot(lam[norm_mask], spec[norm_mask], 'o-', ms=3., c=CBcols[1])
	axes[0].plot(lam, filt, c='black')
	axes[0].set_ylim(0.95*np.min(spec), 1.05*np.max(spec))
	axins = axes[0].inset_axes([0.05, 0.1, 0.5, 0.5])
	axins.plot(lam, spec, c=CBcols[0])
	axins.plot(lam[norm_mask], spec[norm_mask], 'o-', ms=3., c=CBcols[1])
	axins.plot(lam, filt, c='black')
	axins.set_xlim(np.mean(lam)+3e-3, np.mean(lam)+3e-3+3e-3)
	axins.set_ylim(0.90*np.median(spec[norm_mask]), 1.05*np.median(spec[norm_mask]))
	axins.set_xticks(ticks=[])
	axins.set_yticks(ticks=[])
	mark_inset(axes[0], axins, loc1=2, loc2=1, fc='None', ec='k', zorder=10)
	axes[0].set_ylabel('Spectra', fontsize=26)
	axes[0].tick_params(axis='x', labelsize=26)
	axes[0].tick_params(axis='y', labelsize=26)

	axes[1].plot(lam, spec/filt, c=CBcols[0], label='Spectra')
	axes[1].axhline(1.0, ls='dashed', c='grey')
	axins = axes[1].inset_axes([0.05, 0.1, 0.5, 0.5])
	axins.plot(lam, spec/filt)
	axins.axhline(1.0, ls='dashed', c='grey')
	axins.set_xlim(np.mean(lam)+3e-3, np.mean(lam)+3e-3+3e-3)
	axins.set_ylim(0.9,1.05)
	axins.set_xticks(ticks=[])
	axins.set_yticks(ticks=[])
	mark_inset(axes[1], axins, loc1=2, loc2=1, fc='None', ec='k', zorder=10)
	axes[1].set_xlabel('Wavelength ('+r'$\mu$'+'m)', fontsize=26)
	axes[1].set_ylabel('Spectra', fontsize=26)
	axes[1].tick_params(axis='x', labelsize=26)
	axes[1].tick_params(axis='y', labelsize=26)

	plt.subplots_adjust(hspace=0)

	if savefig:		
		plt.savefig(dirname+filename+'.pdf')
		plt.savefig(dirname+filename+'.jpeg', dpi=300)

	#plt.show()
	plt.close()

def plot_atm_prof(profiles, profiles_names, plot_leg=False, savefig=False, dirname=save_directory, filename='atm_profile'):
	linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
	markerstyles = ['s', '*']
	fig = plt.figure(figsize=(12,10))
	#fig = plt.figure()
	host = plt.subplot2grid((1,1), (0,0), rowspan=1, colspan=1)

	part1 = host.twiny()
	part2 = host.twiny()
	part3 = host.twiny()

	part2.spines.top.set_position(('axes', 1.10))
	part3.spines.top.set_position(('axes', 1.20))

	handles = []
	for i, profile in enumerate(profiles):
		height = profile[0]	
		pressure = profile[1]
		temperature = profile[2]
		X_CO2 = profile[3]
		X_H2O = profile[4]
		p1, = host.plot(pressure, height, ls=linestyles[i], c=CBcols[0], label='Pressure', zorder=i)
		p2, = part1.plot(temperature, height, ls=linestyles[i], c=CBcols[1], label='Temperature', zorder=i)
		p3, = part2.plot(X_CO2, height, ls=linestyles[i], c=CBcols[2], label='X_CO2', zorder=i)
		p4, = part3.plot(X_H2O, height, ls=linestyles[i], c=CBcols[3], label='X_H2O', zorder=i)

	for i in range(len(profiles_names)):
		handles.append(Line2D([0], [0], ls=linestyles[i], lw=3., c='grey'))

	host.set_ylabel('Height (m)', fontsize=16)
	host.set_xlabel('Pressure (hPa)', color=p1.get_color(), fontsize=16)
	part1.set_xlabel('Temperature (K)', color=p2.get_color(), fontsize=16)
	part2.set_xlabel('X_CO2 (ppmv)', color=p3.get_color(), fontsize=16)
	part3.set_xlabel('X_H2O (ppmv)', color=p4.get_color(), fontsize=16)

	if plot_leg:
		plt.legend(loc='best', handles=handles, labels=profiles_names, fontsize=16).set_zorder(20)

	#host.set_ylim(2000,7000)
	host.tick_params(axis='y', labelsize=16)
	host.tick_params(axis='x', labelsize=16)
	part1.tick_params(axis='x', labelsize=16)
	part2.tick_params(axis='x', labelsize=16)
	part3.tick_params(axis='x', labelsize=16)
	host.invert_xaxis()
	part1.invert_xaxis()
	part2.invert_xaxis()
	part3.invert_xaxis()

	if savefig:
		fig.savefig(dirname+filename+'.pdf', bbox_inches='tight')
		fig.savefig(dirname+filename+'.jpeg', dpi=350, bbox_inches='tight')

	#plt.show()
	plt.close()

def plot_all_atm_profs(profiles, profiles_names, plot_leg=False, savefig=False, dirname=save_directory, filename='all_atm_profs'):
	linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
	markerstyles = ['s', '*']
	cols = [CBcols[0],CBcols[1],CBcols[2],CBcols[3],CBcols[5],CBcols[4]]

	fig = plt.figure(figsize=(12,10))
	host = plt.subplot2grid((1,1), (0,0), rowspan=1, colspan=1)

	part1 = host.twiny()
	part2 = host.twiny()
	part3 = host.twiny()

	part2.spines.top.set_position(('axes', 1.12))
	part3.spines.top.set_position(('axes', 1.24))

	handles = []
	for i, profile_list in enumerate(profiles):
		for j in range(len(profile_list)):
			height = 1e-3*profile_list[0]	
			pressure = profile_list[1]
			temperature = profile_list[2]
			X_CO2 = profile_list[3]['CO2']
			X_H2O = profile_list[3]['H2O']
			p1, = host.plot(pressure, height, ls=linestyles[0], c=CBcols[0], label='Pressure', alpha=0.15)
			p2, = part1.plot(temperature, height, ls=linestyles[0], c=CBcols[1], label='Temperature', alpha=0.15)
			p3, = part2.plot(X_CO2, height, ls=linestyles[0], c=CBcols[2], label='X_CO2', alpha=0.15)
			p4, = part3.plot(X_H2O, height, ls=linestyles[0], c=CBcols[3], label='X_H2O', alpha=0.15)
		handles.append(Line2D([0], [0], ls=linestyles[0], lw=3., c='gray'))

	host.set_ylabel('Height (km)', fontsize=28)
	host.set_xlabel('Pressure (hPa)', color=p1.get_color(), fontsize=28)
	part1.set_xlabel('Temperature (K)', color=p2.get_color(), fontsize=28)
	part2.set_xlabel('X_CO2 (ppmv)', color=p3.get_color(), fontsize=28)
	part3.set_xlabel('X_H2O (ppmv)', color=p4.get_color(), fontsize=28)

	#host.set_ylim(2000,10000)
	host.tick_params(axis='y', labelsize=28)
	host.tick_params(axis='x', labelsize=28)
	part1.tick_params(axis='x', labelsize=28)
	part2.tick_params(axis='x', labelsize=28)
	part3.tick_params(axis='x', labelsize=28)
	host.invert_xaxis()
	part1.invert_xaxis()
	part2.invert_xaxis()
	part3.invert_xaxis()

	if plot_leg:
		plt.legend(loc='best', handles=handles, labels=profiles_names, fontsize=28).set_zorder(20)

	if savefig:
		fig.savefig(dirname+filename+'.pdf', bbox_inches='tight')
		fig.savefig(dirname+filename+'.jpeg', dpi=350, bbox_inches='tight')

	#plt.show()
	plt.close()

def plot_molecs_atm_prof(profiles, profiles_names, plot_leg=False, savefig=False, dirname=save_directory, filename='molecs_atm_profile'):
	linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
	markerstyles = ['s', '*']
	fig = plt.figure(figsize=(12,10))
	host = plt.subplot2grid((1,1), (0,0), rowspan=1, colspan=1)

	part1 = host.twiny()
	part2 = host.twiny()
	part3 = host.twiny()
	part4 = host.twiny()

	part2.spines.top.set_position(('axes', 1.10))
	part3.spines.top.set_position(('axes', 1.20))
	part4.spines.top.set_position(('axes', 1.30))

	handles = []
	for i, profile in enumerate(profiles):
		height = profile[0]	
		X_CO2 = profile[1]
		X_CH4 = profile[2]
		X_H2O = profile[3]
		X_O2 = profile[4]
		X_N2 = profile[5]
		p1, = host.plot(X_CO2, height, ls=linestyles[i], c=CBcols[0], label='X_CO2', zorder=i)
		p2, = part1.plot(X_CH4, height, ls=linestyles[i], c=CBcols[1], label='X_CH4', zorder=i)
		p3, = part2.plot(X_H2O, height, ls=linestyles[i], c=CBcols[2], label='X_H2O', zorder=i)
		p4, = part3.plot(X_O2, height, ls=linestyles[i], c=CBcols[3], label='X_O2', zorder=i)
		p5, = part4.plot(X_N2, height, ls=linestyles[i], c=CBcols[4], label='X_N2', zorder=i)

	for i in range(len(profiles_names)):
		handles.append(Line2D([0], [0], ls=linestyles[i], lw=3., c='grey'))

	host.set_ylabel('Height (m)', fontsize=16)
	host.set_xlabel('X_CO2 (ppmv)', color=p1.get_color(), fontsize=16)
	part1.set_xlabel('X_CH4 (ppmv)', color=p2.get_color(), fontsize=16)
	part2.set_xlabel('X_H2O (ppmv)', color=p3.get_color(), fontsize=16)
	part3.set_xlabel('X_O2 (ppmv)', color=p4.get_color(), fontsize=16)
	part4.set_xlabel('X_N2 (ppmv)', color=p5.get_color(), fontsize=16)

	if plot_leg:
		plt.legend(loc='best', handles=handles, labels=profiles_names, fontsize=16).set_zorder(20)

	#host.set_ylim(2000,7000)
	host.tick_params(axis='y', labelsize=16)
	host.tick_params(axis='x', labelsize=16)
	part1.tick_params(axis='x', labelsize=16)
	part2.tick_params(axis='x', labelsize=16)
	part3.tick_params(axis='x', labelsize=16)
	part4.tick_params(axis='x', labelsize=16)
	host.invert_xaxis()
	part1.invert_xaxis()
	part2.invert_xaxis()
	part3.invert_xaxis()
	part4.invert_xaxis()

	if savefig:
		fig.savefig(dirname+filename+'.pdf', bbox_inches='tight')
		fig.savefig(dirname+filename+'.jpeg', dpi=350, bbox_inches='tight')

	#plt.show()
	plt.close()

def plot_tel_lines(lam, specs, labels=['Obs. spectra', 'Model spectra'], zoom=True, title=None, xlim=None, ticks=None, cols=CBcols, savefig=False, dirname=save_directory, filename='model_spectra'):
	#fig, axes = plt.subplots(2, figsize=(10, 8), sharex=True)
	fig = plt.figure(figsize=(10,8))
	ax1 = plt.subplot2grid((3,1), (0,0), rowspan=2, colspan=1)
	ax2 = plt.subplot2grid((3,1), (2,0), rowspan=1, colspan=1)#, sharex=ax1)

	if title:
		plt.title(title, fontsize=16, fontweight='bold')

	spec_obs = specs[0]
	ax1.plot(lam, spec_obs, lw=2., label=labels[0], c=cols[0], zorder=5)
	ax1.axhline(1, ls='dashed', c='black', zorder=6)

	if zoom:
		axins = ax1.inset_axes([0.05, 0.1, 0.5, 0.5])
		axins.plot(lam, spec_obs, c=cols[0])
		axins.axhline(1, ls='dashed', c='black')
		axins.set_xlim(xlim[0], xlim[1])
		axins.set_ylim(0.95*np.min(spec_obs[(lam > xlim[0]) & (lam < xlim[1])]), 1.03)
		axins.set_xticks(ticks=[])
		axins.set_yticks(ticks=[])
		mark_inset(ax1, axins, loc1=2, loc2=1, fc='None', ec='k', zorder=10)
	
	#ax1.tick_params(axis='x', labelsize=16)
	ax1.tick_params(axis='y', labelsize=16)
	ax1.set_xticks([])
	ax1.set_ylabel('Transmission', fontsize=16)
	if not zoom and xlim:
		ax1.set_xlim(xlim[0], xlim[1])
		ax2.set_xlim(xlim[0], xlim[1])
	if ticks:
		ax2.set_xticks(ticks, labels=ticks)

	ax2.axhline(0, ls='dashed', lw=2., c='black')
	#ax2.set_ylim(-0.06, 0.06)
	#ax1.set_ylim(0.90, 1.005)
	#ax2.set_ylim(-0.03, 0.004)
	ax2.tick_params(axis='x', labelsize=16)
	ax2.tick_params(axis='y', labelsize=16)
	ax2.set_ylabel('(O-C)', fontsize=16)
	ax2.set_xlabel(r'$Wavelength\ (\mu m)$', fontsize=16)
	
	for i in range(1,len(specs)):
		ax1.plot(lam, specs[i], lw=2., label=labels[i], c=cols[i], zorder=5)
		if zoom:
			axins.plot(lam, specs[i], c=cols[i])
		ax2.plot(lam, spec_obs - specs[i], 'o', ms=2., c=cols[i])

	if zoom:
		ax1.legend(loc='lower right', fontsize=12).set_zorder(15)
	else:
		ax1.legend(loc='best', fontsize=12).set_zorder(15)
	plt.subplots_adjust(hspace=0)

	if savefig:		
		plt.savefig(dirname+filename+'.pdf')
		plt.savefig(dirname+filename+'.jpeg', dpi=300)

	#plt.show()
	plt.close()

def make_autocorr_plot(check_each_n_steps, savefig=True, dirname=save_directory, filename='autocorr'):
	chaincor, acortime = np.loadtxt(dirname+filename+'.txt', unpack=True)

	plt.figure(figsize=(12,8))
	plt.plot(chaincor, acortime)
	plt.plot(np.arange(check_each_n_steps,np.max(chaincor)), np.arange(check_each_n_steps,np.max(chaincor))/100, ls='dashed', c='grey')
	plt.xlim(0,np.max(chaincor))
	plt.ylim(0,1.1*np.max(acortime))
	plt.xlabel('Chain number', fontsize=16)
	plt.ylabel('Auto-correlation time', fontsize=16)
	
	if savefig:		
		plt.savefig(dirname+filename+'.pdf')
		plt.savefig(dirname+filename+'.jpeg', dpi=300)

	#plt.show()
	plt.close()

def make_chain_plot(samples, ndim, labels, savefig=True, dirname=save_directory, filename='chain'):
	fig, axes = plt.subplots(ndim, figsize=(14, 8), sharex=True)
	for i in range(ndim):
		if ndim == 1:
			ax = axes
		else:
			ax = axes[i]
		ax.plot(samples[:, :, i], "k", alpha=0.3)
		ax.set_xlim(0, len(samples))
		ax.set_ylabel(labels[i], fontsize=16)
		ax.yaxis.set_label_coords(-0.1+np.heaviside((-1)**(i),0)*0.03, 0.5)

	if ndim == 1:
		axes.set_xlabel("Step number", fontsize=16)
	else:
		axes[-1].set_xlabel("Step number", fontsize=16)
	plt.subplots_adjust(hspace=0.6)
	
	if savefig:		
		plt.savefig(dirname+filename+'.pdf')
		plt.savefig(dirname+filename+'.jpeg', dpi=300)

	#plt.show()
	plt.close()

def make_corner_plot(ndim, flat_samples, logpflat, labels, truths, savefig=True, dirname=save_directory, filename='corner'):
	fig = corner.corner(flat_samples, labels=labels, truths=truths)

	# Extract the axes
	axes = np.array(fig.axes).reshape((ndim, ndim))

	mcmc_median = np.median(flat_samples, axis=0)
	mcmc_mean = np.mean(flat_samples, axis=0)
	mcmc_std = np.std(flat_samples, axis=0)
	highest_prob = flat_samples[np.argmax(logpflat),:]

	# Loop over the diagonal
	for i in range(ndim):
		ax = axes[i, i]
		ax.set_title(f"{mcmc_median[i]:.3f} +- {mcmc_std[i]:.3f} ")
		ax.axvline(mcmc_median[i], color=CBcols[2])
		ax.axvline(highest_prob[i], color=CBcols[7])

	# Loop over the histograms
	for yi in range(ndim):
		for xi in range(yi):
			ax = axes[yi, xi]
			ax.axvline(mcmc_median[xi], color=CBcols[2])
			ax.axhline(mcmc_median[yi], color=CBcols[2])
			ax.plot(mcmc_median[xi], mcmc_median[yi], color=CBcols[2], marker='s')
			ax.axvline(highest_prob[xi], color=CBcols[7])
			ax.axhline(highest_prob[yi], color=CBcols[7])
			ax.plot(highest_prob[xi], highest_prob[yi], color=CBcols[7], marker='s')

	if savefig:		
		plt.savefig(dirname+filename+'.pdf')
		plt.savefig(dirname+filename+'.jpeg', dpi=300)

	#plt.show()
	plt.close()

def plot_par_evol(times, u_times, Xs, u_Xs, legend_labels, xlim, ylim, ticks, tick_labels, xlabel, ylabel, markerstyles, markersizes, linestyles, zorders, cols=CBcols, xscale='linear', yscale='linear', legend=False, leg_loc='best', nleg_cols=1, savefig=False, dirname=save_directory, filename='abundance_evolution'):

	plt.figure()
	host = plt.subplot2grid((1,1), (0,0), rowspan=1, colspan=1)

	lgs = []
	for i in range(len(Xs)):
		p1 = host.errorbar(times[i], Xs[i], yerr=u_Xs[i], xerr=u_times[i], fmt=markerstyles[i], ms=markersizes[i], ls=linestyles[i], c=cols[i], label=legend_labels[i], zorder=zorders[i])
		if not legend_labels[i] == None:
			lgs.append(p1)

	if xscale == 'log':
		plt.xscale('log')
	else:
		plt.xlim(xlim[0], xlim[1])
	if yscale == 'log':
		plt.yscale('log')
	else:
		plt.ylim(ylim[0], ylim[1])
	
	plt.xticks(ticks, labels=tick_labels, rotation=-60)
	plt.tick_params(axis='x', labelsize=14)
	plt.tick_params(axis='y', labelsize=14)

	plt.xlabel(xlabel, fontsize=14)
	plt.ylabel(ylabel, fontsize=14)

	if legend:
		host.legend(loc=leg_loc, handles=lgs, fontsize=10, framealpha=0.5, ncols=nleg_cols).set_zorder(20)

	if savefig:		
		plt.savefig(dirname+filename+'.pdf', bbox_inches='tight')
		plt.savefig(dirname+filename+'.jpeg', dpi=350, bbox_inches='tight')

	plt.show()
	plt.close()

def plot_par_evol_with_res(times, u_times, Xs, u_Xs, truth_vals, legend_labels, xlim, ylim, ticks, tick_labels, xlabel, ylabel, markerstyles, markersizes, linestyles, zorders, cols=CBcols, xscale='linear', yscale='linear', legend=False, leg_loc='best', nleg_cols=1, savefig=False, dirname=save_directory, filename='abundance_evolution'):

	fig = plt.figure(figsize=(20,8))
	ax1 = plt.subplot2grid((3,1), (0,0), rowspan=2, colspan=1)
	ax2 = plt.subplot2grid((3,1), (2,0), rowspan=1, colspan=1)

	lgs = []
	for i in range(len(Xs)):
		ax1.plot(times[i], truth_vals, ls='dashed', c='black', zorder=zorders[i])
		p1 = ax1.errorbar(times[i], Xs[i], yerr=u_Xs[i], xerr=u_times[i], fmt=markerstyles[i], ms=markersizes[i], ls=linestyles[i], c=cols[i], label=legend_labels[i], zorder=zorders[i])
		ax2.errorbar(times[i], Xs[i]-truth_vals, yerr=u_Xs[i], xerr=u_times[i], fmt=markerstyles[i], ms=markersizes[i], ls=linestyles[i], c=cols[i], label=legend_labels[i], zorder=zorders[i])
		if not legend_labels[i] == None:
			lgs.append(p1)

	ax2.axhline(0, ls='dashed', lw=2., c='black', zorder=7)
	
	ax1.tick_params(axis='y', labelsize=28)
	ax1.set_xticks([])

	if xscale == 'log':
		ax1.set_xscale('log')
		ax2.set_xscale('log')
	else:
		plt.xlim(xlim[0], xlim[1])
	if yscale == 'log':
		ax1.set_yscale('log')
	else:
		plt.ylim(ylim[0], ylim[1])
	
	ax1.set_xticks([])
	ax1.tick_params(axis='y', labelsize=28)
	ax1.set_yticks(ticks, labels=tick_labels)
	ax1.set_ylabel(ylabel, fontsize=28)
	
	ax2.tick_params(axis='x', labelsize=28)
	ax2.tick_params(axis='y', labelsize=28)
	ax2.set_xlabel(xlabel, fontsize=28)
	ax2.set_xticks(ticks, labels=tick_labels)
	#ax2.set_ylim(-3,3)
	ax2.set_ylabel('(O-C)', fontsize=28)
	ax2.set_xlabel(xlabel, fontsize=28)

	if legend:
		ax1.legend(loc=leg_loc, handles=lgs, fontsize=25, framealpha=0.5, ncols=nleg_cols).set_zorder(20)

	plt.subplots_adjust(hspace=0)

	if savefig:		
		plt.savefig(dirname+filename+'.pdf', bbox_inches='tight')
		plt.savefig(dirname+filename+'.jpeg', dpi=350, bbox_inches='tight')

	plt.show()
	plt.close()

def plot_BIG_par_evol(times, u_times, Xs, u_Xs, legend_labels, xlim, ylim, ticks, tick_labels, xlabel, ylabel, markerstyles, markersizes, linestyles, zorders, cols=CBcols, xscale='linear', yscale='linear', legend=False, leg_loc='best', nleg_cols=1, savefig=False, dirname=save_directory, filename='abundance_evolution'):

	n_subplots = len(Xs)
	fig = plt.figure()#figsize=(20,5*n_subplots))
	axes = []
	for n in range(n_subplots):
		axes.append(plt.subplot2grid((n_subplots*n_subplots,1), (n_subplots*n,0), rowspan=n_subplots, colspan=1))

	for i, ax in enumerate(axes):
		lgs = []
		for j in range(len(Xs[i])):
			p1 = ax.errorbar(times[i][j], Xs[i][j], yerr=u_Xs[i][j], xerr=u_times[i][j], fmt=markerstyles[i][j], ms=markersizes[i][j], ls=linestyles[i][j], c=cols[i][j], label=legend_labels[j], zorder=zorders[i][j])
			if not legend_labels[j] == None:
				lgs.append(p1)

		if xscale == 'log':
			ax.set_xscale('log')
		else:
			ax.set_xlim(xlim[0], xlim[1])
		if yscale == 'log':
			ax.set_yscale('log')
		else:
			ax.set_ylim(ylim[i][0], ylim[i][1])
		
		ax.set_xticks([])
		ax.tick_params(axis='y', labelsize=14)

		ax.set_ylabel(ylabel[i], fontsize=14)

		if legend[i]:
			axes[0].legend(loc=leg_loc[i], handles=lgs, fontsize=12, framealpha=0.5, ncols=nleg_cols).set_zorder(20)

	ax.set_xticks(ticks, labels=tick_labels)
	ax.tick_params(axis='x', labelsize=14)
	ax.set_xlabel(xlabel, fontsize=14)

	fig.align_ylabels()
	plt.subplots_adjust(hspace=0)

	if savefig:		
		plt.savefig(dirname+filename+'.pdf', bbox_inches='tight')
		plt.savefig(dirname+filename+'.jpeg', dpi=350, bbox_inches='tight')

	plt.show()
	plt.close()

def plot_all_spec_with_mods(data, labels=['Obs. spectra', 'Model spectra', 'New model spectra'], xlim=None, xx=None, cols=CBcols, savefig=False, dirname=save_directory, filename='all_specs_and_mods'):
	#fig, axes = plt.subplots(2, figsize=(10, 8), sharex=True)
	fig = plt.figure()#figsize=(12,8))
	ax1 = plt.subplot2grid((3,1), (0,0), rowspan=2, colspan=1)
	ax2 = plt.subplot2grid((3,1), (2,0), rowspan=1, colspan=1)#, sharex=ax1)

	ax1.axhline(1, ls='dashed', c='black', zorder=7)
	ax2.axhline(0, ls='dashed', lw=2., c='black', zorder=7)
	if xx:
		ax1.axvline(xx[0], lw=2, c='k', zorder=7)
		ax1.axvline(xx[1], lw=2, c='k', zorder=7)

	RMSEs = np.zeros(len(data))
	RMSEs_new = np.zeros(len(data))
	for i in range(len(data)):
		lam = data[i][0]
		spec_obs = data[i][1]
		spec_mod = data[i][2]
		spec_mod_new = data[i][3]
		ax1.plot(lam, spec_obs, lw=2., c=cols[0], alpha=0.3, zorder=5)
		ax1.plot(lam, spec_mod, lw=2., c=cols[1], alpha=0.3, zorder=6)
		ax1.plot(lam, spec_mod_new, lw=2., c=cols[2], alpha=0.3, zorder=6)
		ax2.plot(lam, spec_obs - spec_mod, 'o', ms=2., c=cols[1], alpha=0.3)
		ax2.plot(lam, spec_obs - spec_mod_new, 'o', ms=2., c=cols[2], alpha=0.3)
		RMSEs[i] = np.sqrt(np.mean((spec_obs - spec_mod)**2))
		RMSEs_new[i] = np.sqrt(np.mean((spec_obs - spec_mod_new)**2))

	handles = []
	for i in range(len(labels)):
		handles.append(Line2D([0], [0], ls="solid", c=cols[i], label=labels[i]))
	
	ax1.tick_params(axis='y', labelsize=14)
	ax1.set_xticks([])
	ax1.set_ylabel('Transmission', fontsize=14)
	if xlim:
		ax1.set_xlim(xlim[0], xlim[1])
		ax2.set_xlim(xlim[0], xlim[1])
	ax1.set_ylim(-0.01,1.2)
	#ax2.set_ylim(-0.15,0.15)

	RMSEs_mean = np.mean(RMSEs)
	ax2.text(0.05, 0.8, r'$\overline{RMSE}=$'+f'{RMSEs_mean:.2f}', fontsize=12, bbox=dict(facecolor='none', edgecolor='grey', boxstyle='round'), transform = ax2.transAxes)
	RMSEs_new_mean = np.mean(RMSEs_new)
	ax2.text(0.45, 0.8, r'$\overline{RMSE}=$'+f'{RMSEs_new_mean:.2f}', fontsize=12, bbox=dict(facecolor='none', edgecolor='grey', boxstyle='round'), transform = ax2.transAxes)

	ax2.tick_params(axis='x', labelsize=14)
	ax2.tick_params(axis='y', labelsize=14)
	ax2.set_ylabel('Residuals', fontsize=14)
	ax2.set_xlabel('Wavelength ('+r'$\mu$'+'m)', fontsize=14)

	ax1.legend(loc='lower right', handles=handles, fontsize=12).set_zorder(15)

	fig.align_ylabels()
	plt.tight_layout()
	plt.subplots_adjust(hspace=0)

	if savefig:		
		#plt.savefig(dirname+filename+'.pdf')
		plt.savefig(dirname+filename+'.jpeg', dpi=300)

	plt.show()
	plt.close()

def plot_cubes(lam, flux, orders, orders_to_plot, title, vmin=None, vmax=None, savefig=False, dirname=save_directory, filename='flux_cubes'):
	norders, nobs, nwaves = np.shape(lam)
	norders = len(orders_to_plot)
	x_arr = np.arange(nwaves+1,step=int(nwaves/6)-1, dtype=int)
	fig, axes = plt.subplots(norders, figsize=(10,int(norders*2)))
	mid = (fig.subplotpars.right + fig.subplotpars.left)/2
	fig.suptitle(title, fontsize=22, x=mid)
	for i in range(norders):
		ord_idx = orders.index(orders_to_plot[i])
		axes[i].set_title(f'Order {orders_to_plot[i]}')
		if vmin and vmax:
			axes[i].imshow(flux[ord_idx,:,:], origin='lower', aspect='auto', interpolation='None', vmin=vmin, vmax=vmax)
		else:
			axes[i].imshow(flux[ord_idx,:,:], origin='lower', aspect='auto', interpolation='None')

		axes[i].set(xticks=x_arr, xticklabels=[f'{lam[ord_idx,0,x]*1e6:.3f}' for x in x_arr])
	fig.supxlabel('Wavelength ('+r'$\mu$'+'m)')
	fig.supylabel('Obs. nr.')

	plt.tight_layout()
	#if len(orders_to_plot) > 7:
	#	plt.subplots_adjust(bottom=0.025, top=0.965)
	plt.subplots_adjust(bottom=0.035, top=0.95)

	if savefig:		
		plt.savefig(dirname+filename+'.pdf')
		plt.savefig(dirname+filename+'.jpeg', dpi=300)

	#plt.show()
	plt.close()

def plot_cube_steps(lam, cube_steps, vmins, vmaxs, y_phase=[], savefig=False, dirname=save_directory, filename='cube_steps'):
	alphabet = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)', '(m)', '(n)', '(o)', '(p)']
	nsteps = len(cube_steps)
	nwaves = len(lam)
	x_arr = np.arange(nwaves+1,step=int(nwaves/6)-1, dtype=int)
	fig, axes = plt.subplots(nsteps)
	#fig, axes = plt.subplots(nsteps, figsize=(10,int(nsteps*2)))
	mid = (fig.subplotpars.right + fig.subplotpars.left)/2
	for i in range(nsteps):
		im = axes[i].imshow(cube_steps[i], origin='lower', aspect='auto', interpolation='None', vmin=vmins[i], vmax=vmaxs[i])
		axes[i].set(xticks=[])
		fig.colorbar(im, aspect=5, pad=0.01)

		## Checking if the left y-axis is observation number or phase and adding tick labels accordingly
		if len(y_phase) > 0:
			y_arr = np.arange(len(y_phase),step=int(len(y_phase)/3), dtype=int)
			axes[i].set(yticks=y_arr, yticklabels=[f'{y_phase[y]:.3f}' for y in y_arr])
		else:
			nobs = len(cube_steps[i])
			y_arr = np.arange(nobs,step=int(nobs/3), dtype=int)
			axes[i].set(yticks=y_arr, yticklabels=[f'{y+1:.0f}' for y in y_arr])
		axes[i].text(0.02, 0.70, alphabet[i], transform=axes[i].transAxes)

	axes[i].set(xticks=x_arr, xticklabels=[f'{lam[x]*1e6:.3f}' for x in x_arr])
	fig.supxlabel('Wavelength ('+r'$\mu$'+'m)')
	if len(y_phase) > 0:
		fig.supylabel('Orbital phase')
	else:
		fig.supylabel('Obs. nr.')

	plt.tight_layout()
	plt.subplots_adjust(bottom=0.12, left=0.10, hspace=0)

	if savefig:		
		plt.savefig(dirname+filename+'.pdf')
		plt.savefig(dirname+filename+'.jpeg', dpi=300)

	#plt.show()
	plt.close()

def plot_squares(square, extent, xlabel, ylabel, origin='lower', aspect='auto', interpolation='None', savefig=False, dirname=save_directory, filename='flux_squares'):
	fig = plt.figure()
	im = plt.imshow(square, origin=origin, aspect=aspect, interpolation=interpolation, extent=extent)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.xlim(extent[0],extent[1])
	plt.ylim(extent[2],extent[3])
	
	if savefig:		
		plt.savefig(dirname+filename+'.pdf')
		plt.savefig(dirname+filename+'.jpeg', dpi=300)

	#plt.show()
	plt.close()

def plot_SNR_square(SNR, extent, xlabel, ylabel, vRest, kpVec, l1=[], l2=[], SNR_max=None, vRest_max_idx=None, kpVec_max_idx=None, SNR_inj_sig_pos=None, vRest_inj_sig_pos_idx=None, kpVec_inj_sig_pos_idx=None, origin='lower', aspect='auto', interpolation='None', savefig=False, dirname=save_directory, filename='SNR_square'):
	fig = plt.figure()
	im = plt.imshow(SNR, origin=origin, aspect=aspect, interpolation=interpolation, extent=extent)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	if len(l1) > 0:
		plt.plot(vRest, l1, ls='dashed', c='black')
	if len(l2) > 0:
		plt.plot(vRest, l2, ls='dashed', c='black')
	plt.xlim(extent[0],extent[1])
	plt.ylim(extent[2],extent[3])

	if SNR_max:
		## Plotting the maximum at the exact location of the injected signal
		plabel = f'SNR={SNR_inj_sig_pos:.2f}\nVrest={vRest[vRest_inj_sig_pos_idx]:.2f}\ km/s\nKp={kpVec[kpVec_inj_sig_pos_idx]:.2f}\ km/s'
		plt.plot(vRest[vRest_inj_sig_pos_idx], kpVec[kpVec_inj_sig_pos_idx], 's', c='red', label=plabel)
		## Plotting the maximum within the region delimited by l1 and l2
		plabel = f'SNR={SNR_max:.2f}\nVrest={vRest[vRest_max_idx]:.2f}\ km/s\nKp={kpVec[kpVec_max_idx]:.2f}\ km/s'
		plt.plot(vRest[vRest_max_idx], kpVec[kpVec_max_idx], 'o', c='gray', label=plabel)
		## Plotting the global maximum
		plabel = f'SNR={np.max(SNR):.2f}\nVrest={vRest[np.where(SNR == np.max(SNR))[1][0]]:.2f}\ km/s\nKp={kpVec[np.where(SNR == np.max(SNR))[0][0]]:.2f}\ km/s'
		plt.plot(vRest[np.where(SNR == np.max(SNR))[1][0]], kpVec[np.where(SNR == np.max(SNR))[0][0]], '*', c='black', label=plabel)
		plt.legend(loc='upper left')
	cbar = fig.colorbar(im)
	cbar.set_label('SNR')
	
	if savefig:		
		plt.savefig(dirname+filename+'.pdf')
		plt.savefig(dirname+filename+'.jpeg', dpi=300)

	#plt.show()
	plt.close()