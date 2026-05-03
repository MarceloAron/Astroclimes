# =====================================================================================
# Basic packages
# =====================================================================================
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#sys.path.append('/home/marceloaron/MarceloAron/Masters')

# =====================================================================================
# Scripts
# =====================================================================================

rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)

## Color-blind friendly colors: blue, orange, green, pink, brown, purple, gray, red and yellow
CBcols= ['#377EB8', '#FF7F00', '#4DAF4A',
		 '#F781BF', '#A65628', '#984EA3',
		 '#999999', '#E41A1C', '#DEDE00']

## Added extra: blue, light blue, orange, light orange, green, light freen, red, light red,
## purple, light purple, brown, light brown, pink, light pink, gray, light gray, beige-yellow,
## light beige-yellow, turquoise, light turquoise
CBcols2 = [	'#1F77B4', '#AEC7E8', '#FF7F0E', 
			'#FFBB78', '#2CA02C', '#98DF8A', 
			'#D62728', '#FF9896', '#9467BD', 
			'#C5B0D5', '#8C564B', '#C49C94', 
			'#E377C2', '#F7B6D2', '#7F7F7F', 
			'#C7C7C7', '#BCBD22', '#DBDB8D', 
			'#17BECF', '#9EDAE5']

kB = 1.380649*1e-23	# Boltzmann constant, in m^2 kg s^-2 K^-1
R = 8.31446 		# Gas constant, in J mol^-1 K^-1
N_A = 6.022*1e23	# Avogadro constant, in molecules/mole
c = 299792458		# Speed of light, in m/s
DU = 2.69*1e20		# Dobson unit, in molecules/m^2
M_air = 0.02896546 	# Molar mass of dry air, in kg/mol
M_CO2 = 0.0440095	# Molar mass of CO2, in kg/mol
M_CH4 = 0.0160425	# Molar mass of CH4, in kg/mol
M_H2O = 0.0180153	# Molar mass of H2O, in kg/mol

Kp = 110 			# Planet semi-amplitude in km/s, from Brogi et al. (2012), also consistent with Lockwood et al. (2014; 111 km/s), they also report Vsys = -16.4 km/s
Vsys = -16.4 		# Systemic velocity for Tau Bootis in km/s, from Brogi et al. (2012). There is also -16.94 km/s from Gaia DR2 or -16.03 km/s from Nidever et al. (2002) - not sure why I chose this one, maybe because it is the most cited paper?

home_directory_R2D2 = '/media/marceloaron/New Volume/'
home_directory_Mac = '/Users/marceloaronkeniger/'
home_directory_hebei = '/home/physics/phrgmq/'
home_directory = home_directory_R2D2

## Unpacking the results from the MCMC
main_directory = home_directory+'PhD/thesis_work/telluric_removal/TauBoo/2018_03_26_Kp_110/'

## Inj. sign. strength at which MCMC runs started to converge (for Astroclimes and each PCA case)
#mcmc_converge_inj_str = [3,3.25,3.75,3.75,3.75,4,5] # For 2018_03_26_Kp_110, dlt=0.7
#mcmc_converge_inj_str = [3,3.50,3.50,3.50,3.50,4.25,4.75] # For 2018_03_26_Kp_110, dlt=0.2, with post-detrending mask
mcmc_converge_inj_str = [3,3.25,3.50,3.50,3.75,4.25,5,2] # For 2018_03_26_Kp_110, last one is Molecfit
#mcmc_converge_inj_str = [2,3.5,4,4,4,4.75,6,2] 			 # For 2018_03_26_Kp_80, last one is Molecfit
#mcmc_converge_inj_str = [2,4.50,7,7,9,9,11,2]	   			 # For 2018_03_26_Kp_50, last one is Molecfit (Nc_PCA = 10 did not converge)
#mcmc_converge_inj_str = [2,8,3,3,2,2,2] 				# For 2019_03_15

filename_MCMC_results_AC = main_directory+'results_inj_rec_mcmc_ac.txt'
inj_str_MCMC_AC, sf_AC, u_sf_AC, log_sf_AC, u_log_sf_AC, Kp_AC, u_Kp_AC, Vsys_AC, u_Vsys_AC = np.loadtxt(filename_MCMC_results_AC, usecols=(0,1,2,3,4,5,6,7,8), unpack=True)
rr_MCMC_AC = inj_str_MCMC_AC >= mcmc_converge_inj_str[0] # Skipping rows from results that didn't converge
inj_str_MCMC_AC, sf_AC, u_sf_AC, log_sf_AC, u_log_sf_AC, Kp_AC, u_Kp_AC, Vsys_AC, u_Vsys_AC = inj_str_MCMC_AC[rr_MCMC_AC], sf_AC[rr_MCMC_AC], u_sf_AC[rr_MCMC_AC], log_sf_AC[rr_MCMC_AC], u_log_sf_AC[rr_MCMC_AC], Kp_AC[rr_MCMC_AC], u_Kp_AC[rr_MCMC_AC], Vsys_AC[rr_MCMC_AC], u_Vsys_AC[rr_MCMC_AC]

filename_MCMC_results_MOLECFIT = main_directory+'results_inj_rec_mcmc_molecfit.txt'
inj_str_MCMC_MOLECFIT, sf_MOLECFIT, u_sf_MOLECFIT, log_sf_MOLECFIT, u_log_sf_MOLECFIT, Kp_MOLECFIT, u_Kp_MOLECFIT, Vsys_MOLECFIT, u_Vsys_MOLECFIT = np.loadtxt(filename_MCMC_results_MOLECFIT, usecols=(0,1,2,3,4,5,6,7,8), unpack=True)
rr_MCMC_MOLECFIT = inj_str_MCMC_MOLECFIT >= mcmc_converge_inj_str[7] # Skipping rows from results that didn't converge
inj_str_MCMC_MOLECFIT, sf_MOLECFIT, u_sf_MOLECFIT, log_sf_MOLECFIT, u_log_sf_MOLECFIT, Kp_MOLECFIT, u_Kp_MOLECFIT, Vsys_MOLECFIT, u_Vsys_MOLECFIT = inj_str_MCMC_MOLECFIT[rr_MCMC_MOLECFIT], sf_MOLECFIT[rr_MCMC_MOLECFIT], u_sf_MOLECFIT[rr_MCMC_MOLECFIT], log_sf_MOLECFIT[rr_MCMC_MOLECFIT], u_log_sf_MOLECFIT[rr_MCMC_MOLECFIT], Kp_MOLECFIT[rr_MCMC_MOLECFIT], u_Kp_MOLECFIT[rr_MCMC_MOLECFIT], Vsys_MOLECFIT[rr_MCMC_MOLECFIT], u_Vsys_MOLECFIT[rr_MCMC_MOLECFIT]

filename_MCMC_results_PCA_nc_2 = main_directory+'results_inj_rec_mcmc_pca_nc_2.txt'
inj_str_MCMC_PCA_nc_2, nc_PCA_nc_2, sf_PCA_nc_2, u_sf_PCA_nc_2, log_sf_PCA_nc_2, u_log_sf_PCA_nc_2, Kp_PCA_nc_2, u_Kp_PCA_nc_2, Vsys_PCA_nc_2, u_Vsys_PCA_nc_2 = np.loadtxt(filename_MCMC_results_PCA_nc_2, usecols=(0,1,2,3,4,5,6,7,8,9), unpack=True)
rr_MCMC_PCA_nc_2 = inj_str_MCMC_PCA_nc_2 >= mcmc_converge_inj_str[1] # Skipping rows from results that didn't converge
inj_str_MCMC_PCA_nc_2, nc_PCA_nc_2, sf_PCA_nc_2, u_sf_PCA_nc_2, log_sf_PCA_nc_2, u_log_sf_PCA_nc_2, Kp_PCA_nc_2, u_Kp_PCA_nc_2, Vsys_PCA_nc_2, u_Vsys_PCA_nc_2 = inj_str_MCMC_PCA_nc_2[rr_MCMC_PCA_nc_2], nc_PCA_nc_2[rr_MCMC_PCA_nc_2], sf_PCA_nc_2[rr_MCMC_PCA_nc_2], u_sf_PCA_nc_2[rr_MCMC_PCA_nc_2], log_sf_PCA_nc_2[rr_MCMC_PCA_nc_2], u_log_sf_PCA_nc_2[rr_MCMC_PCA_nc_2], Kp_PCA_nc_2[rr_MCMC_PCA_nc_2], u_Kp_PCA_nc_2[rr_MCMC_PCA_nc_2], Vsys_PCA_nc_2[rr_MCMC_PCA_nc_2], u_Vsys_PCA_nc_2[rr_MCMC_PCA_nc_2]

filename_MCMC_results_PCA_nc_3 = main_directory+'results_inj_rec_mcmc_pca_nc_3.txt'
inj_str_MCMC_PCA_nc_3, nc_PCA_nc_3, sf_PCA_nc_3, u_sf_PCA_nc_3, log_sf_PCA_nc_3, u_log_sf_PCA_nc_3, Kp_PCA_nc_3, u_Kp_PCA_nc_3, Vsys_PCA_nc_3, u_Vsys_PCA_nc_3 = np.loadtxt(filename_MCMC_results_PCA_nc_3, usecols=(0,1,2,3,4,5,6,7,8,9), unpack=True)
rr_MCMC_PCA_nc_3 = inj_str_MCMC_PCA_nc_3 >= mcmc_converge_inj_str[2] # Skipping rows from results that didn't converge
inj_str_MCMC_PCA_nc_3, nc_PCA_nc_3, sf_PCA_nc_3, u_sf_PCA_nc_3, log_sf_PCA_nc_3, u_log_sf_PCA_nc_3, Kp_PCA_nc_3, u_Kp_PCA_nc_3, Vsys_PCA_nc_3, u_Vsys_PCA_nc_3 = inj_str_MCMC_PCA_nc_3[rr_MCMC_PCA_nc_3], nc_PCA_nc_3[rr_MCMC_PCA_nc_3], sf_PCA_nc_3[rr_MCMC_PCA_nc_3], u_sf_PCA_nc_3[rr_MCMC_PCA_nc_3], log_sf_PCA_nc_3[rr_MCMC_PCA_nc_3], u_log_sf_PCA_nc_3[rr_MCMC_PCA_nc_3], Kp_PCA_nc_3[rr_MCMC_PCA_nc_3], u_Kp_PCA_nc_3[rr_MCMC_PCA_nc_3], Vsys_PCA_nc_3[rr_MCMC_PCA_nc_3], u_Vsys_PCA_nc_3[rr_MCMC_PCA_nc_3]

filename_MCMC_results_PCA_nc_4 = main_directory+'results_inj_rec_mcmc_pca_nc_4.txt'
inj_str_MCMC_PCA_nc_4, nc_PCA_nc_4, sf_PCA_nc_4, u_sf_PCA_nc_4, log_sf_PCA_nc_4, u_log_sf_PCA_nc_4, Kp_PCA_nc_4, u_Kp_PCA_nc_4, Vsys_PCA_nc_4, u_Vsys_PCA_nc_4 = np.loadtxt(filename_MCMC_results_PCA_nc_4, usecols=(0,1,2,3,4,5,6,7,8,9), unpack=True)
rr_MCMC_PCA_nc_4 = inj_str_MCMC_PCA_nc_4 >= mcmc_converge_inj_str[3] # Skipping rows from results that didn't converge
inj_str_MCMC_PCA_nc_4, nc_PCA_nc_4, sf_PCA_nc_4, u_sf_PCA_nc_4, log_sf_PCA_nc_4, u_log_sf_PCA_nc_4, Kp_PCA_nc_4, u_Kp_PCA_nc_4, Vsys_PCA_nc_4, u_Vsys_PCA_nc_4 = inj_str_MCMC_PCA_nc_4[rr_MCMC_PCA_nc_4], nc_PCA_nc_4[rr_MCMC_PCA_nc_4], sf_PCA_nc_4[rr_MCMC_PCA_nc_4], u_sf_PCA_nc_4[rr_MCMC_PCA_nc_4], log_sf_PCA_nc_4[rr_MCMC_PCA_nc_4], u_log_sf_PCA_nc_4[rr_MCMC_PCA_nc_4], Kp_PCA_nc_4[rr_MCMC_PCA_nc_4], u_Kp_PCA_nc_4[rr_MCMC_PCA_nc_4], Vsys_PCA_nc_4[rr_MCMC_PCA_nc_4], u_Vsys_PCA_nc_4[rr_MCMC_PCA_nc_4]

filename_MCMC_results_PCA_nc_5 = main_directory+'results_inj_rec_mcmc_pca_nc_5.txt'
inj_str_MCMC_PCA_nc_5, nc_PCA_nc_5, sf_PCA_nc_5, u_sf_PCA_nc_5, log_sf_PCA_nc_5, u_log_sf_PCA_nc_5, Kp_PCA_nc_5, u_Kp_PCA_nc_5, Vsys_PCA_nc_5, u_Vsys_PCA_nc_5 = np.loadtxt(filename_MCMC_results_PCA_nc_5, usecols=(0,1,2,3,4,5,6,7,8,9), unpack=True)
rr_MCMC_PCA_nc_5 = inj_str_MCMC_PCA_nc_5 >= mcmc_converge_inj_str[4] # Skipping rows from results that didn't converge
inj_str_MCMC_PCA_nc_5, nc_PCA_nc_5, sf_PCA_nc_5, u_sf_PCA_nc_5, log_sf_PCA_nc_5, u_log_sf_PCA_nc_5, Kp_PCA_nc_5, u_Kp_PCA_nc_5, Vsys_PCA_nc_5, u_Vsys_PCA_nc_5 = inj_str_MCMC_PCA_nc_5[rr_MCMC_PCA_nc_5], nc_PCA_nc_5[rr_MCMC_PCA_nc_5], sf_PCA_nc_5[rr_MCMC_PCA_nc_5], u_sf_PCA_nc_5[rr_MCMC_PCA_nc_5], log_sf_PCA_nc_5[rr_MCMC_PCA_nc_5], u_log_sf_PCA_nc_5[rr_MCMC_PCA_nc_5], Kp_PCA_nc_5[rr_MCMC_PCA_nc_5], u_Kp_PCA_nc_5[rr_MCMC_PCA_nc_5], Vsys_PCA_nc_5[rr_MCMC_PCA_nc_5], u_Vsys_PCA_nc_5[rr_MCMC_PCA_nc_5]

filename_MCMC_results_PCA_nc_6 = main_directory+'results_inj_rec_mcmc_pca_nc_6.txt'
inj_str_MCMC_PCA_nc_6, nc_PCA_nc_6, sf_PCA_nc_6, u_sf_PCA_nc_6, log_sf_PCA_nc_6, u_log_sf_PCA_nc_6, Kp_PCA_nc_6, u_Kp_PCA_nc_6, Vsys_PCA_nc_6, u_Vsys_PCA_nc_6 = np.loadtxt(filename_MCMC_results_PCA_nc_6, usecols=(0,1,2,3,4,5,6,7,8,9), unpack=True)
rr_MCMC_PCA_nc_6 = inj_str_MCMC_PCA_nc_6 >= mcmc_converge_inj_str[5] # Skipping rows from results that didn't converge
inj_str_MCMC_PCA_nc_6, nc_PCA_nc_6, sf_PCA_nc_6, u_sf_PCA_nc_6, log_sf_PCA_nc_6, u_log_sf_PCA_nc_6, Kp_PCA_nc_6, u_Kp_PCA_nc_6, Vsys_PCA_nc_6, u_Vsys_PCA_nc_6 = inj_str_MCMC_PCA_nc_6[rr_MCMC_PCA_nc_6], nc_PCA_nc_6[rr_MCMC_PCA_nc_6], sf_PCA_nc_6[rr_MCMC_PCA_nc_6], u_sf_PCA_nc_6[rr_MCMC_PCA_nc_6], log_sf_PCA_nc_6[rr_MCMC_PCA_nc_6], u_log_sf_PCA_nc_6[rr_MCMC_PCA_nc_6], Kp_PCA_nc_6[rr_MCMC_PCA_nc_6], u_Kp_PCA_nc_6[rr_MCMC_PCA_nc_6], Vsys_PCA_nc_6[rr_MCMC_PCA_nc_6], u_Vsys_PCA_nc_6[rr_MCMC_PCA_nc_6]

filename_MCMC_results_PCA_nc_10 = main_directory+'results_inj_rec_mcmc_pca_nc_10.txt'
inj_str_MCMC_PCA_nc_10, nc_PCA_nc_10, sf_PCA_nc_10, u_sf_PCA_nc_10, log_sf_PCA_nc_10, u_log_sf_PCA_nc_10, Kp_PCA_nc_10, u_Kp_PCA_nc_10, Vsys_PCA_nc_10, u_Vsys_PCA_nc_10 = np.loadtxt(filename_MCMC_results_PCA_nc_10, usecols=(0,1,2,3,4,5,6,7,8,9), unpack=True)
rr_MCMC_PCA_nc_10 = inj_str_MCMC_PCA_nc_10 >= mcmc_converge_inj_str[6] # Skipping rows from results that didn't converge
inj_str_MCMC_PCA_nc_10, nc_PCA_nc_10, sf_PCA_nc_10, u_sf_PCA_nc_10, log_sf_PCA_nc_10, u_log_sf_PCA_nc_10, Kp_PCA_nc_10, u_Kp_PCA_nc_10, Vsys_PCA_nc_10, u_Vsys_PCA_nc_10 = inj_str_MCMC_PCA_nc_10[rr_MCMC_PCA_nc_10], nc_PCA_nc_10[rr_MCMC_PCA_nc_10], sf_PCA_nc_10[rr_MCMC_PCA_nc_10], u_sf_PCA_nc_10[rr_MCMC_PCA_nc_10], log_sf_PCA_nc_10[rr_MCMC_PCA_nc_10], u_log_sf_PCA_nc_10[rr_MCMC_PCA_nc_10], Kp_PCA_nc_10[rr_MCMC_PCA_nc_10], u_Kp_PCA_nc_10[rr_MCMC_PCA_nc_10], Vsys_PCA_nc_10[rr_MCMC_PCA_nc_10], u_Vsys_PCA_nc_10[rr_MCMC_PCA_nc_10]

## Unpacking the SNR values
filename_SNR_results_AC = main_directory+'SNR_logs_ac.txt'
inj_str_SNR_AC, SNR_inj_pos_AC = np.loadtxt(filename_SNR_results_AC, usecols=(0,1), unpack=True)

filename_SNR_results_MOLECFIT = main_directory+'SNR_logs_molecfit.txt'
inj_str_SNR_MOLECFIT, SNR_inj_pos_MOLECFIT = np.loadtxt(filename_SNR_results_MOLECFIT, usecols=(0,1), unpack=True)

filename_SNR_results_PCA_nc_2 = main_directory+'SNR_logs_pca_nc_2.txt'
inj_str_SNR_PCA_nc_2, SNR_inj_pos_PCA_nc_2 = np.loadtxt(filename_SNR_results_PCA_nc_2, usecols=(0,2), unpack=True)

filename_SNR_results_PCA_nc_3 = main_directory+'SNR_logs_pca_nc_3.txt'
inj_str_SNR_PCA_nc_3, SNR_inj_pos_PCA_nc_3 = np.loadtxt(filename_SNR_results_PCA_nc_3, usecols=(0,2), unpack=True)

filename_SNR_results_PCA_nc_4 = main_directory+'SNR_logs_pca_nc_4.txt'
inj_str_SNR_PCA_nc_4, SNR_inj_pos_PCA_nc_4 = np.loadtxt(filename_SNR_results_PCA_nc_4, usecols=(0,2), unpack=True)

filename_SNR_results_PCA_nc_5 = main_directory+'SNR_logs_pca_nc_5.txt'
inj_str_SNR_PCA_nc_5, SNR_inj_pos_PCA_nc_5 = np.loadtxt(filename_SNR_results_PCA_nc_5, usecols=(0,2), unpack=True)

filename_SNR_results_PCA_nc_6 = main_directory+'SNR_logs_pca_nc_6.txt'
inj_str_SNR_PCA_nc_6, SNR_inj_pos_PCA_nc_6 = np.loadtxt(filename_SNR_results_PCA_nc_6, usecols=(0,2), unpack=True)

filename_SNR_results_PCA_nc_10 = main_directory+'SNR_logs_pca_nc_10.txt'
inj_str_SNR_PCA_nc_10, SNR_inj_pos_PCA_nc_10 = np.loadtxt(filename_SNR_results_PCA_nc_10, usecols=(0,2), unpack=True)

#filename_SNR_results_PCA_nc_12 = main_directory+'SNR_logs_pca_nc_12.txt'
#inj_str_SNR_PCA_nc_12, SNR_inj_pos_PCA_nc_12 = np.loadtxt(filename_SNR_results_PCA_nc_12, usecols=(0,2), unpack=True)

#filename_SNR_results_PCA_nc_13 = main_directory+'SNR_logs_pca_nc_13.txt'
#inj_str_SNR_PCA_nc_13, SNR_inj_pos_PCA_nc_13 = np.loadtxt(filename_SNR_results_PCA_nc_13, usecols=(0,2), unpack=True)

#filename_SNR_results_PCA_nc_14 = main_directory+'SNR_logs_pca_nc_14.txt'
#inj_str_SNR_PCA_nc_14, SNR_inj_pos_PCA_nc_14 = np.loadtxt(filename_SNR_results_PCA_nc_14, usecols=(0,2), unpack=True)

#filename_SNR_results_PCA_nc_15 = main_directory+'SNR_logs_pca_nc_15.txt'
#inj_str_SNR_PCA_nc_15, SNR_inj_pos_PCA_nc_15 = np.loadtxt(filename_SNR_results_PCA_nc_15, usecols=(0,2), unpack=True)

#filename_SNR_results_PCA_nc_16 = main_directory+'SNR_logs_pca_nc_16.txt'
#inj_str_SNR_PCA_nc_16, SNR_inj_pos_PCA_nc_16 = np.loadtxt(filename_SNR_results_PCA_nc_16, usecols=(0,2), unpack=True)

## Printing the average signal retention for each case
print(f"{' ':15s}\t\t {'Delta':10s} \t {'Kp':10s} \t {'Vsys':10s}")
print(f"{'Astroclimes:':15s} {np.mean(sf_AC[SNR_inj_pos_AC[rr_MCMC_AC] > 4]):10.2f} \t {np.mean(Kp_AC[SNR_inj_pos_AC[rr_MCMC_AC] > 4]):10.2f} \t {np.mean(Vsys_AC[SNR_inj_pos_AC[rr_MCMC_AC] > 4]):10.2f}")
print(f"{'Molecfit:':15s} {np.mean(sf_MOLECFIT[SNR_inj_pos_MOLECFIT[rr_MCMC_MOLECFIT] > 4]):10.2f} \t {np.mean(Kp_MOLECFIT[SNR_inj_pos_MOLECFIT[rr_MCMC_MOLECFIT] > 4]):10.2f} \t {np.mean(Vsys_MOLECFIT[SNR_inj_pos_MOLECFIT[rr_MCMC_MOLECFIT] > 4]):10.2f}")
print(f"{'PCA N_c = 2:':15s} {np.mean(sf_PCA_nc_2[SNR_inj_pos_PCA_nc_2[rr_MCMC_PCA_nc_2] > 4]):10.2f} \t {np.mean(Kp_PCA_nc_2[SNR_inj_pos_PCA_nc_2[rr_MCMC_PCA_nc_2] > 4]):10.2f} \t {np.mean(Vsys_PCA_nc_2[SNR_inj_pos_PCA_nc_2[rr_MCMC_PCA_nc_2] > 4]):10.2f}")
print(f"{'PCA N_c = 3:':15s} {np.mean(sf_PCA_nc_3[SNR_inj_pos_PCA_nc_3[rr_MCMC_PCA_nc_3] > 4]):10.2f} \t {np.mean(Kp_PCA_nc_3[SNR_inj_pos_PCA_nc_3[rr_MCMC_PCA_nc_3] > 4]):10.2f} \t {np.mean(Vsys_PCA_nc_3[SNR_inj_pos_PCA_nc_3[rr_MCMC_PCA_nc_3] > 4]):10.2f}")
print(f"{'PCA N_c = 4:':15s} {np.mean(sf_PCA_nc_4[SNR_inj_pos_PCA_nc_4[rr_MCMC_PCA_nc_4] > 4]):10.2f} \t {np.mean(Kp_PCA_nc_4[SNR_inj_pos_PCA_nc_4[rr_MCMC_PCA_nc_4] > 4]):10.2f} \t {np.mean(Vsys_PCA_nc_4[SNR_inj_pos_PCA_nc_4[rr_MCMC_PCA_nc_4] > 4]):10.2f}")
print(f"{'PCA N_c = 5:':15s} {np.mean(sf_PCA_nc_5[SNR_inj_pos_PCA_nc_5[rr_MCMC_PCA_nc_5] > 4]):10.2f} \t {np.mean(Kp_PCA_nc_5[SNR_inj_pos_PCA_nc_5[rr_MCMC_PCA_nc_5] > 4]):10.2f} \t {np.mean(Vsys_PCA_nc_5[SNR_inj_pos_PCA_nc_5[rr_MCMC_PCA_nc_5] > 4]):10.2f}")
print(f"{'PCA N_c = 6:':15s} {np.mean(sf_PCA_nc_6[SNR_inj_pos_PCA_nc_6[rr_MCMC_PCA_nc_6] > 4]):10.2f} \t {np.mean(Kp_PCA_nc_6[SNR_inj_pos_PCA_nc_6[rr_MCMC_PCA_nc_6] > 4]):10.2f} \t {np.mean(Vsys_PCA_nc_6[SNR_inj_pos_PCA_nc_6[rr_MCMC_PCA_nc_6] > 4]):10.2f}")
print(f"{'PCA N_c = 10:':15s} {np.mean(sf_PCA_nc_10[SNR_inj_pos_PCA_nc_10[rr_MCMC_PCA_nc_10] > 4]):10.2f} \t {np.mean(Kp_PCA_nc_10[SNR_inj_pos_PCA_nc_10[rr_MCMC_PCA_nc_10] > 4]):10.2f} \t {np.mean(Vsys_PCA_nc_10[SNR_inj_pos_PCA_nc_10[rr_MCMC_PCA_nc_10] > 4]):10.2f}")

'''
## Plotting the SNR versus injected planetary signal strength
fig = plt.figure()
plt.plot(inj_str_SNR_AC, SNR_inj_pos_AC, 'o-', c='black', label='Astroclimes', zorder=10)
plt.plot(inj_str_SNR_MOLECFIT, SNR_inj_pos_MOLECFIT, 'o-', c='lightgray', label='Molecfit', zorder=9.5)
plt.plot(inj_str_SNR_PCA_nc_2, SNR_inj_pos_PCA_nc_2, 'o-', c=CBcols[0], label='PCA '+r'$N_c = 2$', zorder=9)
plt.plot(inj_str_SNR_PCA_nc_3, SNR_inj_pos_PCA_nc_3, 'o-', c=CBcols[1], label='PCA '+r'$N_c = 3$', zorder=8)
plt.plot(inj_str_SNR_PCA_nc_4, SNR_inj_pos_PCA_nc_4, 'o-', c=CBcols[2], label='PCA '+r'$N_c = 4$', zorder=7)
plt.plot(inj_str_SNR_PCA_nc_5, SNR_inj_pos_PCA_nc_5, 'o-', c=CBcols[3], label='PCA '+r'$N_c = 5$', zorder=6)
plt.plot(inj_str_SNR_PCA_nc_6, SNR_inj_pos_PCA_nc_6, 'o-', c=CBcols[7], label='PCA '+r'$N_c = 6$', zorder=5)
plt.plot(inj_str_SNR_PCA_nc_10, SNR_inj_pos_PCA_nc_10, 'o-', c=CBcols[8], label='PCA '+r'$N_c = 10$', zorder=4)
plt.xlabel('Injected planet signal multiplier', fontsize=14)
plt.ylabel('SNR', fontsize=14)
plt.tick_params(axis='x', labelsize=14)
plt.tick_params(axis='y', labelsize=14)
plt.axhline(4, ls='dashed', c='gray')
plt.legend(loc='best', framealpha=0.5)
plt.tight_layout()
plt.savefig(main_directory+'SNR.pdf')
plt.savefig(main_directory+'SNR.jpeg', dpi=300)
plt.show()
plt.close()

## Plotting the MCMC scaling factor versus injected planetary signal strength
plt.figure()
plt.errorbar(inj_str_MCMC_AC, sf_AC, yerr=u_sf_AC, fmt='o-', c='black',label= 'Astroclimes', zorder=10)
plt.errorbar(inj_str_MCMC_MOLECFIT, sf_MOLECFIT, yerr=u_sf_MOLECFIT, fmt='o-', c='lightgray', label='Molecfit', zorder=9.5)
plt.errorbar(inj_str_MCMC_PCA_nc_2, sf_PCA_nc_2, yerr=u_sf_PCA_nc_2, fmt='o-', c=CBcols[0], label='PCA '+r'$N_c = 2$', zorder=9)
plt.errorbar(inj_str_MCMC_PCA_nc_3, sf_PCA_nc_3, yerr=u_sf_PCA_nc_3, fmt='o-', c=CBcols[1], label='PCA '+r'$N_c = 3$', zorder=8)
plt.errorbar(inj_str_MCMC_PCA_nc_4, sf_PCA_nc_4, yerr=u_sf_PCA_nc_4, fmt='o-', c=CBcols[2], label='PCA '+r'$N_c = 4$', zorder=7)
plt.errorbar(inj_str_MCMC_PCA_nc_5, sf_PCA_nc_5, yerr=u_sf_PCA_nc_5, fmt='o-', c=CBcols[3], label='PCA '+r'$N_c = 5$', zorder=6)
plt.errorbar(inj_str_MCMC_PCA_nc_6, sf_PCA_nc_6, yerr=u_sf_PCA_nc_6, fmt='o-', c=CBcols[7], label='PCA '+r'$N_c = 6$', zorder=5)
plt.errorbar(inj_str_MCMC_PCA_nc_10, sf_PCA_nc_10, yerr=u_sf_PCA_nc_10, fmt='o-', c=CBcols[8], label='PCA '+r'$N_c = 10$', zorder=4)
plt.xlabel('Injected planet signal multiplier', fontsize=14)
plt.ylabel('Scaling factor '+r'$\delta$', fontsize=14)
plt.tick_params(axis='x', labelsize=14)
plt.tick_params(axis='y', labelsize=14)
plt.legend(loc=(0.04,0.518), fontsize=10, ncols=3, framealpha=0.5) 	# For 2018/03/26, Kp = 110 km/s (with dlt=0.7 as well)
#plt.legend(loc=(0.04,0.418), fontsize=10, ncols=3, framealpha=0.5)		
#plt.legend(loc='upper right', fontsize=10, ncols=3, framealpha=0.5)		# For 2019/03/15, Kp = 110 km/s
plt.tight_layout()
plt.savefig(main_directory+'MCMC_sf.pdf')
plt.savefig(main_directory+'MCMC_sf.jpeg', dpi=300)
plt.show()
plt.close()

## Plotting the MCMC Vsys versus injected planetary signal strength
fig = plt.figure()
ax1 = plt.subplot2grid((3,1), (0,0), rowspan=2, colspan=1)
ax2 = plt.subplot2grid((3,1), (2,0), rowspan=1, colspan=1)

ax1.axhline(Vsys, ls='dashed', c='gray')
ax1.errorbar(inj_str_MCMC_AC, Vsys_AC, yerr=u_Vsys_AC, fmt='o-', c='black', label='Astroclimes', zorder=10)
ax1.errorbar(inj_str_MCMC_MOLECFIT, Vsys_MOLECFIT, yerr=u_Vsys_MOLECFIT, fmt='o-', c='lightgray', label='Molecfit', zorder=9.5)
ax1.errorbar(inj_str_MCMC_PCA_nc_2, Vsys_PCA_nc_2, yerr=u_Vsys_PCA_nc_2, fmt='o-', c=CBcols[0], label='PCA '+r'$N_c = 2$', zorder=9)
ax1.errorbar(inj_str_MCMC_PCA_nc_3, Vsys_PCA_nc_3, yerr=u_Vsys_PCA_nc_3, fmt='o-', c=CBcols[1], label='PCA '+r'$N_c = 3$', zorder=8)
ax1.errorbar(inj_str_MCMC_PCA_nc_4, Vsys_PCA_nc_4, yerr=u_Vsys_PCA_nc_4, fmt='o-', c=CBcols[2], label='PCA '+r'$N_c = 4$', zorder=7)
ax1.errorbar(inj_str_MCMC_PCA_nc_5, Vsys_PCA_nc_5, yerr=u_Vsys_PCA_nc_5, fmt='o-', c=CBcols[3], label='PCA '+r'$N_c = 5$', zorder=6)
ax1.errorbar(inj_str_MCMC_PCA_nc_6, Vsys_PCA_nc_6, yerr=u_Vsys_PCA_nc_6, fmt='o-', c=CBcols[7], label='PCA '+r'$N_c = 6$', zorder=5)
ax1.errorbar(inj_str_MCMC_PCA_nc_10, Vsys_PCA_nc_10, yerr=u_Vsys_PCA_nc_10, fmt='o-', c=CBcols[8], label='PCA '+r'$N_c = 10$', zorder=4)
ax1.set_xticks([])
ax1.set_ylabel(r'$V_{sys}$'+' (km/s)', fontsize=14)
ax1.tick_params(axis='y', labelsize=14)
ax1.legend(loc='lower right', ncols=2, framealpha=0.5) # For 2018/03/26, Kp = 110 km/s (with dlt=0.7 as well) and 2019/03/15, Kp = 110 km/s
#ax1.legend(loc='upper right', ncols=2, framealpha=0.5) # For 2018/03/26, Kp = 80 km/s
#ax1.legend(loc='upper left', ncols=2, framealpha=0.5) # For 2018/03/26, Kp = 80 km/s

ax2.errorbar(inj_str_MCMC_AC, Vsys - Vsys_AC, yerr=u_Vsys_AC, fmt='o-', c='black', label='Astroclimes', zorder=10)
ax2.errorbar(inj_str_MCMC_MOLECFIT, Vsys - Vsys_MOLECFIT, yerr=u_Vsys_MOLECFIT, fmt='o-', c='lightgray', label='Molecfit', zorder=9.5)
ax2.errorbar(inj_str_MCMC_PCA_nc_2, Vsys - Vsys_PCA_nc_2, yerr=u_Vsys_PCA_nc_2, fmt='o-', c=CBcols[0], label='PCA '+r'$N_c = 2$', zorder=9)
ax2.errorbar(inj_str_MCMC_PCA_nc_3, Vsys - Vsys_PCA_nc_3, yerr=u_Vsys_PCA_nc_3, fmt='o-', c=CBcols[1], label='PCA '+r'$N_c = 3$', zorder=8)
ax2.errorbar(inj_str_MCMC_PCA_nc_4, Vsys - Vsys_PCA_nc_4, yerr=u_Vsys_PCA_nc_4, fmt='o-', c=CBcols[2], label='PCA '+r'$N_c = 4$', zorder=7)
ax2.errorbar(inj_str_MCMC_PCA_nc_5, Vsys - Vsys_PCA_nc_5, yerr=u_Vsys_PCA_nc_5, fmt='o-', c=CBcols[3], label='PCA '+r'$N_c = 5$', zorder=6)
ax2.errorbar(inj_str_MCMC_PCA_nc_6, Vsys - Vsys_PCA_nc_6, yerr=u_Vsys_PCA_nc_6, fmt='o-', c=CBcols[7], label='PCA '+r'$N_c = 6$', zorder=5)
ax2.errorbar(inj_str_MCMC_PCA_nc_10, Vsys - Vsys_PCA_nc_10, yerr=u_Vsys_PCA_nc_10, fmt='o-', c=CBcols[8], label='PCA '+r'$N_c = 10$', zorder=4)
ax2.axhline(0, ls='dashed', c='gray')
ax2.set_ylabel('Residuals (km/s)', fontsize=14)
ax2.set_xlabel('Injected planet signal multiplier', fontsize=14)
ax2.tick_params(axis='x', labelsize=14)
ax2.tick_params(axis='y', labelsize=14)

fig.align_ylabels()
plt.tight_layout()
plt.subplots_adjust(hspace=0)

plt.savefig(main_directory+'MCMC_Vsys.pdf')
plt.savefig(main_directory+'MCMC_Vsys.jpeg', dpi=300)

plt.show()
plt.close()

## Plotting the MCMC Kp versus injected planetary signal strength
fig = plt.figure()
ax1 = plt.subplot2grid((3,1), (0,0), rowspan=2, colspan=1)
ax2 = plt.subplot2grid((3,1), (2,0), rowspan=1, colspan=1)

ax1.axhline(Kp, ls='dashed', c='gray')
ax1.errorbar(inj_str_MCMC_AC, Kp_AC, yerr=u_Kp_AC, fmt='o-', c='black', label='Astroclimes', zorder=10)
ax1.errorbar(inj_str_MCMC_MOLECFIT, Kp_MOLECFIT, yerr=u_Kp_MOLECFIT, fmt='o-', c='lightgray', label='Molecfit', zorder=9.5)
ax1.errorbar(inj_str_MCMC_PCA_nc_2, Kp_PCA_nc_2, yerr=u_Kp_PCA_nc_2, fmt='o-', c=CBcols[0], label='PCA '+r'$N_c = 2$', zorder=9)
ax1.errorbar(inj_str_MCMC_PCA_nc_3, Kp_PCA_nc_3, yerr=u_Kp_PCA_nc_3, fmt='o-', c=CBcols[1], label='PCA '+r'$N_c = 3$', zorder=8)
ax1.errorbar(inj_str_MCMC_PCA_nc_4, Kp_PCA_nc_4, yerr=u_Kp_PCA_nc_4, fmt='o-', c=CBcols[2], label='PCA '+r'$N_c = 4$', zorder=7)
ax1.errorbar(inj_str_MCMC_PCA_nc_5, Kp_PCA_nc_5, yerr=u_Kp_PCA_nc_5, fmt='o-', c=CBcols[3], label='PCA '+r'$N_c = 5$', zorder=6)
ax1.errorbar(inj_str_MCMC_PCA_nc_6, Kp_PCA_nc_6, yerr=u_Kp_PCA_nc_6, fmt='o-', c=CBcols[7], label='PCA '+r'$N_c = 6$', zorder=5)
ax1.errorbar(inj_str_MCMC_PCA_nc_10, Kp_PCA_nc_10, yerr=u_Kp_PCA_nc_10, fmt='o-', c=CBcols[8], label='PCA '+r'$N_c = 10$', zorder=4)
ax1.set_xticks([])
ax1.set_ylabel(r'$K_p$'+' (km/s)', fontsize=14)
ax1.tick_params(axis='y', labelsize=14)
#ax1.legend(loc='upper right', ncols=2, framealpha=0.5) # For 2018/03/26, Kp = 110, 80 km/s and 2019/03/15, Kp = 110 km/s
#ax1.legend(loc='upper left', ncols=2, framealpha=0.5) # For 2018/03/26, Kp = 50 km/s
ax1.legend(loc='lower right', ncols=2, framealpha=0.5) # For 2018/03/26, Kp = 110, dlt = 0.7

ax2.errorbar(inj_str_MCMC_AC, Kp - Kp_AC, yerr=u_Kp_AC, fmt='o-', c='black', label='Astroclimes', zorder=10)
ax2.errorbar(inj_str_MCMC_MOLECFIT, Kp - Kp_MOLECFIT, yerr=u_Kp_MOLECFIT, fmt='o-', c='lightgray', label='Molecfit', zorder=9.5)
ax2.errorbar(inj_str_MCMC_PCA_nc_2, Kp - Kp_PCA_nc_2, yerr=u_Kp_PCA_nc_2, fmt='o-', c=CBcols[0], label='PCA '+r'$N_c = 2$', zorder=9)
ax2.errorbar(inj_str_MCMC_PCA_nc_3, Kp - Kp_PCA_nc_3, yerr=u_Kp_PCA_nc_3, fmt='o-', c=CBcols[1], label='PCA '+r'$N_c = 3$', zorder=8)
ax2.errorbar(inj_str_MCMC_PCA_nc_4, Kp - Kp_PCA_nc_4, yerr=u_Kp_PCA_nc_4, fmt='o-', c=CBcols[2], label='PCA '+r'$N_c = 4$', zorder=7)
ax2.errorbar(inj_str_MCMC_PCA_nc_5, Kp - Kp_PCA_nc_5, yerr=u_Kp_PCA_nc_5, fmt='o-', c=CBcols[3], label='PCA '+r'$N_c = 5$', zorder=6)
ax2.errorbar(inj_str_MCMC_PCA_nc_6, Kp - Kp_PCA_nc_6, yerr=u_Kp_PCA_nc_6, fmt='o-', c=CBcols[7], label='PCA '+r'$N_c = 6$', zorder=5)
ax2.errorbar(inj_str_MCMC_PCA_nc_10, Kp - Kp_PCA_nc_10, yerr=u_Kp_PCA_nc_10, fmt='o-', c=CBcols[8], label='PCA '+r'$N_c = 10$', zorder=4)
ax2.axhline(0, ls='dashed', c='gray')
ax2.set_ylabel('Residuals (km/s)', fontsize=14)
ax2.set_xlabel('Injected planet signal multiplier', fontsize=14)
ax2.tick_params(axis='x', labelsize=14)
ax2.tick_params(axis='y', labelsize=14)

fig.align_ylabels()
plt.tight_layout()
plt.subplots_adjust(hspace=0)

plt.savefig(main_directory+'MCMC_Kp.pdf')
plt.savefig(main_directory+'MCMC_Kp.jpeg', dpi=300)

plt.show()
plt.close()
'''
## Plotting all of the above plots on the same plot (2 column option)
fig = plt.figure(figsize=(10,7))
ax1 = plt.subplot2grid((4,2), (0,0), rowspan=2, colspan=1)
ax2 = plt.subplot2grid((4,2), (0,1), rowspan=2, colspan=1)
ax3 = plt.subplot2grid((4,2), (2,0), rowspan=1, colspan=1)
ax4 = plt.subplot2grid((4,2), (3,0), rowspan=1, colspan=1)
ax5 = plt.subplot2grid((4,2), (2,1), rowspan=1, colspan=1)
ax6 = plt.subplot2grid((4,2), (3,1), rowspan=1, colspan=1)

## SNR VS inj. signal strength
ax1.plot(inj_str_SNR_AC, SNR_inj_pos_AC, 'o-', c='black', label='Astroclimes', zorder=10)
ax1.plot(inj_str_SNR_MOLECFIT, SNR_inj_pos_MOLECFIT, 'o-', c='lightgray', label='Molecfit', zorder=9.5)
ax1.plot(inj_str_SNR_PCA_nc_2, SNR_inj_pos_PCA_nc_2, 'o-', c=CBcols[0], label='PCA '+r'$N_c = 2$', zorder=9)
ax1.plot(inj_str_SNR_PCA_nc_3, SNR_inj_pos_PCA_nc_3, 'o-', c=CBcols[1], label='PCA '+r'$N_c = 3$', zorder=8)
ax1.plot(inj_str_SNR_PCA_nc_4, SNR_inj_pos_PCA_nc_4, 'o-', c=CBcols[2], label='PCA '+r'$N_c = 4$', zorder=7)
ax1.plot(inj_str_SNR_PCA_nc_5, SNR_inj_pos_PCA_nc_5, 'o-', c=CBcols[3], label='PCA '+r'$N_c = 5$', zorder=6)
ax1.plot(inj_str_SNR_PCA_nc_6, SNR_inj_pos_PCA_nc_6, 'o-', c=CBcols[7], label='PCA '+r'$N_c = 6$', zorder=5)
ax1.plot(inj_str_SNR_PCA_nc_10, SNR_inj_pos_PCA_nc_10, 'o-', c=CBcols[8], label='PCA '+r'$N_c = 10$', zorder=4)
ax1.set_xlim(0.5,10.5)
ax1.set_ylim(0,11) 		# For 2018/03/26, Kp = 50 km/s
#ax1.set_ylim(0,13) 	# For 2018/03/26, Kp = 80 km/s
#ax1.set_ylim(-1,13) 	# For 2018/03/26, Kp = 110 km/s
#ax1.set_ylim(-1,15) 	# For 2018/03/26, Kp = 110 km/s, dlt = 0.7
#ax1.set_ylim(-3.5,22) 	# For 2019/03/15, Kp = 110 km/s
ax1.set_ylabel('SNR', fontsize=14)
ax1.set_xticks([])
#ax1.tick_params(axis='x', direction='in', pad=-15, labelsize=14)
ax1.tick_params(axis='y', labelsize=14)
ax1.axhline(4, ls='dashed', c='gray')
ax1.legend(loc='best', framealpha=0.5)
ax1.text(0.85,0.05, '(a)', transform = ax1.transAxes, fontsize=18)

## MCMC scaling factor VS inj. signal strength
ax2.errorbar(inj_str_MCMC_AC, sf_AC, yerr=u_sf_AC, fmt='o-', c='black', label='Astroclimes', zorder=10)
ax2.errorbar(inj_str_MCMC_MOLECFIT, sf_MOLECFIT, yerr=u_sf_MOLECFIT, fmt='o-', c='lightgray', label='Molecfit', zorder=9.5)
ax2.errorbar(inj_str_MCMC_PCA_nc_2, sf_PCA_nc_2, yerr=u_sf_PCA_nc_2, fmt='o-', c=CBcols[0], label='PCA '+r'$N_c = 2$', zorder=9)
ax2.errorbar(inj_str_MCMC_PCA_nc_3, sf_PCA_nc_3, yerr=u_sf_PCA_nc_3, fmt='o-', c=CBcols[1], label='PCA '+r'$N_c = 3$', zorder=8)
ax2.errorbar(inj_str_MCMC_PCA_nc_4, sf_PCA_nc_4, yerr=u_sf_PCA_nc_4, fmt='o-', c=CBcols[2], label='PCA '+r'$N_c = 4$', zorder=7)
ax2.errorbar(inj_str_MCMC_PCA_nc_5, sf_PCA_nc_5, yerr=u_sf_PCA_nc_5, fmt='o-', c=CBcols[3], label='PCA '+r'$N_c = 5$', zorder=6)
ax2.errorbar(inj_str_MCMC_PCA_nc_6, sf_PCA_nc_6, yerr=u_sf_PCA_nc_6, fmt='o-', c=CBcols[7], label='PCA '+r'$N_c = 6$', zorder=5)
ax2.errorbar(inj_str_MCMC_PCA_nc_10, sf_PCA_nc_10, yerr=u_sf_PCA_nc_10, fmt='o-', c=CBcols[8], label='PCA '+r'$N_c = 10$', zorder=4)
ax2.set_ylabel('Scaling factor '+r'$\delta$', fontsize=14)
#ax2.set_xlim(0.5,10.5)
ax2.set_ylim(-0.1,1.3)	# For 2018/03/26, Kp = 50 km/s
ax2.set_xticks([])
ax2.tick_params(axis='y', labelsize=14)
ax2.text(0.85,0.05, '(b)', transform = ax2.transAxes, fontsize=18)

## MCMC Vsys VS inj. signal strength
ax3.axhline(Vsys, ls='dashed', c='gray')
ax3.errorbar(inj_str_MCMC_AC, Vsys_AC, yerr=u_Vsys_AC, fmt='o-', c='black', label='Astroclimes', zorder=10)
ax3.errorbar(inj_str_MCMC_MOLECFIT, Vsys_MOLECFIT, yerr=u_Vsys_MOLECFIT, fmt='o-', c='lightgray', label='Molecfit', zorder=9.5)
ax3.errorbar(inj_str_MCMC_PCA_nc_2, Vsys_PCA_nc_2, yerr=u_Vsys_PCA_nc_2, fmt='o-', c=CBcols[0], label='PCA '+r'$N_c = 2$', zorder=9)
ax3.errorbar(inj_str_MCMC_PCA_nc_3, Vsys_PCA_nc_3, yerr=u_Vsys_PCA_nc_3, fmt='o-', c=CBcols[1], label='PCA '+r'$N_c = 3$', zorder=8)
ax3.errorbar(inj_str_MCMC_PCA_nc_4, Vsys_PCA_nc_4, yerr=u_Vsys_PCA_nc_4, fmt='o-', c=CBcols[2], label='PCA '+r'$N_c = 4$', zorder=7)
ax3.errorbar(inj_str_MCMC_PCA_nc_5, Vsys_PCA_nc_5, yerr=u_Vsys_PCA_nc_5, fmt='o-', c=CBcols[3], label='PCA '+r'$N_c = 5$', zorder=6)
ax3.errorbar(inj_str_MCMC_PCA_nc_6, Vsys_PCA_nc_6, yerr=u_Vsys_PCA_nc_6, fmt='o-', c=CBcols[7], label='PCA '+r'$N_c = 6$', zorder=5)
ax3.errorbar(inj_str_MCMC_PCA_nc_10, Vsys_PCA_nc_10, yerr=u_Vsys_PCA_nc_10, fmt='o-', c=CBcols[8], label='PCA '+r'$N_c = 10$', zorder=4)
ax3.set_xticks([])
ax3.set_ylabel(r'$V_{sys}$'+' (km/s)', fontsize=14)
ax3.tick_params(axis='y', labelsize=14)
#ax3.text(0.85,0.05, '(c)', transform = ax3.transAxes, fontsize=18)
ax3.spines['bottom'].set_linestyle((0,(3,5,1,5)))
ax3.set_xlim(0.5,10.5)

ax4.errorbar(inj_str_MCMC_AC, Vsys - Vsys_AC, yerr=u_Vsys_AC, fmt='o-', c='black', label='Astroclimes', zorder=10)
ax4.errorbar(inj_str_MCMC_MOLECFIT, Vsys - Vsys_MOLECFIT, yerr=u_Vsys_MOLECFIT, fmt='o-', c='lightgray', label='Molecfit', zorder=9.5)
ax4.errorbar(inj_str_MCMC_PCA_nc_2, Vsys - Vsys_PCA_nc_2, yerr=u_Vsys_PCA_nc_2, fmt='o-', c=CBcols[0], label='PCA '+r'$N_c = 2$', zorder=9)
ax4.errorbar(inj_str_MCMC_PCA_nc_3, Vsys - Vsys_PCA_nc_3, yerr=u_Vsys_PCA_nc_3, fmt='o-', c=CBcols[1], label='PCA '+r'$N_c = 3$', zorder=8)
ax4.errorbar(inj_str_MCMC_PCA_nc_4, Vsys - Vsys_PCA_nc_4, yerr=u_Vsys_PCA_nc_4, fmt='o-', c=CBcols[2], label='PCA '+r'$N_c = 4$', zorder=7)
ax4.errorbar(inj_str_MCMC_PCA_nc_5, Vsys - Vsys_PCA_nc_5, yerr=u_Vsys_PCA_nc_5, fmt='o-', c=CBcols[3], label='PCA '+r'$N_c = 5$', zorder=6)
ax4.errorbar(inj_str_MCMC_PCA_nc_6, Vsys - Vsys_PCA_nc_6, yerr=u_Vsys_PCA_nc_6, fmt='o-', c=CBcols[7], label='PCA '+r'$N_c = 6$', zorder=5)
ax4.errorbar(inj_str_MCMC_PCA_nc_10, Vsys - Vsys_PCA_nc_10, yerr=u_Vsys_PCA_nc_10, fmt='o-', c=CBcols[8], label='PCA '+r'$N_c = 10$', zorder=4)
ax4.axhline(0, ls='dashed', c='gray')
ax4.set_ylabel('Residuals (km/s)', fontsize=14)
ax4.set_xlabel('Injected planet signal multiplier', fontsize=14)
ax4.tick_params(axis='x', labelsize=14)
ax4.tick_params(axis='y', labelsize=14)
ax4.text(0.85,0.05, '(c)', transform = ax4.transAxes, fontsize=18)
ax4.spines['top'].set_linestyle((0,(3,5,1,5)))
ax4.set_xlim(0.5,10.5)

## MCMC Kp VS inj. signal strength
ax5.axhline(Kp, ls='dashed', c='gray')
ax5.errorbar(inj_str_MCMC_AC, Kp_AC, yerr=u_Kp_AC, fmt='o-', c='black', label='Astroclimes', zorder=10)
ax5.errorbar(inj_str_MCMC_MOLECFIT, Kp_MOLECFIT, yerr=u_Kp_MOLECFIT, fmt='o-', c='lightgray', label='Molecfit', zorder=9.5)
ax5.errorbar(inj_str_MCMC_PCA_nc_2, Kp_PCA_nc_2, yerr=u_Kp_PCA_nc_2, fmt='o-', c=CBcols[0], label='PCA '+r'$N_c = 2$', zorder=9)
ax5.errorbar(inj_str_MCMC_PCA_nc_3, Kp_PCA_nc_3, yerr=u_Kp_PCA_nc_3, fmt='o-', c=CBcols[1], label='PCA '+r'$N_c = 3$', zorder=8)
ax5.errorbar(inj_str_MCMC_PCA_nc_4, Kp_PCA_nc_4, yerr=u_Kp_PCA_nc_4, fmt='o-', c=CBcols[2], label='PCA '+r'$N_c = 4$', zorder=7)
ax5.errorbar(inj_str_MCMC_PCA_nc_5, Kp_PCA_nc_5, yerr=u_Kp_PCA_nc_5, fmt='o-', c=CBcols[3], label='PCA '+r'$N_c = 5$', zorder=6)
ax5.errorbar(inj_str_MCMC_PCA_nc_6, Kp_PCA_nc_6, yerr=u_Kp_PCA_nc_6, fmt='o-', c=CBcols[7], label='PCA '+r'$N_c = 6$', zorder=5)
ax5.errorbar(inj_str_MCMC_PCA_nc_10, Kp_PCA_nc_10, yerr=u_Kp_PCA_nc_10, fmt='o-', c=CBcols[8], label='PCA '+r'$N_c = 10$', zorder=4)
ax5.set_xticks([])
ax5.set_ylabel(r'$K_p$'+' (km/s)', fontsize=14)
ax5.tick_params(axis='y', labelsize=14)
#ax5.text(0.85,0.05, '(d)', transform = ax5.transAxes, fontsize=18)
ax5.spines['bottom'].set_linestyle((0,(3,5,1,5)))
#ax5.set_xlim(0.5,10.5)

ax6.errorbar(inj_str_MCMC_AC, Kp - Kp_AC, yerr=u_Kp_AC, fmt='o-', c='black', label='Astroclimes', zorder=10)
ax6.errorbar(inj_str_MCMC_MOLECFIT, Kp - Kp_MOLECFIT, yerr=u_Kp_MOLECFIT, fmt='o-', c='lightgray', label='Molecfit', zorder=9.5)
ax6.errorbar(inj_str_MCMC_PCA_nc_2, Kp - Kp_PCA_nc_2, yerr=u_Kp_PCA_nc_2, fmt='o-', c=CBcols[0], label='PCA '+r'$N_c = 2$', zorder=9)
ax6.errorbar(inj_str_MCMC_PCA_nc_3, Kp - Kp_PCA_nc_3, yerr=u_Kp_PCA_nc_3, fmt='o-', c=CBcols[1], label='PCA '+r'$N_c = 3$', zorder=8)
ax6.errorbar(inj_str_MCMC_PCA_nc_4, Kp - Kp_PCA_nc_4, yerr=u_Kp_PCA_nc_4, fmt='o-', c=CBcols[2], label='PCA '+r'$N_c = 4$', zorder=7)
ax6.errorbar(inj_str_MCMC_PCA_nc_5, Kp - Kp_PCA_nc_5, yerr=u_Kp_PCA_nc_5, fmt='o-', c=CBcols[3], label='PCA '+r'$N_c = 5$', zorder=6)
ax6.errorbar(inj_str_MCMC_PCA_nc_6, Kp - Kp_PCA_nc_6, yerr=u_Kp_PCA_nc_6, fmt='o-', c=CBcols[7], label='PCA '+r'$N_c = 6$', zorder=5)
ax6.errorbar(inj_str_MCMC_PCA_nc_10, Kp - Kp_PCA_nc_10, yerr=u_Kp_PCA_nc_10, fmt='o-', c=CBcols[8], label='PCA '+r'$N_c = 10$', zorder=4)
ax6.axhline(0, ls='dashed', c='gray')
ax6.set_ylabel('Residuals (km/s)', fontsize=14)
ax6.set_xlabel('Injected planet signal multiplier', fontsize=14)
ax6.tick_params(axis='x', labelsize=14)
ax6.tick_params(axis='y', labelsize=14)
ax6.text(0.85,0.05, '(d)', transform = ax6.transAxes, fontsize=18)
ax6.spines['top'].set_linestyle((0,(3,5,1,5)))
#ax6.set_xlim(0.5,10.5)

fig.align_ylabels()
plt.tight_layout()
plt.subplots_adjust(hspace=0)

plt.savefig(main_directory+'all_plots_2columns.pdf', bbox_inches='tight')
plt.savefig(main_directory+'all_plots_2columns.jpeg', dpi=300, bbox_inches='tight')

plt.show()
plt.close()
'''
## Unpacking the SNR values from the Delta_CCF
filename_SNR_results_PCA_nc_2_DCCF = main_directory+'SNR_logs_pca_DeltaCCF_nc_2.txt'
inj_str_SNR_PCA_nc_2_DCCF, SNR_inj_pos_PCA_nc_2_DCCF = np.loadtxt(filename_SNR_results_PCA_nc_2_DCCF, usecols=(0,2), unpack=True)

filename_SNR_results_PCA_nc_3_DCCF = main_directory+'SNR_logs_pca_DeltaCCF_nc_3.txt'
inj_str_SNR_PCA_nc_3_DCCF, SNR_inj_pos_PCA_nc_3_DCCF = np.loadtxt(filename_SNR_results_PCA_nc_3_DCCF, usecols=(0,2), unpack=True)

filename_SNR_results_PCA_nc_4_DCCF = main_directory+'SNR_logs_pca_DeltaCCF_nc_4.txt'
inj_str_SNR_PCA_nc_4_DCCF, SNR_inj_pos_PCA_nc_4_DCCF = np.loadtxt(filename_SNR_results_PCA_nc_4_DCCF, usecols=(0,2), unpack=True)

filename_SNR_results_PCA_nc_5_DCCF = main_directory+'SNR_logs_pca_DeltaCCF_nc_5.txt'
inj_str_SNR_PCA_nc_5_DCCF, SNR_inj_pos_PCA_nc_5_DCCF = np.loadtxt(filename_SNR_results_PCA_nc_5_DCCF, usecols=(0,2), unpack=True)

filename_SNR_results_PCA_nc_6_DCCF = main_directory+'SNR_logs_pca_DeltaCCF_nc_6.txt'
inj_str_SNR_PCA_nc_6_DCCF, SNR_inj_pos_PCA_nc_6_DCCF = np.loadtxt(filename_SNR_results_PCA_nc_6_DCCF, usecols=(0,2), unpack=True)

filename_SNR_results_PCA_nc_10_DCCF = main_directory+'SNR_logs_pca_DeltaCCF_nc_10.txt'
inj_str_SNR_PCA_nc_10_DCCF, SNR_inj_pos_PCA_nc_10_DCCF = np.loadtxt(filename_SNR_results_PCA_nc_10_DCCF, usecols=(0,2), unpack=True)

#filename_SNR_results_PCA_nc_12_DCCF = main_directory+'SNR_logs_pca_DeltaCCF_nc_12.txt'
#inj_str_SNR_PCA_nc_12_DCCF, SNR_inj_pos_PCA_nc_12_DCCF = np.loadtxt(filename_SNR_results_PCA_nc_12_DCCF, usecols=(0,2), unpack=True)

#filename_SNR_results_PCA_nc_13_DCCF = main_directory+'SNR_logs_pca_DeltaCCF_nc_13.txt'
#inj_str_SNR_PCA_nc_13_DCCF, SNR_inj_pos_PCA_nc_13_DCCF = np.loadtxt(filename_SNR_results_PCA_nc_13_DCCF, usecols=(0,2), unpack=True)

#filename_SNR_results_PCA_nc_14_DCCF = main_directory+'SNR_logs_pca_DeltaCCF_nc_14.txt'
#inj_str_SNR_PCA_nc_14_DCCF, SNR_inj_pos_PCA_nc_14_DCCF = np.loadtxt(filename_SNR_results_PCA_nc_14_DCCF, usecols=(0,2), unpack=True)

#filename_SNR_results_PCA_nc_15_DCCF = main_directory+'SNR_logs_pca_DeltaCCF_nc_15.txt'
#inj_str_SNR_PCA_nc_15_DCCF, SNR_inj_pos_PCA_nc_15_DCCF = np.loadtxt(filename_SNR_results_PCA_nc_15_DCCF, usecols=(0,2), unpack=True)

#filename_SNR_results_PCA_nc_16_DCCF = main_directory+'SNR_logs_pca_DeltaCCF_nc_16.txt'
#inj_str_SNR_PCA_nc_16_DCCF, SNR_inj_pos_PCA_nc_16_DCCF = np.loadtxt(filename_SNR_results_PCA_nc_16_DCCF, usecols=(0,2), unpack=True)

## Creating new variables to plot the effect on the SNR/planet signal strength as a function of PCA components for a given injected planet signal strength
n_PCA_comps = np.array([2, 3, 4, 5, 6, 10])
#n_PCA_comps = np.array([2, 3, 4, 5, 6, 10, 12, 13, 14, 15, 16])
#inj_sig_str = 10
for inj_sig_str in [1,2,3,3.25,3.50,3.75,4,4.25,4.50,4.75,5,6,7,8,9,10]:
	SNR_per_PCA_comp = np.array([SNR_inj_pos_PCA_nc_2[np.where(inj_str_SNR_PCA_nc_2 == inj_sig_str)], SNR_inj_pos_PCA_nc_3[np.where(inj_str_SNR_PCA_nc_3 == inj_sig_str)],
								SNR_inj_pos_PCA_nc_4[np.where(inj_str_SNR_PCA_nc_4 == inj_sig_str)], SNR_inj_pos_PCA_nc_5[np.where(inj_str_SNR_PCA_nc_5 == inj_sig_str)],
								SNR_inj_pos_PCA_nc_6[np.where(inj_str_SNR_PCA_nc_6 == inj_sig_str)], SNR_inj_pos_PCA_nc_10[np.where(inj_str_SNR_PCA_nc_10 == inj_sig_str)]])
#								SNR_inj_pos_PCA_nc_12[np.where(inj_str_SNR_PCA_nc_12 == inj_sig_str)], SNR_inj_pos_PCA_nc_13[np.where(inj_str_SNR_PCA_nc_13 == inj_sig_str)],
#								SNR_inj_pos_PCA_nc_14[np.where(inj_str_SNR_PCA_nc_14 == inj_sig_str)], SNR_inj_pos_PCA_nc_15[np.where(inj_str_SNR_PCA_nc_15 == inj_sig_str)],
#								SNR_inj_pos_PCA_nc_16[np.where(inj_str_SNR_PCA_nc_16 == inj_sig_str)]])

	SNR_DCCF_per_PCA_comp = np.array([	SNR_inj_pos_PCA_nc_2_DCCF[np.where(inj_str_SNR_PCA_nc_2_DCCF == inj_sig_str)], SNR_inj_pos_PCA_nc_3_DCCF[np.where(inj_str_SNR_PCA_nc_3_DCCF == inj_sig_str)],
										SNR_inj_pos_PCA_nc_4_DCCF[np.where(inj_str_SNR_PCA_nc_4_DCCF == inj_sig_str)], SNR_inj_pos_PCA_nc_5_DCCF[np.where(inj_str_SNR_PCA_nc_5_DCCF == inj_sig_str)],
										SNR_inj_pos_PCA_nc_6_DCCF[np.where(inj_str_SNR_PCA_nc_6_DCCF == inj_sig_str)], SNR_inj_pos_PCA_nc_10_DCCF[np.where(inj_str_SNR_PCA_nc_10_DCCF == inj_sig_str)]])
#										SNR_inj_pos_PCA_nc_12_DCCF[np.where(inj_str_SNR_PCA_nc_12_DCCF == inj_sig_str)], SNR_inj_pos_PCA_nc_13_DCCF[np.where(inj_str_SNR_PCA_nc_13_DCCF == inj_sig_str)],
#										SNR_inj_pos_PCA_nc_14_DCCF[np.where(inj_str_SNR_PCA_nc_14_DCCF == inj_sig_str)], SNR_inj_pos_PCA_nc_15_DCCF[np.where(inj_str_SNR_PCA_nc_15_DCCF == inj_sig_str)],
#										SNR_inj_pos_PCA_nc_16_DCCF[np.where(inj_str_SNR_PCA_nc_16_DCCF == inj_sig_str)]])

	#sf_per_PCA_comp = np.array([sf_PCA_nc_2[np.where(inj_str_MCMC_PCA_nc_2 == inj_sig_str)], sf_PCA_nc_3[np.where(inj_str_MCMC_PCA_nc_3 == inj_sig_str)],
	#							sf_PCA_nc_4[np.where(inj_str_MCMC_PCA_nc_4 == inj_sig_str)], sf_PCA_nc_5[np.where(inj_str_MCMC_PCA_nc_5 == inj_sig_str)],
	#							sf_PCA_nc_6[np.where(inj_str_MCMC_PCA_nc_6 == inj_sig_str)], sf_PCA_nc_10[np.where(inj_str_MCMC_PCA_nc_10 == inj_sig_str)]])

	plt.figure()
	plt.plot(n_PCA_comps, SNR_per_PCA_comp, 'o-', c='black', label='CCF'+r'$_{inj}$')
	plt.plot(n_PCA_comps, SNR_DCCF_per_PCA_comp, 'o-', c='red', label=r'$\Delta$'+'CCF')
	plt.plot(n_PCA_comps[np.argmax(SNR_per_PCA_comp)], np.max(SNR_per_PCA_comp), '*', ms=12., c='black')
	plt.plot(n_PCA_comps[np.argmax(SNR_DCCF_per_PCA_comp)], np.max(SNR_DCCF_per_PCA_comp), '*', ms=12., c='red')
	plt.ylabel('SNR', fontsize=14)
	plt.xlabel('PCA components', fontsize=14)
	plt.xticks([2,3,4,5,6,10])
	#plt.xticks([2,3,4,5,6,10,12,13,14,15,16])
	plt.tick_params(axis='x', labelsize=14)
	plt.tick_params(axis='y', labelsize=14)
	plt.legend(loc='best', framealpha=0.5)
	plt.tight_layout()
	plt.savefig(main_directory+f'SNR_per_PCA_comp_{inj_sig_str}x.pdf', bbox_inches='tight')
	plt.savefig(main_directory+f'SNR_per_PCA_comp_{inj_sig_str}x.jpeg', dpi=300, bbox_inches='tight')
	plt.show()
	plt.close()

#plt.figure()
#plt.plot(n_PCA_comps, sf_per_PCA_comp, 'o-', c='black')
#plt.ylabel('Scaling factor '+r'$\delta$', fontsize=14)
#plt.xlabel('PCA components', fontsize=14)
#plt.tick_params(axis='x', labelsize=14)
#plt.tick_params(axis='y', labelsize=14)
#plt.tight_layout()
#plt.savefig(main_directory+f'sf_per_PCA_comp_{inj_sig_str}x.pdf')
#plt.savefig(main_directory+f'sf_per_PCA_comp_{inj_sig_str}x.jpeg', dpi=300)
#plt.show()
#plt.close()
'''
print('This is the end')