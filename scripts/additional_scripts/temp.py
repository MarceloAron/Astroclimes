'''
Script to combine the individual files containing the SNR and MCMC results into just one file (if, for example, the algorithm was run with SLURM)
Simply change the name of the relevant files (you may need to change the indexes in te readline() statements, depending on how your files were generated)
It can also be used to check the slurm.out files to see if all of the runs finished properly
'''

# =====================================================================================
# Basic packages
# =====================================================================================
import glob
from subprocess import call

#for nc_PCA in [2,3,4,5,6,10]:
#file_write = '/media/marceloaron/New Volume/PhD/thesis_work/MCMC_results/tauboo/2018_03_26/molecfit/results_mcmc.txt'
#files_read_common_name = '/media/marceloaron/New Volume/PhD/thesis_work/MCMC_results/tauboo/2018_03_26/molecfit/results_mcmc_'
file_write = f'/media/marceloaron/New Volume/PhD/thesis_work/telluric_removal/TauBoo/2018_03_26_Kp_50/molecfit/results_inj_rec_mcmc_molecfit.txt'
files_read_common_name = f'/media/marceloaron/New Volume/PhD/thesis_work/telluric_removal/TauBoo/2018_03_26_Kp_50/molecfit/results_inj_rec_mcmc_molecfit_sf_'
files_read_suffix = '.txt'

wtxt = open(file_write, 'w') 
rtxt = open(glob.glob(files_read_common_name+'*')[0], 'r')
text = rtxt.readlines()[0] 
wtxt.write(text) 
#for i in range(len(glob.glob(files_read_common_name+'*'))):
for i in [1,2,3,3.25,3.50,3.75,4,4.25,4.50,4.75,5,6,7,8,9,10]:
	try:
		rtxt = open(files_read_common_name+str(i).replace('.','_')+files_read_suffix, 'r')
		text = rtxt.readlines()[2] 
	except IndexError:
		rtxt = open(files_read_common_name+str(i).replace('.','_')+files_read_suffix, 'r')
		text = rtxt.readlines()[1] 
	wtxt.write(text) 
wtxt.close()


#for nc_PCA in [2,3,4,5,6,10]:
file_write = f'/media/marceloaron/New Volume/PhD/thesis_work/telluric_removal/TauBoo/2018_03_26_Kp_50/molecfit/SNR_logs_molecfit.txt'
files_read_common_name = f'/media/marceloaron/New Volume/PhD/thesis_work/telluric_removal/TauBoo/2018_03_26_Kp_50/molecfit/SNR_logs_molecfit_sf_'
files_read_suffix = '.txt'

wtxt = open(file_write, 'w') 
#rtxt = open(files_read_common_name+'1'+files_read_suffix, 'r')
rtxt = open(f'/media/marceloaron/New Volume/PhD/thesis_work/telluric_removal/TauBoo/2018_03_26_Kp_50/molecfit/SNR_logs_molecfit_sf_1.txt', 'r')
text = rtxt.readlines()
wtxt.write(text[0]) 
wtxt.write(text[1]) 
wtxt.write(text[2]) 
#for i in range(len(glob.glob(files_read_common_name+'*'))):
for i in [1,2,3,3.25,3.50,3.75,4,4.25,4.50,4.75,5,6,7,8,9,10]:
	rtxt = open(files_read_common_name+str(i).replace('.','_')+files_read_suffix, 'r')
	text = rtxt.readlines()[3] 
	wtxt.write(text) 
wtxt.close()

#files_slurm = glob.glob('/media/marceloaron/New Volume/PhD/thesis_work/telluric_removal/TauBoo/2018_03_26/molecfit/slurm_logs/*.out')
#files_slurm.sort()
#for f in files_slurm:
#	print(f)
#	call(['cat', f])