'''
Script to combine the individual files containing the SNR and MCMC results into just one file (if, for example, the algorithm was run with SLURM)
Simply change the name of the relevant files (you may need to change the indexes in te readline() statements, depending on how your files were generated)
'''
# =====================================================================================
# Basic packages
# =====================================================================================
import glob

file_write = '/media/marceloaron/New Volume/PhD/thesis_work/telluric_removal/TauBoo/2018_03_26/rerun_test/results_telrem_mcmc_pca_nc_10.txt'
files_read_common_name = '/media/marceloaron/New Volume/PhD/thesis_work/telluric_removal/TauBoo/2018_03_26/rerun_test/results_telrem_mcmc_pca_nc_10_sf_'
files_read_suffix = '.txt'

wtxt = open(file_write, 'w') 
rtxt = open(glob.glob(files_read_common_name+'*')[0], 'r')
text = rtxt.readlines()[0] 
wtxt.write(text) 
#for i in range(len(glob.glob(files_read_common_name+'*'))):
for i in [1,2,3,3.25,3.50,3.75,4,4.25,4.50,4.75,5,6,7,8,9,10]:
	rtxt = open(files_read_common_name+f'{i}'+files_read_suffix, 'r')
	text = rtxt.readlines()[1] 
	wtxt.write(text) 
wtxt.close()

file_write = '/media/marceloaron/New Volume/PhD/thesis_work/telluric_removal/TauBoo/2018_03_26/rerun_test/SNR_logs_pca_nc_10.txt'
files_read_common_name = '/media/marceloaron/New Volume/PhD/thesis_work/telluric_removal/TauBoo/2018_03_26/rerun_test/SNR_logs_pca_nc_10_sf_'
files_read_suffix = '.txt'

wtxt = open(file_write, 'w') 
rtxt = open(glob.glob(files_read_common_name+'*')[0], 'r')
text = rtxt.readlines()
#wtxt.write(text[0]) 
#wtxt.write(text[1]) 
#wtxt.write(text[2]) 
#for i in range(len(glob.glob(files_read_common_name+'*'))):
for i in [1,2,3,3.25,3.50,3.75,4,4.25,4.50,4.75,5,6,7,8,9,10]:
	rtxt = open(files_read_common_name+f'{i}'+files_read_suffix, 'r')
	text = rtxt.readlines()[0] 
	wtxt.write(text) 
wtxt.close()