import numpy as np
from subprocess import call
import os
from datetime import datetime, timedelta
import glob

'''
This script can be used to create and download a list of the GEOS auxiliary files needed to generate atmospheric profiles with ginput
'''

## Specify the "home" directory
home_directory = '/home/marceloaron/MarceloAron/PhD/thesis_work/'

## Parent directory where subdirectories will be created
## Because of how Ginput works, the GEOS-FP files must be organised in a parent directory GEOSFP/ 
## with the files divided in the subdirectories Nv/ and Nx/ depending on their type
parent_directory = home_directory+'Astroclimes/data/GEOSFP/'

## Creating parent directory, if it doesn't already exist
if not os.path.isdir(parent_directory):
	call(['mkdir', parent_directory])

## Creating Nv subdirectory, if it doesn't already exist
if not os.path.isdir(parent_directory+'Nv/'):
	call(['mkdir', parent_directory+'Nv/'])

## Creating Nx subdirectory, if it doesn't already exist
if not os.path.isdir(parent_directory+'Nx/'):
	call(['mkdir', parent_directory+'Nx/'])

def create_GEOS_FP_file_list(dates):
	'''
	Creating a list of URLs to download the relevant auxiliary files to generate Ginput (GGG2020) atmospheric profiles
	This function takes a list of dates and generates a list containing only the URLs for the two auxiliary files closest to it
	'''

	## GEOS-FP files are available every 3 hours
	GEOS_FP_hours = np.array([0, 3, 6, 9, 12, 15, 18, 21])
	txt_asm_3D = open(parent_directory+'Nv/getFP_asm_Nv.dat', 'w')
	txt_chm_3D = open(parent_directory+'Nv/getFP_chm_Nv.dat', 'w')
	txt_asm_2D = open(parent_directory+'Nx/getFP_asm_Nx.dat', 'w')
	selected_GEOS_FP_dates = np.empty(0)
	for date in dates:
		year = int(date[:4])
		month = int(date[4:6])
		day = int(date[6:8])
		hour = int(date[9:11])
		minute = int(date[12:14])
		second = int(date[15:17])
		## Getting the current date in datetime format
		obs_date = datetime(year, month, day, hour, minute, second)
		## Creating a list of all the possible GEOS-FP timestamps for that day, plus the last timestamp of the previous day and the first timestamp of the next day
		GEOS_FP_dates = np.array([datetime(year, month, day, 21, 0, 0)-timedelta(days=1)]+[datetime(year, month, day, G_h, 0, 0) for G_h in GEOS_FP_hours]+[datetime(year, month, day, 0, 0, 0)+timedelta(days=1)])
		## Calculating the difference (in seconds) between the GEOS-FP timestamps and the actual time of observation
		time_deltas_secs = np.array([abs((GEOS_FP_dates[i] - obs_date).total_seconds()) for i in range(len(GEOS_FP_dates))])
		## Selecting only the two closest timestamps to the time of observation
		selected_GEOS_FP_dates = np.concatenate((selected_GEOS_FP_dates, np.sort(GEOS_FP_dates[np.argsort(time_deltas_secs)][:2])))

	## Making sure there are no duplicate timestamps to avoid downloading the same file twice
	selected_GEOS_FP_dates = np.unique(selected_GEOS_FP_dates)
	for sd in selected_GEOS_FP_dates:
		txt_asm_3D.write(f"https://portal.nccs.nasa.gov/datashare/gmao_ops/pub/fp/das/Y{sd.year}/M{sd.month:0>2}/D{sd.day:0>2}/GEOS.fp.asm.inst3_3d_asm_Nv.{sd.year}{sd.month:0>2}{sd.day:0>2}_{sd.hour:0>2}{sd.minute:0>2}.V01.nc4\n")
		txt_chm_3D.write(f"https://portal.nccs.nasa.gov/datashare/gmao_ops/pub/fp/das/Y{sd.year}/M{sd.month:0>2}/D{sd.day:0>2}/GEOS.fp.asm.inst3_3d_chm_Nv.{sd.year}{sd.month:0>2}{sd.day:0>2}_{sd.hour:0>2}{sd.minute:0>2}.V01.nc4\n")
		txt_asm_2D.write(f"https://portal.nccs.nasa.gov/datashare/gmao_ops/pub/fp/das/Y{sd.year}/M{sd.month:0>2}/D{sd.day:0>2}/GEOS.fp.asm.inst3_2d_asm_Nx.{sd.year}{sd.month:0>2}{sd.day:0>2}_{sd.hour:0>2}{sd.minute:0>2}.V01.nc4\n")

## Getting the list of dates for which GEOS-FP files are needed, from the observational spectra
spec_file_list = glob.glob(home_directory+'Astroclimes/data/CARMENES/*-nir_A.fits') 
spec_file_list.sort()
spec_dates = np.array([f.split('/')[-1].split('-')[1] for f in spec_file_list])
create_GEOS_FP_file_list(spec_dates)

## Downloading the GEOS-FP files from the URLs list just created
download_asm_3D = True
if download_asm_3D:
	GEOS_FP_file_asm_3D = np.loadtxt(parent_directory+'Nv/getFP_asm_Nv.dat', dtype=str)
	for f in GEOS_FP_file_asm_3D:
		call(['wget', f, '-P', parent_directory+'Nv/'])

download_chm_3D = True
if download_chm_3D:
	GEOS_FP_file_chm_3D = np.loadtxt(parent_directory+'Nv/getFP_chm_Nv.dat', dtype=str)
	for f in GEOS_FP_file_chm_3D:
		call(['wget', f, '-P', parent_directory+'Nv/'])

download_asm_2D = True
if download_asm_2D:
	GEOS_FP_file_asm_2D = np.loadtxt(parent_directory+'Nx/getFP_asm_Nx.dat', dtype=str)
	for f in GEOS_FP_file_asm_2D:
		call(['wget', f, '-P', parent_directory+'Nx/'])