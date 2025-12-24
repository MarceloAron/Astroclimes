import numpy as np
from subprocess import call
from datetime import datetime, timedelta
import glob

'''
This script can be used to create and download a list of auxiliary files as well as to generate atmospheric profiles using Ginput
'''

## Specify the home directory
home_directory = '/home/marceloaron/MarceloAron/PhD/thesis_work/'

def GEOSdate_2_datetime(date):
	'''
	Converting from  the date format needed to run Ginput to datetime format
	'''
	year = int(date[:4])
	month = int(date[4:6])
	day = int(date[6:8])
	hour = int(date[9:])
	
	return datetime(year, month, day, hour)

def datetime_2_GEOSdate(date):
	'''
	Converting from datetime format to the date format needed to run Ginput
	'''
	return f'{date.year:4d}{date.month:02d}{date.day:02d}_{date.hour:02d}'

## Here we will get the dates from the list of GEOSFP files
path_GEOS = home_directory+'Astroclimes/data/GEOSFP/'
file_list = glob.glob(path_GEOS+'Nv/*asm_Nv*.nc4')
file_list.sort()
dates = np.array([f.split('/')[-1].split('.')[-3][:-2] for f in file_list])

save_path_mod = home_directory+'Astroclimes/atmosphere_profiles/GGG2020/'
site = 'al'
mode = 'fp-eta'
product = 'fp'
map_format = 'txt'

for i in range(len(dates)):
	start_date = dates[i]
	end_date = datetime_2_GEOSdate(GEOSdate_2_datetime(dates[i]) + timedelta(hours=3))
	print(f'Progress:{i+1:d}/{len(dates):d}',f'Date range: {start_date}-{end_date}')
	call(['./run_ginput.py', 'tccon-mod', f'{start_date}-{end_date}', path_GEOS, '--site', site, '--include-chem', '--mode', mode, '-s', save_path_mod])
	call(['./run_ginput.py', 'vmr', f'{start_date}-{end_date}', save_path_mod+f'fp/{site}/vertical/', '--site', site, '--product', product, '-s', save_path_mod])
	call(['./run_ginput.py', 'map', '-r', save_path_mod, '-s', save_path_mod+f'fp/{site}/maps-vertical/', '--product', product, '--site', site, '-f', map_format, f'{start_date}-{end_date}'])
