#!/usr/local/bin/python
"""  This script generates temp.th and salt.th for a selfe run.

Db22 uses values from Beaver Army and Longview.
Db26 uses values from Bonneville, Newberg, Longview, and/or Warrendale
"""
import sys
import os
import time
import datetime
import linecache
import numpy 

TH_PATH = "/home/workspace/ccalmr/elcirc/inputs/"
DATA_PATH  = "/home/workspace/ccalmr/data/external/river_data/verified/"
DT = 90
file_num = 0

def __readData(start_jd, n_days, DATA_PATH, file, run_times):
	global file_num
	file_num += 1
	#log = open("log_%i.txt" % file_num, "wb")
# If file is newberg_discharge it is update every 15 min, otherwise 30
	if (file.find("newberg_discharge") != -1):
		FREQ = 1/48.0
		DIV  = 48
	elif (file.find("beaverarmy_discharge") != -1):
		FREQ = 1/48.0
		DIV  = 48
	elif (file.find("longview") != -1):
		FREQ = 1/240.0
		DIV  = 240
	else:
		FREQ = 1/24.0
		DIV  = 24
	_data = []
	_time = []
# Find the first day number in the file.  Then find the days of interest for
# the run in the time frame of the file date scheme.  I think this is corie
# time.
	file_day_1 = linecache.getline("%s/%s" % (DATA_PATH,file), 1)
	if file_day_1 == '':
		print("Error: Unable to open %s/%s" % (DATA_PATH,file))
		raise IOError
	file_day_1 = file_day_1.split()
	file_day_1 = int(float(file_day_1[0]))
	file_start_day = file_day_1 + (start_jd-1)
	file_end_day   = file_start_day + n_days
	end_jd = start_jd + n_days
	#log.write("%s\n" % file)
	#log.write("First day in file: %i\n" 
	#          "File start day: %i\n" 
	#		  "File end day: %i\n" 
	#		  "JD start day: %i\n"
	#		  "JD end day: %i\n" % (file_day_1, file_start_day, file_end_day, start_jd, end_jd))
# Check if the predicted line number is the correct day.  Only need to check
# because data is missing in some files.
	try:
		if start_jd == 1:
			line_num = 1
		else:
			line_num = (start_jd-1)*DIV
		_input = linecache.getline("%s/%s" % (DATA_PATH,file), line_num)	
		#log.write("Line 1: %s" % _input)
# Check to make sure we are not at the end of a file.  If we are, go back and
# find a readable line.
		if not _input:
			while True:
				line_num -= 5
				_input = linecache.getline("%s/%s" % (DATA_PATH,file), line_num)
				#log.write("Searching: %s" % _input)
				if _input:
					break;
		_input = _input.split()
# If the date just read in from the read is later than the date we need, it 
# means data is missing and we need to find the right day.
# It will probably be close so linearly searching backwards is fine
		if float(_input[0]) > file_start_day:
			line_num = line_num-1
			while True:
				_input = linecache.getline("%s/%s" % (DATA_PATH,file), line_num)	
				#log.write("Searching back: %s" % _input)
				if _input == '':
					exit("Error: Unable to read from %s/%s line %i during back "
					     "find." % (DATA_PATH,file,line_num))
				_input = _input.split()
				_input = float(_input[0])
				if _input > file_start_day:
					line_num -= 1
				elif _input < file_start_day:
					line_num += 1
					break
				else:
					break
		elif float(_input[0]) < file_start_day:
			line_num = line_num+1
			while True:
				_input = linecache.getline("%s/%s" % (DATA_PATH,file), line_num)	
				#log.write("Searching forward: %s" % _input)
				if _input == '':
					exit("Error: Unable to read from %s/%s line %i during " 
						 "forward find." % (DATA_PATH,file,line_num))
				_input = _input.split()
				_input = float(_input[0])
				if _input < file_start_day:
					line_num += 1
				elif _input > file_start_day:
					line_num -= 1 
					break
				else:
					break
		else:
			pass
# Now read in values. Number of values to read in should be # of data per day *
# number of days.  Need to make sure that the file_day is not > than the end 
# file day in case data is missing.
		for j in range(DIV*n_days):
			_input = linecache.getline("%s/%s" % (DATA_PATH,file), line_num)	
			#log.write("Line used: %s" % _input)
			if _input == '':
				exit("Error: Unable to read data from %s/%s line %i" % 
				     (DATA_PATH,file,line_num))
			_input = _input.split()
			if float(_input[0]) < file_end_day:
# Append time in seconds with the run start day as 0.
				_time.append((float(_input[0])-file_start_day)*86400)
				_data.append(float(_input[1]))
			else:
				break
			line_num += 1
# Return the file data linearly interpolated to run_times, the run times based
# on dt, for the run.
		return  numpy.interp(run_times, _time, _data)
		linecache.clearcache()
	except IOError:
		exit("Problem reading %s/%s, check if data is available." %			  
			(DATA_PATH, file))
	#log.close()

def genFluxTemp(startDate, endDate, grid):

# Get day numbers in "Julian" days and run_times
	start_jd  = startDate[7]
	end_jd    = endDate[7] 
	n_days    = end_jd - start_jd + 1
	run_times = range(90, (n_days*86400)+1, 90)
	
	if(grid!=15 and grid!=22 and grid!=26):
		sys.exit("Error: wrong grid specified in genFluxTemp")
# Grid 15 (hi-res) ends at Beavery Army... so just treat it the same as 22
	if(grid==15):
		grid = 22

	if(grid==22):
# Beaver Army data for flux
# Longview used for 2008,2009, and 2010 because Beaver Army data is not available
# in /home/workspace/ccalmr/data/external/river_data/verified
		if startDate[0] > 2007:
			beaver_temp = __readData(start_jd, n_days, DATA_PATH,				   
				"longview_temperature.%i" % startDate[0], run_times)
		else: 
			beaver_temp = __readData(start_jd, n_days, DATA_PATH,				   
				"beaveryarmy_temperature.%i" % startDate[0], run_times)
		beaver_flux = __readData(start_jd, n_days, DATA_PATH,				   
			"beaverarmy_discharge.%i" % startDate[0], run_times)

# Write the files
# -flux ensures that the flux is negative to indicate flow into the domain.
		try:
			fluxFile = open("flux.th", "wb")
			for i in range(len(run_times)):
				fluxFile.write("%i\t%f\n" % (run_times[i],-beaver_flux[i]))
		except IOError:
			exit("Error: problem opening flux.th file")
		finally:
			fluxFile.close()

		try:
			tempFile = open("temp.th", "wb")
			for i in range(len(run_times)):
				tempFile.write("%i\t%f\n" % (run_times[i], beaver_temp[i]))
		except IOError:
			exit("Error: problem opening temp.th file")
		finally:
			tempFile.close()

# DB26, multiple for-loops used just for clarity in reading different files
	else:
# Willamette temperature and flux from Newberg
# Unless the year is 2010 and then the values are from Morrison gage height.
# Gage height values are converted to discharge by using the transformation.
# d(h) = 307.1*h^0.9327 
# that I calculated by using the
# discharge and gage height from Morrison for 2009.
# Gage height in ft. and output is in m^3s^-1.
		willam_temp = __readData(start_jd, n_days, DATA_PATH,				   
			"newberg_temperature.%i" % startDate[0], run_times)
		if startDate[0] == 2010:
			willam_flux = __readData(start_jd, n_days, DATA_PATH,
				"morrison_stage.%i" % startDate[0], run_times)		
			for i in range(len(willam_flux)):
				h = willam_flux[i]   
				willam_flux[i] = 307.1*h**0.9327
		else:
			willam_flux = __readData(start_jd, n_days, DATA_PATH,				   
				"newberg_discharge.%i" % startDate[0], run_times)

# Columbia discharge from Bonneville 
		bonne_flux = __readData(start_jd, n_days, DATA_PATH,				  
			"columbia_discharge.%i" % startDate[0], run_times)

# Temp data for Winter is from Warrendale because of close correlation with 
# Bonneville. Temp data for Summer is from Longview
		if(start_jd<=90 or start_jd>=270):
			bonne_temp = __readData(start_jd, n_days, DATA_PATH,			  
				"warrendale_temperature.%i" % startDate[0], run_times)
		else:
			bonne_temp = __readData(start_jd, n_days, DATA_PATH,		
				"longview_temperature.%i" % startDate[0], run_times)

# Write the files
# Flux from rivir must be negative to indicate flow into domain.
# Some flux data is positive, some negative, so it's enforced with -abs(n).
		try:
			fluxFile = open("flux.th", "wb")
			for i in range(len(run_times)):
				fluxFile.write("%i\t%f\t%f\n" % (run_times[i], 
	            	-bonne_flux[i], -willam_flux[i]))
		except IOError:
			exit("Error: problem opening flux.th file")
		finally:
			fluxFile.close()

		try:
			tempFile = open("temp.th", "wb")
			for i in range(len(run_times)):
				tempFile.write("%i\t%f\t%f\n" % (run_times[i], bonne_temp[i], 
					willam_temp[i]))
		except IOError:
			exit("Error: problem opening temp.th file")
		finally:
			tempFile.close()
			

if __name__ == '__main__':
# grab start and number of days 
	if len(sys.argv) != 4:
		sys.exit("usage: %s [begin_date mm-dd-yyyy] [number of days in run] "
				 "[grid (22 or 26)]" % sys.argv[0])
	else:
		start = sys.argv[1]
		nDays = int(sys.argv[2])
		grid  = int(sys.argv[3])

# deal with start and ending days
	(month, day, year) = start.split('-')
	startDay   = int(day)
	startMonth = int(month)
	startYear  = int(year) 

# Python style dates for script.
	startDate = datetime.date(startYear, startMonth, startDay)
	endDate   = startDate + datetime.timedelta(days=(nDays-1))
	startDate = startDate.timetuple()
	endDate   = endDate.timetuple()

	genFluxTemp(startDate, endDate, grid)
