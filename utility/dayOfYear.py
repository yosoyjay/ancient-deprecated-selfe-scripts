#!/usr/local/bin/python
""" This script prints the day of year given a date as input.
"""
import sys
import datetime

def dayOfYear(date):
	print date[7]

if __name__ == "__main__":
	if len(sys.argv) != 2:
		sys.exit("Usage %s [date mm-dd-yyyy]" % sys.argv[0])
	else:
		start = sys.argv[1]

# deal with start and ending days
	(month, day, year) = start.split('-')
	startDay   = int(day)
	startMonth = int(month)
	startYear  = int(year) 

# Python style dates for script.
	startDate = datetime.date(startYear, startMonth, startDay)
	startDate = startDate.timetuple()

	dayOfYear(startDate)	
