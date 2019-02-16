import urllib
import cmop.db as db
import datetime

if len(sys.argv) != 6:
  sys.exit("usage: python quick_lobo19.py [start yr] [start day] [end yr] [end day] [out file]")

# comma delimited?
DELIM = "\t"

# set up start and end time
endtime = datetime.datetime.now()
delta = datetime.timedelta(hours=24.0)
starttime = endtime - delta
mindate = starttime.strftime("%Y%m%d")
maxdate = endtime.strftime("%Y%m%d")

node = '0019'
varlist = 'salinity,temperature'

url = "http://yaquina.satlantic.com/cgi-data/nph-data.cgi?x=date&y=%s&min_date=%s&max_date=%s&node=%s&data_format=text" % (varlist,mindate,maxdate,node)

file_obj = urllib.urlopen(url)

lines = file_obj.readlines()
for line in lines[3:]:
  line = line.strip('\n') # remove newline
  tuple = line.split(DELIM)
  tuple[0] += ' PST' # time is in PST, convert to explicit PST
  tuple.insert(0,deploymentid)
  tuple.append('PD0')
  linetest = 0
  for i in range(2,len(tuple)-1): # first two tuple values are time and deploymentid, last tuple value is process level
    if tuple[i] == '':
       tuple[i] = 'NULL'
    else:
      if cols[i-2] <> 'voltage':
        linetest += 1
  # convert time format, and do other processing here
  if linetest > 0:
   # print ins % zip(*(zip(tuple)))[0]
    conn.execCommand(ins % zip(*(zip(tuple)))[0])
  else:
    print 'line %s has no valid data' % line
