import gamin, threading
import os, os.path, sys, tempfile
import traceback, time, datetime, dateutil.parser
import Queue
import cmop
import cmop.db
import re
import signal
import cmop.seabird
from cmop import *


QSIZE = 2000
killall = Queue.Queue(1)

def Die(msg=None):
  if not msg:
    msg = "Unhandled Exception: %s" % ("".join(traceback.format_exception(*sys.exc_info())),)
  cmop.error(msg)
  killall.put(sys.exc_info())
  raise

def ErrorRaised():
  return killall.full()

class config:
  apppath = '.'

def logerror():
  pass
  #raise
  #traceback.print_exc(file=stdout)

class Source:
  '''a wrapper for a threadsafe queue'''
  def __init__(self, maxsize=QSIZE):
    self.Q = Queue.Queue(maxsize)

  def qsize(self): 
    return self.Q.qsize()

  def Dequeue(self, wait=True):
    return self.Q.get(wait)

  def Enqueue(self, item):
    self.Q.put(item, True)

  def start(self):
    # catch recursive call to start a pipeline
    pass

# A Scanner converts a stream of file 
# append operations to a stream of data values
# Each scanner subclass corresponds to a different 
# data type/file format
class Scanner(Source):
  def __init__(self, maxsize=QSIZE):
    # TODO: make this a persistent dictionary
    self.state = {}
    Source.__init__(self, maxsize)

  def Exception(self, e):
    self.Enqueue(e)

  def OnCreate(self, filename):
    pass

  def OnChange(self, filename):
    pass

  def OnExists(self, filename):
    pass

# Each scanner must implement three methods:
# "match", "OnCreate" and "OnChange." 
# "match" takes a filename and returns a boolean 
# indicating whether or not it knows how to 
# process the file.  
# "OnCreate" and "OnChange" each take a filename and 
# extracts data from the associated file, placing parsed 
# data items on the output queue.
# "OnCreate" corresponds to a file creation event, and 
# "OnChange" corresponds to a file update event

# A TestScanner that matches all files and
# places blocks of uninterpreted data on the 
# output queue

class TestScanner(Scanner):
  def match(self, filename):
    return True

  def OnChange(self, filename):
    pos = self.state.setdefault(filename,0)
    print "TestScanner: ", filename, pos

    if os.path.isdir(filename): return 

    f = file(filename)
    f.seek(pos)
    data = f.read()
    self.state[filename] = pos + len(data)
    lines = data.split()
    f.close()
    for l in lines:
      self.Enqueue(l)

  def OnCreate(self, filename):
    self.Enqueue(filename)

  def OnExists(self, filename):
    self.OnCreate(filename)
    self.OnChange(filename)

def dospath(path):
   return path.replace('/', '\\')

def dosabspath(path, drive='C'):
   ap = os.path.abspath(path)
   return "%s:%s" % (drive, dospath)

class PassThroughScanner(Scanner):
  def __init__(self, condition, maxsize=QSIZE):
    self.condition = condition
    Scanner.__init__(self, maxsize)

  def match(self, filename):
    return self.condition(filename)

  def OnCreate(self, filename):
    self.Enqueue(filename)

  def OnChange(self, filename):   
    pass
    #self.Enqueue(filename)

  def OnExists(self, filename):   
    self.Enqueue(filename)

class RecordStream(Scanner):
  '''scanner for files consisting of 
a header followed by a stream of records.
Implement 
ReadHeader(fileobj) and 
ReadBlock(header, fileobj)

ReadHeader accepts an open file and returns the header,
leaving the file pointer just after the header (if any).

ReadBlock reads the remainder of the file and yields 
records. (ReadBlock is a generator)

constructor accepts a function condition :: filename -> Bool
that determines which files will be scanned
'''

  def __init__(self, condition, maxsize=QSIZE):
    self.condition = condition
    Scanner.__init__(self, maxsize)

  def match(self, filename):
    return self.condition(filename)
  
  def ReadBlock(self, f):
    raise TypeError("ReadBlock not implemented; class abstract")

  def ReadHeader(self, f):
    raise TypeError("ReadHeader not implemented; class abstract")

  def Scan(self, filename):
    f = file(filename)
    hdr = self.ReadHeader(f)
    pos, item = self.state.get(filename, (0,0))

    cmop.debug("Scanning %s as RecordStream" % (filename,), 8)

    if pos == 0:
      pos = f.tell()

    opos, oitem = pos, item

    f.seek(pos)
    for x in self.ReadBlock(hdr, f):
      x['row'] = item
      item += 1
      # add a column for filename
      x['file'] = filename
      self.Enqueue(x)
  
    pos = f.tell()
    if item != oitem:
      cmop.info("Scanned %s records (%s bytes) from %s" % (item-oitem, pos-opos, filename))

    self.state[filename] = (pos, item)

  def OnCreate(self, filename):
    self.Scan(filename)

  def OnChange(self, filename):   
    self.Scan(filename)

  def OnExists(self, filename):   
    self.Scan(filename)

class DelimitedRecordStreamWithHeader(RecordStream):
  '''Parse a record stream consisting of a one-line
delimited header followed by N one-line delimited records'''

  def __init__(self, condition, delimiter=",", maxsize=QSIZE):
    self.delimiter = delimiter
    RecordStream.__init__(self, condition, maxsize)

  def ReadHeader(self, f):
    hdr = f.readline()
    return hdr.split(self.delimiter)

  def ReadBlock(self, header, f):
    for line in f:
      t = tuple(line.split(self.delimiter))
      assert(len(t) == len(header))
      d = dict(zip(header,t))
      yield d

class DelimitedRecordStreamNoHeader(DelimitedRecordStreamWithHeader):
  '''Parse a record stream consisting of N 
one-line delimited records.  Constructor requires
a list of field names.'''

  def __init__(self, condition, header, delimiter=",", maxsize=QSIZE): 
    DelimitedRecordStreamWithHeader.__init__(self, condition, delimiter, maxsize=QSIZE)
    self.header = header
   
  def ReadHeader(self, f):
    return self.header

class Filter(Source, threading.Thread):
  ''' Consume data from a queue in a separate thread '''

  def __init__(self, source, maxsize=QSIZE):
    ''' save the input queue'''
    Source.__init__(self, maxsize)
    threading.Thread.__init__(self)
    self.setDaemon(1)
    self.source = source 

  def Sizes(self):
    return self.source.qsize(), self.qsize()

  def start(self):
    # recursively start my source
    self.source.start()
    threading.Thread.start(self)

  def RunProcess(self, data):
   
    if ErrorRaised(): sys.exit()

    # unnest a sequence of results
    # unless the result is a string
    # also capture some timing information
    start = time.time()
    try: 
      for r in self.process(data):
        self.Enqueue(r)
    except KeyboardInterrupt:
      Die()
    except: 
      msg = "Ingest Error: %s" % ("".join(traceback.format_exception(*sys.exc_info())),)
      cmop.error(msg)
      #Die()  # signals all threads to stop.  useful for debugging

    end = time.time()

    msg = str(self.__class__.__name__)
    msg += ": %s seconds, %s inbound, %s outbound"
    i, o = self.Sizes()
    msg = msg % (end - start, i, o)
    profile(msg)

  def run(self):
    ''' overload the run method of threading.Thread'''
    while True:
      data = self.source.Dequeue()
      self.RunProcess(data)

class CoalesceFilter(Filter):
  ''' 
Filter that groups redundant messages
found in its input queue.  Useful for reducing 
wasted work if messages are idempotent and each 
takes a lot of time to process.

Use the sleeptime argument to force a delay between calls
'''

  def __init__(self, source, sleeptime=30, maxsize=QSIZE):
    Filter.__init__(self, source, maxsize)
    self.sleeptime=sleeptime

  def run(self):

    uniquemessages = set([])

    while True:
      try:
        # consume an item, but throw an exception if queue is empty
        data = self.source.Dequeue(False)
        uniquemessages.add(data)

      except Queue.Empty:
        cmop.debug("Coalesced %s messages" % (len(uniquemessages),), 8)
        for m in uniquemessages:
          self.RunProcess(m)

        uniquemessages.clear()
        time.sleep(self.sleeptime)

  def process(self, data): 
    # do nothing by default; subclasses may override
    return [data]

class FilePunctuatedRecordStream(RecordStream):

  def MakeGroup(self, filename):
     data = {}
     data['group'] = filename
     self.Enqueue(data)
     return data

  def OnCreate(self, filename):
     self.MakeGroup(filename)
     RecordStream.OnCreate(self, filename)

  def OnExists(self, filename):
     self.MakeGroup(filename)
     RecordStream.OnExists(self, filename)

class CMOPMetadata(Filter):
  def __init__(self, vessel, cruise, instrument, source, instrumenttype='CTD', interval=0.25, maxsize=QSIZE):
    self.cruise = cruise
    self.vessel = vessel
    self.instrument = instrument
    self.instrumenttype = instrumenttype
    self.interval = interval
    Filter.__init__(self, source, maxsize)

  def AttachMetadata(self, data):
    data['cruise'] = self.cruise
    data['vessel'] = self.vessel

class ForerunnerCastScanner(DelimitedRecordStreamNoHeader, FilePunctuatedRecordStream):
  # mixed-in class 
  OnCreate = FilePunctuatedRecordStream.OnCreate
  OnExists = FilePunctuatedRecordStream.OnExists

class DropFilter(Filter):
  def process(self, data):
    return []

class WecomaTSGCleaner(CMOPMetadata):

  def GetLatLon(self, tup, data): 
      if data.has_key('pcode_lat_deg') and data.has_key('pcode_long_deg') \
         and data.has_key('pcode_lat_min') and data.has_key('pcode_long_min'):
         latdeg = data['pcode_lat_deg']
         londeg = data['pcode_long_deg']
         latmin = data['pcode_lat_min']
         lonmin = data['pcode_long_min']
         lat = float(latdeg) + float(latmin)/60
         lon = float(londeg) + float(lonmin)/60

      elif data.has_key('pcode_decimal_latitude') and data.has_key('pcode_decimal_longitude'):
         lat, lon = data['pcode_decimal_latitude'], data['pcode_decimal_longitude']
         lat, lon = float(lat), float(lon)
 
      else: 
         raise ValueError("No suitable Lat/Lon values found (%s)" % (data,))

      if lat >= 80 or lat <= 40 or lon > -100 or lon < -130:
        raise ValueError("Lat/Lon out of range (%s, %s)" % (lat, lon))

      tup['latitude']  = lat        
      tup['longitude']  = lon
 
  def MakeTime(self, yearday, hour, min, sec=0):
     year = time.localtime()[0]
     return "('12/31/%s'::timestamp + '%s days'::interval + '%s:%s:%s'::time || 'UTC')::timestamptz" % (int(year)-1, int(yearday), int(hour), int(min), int(sec))

  def GetTime(self, tup, data): 
     if data.has_key('truetime'):
       timerec = data['truetime'].split(':')
       timerec[3] = re.sub("[^0-9]", "", timerec[3])
       time = self.MakeTime(*tuple(timerec))
     elif data.has_key('truetime_day') \
      and data.has_key('truetime_hour') \
      and data.has_key('truetime_minute'):
       day = data['truetime_day']
       hour = data['truetime_hour']
       mint = data['truetime_minute']
       time = self.MakeTime(day, hour, mint)
     else:
       raise ValueError("No suitable time atteibutes found: %s" % (data,))

     tup['time'] = time

  def process(self, data):
    try:
      tuple = {}

      self.GetLatLon(tuple, data)
      self.GetTime(tuple, data)

      tuple['table'] = "cruise.tsg"
      tuple['vessel'] = self.vessel
      tuple['cruise'] = self.cruise
      tuple['instrument'] = self.instrument
      tuple['instrumenttype'] = self.instrumenttype
      tuple['salinity'] = float(data['computed_salinity_flothru'])
      #tuple['temperature'] = float(data['water_temp_seabird_flothru'])
      tuple['temperature'] = float(data['surface_water_temp_seabird'])
      tuple['conductivity'] = float(data['conductivity_seabird_flothru'])
      tuple['winddirection'] = float(data.get('wind_heading_ultrasonic_true', -999999))
      tuple['windspeed'] = float(data.get('wind_speed_ultrasonic_true(knots)', -999999))
      tuple['atmosphericpressure'] = float(data['barometer'])
      tuple['atmospherictemperature'] = float(data['air_temp_rmyoung_doghouse'])

      for k, v in tuple.items():
        if str(v) == 'nan':
          tuple[k] = None

      cmop.debug("WecomaTSGCleaner generated a tuple: %s" % (tuple,), 8)
      return [tuple]

    except:
      msg = traceback.format_exc()
      s = ",".join(["(%s=%s)" % (k,v) for k,v in data.iteritems()])
      cmop.info("%s : Skipping bad TSG record: %s" % (msg,s))
      return []

class BarnesTSGCleaner(Filter):

  def __init__(self, cruise, filter, instrument='barnestsg', vessel='Barnes'):
    Filter.__init__(self, filter)
    self.cruise = cruise
    self.instrument = instrument
    self.instrumenttype = 'CTD'
    self.vessel = vessel
    self.table = 'cruise.tsg'

  def process(self, data):
    try:
      tuple = {}
      lat, lon = data['GPS-LAT'], data['GPS-LON']
      tuple['latitude'] = float("%s" % (float(lat[0:2]) + float(lat[2:-1])/60,))
      tuple['longitude'] = float("%s" % (float(lon[0:3]) + float(lon[3:-1])/60,))
      tuple['time'] = "%s %s" % (data['Date'], data['Time'])
      tuple['salinity'] = float(data['SBE21-Sal'])
      tuple['temperature'] = float(data['SBE21-Temp-ITS-90'])
      tuple['conductivity'] = float(data['SBE21-Cond'])
      tuple['cruise'] = self.cruise
      tuple['instrument'] = self.instrument
      tuple['instrumenttype'] = self.instrumenttype
      tuple['vessel'] = self.vessel
      tuple['table'] = self.table
      return [tuple]
    except:
      s = ",".join(["(%s=%s)" % (k,v) for k,v in data.iteritems()])
      cmop.error("Bad TSG data (%s): %s" % ("".join(traceback.format_exception(*sys.exc_info())), s))
      return []

class DBLookup(Filter):
  '''Base class for filter using a database connection'''
  def __init__(self, filter, host=cmop.config.hostname, user=cmop.config.user, password=cmop.config.password, maxsize=QSIZE):
    Filter.__init__(self, filter, maxsize)
    self.db = cmop.db.DB(host=host, user=user, password=password)

class ForerunnerObservation(Filter):
  def __init__(self, cruise, source):
    self.cruise = cruise
    self.vessel = 'Forerunner'
    Filter.__init__(self, source)

class ForerunnerTSG(ForerunnerObservation):
  def process(self, data):
    pass

class CastCleaner(CMOPMetadata):

  def makeCast(self, data):
     raise TypeError("not implemented; class abstract")

  def makeObservation(self, data):
     raise TypeError("not implemented; class abstract")

  def process(self, data):
    if data.has_key('group'):
      return [self.makeCast(data)]
    else:
      return [self.makeObservation(data)]

class KnownCast(CoalesceFilter, DBLookup):
  ''' filter out known casts'''
  def __init__(self, source, groupseconds=5, maxsize=QSIZE, host=cmop.config.hostname):
    DBLookup.__init__(self, source, host) 
    CoalesceFilter.__init__(self, source, groupseconds, maxsize) 
 
  def process(self, fname):
    base, name = os.path.split(fname)
    root, ext = os.path.splitext(name)
    sel = '''SELECT * FROM cruise.ctdcast WHERE file like '%%%s%%' ''' % (root,)
    rs = self.db.execQuery(sel)
    if rs:
      return []
    else:
      return [fname]


class SeabirdCastCleaner(CoalesceFilter, CastCleaner):
  ''' 
Construct a Seabird Cast object and extract the relevant data.
''' 

  def __init__(self,  vessel, cruise, instrument, source, mapping, instrumenttype='CTD', interval=-1, groupseconds=5, timezone=time.tzname[1]):
    self.mapping = mapping
    CastCleaner.__init__(self,  vessel, cruise, instrument, source, instrumenttype, interval)
    CoalesceFilter.__init__(self, source, groupseconds, maxsize=20000)
    self.timezone=timezone

  def process(self, fname):
    cmop.debug("Seabird file found: %s" % (fname,))
    cast = seabird.SeabirdCast(fname) 
    cast.GetConverted()
    cast.GetBottle()
     
    # uncomment for speedier testing
    #f = file('sample.cnv')
    #cast.datcnvoutput = f.read()
    cast.Parse()
    cast.ParseBottle()
    record  = self.makeCast(cast)
    
    cmop.debug("Parsed Seabird File %s" % (fname,))
   # yield {"sqlcommand":"begin"}
    try:
      yield record
      for o in cast.IterateBottle(self.mapping['depth']):
        yield self.makeBottle(record,o)
      for o in cast.Iterate():
        yield self.makeObservation(record, o)
    
    #  yield {"sqlcommand":"commit"}
    except: 
      pass
     # yield {"sqlcommand":"rollback"}
    
    cmop.info("Processed %s observations for file %s" % (len(cast.rows),fname))
 
  def GetInterval(self, seacast):
    '''Compute the interval from a seabrd cast and an observation record, or return a default value'''
    if hasattr(seacast, 'interval'):
      # seabird emits "units: value"
      interval = seacast.interval.split(':')
      if interval[0] == 'seconds':
        return float(interval[1])
      #elif milliseconds, etc?

    return self.interval
   
  def makeBottle(self, record, b):
    '''
Consumes bottle tuples as generated by SeabirdCast objects
(see cmop/seabird.py)
'''
    bottle = b
    self.AttachMetadata(bottle)
    bottlefields = ('bottle','firesequence','time','startscan','endscan','depth')     

    bottle['castid'] = record['castid']

    bottle['time'] = "'%s %s'::timestamptz" % (bottle['time'], self.timezone) 
    bottle['table'] = 'cruise.sample'
    return bottle
   
  def makeCast(self, data):
    '''
Consumes cast tuples as generated by SeabirdCast objects
(see cmop/seabird.py)
'''
    cast = {}
    self.AttachMetadata(cast)

    seacast = data

    self.GetTime(seacast, cast)

    self.GetLocation(seacast, cast)

    root, fname = os.path.split(seacast.dat)
    base, ext = os.path.splitext(fname)
    cast['castid'] = base

    cast['file'] = seacast.dat

    cast['instrument'] = self.instrument
    cast['instrumenttype'] = self.instrumenttype
    cast['interval'] = self.GetInterval(seacast)

    if hasattr(seacast, 'station'):
      cast['site'] = seacast.station

    if hasattr(seacast, 'beddepth'):
      cast['beddepth'] = seacast.beddepth

    cast['table'] = 'cruise.ctdcast'
    return cast

  def GetObservationTime(self, castrecord, o):

    timecalc = "'%s'::timestamptz + '%s seconds'"
    if o.has_key('time'):
      elapsed = float(o['time'])
    elif o.has_key('scan: Scan Count') and castrecord.has_key('interval'):
      elapsed = int(o['scan: Scan Count']) * float(castrecord['interval'])
    else:
      raise ValueError("Cannot compute time for observation.  Require either elpased time in seconds or scan number and interval. Cast: %s \n Observation: \n %s" % (castrecord, o))

    return timecalc % (castrecord['time'], elapsed)

  def makeObservation(self, castrecord, o):
    obs = {}

    obs['time'] = self.GetObservationTime(castrecord, o)
 
    obs['castid'] = castrecord['castid']

    def setval(n1, n2, conv=float):
      
      matches = [k for k in o if n2.lower() in k.lower()]
      
      if matches:
        v = o.get(matches[0], None)
        if v: 
          obs[n1] = conv(v)
     
    for k,v in self.mapping.items():
      setval(k,v)
    
    obs['table'] = 'cruise.castobservation'

    self.AttachMetadata(obs)
    return obs

  def GetLocation(self, seacast, cast):
    if hasattr(seacast,'lat') and hasattr(seacast, 'lon'):
      cast['latitude'] = cmop.DM2DD(seacast.lat)
      cast['longitude'] = cmop.DM2DD(seacast.lon)
    else:
      # Need to look up most recent position from TSG stream
      # Doing this with the database for now;
      # assumes TSG stream is functioning and current
      pass     

  def GetTime(self, seacast, cast):
    ''' 
Deduce the time of the cast. 
Set the 'time' key on the cast dictionary to an 
expression that can be evaluated by Postgresql
'''

    if hasattr(seacast, 'NMEAtime'):
         # System UpLoad Time is in local time and comes
         # from the system clock of the machine running
         # seasave
         # NMEA time comes from the device, buit has no date 
         # information.  
         # must assume a bound on clock skew to guess the right date
         test_time1 = dateutil.parser.parse(seacast.NMEAtime,default=datetime.datetime(year=1,month=2,day=3))
         test_time5 = dateutil.parser.parse(seacast.NMEAtime,default=datetime.datetime(year=5,month=1,day=1))
	 test = test_time1-test_time5 
   	 cmop.debug('NMEAtime %s' % seacast.NMEAtime)
         if test == datetime.timedelta(0):  # NMEAtime parses as a full date/time
            cast['time'] = "%s %s" % (seacast.NMEAtime, 'UTC')
            cmop.debug('full NMEA %s' % cast['time'])
            return 
         if hasattr(seacast, 'start_time'):
            cmop.debug('start time %s' % seacast.start_time)
            start_time = dateutil.parser.parse(seacast.start_time)
            NMEAtime = dateutil.parser.parse(seacast.NMEAtime)
            comb_time = datetime.datetime.combine(start_time.date(),NMEAtime.time())
            if start_time > comb_time + datetime.timedelta(0,43200):
		comb_time += datetime.timedelta(1)
            elif  comb_time > start_time + datetime.timedelta(0,43200):
                comb_time -= datetime.timedelta(1)
            cast['time'] = "%s %s" % (comb_time, 'UTC')
            cmop.debug('NMEA start %s' % cast['time'])
            return 
         if hasattr(seacast, 'uploadtime'):
            cmop.debug('upload time %s' % seacast.uploadtime)
            uploadtime = dateutil.parser.parse(seacast.uploadtime)
            NMEAtime = dateutil.parser.parse(seacast.NMEAtime)
            comb_time = datetime.datetime.combine(uploadtime.date(),NMEAtime.time())
            if uploadtime > comb_time + datetime.timedelta(0,43200):
		comb_time += datetime.timedelta(1)
            elif  comb_time > uploadtime + datetime.timedelta(0,43200):
                comb_time -= datetime.timedelta(1)
            cast['time'] = "%s %s" % (comb_time, 'UTC')
            cmop.debug('NMEA upload %s' % cast['time'])
            return 
    if hasattr(seacast, 'GMTtime') and hasattr(seacast, 'GMTdate'):
         time = seacast.GMTtime
         time.replace(' ','')
         cmop.debug('trying to merge %s and %s' % (seacast.GMTdate,seacast.GMTtime))
         if re.search('^\d{4}',seacast.GMTtime.strip()):
           time = datetime.time(int(time[0:2]), int(time[2:4]))
           date = seacast.GMTdate
           try:
              date = int(date)
              date = datetime.date.fromordinal(date)
           except:
              date = dateutil.parser.parse(date)
           try: 
              cmop.debug('preparing to merge %s and %s' % (date,time))
              date = datetime.datetime.combine(date,time)
              cmop.debug('converted %s to %s' % (seacast.GMTdate,date))
              if hasattr(seacast, 'start_time'):
                cmop.debug('start time %s' % seacast.start_time)
                altdate = dateutil.parser.parse(seacast.start_time)
              elif hasattr(seacast, 'uploadtime'):
                cmop.debug('upload time %s' % seacast.uploadtime)
                altdate = dateutil.parser.parse(seacast.uploadtime)
              date = date.replace(year=altdate.year)
              if altdate > date + datetime.timedelta(100):
	        date.replace(year=date.year-1)
              elif date > altdate + datetime.timedelta(100):
	        date.replace(year=date.year+1)
              cast['time'] = "%s %s" % (date, 'UTC')
              cmop.debug('GMT time %s' % cast['time'])
              return 
           except:
              raise TypeError('could not convert %s' % seacast.GMTdate)
    if hasattr(seacast, 'uploadtime'):
         cast['time'] = "%s %s" % (seacast.uploadtime, self.timezone)
         cmop.debug('upload time %s' % cast['time'])
         return
    if hasattr(seacast, 'start_time'):
         cast['time'] = "%s %s" % (seacast.start_time, self.timezone)
         cmop.debug('start time %s' % cast['time'])
         return
    raise ValueError("Cannot calculate cast date: Need either GMTdate and GMTtime, or System Upload Time and NMEA Time.  \n%s" % (seacast.__dict__.keys(),))
 
class Select(Filter):
  def __init__(self, condition, source):
    Filter.__init__(self, source)
    self.condition = condition

  def process(self, data):
    if self.condition(data):
      return [data]
    else:
      return []

class Transform(Filter):
  def __init__(self, transform, source):
    Filter.__init__(self, source)
    self.transform = transform

  def process(self, data):
    return self.transform(data)


class JoinGPSStream(DBLookup):
  def __init__(self, gpsinstrument, source,  gpsinstrumenttype= 'CTD', host=cmop.config.hostname):
    DBLookup.__init__(self, source, host=host)
    self.gpsinstrument = gpsinstrument
    self.gpsinstrumenttype = gpsinstrumenttype

  def process(self, data):
    if not data.has_key('table'):
      return [data]
    elif data['table'] != 'cruise.ctdcast':
      return [data]
    else:
      sel = '''
SELECT location
  FROM cruise.tsg t, cruise.deployment d
 WHERE d.vessel = '%s'
   AND d.cruise = '%s' 
   AND d.instrument = '%s' 
   AND d.instrumenttype = '%s' 
   AND d.deploymentid = t.deploymentid
   AND t.time < '%s'
ORDER BY t.time DESC LIMIT 1
''' % (
  data['vessel'], 
  data['cruise'], 
  self.gpsinstrument, 
  self.gpsinstrumenttype,
  data['time']
)
      rs = self.db.execQuery(sel) 
      if len(rs) == 0:
        cmop.warn("No suitable TSG records found with %s" % (sel,))
      else:
        data['location'] = rs[0][0]

      return [data]

class ForerunnerCastCleaner(CastCleaner):
  def __init__(self, cruise, source):
    CastCleaner.__init__(
      self, 
      vessel='Forerunner',
      instrument='castforerunner',
      instrumenttype='CTD',
      interval=0.25,
      source=source
    )

  def LookupLocation(self, path, time):
    loc = file(path)
    lines = loc.readlines()
    for l in lines:
      flds = l.split()
      if flds[3] == time:
        return tuple(flds[0:2])

  def parsePath(self, fn):
    base, name = os.path.split(fn)
    name, ext = os.path.splitext(name)
    sday, stime = name.split('_')
    return base, name, sday, stime

  def getCorieday(self, sday, stime):
    day = float(sday[2:])
    hr, mn = float(stime[0:2]), float(stime[2:4])
    return day + hr/24.0 + mn/60.0/24.0 

  def makeCast(self, data):
    cast = {}
    self.AttachMetadata(cast)

    fn = data['group']
    base, castid, sday, stime = self.parsePath(fn)

    cast['corieday'] = self.getCorieday(sday, stime)
    cast['castid'] = castid

    locpath = os.path.join(base, sday + 'a.LOC')
    lat, lon = self.LookupLocation(locpath, stime[0:4]) 

    cast['latitude'] = float(lat)
    cast['longitude'] = float(lon)

    cast['instrument'] = self.instrument
    cast['instrumenttype'] = self.instrumenttype

    # TODO: get the site id from the user, or the database
    # unclear how to wire this into the application
    cast['site'] = "%3.3f/%2.3f" % (-float(lon), float(lat))

    cast['table'] = 'cruise.ctdcast'

    return cast

  def makeObservation(self, data):
    obs = {}

    fn = data['file']
    base, castid, sday, stime = self.parsePath(fn)
    corieday = self.getCorieday(sday, stime)
    obs['corieday'] = corieday + float(data['row']) * self.interval/86400
    obs['castid'] = castid
    obs['file'] = data['file']
    obs['depth'] = float(data['pressure'])
    obs['salinity'] = float(data['salinity'])
    obs['turbidity'] = float(data['turbidity'])
    obs['temperature'] = float(data['temperature'])

    obs['table'] = 'cruise.castobservation'
    self.AttachMetadata(obs)
    obs.pop('file')
    return obs

from pycodas import codas
import datetime
class CODASScanner(Scanner):

  def __init__(self, maxsize=QSIZE):
    self.expr = re.compile('(.*/)?(.*)(\d\d\d)\.(blk)|(BLK)')
    self.opendbs = {}
    Scanner.__init__(self, maxsize=QSIZE)

  def match(self, filename): 
    if self.expr.match(filename):
      return True
    else:
      return False

  def OnCreate(self, filename):
    dbname, db, id = self.ParseFilename(filename)
    db, enddate = self.OpenDB(dbname)
    self.ExtractConfiguration(db)
    self.ExtractPingsSince(db, enddate)

  def OnChange(self, filename):
    dbname, db, id = self.ParseFilename(filename)
    db, enddate = self.OpenDB(dbname)
    # don't re-extract a configuration
    self.ExtractPingsSince(db, enddate)

  def OnExists(self, filename):
    dbname, db, id = self.ParseFilename(filename)
    db, enddate = self.OpenDB(dbname)
    self.ExtractConfiguration(db)
    self.ExtractPingsSince(db, enddate)

  def ParseFilename(self, filename):
    m = self.expr.match(filename)
    path = m.group(1)
    db = m.group(2)
    id = m.group(3)
    dbpath = os.path.join(path, db)
    return  dbpath, db, id

  def OpenDB(self, dbname):
    # CODAS DB closes itself on delete
    if not dbname in self.opendbs: 
      db = codas.DB(dbname)
      self.opendbs[dbname] = db, db.ymdhms_start

    return self.opendbs[dbname]
  
  def ExtractConfiguration(self, db):
    return []
    # ought to get blanking distance
    # and such out of the configuration
    #for k in profs:
    #  get_structure_def(k)
    #profs = db.get_profiles()
    #print profs
    #sys.exit()

  def GetInstrument(self, db):
    f = db.dbname
    dirs = f.split('/')
    for i, c in enumerate(dirs):
      if c == 'adcpdb':
        return dirs[i-1]
    raise ValueError("Can't extract instrument name form dbname: %s" % (f,))

  def ExtractPingsSince(self, db, since):
    FILL = -32768
 
    now = datetime.datetime.now() 
    profs = db.get_profiles((since, now))
    n = profs['ymdhms'].shape[1]

    instr = self.GetInstrument(db)

    conf = {}

    atime = profs['ymdhms']
    alat = profs['lat_dir'].filled(FILL)
    alon = (profs['lon_dir'] - 360).filled(FILL)
    au = (profs['umeas'] + profs['uship']).filled(FILL)
    av = (profs['vmeas'] + profs['vship']).filled(FILL)
    aw = profs['w'].filled(FILL)
    ae = profs['e'].filled(FILL)
    adepth = profs['depth'].filled(FILL)
 
    ping = {}
    obs = {}
    for i in range(n):
      time = "%s/%s/%s %s:%s:%s" % tuple(atime[:,i])
      depth = adepth[:,i]
      u = au[:,i]
      v = av[:,i]
      w = aw[:,i]
      e = ae[:,i]

      lat = alat[i]
      lon = alon[i]

      depth = adepth[:,i]
      u = au[:,i]
      v = av[:,i]
      w = aw[:,i]
      e = ae[:,i]
      lat = alat[i]
      lon = alon[i]
      ping['table'] = 'cruise.adcpping'
      ping['file'] = db.dbname
      ping['time'] = time
      ping['instrument'] = instr
      ping['ncells'] = len(depth)
      ping['latitude'] = lat
      ping['longitude'] = lon
      self.Enqueue(ping.copy())

      for j in xrange(len(depth)):
        obs['table'] = 'cruise.adcpobservation'
        obs['time'] = time
        obs['depth'] = float(depth[j])
        obs['file'] = db.dbname
        obs['bin'] = j
        obs['latitude'] = lat
        obs['longitude'] = lon
        obs['ncells'] = len(depth)
        obs['v1'] = int(u[j]*1000)
        obs['v2'] = int(v[j]*1000)
        obs['v3'] = int(w[j]*1000)
        obs['v4'] = int(e[j]*1000)
        obs['instrument'] = instr
        self.Enqueue(obs.copy())

     #self.opendbs[dbname] = db, now

class CODASCleaner(CMOPMetadata):

  def process(self,data):
    self.AttachMetadata(data)
    data['instrumenttype'] = self.instrumenttype
    
    return [data]     

class SQLCommand(DBLookup):
  def process(self, data):
    if 'sqlcommand' in data:
      cmd = data['sqlcommand']
      self.db.execCommand(cmd)
      return []
    else:
      return [data]

class TableInsert(DBLookup):
  '''Input is tuples represented as pythons dictionaries.
This Filter loads them into the database specified by the global python config file.
'''
  
  expr = re.compile("\s*('.*')|(\(.*\))\s*")

  def clean(self, data):
    tup = {}
    tup.update(data)
    return tup
 
  def process(self, data):
    tup = {}

    if 'sqlcommand' in data:
      cmd = data['sqlcommand']
      self.db.execCommand(cmd)
      return []

    table = data.pop('table')

    for k,v in data.items():
      tup[k] = db.quote(v)

    tup = self.clean(tup)
    try: 
      yes = self.db.InsertTuple(table, tup.keys(), tup.values())
      cmop.debug("Inserted Tuple: %s" % (tup,), 9)
    except cmop.db.DuplicateError:
      cmop.debug("Duplicate record: %s" % (tup,), 8)


    return [tup]

class CMOPTableInsert(TableInsert):
  '''CMOP-specific tuple cleaning'''
  def clean(self, data):

    tup = TableInsert.clean(self, data)

    # convert lat lon to postgis location
    if tup.has_key('latitude') and tup.has_key('longitude'): 
      lat = tup.pop('latitude')
      lon = tup.pop('longitude')
      if lon > 0: lon = -lon
      tup['location'] = 'SetSRID(MakePoint(%s, %s), 4326)' % (lon, lat) 

    # convert corieday to time
    if tup.has_key('corieday'):
      corieday = tup.pop('corieday')
      tup['time'] = 'daycorie(%s)' % (corieday,)

    return tup

class LookupDeploymentID(DBLookup):
  '''Replace cruise, vessel, instrument, and instrumenttype with deploymentid'''

  def process(self, tup): 
    data = {}
    data.update(tup)
    vessel = data.pop('vessel')
    cruise = data.pop('cruise')
    instrument = data.pop('instrument')
    instrumenttype = data.pop('instrumenttype')

    sql = '''
SELECT deploymentid 
  FROM cruise.deployment 
 WHERE cruise = '%s'
   AND vessel = '%s'
   AND instrument = '%s'
   AND instrumenttype = '%s'
''' % (cruise, vessel, instrument, instrumenttype)
    

    result = self.db.execQuery(sql)
    if result:
      data['deploymentid'] = result[0][0]
      return [data]
    else: 
      ins = '''
INSERT INTO integrated.instrument (instrument, instrumenttype) 
VALUES ('%s', '%s')
''' % (instrument, instrumenttype)

      try:
        self.db.execCommand(ins)
      except db.DuplicateError:
        pass

      ins = '''
INSERT INTO cruise.deployment (vessel, cruise, instrument, instrumenttype, depth)
VALUES ('%s', '%s', '%s', '%s', %s)
''' % (vessel, cruise, instrument, instrumenttype, -1)
      self.db.execCommand(ins)

      result = self.db.execQuery(ins)
      if result: 
        data['deploymentid'] = result[0][0]
        return [data]
      else: 
        raise db.SQLError("Can't find record just inserted! %s" % (ins,))

class LookupADCPConfigurationID(DBLookup):
  '''Replace deploymentid, file with adcpdeploymentid, 
inserting a record if necessary.
'''

  def process(self, tup):
    data = {}
    data.update(tup)
    filename = data.pop('file')
    depid = data.pop('deploymentid')
    ncells = data.pop('ncells')

    sql = '''
SELECT adcpdeploymentid 
  FROM cruise.adcpconfiguration 
 WHERE deploymentid = %s
   AND filename = %s
''' % (depid, db.quote(filename))

    result = self.db.execQuery(sql)
    if result:
      data['adcpdeploymentid'] = result[0][0]
      return [data]

    else:
      ins = '''
INSERT INTO cruise.adcpconfiguration (deploymentid, filename, ncells)
VALUES (%s, %s, %s)
''' % (depid, db.quote(filename), ncells)
      self.db.execCommand(ins)

      result = self.db.execQuery(sql)
      if result:
        data['adcpdeploymentid'] = result[0][0]
        return [data]
      else:
        raise db.SQLError("Can't find record just inserted: %s" % (ins,))

class Log(Filter):
  '''log item to a file in a csv format'''

  def __init__(self, filter, fileprefix=config.apppath):
    Filter.__init__(self, filter)
    filename = os.path.abspath(time.strftime(fileprefix + "%y%m%d_%H%M%S" + ".log"))

    if os.access(filename, os.F_OK):
      self.empty = False
    else:
      self.empty = True

    self.f = file(filename, 'w+')

  def process(self, data):
    if self.empty:
      self.f.write(",".join([str(k) for k in data.keys()]) + "\n")
      self.empty = False
    self.f.write(",".join([str(v) for v in data.values()]) + "\n")
    return [data]

class TestFilter(Filter):
  ''' a test filter that just prints the data'''
  def process(self, data):
    print "DATA (TestFilter): ", data
    return [data]
    self.EnQueue(data)

# Monitor class for watching a directory for changes 
# and dispatching the events to the appropriate scanners
class IngestMonitor:
  def __init__(self):
    self.mon = gamin.WatchMonitor()
    self.scanners = []
    self.dirs = {}

  def RecursiveWatchDirectory(self, directory):
    self.WatchDirectory(directory)
    for root, dirs, files in os.walk(directory):
      for d in dirs:
        self.WatchDirectory(os.path.join(root, d))

  def WatchDirectory(self, directory):
    if not directory in self.dirs:
      directory = os.path.abspath(directory)
      callback = self.MakeDispatch(directory)
      watch = self.mon.watch_directory(directory, callback) 
      self.dirs[directory] = watch

    return self.dirs[directory]

  def AddScanner(self, scanner):
    self.scanners.append(scanner)

  def MakeDispatch(self, directory):
    debug("Making Dispatch function for directory %s" % (directory,), level=5)
    def dispatch(filename, event):
      path = os.path.join(directory, filename)
      debug("Dispatch: %s on '%s'" % (event, path), level=5)
      if os.path.isdir(path):
        self.WatchDirectory(path)
      
      for s in self.scanners:
        try:
          if s.match(path): 
            debug("Match: %s on '%s'" % (event, path), level=8)
            if event == gamin.GAMChanged:
              s.OnChange(path)
            elif event == gamin.GAMCreated:
              s.OnCreate(path)
            elif event == gamin.GAMExists:
              s.OnExists(path)
        except:
          #logerror() 
          cmop.error("Error Dispatching Event: %s" % ("".join(traceback.format_exception(*sys.exc_info()))))
          #Die()

    return dispatch

  def Watch(self):
    time.sleep(1)
    while True:
      self.mon.handle_one_event()

class PollingIngestMonitor(IngestMonitor):
  def __init__(self, delay=3):
    IngestMonitor.__init__(self)
    self.pollingdelay = delay

  def WatchNew(self):
    time.sleep(1)
    currdirs = self.dirs.copy()
    cmop.debug("Watching new directory: %s" % (currdirs,))

    def PhonyEvent(event):
      try:
        for d in currdirs: 
          dispatch = self.MakeDispatch(d)
          for f in os.listdir(d):
            path = os.path.join(d,f)
            dispatch(path, event)
      except:
        cmop.error("Error watching new directories: %s" % ("".join(traceback.format_exception(*sys.exc_info()))))

    while True: 
      PhonyEvent(gamin.GAMExists)
      cmop.debug("sleeping %s seconds" %(self.pollingdelay,))
      time.sleep(self.pollingdelay)


  def Watch(self):
    time.sleep(1)
    currdirs = self.dirs.copy()
    cmop.debug("Watching directory: %s" % (currdirs,))
    def PhonyEvent(event):
      cmop.debug("Triggering Phony Event %s" % (event,), 6)
      for d in currdirs: 
        dispatch = self.MakeDispatch(d)
        try:
          cmop.debug("Checking directory %s" % (d,), 6)
          for f in os.listdir(d):
            cmop.debug("Found file %s" % (os.path.join(d,f),), 6)
            dispatch(os.path.join(d,f), event)
        except:
          cmop.error("Error watching directory: %s" % ("".join(traceback.format_exception(*sys.exc_info()))))
          
    PhonyEvent(gamin.GAMExists)
    while True: 
      #  while self.mon.event_pending():
      if ErrorRaised(): sys.exit()
      #    self.mon.handle_one_event()
      time.sleep(self.pollingdelay)

      PhonyEvent(gamin.GAMChanged)

def testGamin():
    mon = gamin.WatchMonitor()
    def test(filename, event): print filename, event
    watch = mon.watch_directory('.', test)
    time.sleep(1)
    while True:
      mon.handle_one_event()

def getargs():
  if len(sys.argv) < 2: 
    print "Usage: python %s <dir to monitor>" % tuple(sys.argv)
    sys.exit()
  d = sys.argv[1]
  return d

def Barnes():
  d = getargs()

  #t = TestScanner()
  #t = SeabirdScanner()
  def matchELG(filename):
    return ".elg" in filename.lower()

  source = DelimitedRecordStreamWithHeader(condition=matchELG, delimiter=",")
  log = Log(BarnesTSGCleaner('August2007', source), "./tsg")
  dep = LookupDeploymentID(log)
  pipeline = DropFilter(CMOPTableInsert(dep))
  pipeline.start()
  
  m = PollingIngestMonitor(delay=1)
  m.WatchDirectory(d)
  m.AddScanner(source)
  m.Watch()

def Wecoma():
  d = getargs()

  def matchCK(fname):
    return ".ck" in fname.lower()

  source = DelimitedRecordStreamWithHeader(condition=matchCK, delimiter=",", maxsize=2000)
  clean = WecomaTSGCleaner(vessel='Wecoma', cruise='August 2007', instrument='wecomatsg', instrumenttype='CTD', source=source, maxsize=2000)
  
  attach = TestFilter(LookupDeploymentID(clean, host=cmop.config.hostname, maxsize=2000))
  #ins = DropFilter(CMOPTableInsert(attach, host="cdb02", maxsize=2000))
  ins = attach  

  ins.start()

  m = PollingIngestMonitor(delay=1)
  m.WatchDirectory(d)
  m.AddScanner(source)
  m.Watch()

def Forerunner():
  d = getargs()

  def matchCTD(filename):
    return ".ctd" in filename.lower()

  forerunnerhead = [
     "conductivity", 
     "salinity", 
     "temperature", 
     "pressure",  
     "turbidity", 
     "oxygen"
  ]

  source = ForerunnerCastScanner(condition=matchCTD, 
                                 header=forerunnerhead,
                                 delimiter=" ")
  
  cast = ForerunnerCastCleaner(cruise='August2007', source=source)
  dep = LookupDeploymentID(cast)
  log = Log(dep, "./forercast")
  pipeline = DropFilter(CMOPTableInsert(cast))

  pipeline.start()
  
  m = PollingIngestMonitor(delay=1)
  m.WatchDirectory(d)
  m.AddScanner(source)
  m.Watch()

def TestCODAS():
  d = '/home/users/howew/testbed/uhdas/www.soest.hawaii.edu/martech/cruise/2007/w0703a_sanford/ADCP/proc/os75bb/adcpdb/'

  scan = CODASScanner()
  clean = CODASCleaner(vessel='Wecoma', cruise='August2007', instrument='wecoma', instrumenttype='ADP', source=scan)
  pipeline = CMOPTableInsert(clean)

  pipeline.start()
  m = PollingIngestMonitor(delay=1)
  m.WatchDirectory(d)
  m.AddScanner(scan)
  m.Watch()

def TestSeabird():
  d = getargs()

  def matchSeabird(filename):
    return ".dat" in filename.lower()

  scan = PassThroughScanner(matchSeabird)

  cast = SeabirdCastCleaner(
      vessel='Barnes', 
      cruise='August 2007',
      instrument='barnescast', 
      source=scan
  )

  #log = Log(cast, "./barnescast")

  #insert = DropFilter(CMOPTableInsert(TestFilter(cast), host="localhost"))

  pipeline = DropFilter(TestFilter(cast))
  #pipeline = insert

  pipeline.start()
  m = PollingIngestMonitor(delay=1)
  m.WatchDirectory(d)
  m.AddScanner(scan)
  m.Watch()

def TestRecursion():  
  d = getargs()
  scan = TestScanner()
  pipe = TestFilter(scan)
  pipe.start()

  m = PollingIngestMonitor(delay=1)
  m.RecursiveWatchDirectory(d)
  m.AddScanner(scan)
  m.Watch()

def TestSeabirdNoInsert():
  d = getargs()

  def matchSeabird(filename):
    return ".cnv" in filename.lower()

  #receives filenames
  scan =  PassThroughScanner(matchSeabird)

  # filter out casts already in the database
  # (also group redundant messages every 5 seconds)
  #newonly = KnownCast(scan, groupseconds=5, host='cdb02.stccmop.org')
  #pipe = DropFilter(TestFilter(newonly))


  # map CMOP variables to extracted Seabird field names
  # case insensitive, last match overrides previous matches.
  # match is defined as any field containing the value provided,
  # (hmmm...generalize to regular expressions?)
  # Ex:
  # 'fluorescence':'fluorometer' matches 'ECO/APL Fluorometer'
  # and assigns the value to the fluorescence column
  mapping = {
    'depth':'pressure',
    'pressure':'pressure',
    'salinity':'salinity',
    'conductivity':'conductivity',
    'fluorescence':'wetstar',
    'temperature':'temperature',
    'oxygen':'oxygen',
  }

  # (groups redundant messages every 5 seconds)
  cast = SeabirdCastCleaner(
      vessel='Wecoma', 
      cruise='August 2007',
      instrument='wecomacast', 
      mapping=mapping,
      source=scan,
      groupseconds=5
  )

  insert = CMOPTableInsert(cast, host=cmop.config.hostname)
  insert.db.suppressquery = True
  pipeline = DropFilter(insert)


  pipeline.start()
  m = PollingIngestMonitor(delay=1)
  m.WatchDirectory(d)
  m.AddScanner(scan)
  m.Watch()

if __name__ == '__main__':
  #cmop.DebugOn()
  TestSeabirdNoInsert()
  #Barnes()
  #Forerunner()
  #DebugOn()
  #Wecoma()
  #TestRecursion()

#cmop.DebugOn()
