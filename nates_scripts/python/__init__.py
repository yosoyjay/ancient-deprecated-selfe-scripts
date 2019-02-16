import os
import re

import urllib
import popen2

debug_print = False
import ConfigParser
import logging
import traceback
import logging.handlers as handlers

import sys, os.path

path, name = os.path.split(sys.argv[0])
topcaller, ext = os.path.splitext(name)

syslogger = logging.getLogger('cmop')
filelog = handlers.RotatingFileHandler('/tmp/%s.log'  % (topcaller,), maxBytes=1024*1024*5, backupCount=10)
filelog.setLevel(logging.INFO)
frmttr = logging.Formatter('%(levelname)s %(asctime)s: %(message)s')
filelog.setFormatter(frmttr)
syslogger.addHandler(filelog)

consolelog = logging.StreamHandler()
consolelog.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s %(asctime)s: %(message)s')
consolelog.setFormatter(formatter)
syslogger.addHandler(consolelog)

syslogger.setLevel(1)


def DebugOn(level=logging.DEBUG):
  SetLogLevel(level)

def DebugOff():
  SetLogLevel(logging.INFO)

def SetLogLevel(level):
  syslogger.setLevel(level)
  filelog.setLevel(level)
  consolelog.setLevel(level)

def logexception(prefix=""):
  msg = prefix + "".join(traceback.format_exception(*sys.exc_info()))
  error(msg)

def debug(msg, level=logging.DEBUG):
  '''
Log a debugging message.
Second argument is the debugging level.  Usually 1-10.  Defaults to 10 (logging.DEBUG)
'''
  if level == logging.DEBUG:
    syslogger.debug(msg)
  else:
    syslogger.log(level, msg)
  consolelog.flush()

def info(msg):
  syslogger.info(msg)
def critical(msg):
  syslogger.critical(msg)
def error(msg):
  syslogger.error(msg)
def profile(msg):
  '''
Log a profiling message.
Profiling information is rather verbose, so use a lower level than DEBUG
'''
  syslogger.log(5,msg)

class config:
    cfp = ConfigParser.ConfigParser()
    conffile = "/usr/local/cmop/etc/server.ini"
    try:
      info("Reading config file '%s'" % (conffile,))
      cfp.read([conffile])
      hostname = cfp.get("dbserver","dbserver")
      dbname = cfp.get("dbserver","dbname")
      user = cfp.get("dbserver", "dbuser")
      password = cfp.get("dbserver", "dbpasswd")
    except:
      info("Error using config file '%s' ...using defaults" % (conffile,))
      dbname = 'cmop'
      hostname = 'cdb02.stccmop.org'
      user = 'reader'
      password = ''


def DM2DD(lat):
  deg, mn = lat.split()[0:2]
  if "W" in mn or "N" in mn or "E" in mn or "S" in mn:
    mn = mn[:-1]
  return float(deg) + float(mn)/60.0

def LoopCounter(ranges, start, end):
  assert(len(ranges) == len(start))
  for (s, (b,e)) in zip(start, ranges):
    if s < b or s > e:
      msg = "%s is not a valid counter for %s" % (start, ranges)
      raise ValueError(msg)

  sizes = [(s, e-s) for (s,e) in ranges]

  n = len(start)
  counter = list(start[:])

  while tuple(counter) != tuple(end):
    yield counter
    for i in range(1,n+1):
      offset, size = sizes[n-i]
      counter[n-i] = ((counter[n-i] - offset  + 1) % size) + offset
      if counter[n-i] != offset: break

def TimestepIterator(start, end, timestepseconds):
  # week number, day number, timestep number
  modelranges = [(1,53),(1,8),(0,86400/timestepseconds)]
  return LoopCounter(modelranges, start, end)

def UnQuote(s):
  return urllib.unquote(s)

def Quote(s):
  return urllib.quote(s.encode('utf-8'), "/:")

def LockFileName(expr):
  return '/tmp/.%s.lock' % (expr,)

def LockFile(expr):
  f = file(LockFileName(expr), 'w')
  
def Locked(expr):
  return os.path.exists(LockFileName(expr))

def UnlockFile(expr):
  if os.path.exists(LockFileName(expr)):
    os.remove(LockFileName(expr))

def Ping(host, n=2):
  cmd = "ping -W 4 -c %s %s" % (n, host)
  print cmd
  cout, cin = popen2.popen2(cmd)
  response = cout.read()
  m = re.search("(\d+)% packet loss, time (\d+)ms", response)
  if not m:
    return False
  else:
    return m.groups(0)[1]

def AlreadyRunning(expr):
  '''Return count of processes matching expr.  Unix Only.'''
  cmd = "ps -ef | grep %s | wc -l" % (expr,)
  cmd = "ps -ef | grep %s | wc -l" % (expr,)
  cout, cin = popen2.popen2(cmd)
  response = cout.read()
  # one for the shell, one for the grep
  return int(response) - 2
 
   
def WestCoastline(bounds=(-124,33,-120,50)):
  import db
  '''
Returns a list of list of (lon, lat) pairs representing the west coast of the US.
Argument is a tuple (xmin, ymin, xmax, ymax) representing a bounding rectangle.
  '''
  import db
  d = db.DB() 
  sql = '''
SELECT segment, longitude(location), latitude(location) 
  FROM coastlinepoints
 WHERE location
        && setsrid('BOX3D(%s %s, %s %s)'::box3d, 4326)
''' % bounds
  coords = d.execQuery(sql)

  ret = {}
  for s, lon, lat in coords:
    ret.setdefault(s,[]).append((lon, lat))

  return ret.values()
