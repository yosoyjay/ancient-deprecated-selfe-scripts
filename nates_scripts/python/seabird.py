import tempfile
import os.path
import popen2
import cmop
import re

# seabird software included in the python package
# accessing it at runtime is as easy as this, 
# thanks to the amazing setuptools package
from pkg_resources import resource_filename
seabirddir = resource_filename(__name__, 'seabird')

# header values of interest
hdr = {
  'instrumenttype' : '(Sea-Bird SBE \d+)', 
  'version'        : 'Software Version (.+)',
  'ship'           : 'Ship:(.+)',
  'cruise'         : 'Cruise:(.+)',
  'station'        : 'Station.*:(.+)',
  'file'           : 'FileName = (.+)',
  'tempserno'      : 'Temperature SN = (.+)',
  'condserno'      : 'Conductivity SN = (.+)',
  'castnum'        : 'Cast Number:(.+)',
  'station_no'     : 'Station Number: (.+)',
  'GMTdate'        : 'GMT Date \(DD\/MM\/YEAR\): (.+)',
  'GMTtime'        : 'GMT Time \(HHMM\): (.+)',
  'beddepth'       : 'Bottom Depth: (\d+)',
  'lat'            : 'Latitude = (.+) N',
  'lon'            : 'Longitude = (.+) W',
  'uploadtime'     : 'System UpLoad Time = (.+)',
  'NMEAtime'       : 'NMEA UTC \(Time\) = (\d\d:\d\d:\d\d)',
}   

class CONFileError(ValueError):
  pass

def dospath(path):
  return path.replace('/', '\\')

def dosabspath(path, drive='C'):
  ap = os.path.abspath(path)
  return "%s:%s" % (drive, dospath(ap))

def parsepath(f):
  root, name = os.path.split(f)
  base, ext = os.path.splitext(name)
  return root, base, ext

def copyname(d, f):
  ap = os.path.abspath(f)
  base, name = os.path.split(ap)
  path = os.path.join(d, name)
  return ap, path

def link(d, f, newname=None):
  ap = os.path.abspath(f)
  root, name = os.path.split(ap)
  base, ext = os.path.splitext(name)
  if not newname: newname = base
  path = os.path.join(d, "%s%s" % (newname, ext))
  os.symlink(ap, path)
  return path

class SeabirdCast:
  seabirddir = seabirddir
  drive = 'E'

  def __init__(self, dat, con=None, debug=False):
    self.dat = dat
    self.con = con
    self.debug = debug
    self.rows = []

  def GetConverted(self):
    root, fname = os.path.split(self.dat)
    base, ext = os.path.splitext(fname)
    cnvfile = os.path.join(root, "%s.%s" % (base, "cnv"))
    if os.path.exists(cnvfile):
      f = file(cnvfile)
      self.datcnvoutput = f.read()
    else: 
      #try caps?
      print type(cnvfile)
      cnvfile[-3:] = "CNV"
      if os.path.exists(cnvfile):
        f = file(cnvfile)
        self.datcnvoutput = f.read()
      else:
        self.Convert()

  def ConFile(self, dat):
    root, fname = os.path.split(dat)
    base, ext = os.path.splitext(fname)
    confile = os.path.join(root, "%s.%s" % (base, "CON"))
    if os.path.exists(confile):
      return confile
    else: 
      #try caps?
      confile = confile[:-3] + "con"
      if os.path.exists(confile):
        return confile
      else:
        raise CONFileError("No CON file found for DAT file '%s'" % (dat,))

  def makeconfig(self, cfg, d):
    cfgf = file(cfg)
    cfgcontent = cfgf.read()
 
    ab, cfgtemp = copyname(d,cfg)
    cfgftemp = file(cfgtemp, 'w+')
    cfgftemp.write(cfgcontent)
    return cfgtemp

  def GetBottle(self):
    blfile = os.path.splitext(self.dat)[0] + '.bl'
    if os.path.exists(blfile):
      bl = file(blfile)
      self.bloutput = bl.read()
    else:
      #try caps?
      blfile = os.path.splitext(self.dat)[0] + '.BL'
      if os.path.exists(blfile):
        bl = file(blfile)
        self.bloutput = bl.read()
      else: 
        # No converted data found; try and convert it ourselves
	self.bloutput = ''

  def GetConverted(self):
    root, fname = os.path.split(self.dat)
    base, ext = os.path.splitext(fname)
    confile = os.path.join(root, "%s.%s" % (base, "cnv"))
    if os.path.exists(confile):
      datcnv = file(confile)
      self.datcnvoutput = datcnv.read()
    else:
      #try caps?
      confile = confile[:-3] + "CNV"
      if os.path.exists(confile):
        datcnv = file(confile)
        self.datcnvoutput = datcnv.read()
      else: 
        # No converted data found; try and convert it ourselves
        if not self.con: self.con = self.ConFile(dat)
        self.datcnvoutput = self.Convert()

  def Convert(self, cfg="cmop.cfg"):
    '''Use dosemu to run Seabird data conversion.
       datcnv must run in the seabird directory, apparently
       datcnv is slow; prepare for about 260 scans per second
    '''

    dat = self.dat
    con = self.con

    cmop.info("Converting Seabird cast %s, %s" % (dat, con))

    root, base, ext = parsepath(self.dat)

    # file names must be short
    newname = "cast"

    # temp dir   
    d = tempfile.mkdtemp('X', 'X') 
    cmop.debug("tempdir: \n" + d)
  
    # dat file
    datpath = link(d, dat, newname)
    dosdat = dosabspath(datpath, self.drive)

    # con file
    conpath = link(d, con, newname)   
    doscon = dosabspath(conpath, self.drive)

    # cfg file. cfg file path needs to 
    #be relative to seabird directory, inexplicably

    #cfgpath = self.makeconfig(cfg, d)
    #doscfg = dosabspath(cfgpath, self.drive)

    # output file
    out = '%s.cnv' % (newname,)
    outpath = os.path.join(d, out)
    dosout = dosabspath(outpath, self.drive)

    # batch file
    bat = '''
%s:
cd %s 
datcnv.exe -ax -o%s -i%s -c%s -s -e%s
exitemu
'''
    batfile = os.path.join(d, 'seabird.bat')
    f = file(batfile, 'w+')
    dosseabirddir = dosabspath(self.seabirddir, self.drive)
    content = bat % (self.drive, dosseabirddir, dosout, dosdat, doscon, cfg)
    cmop.debug("Batchfile: \n" + content)
    f.write(content)
    f.close()

    # dosemu command
    cmd = 'dosemu -t -quiet -5 -E "%s"' % (batfile,)

    sout, sin = popen2.popen2(cmd)
    response = sout.read()
 
    try:
      f = file(outpath)
      results = f.read()

      os.remove(datpath)
      os.remove(conpath)
      os.remove(outpath)
      os.remove(batfile)
      os.rmdir(d)

      self.datcnvoutput = results
      return results

    except: 
      f = file(os.path.join(d, 'terminal.log'), 'w+')
      f.write(response)
      raise ValueError("Error running datcnv.exe.  See terminal.log in %s" % (d,))
  
  def ParseBottle(self):
    '''Parse the content of a .BL file, expects bottle file to have been read into bloutput attribute.'''
    lines = self.bloutput.split('\n')
    
    self.bottles = []
    for line in lines[2:]:
      line = line.replace('\r', '')
      # data line
      row = re.sub('^\s+', '', line)
      row = re.sub('\s+$', '', row)
      # make sure the line isn't blank
      if row:
        row = [r for r in re.split(',\s+', row)]
        self.bottles.append(row)
   
  def Parse(self):
    '''Parse the content of a datcnv output file. Expects datcnv file to have been read into datcnvoutput attribute.'''

    lines = self.datcnvoutput.split('\n')

    for line in lines:
      line = line.replace('\r', '')

      if re.match('^\*', line):
        # header line
        for key, expr in hdr.items():
          d = re.search(expr, line)
          if d: 
            self.__dict__[key] = d.group(1).strip()
            cmop.debug("HDR: %s=%s" % (key, self.__dict__[key]))
            break

      elif re.match('^\#', line):
        # metadata line

        col = re.search('# (.+) (\d+) = (.+)', line)
        if col:
          cmop.debug("COL: %s" % (col.groups(),))
          self.__dict__.setdefault(col.group(1).strip(), {})
          self.__dict__[col.group(1)][int(col.group(2))] = col.group(3)

        else:
          attr = re.search('# (.+) = (.+)', line)
          if attr:
            cmop.debug("ATTR: %s" % (attr.groups(),))
            self.__dict__[attr.group(1)] = attr.group(2).strip()
          else:
            #unrecognized line
            cmop.debug("unrecognized line format: %s" % (line,))

      else:
        # data line
        row = re.sub('^\s+', '', line)
        row = re.sub('\s+$', '', row)

        # make sure the line isn't blank
        if row:
          row = [r for r in re.split('\s+', row)]
          self.rows.append(row)

  def ValidateNames(self):
    '''Checks the column names for proper formats.
Assumes a datcnv output has already been parsed.'''

    if not hasattr(self, "name"): 
      self.name = {}
    if not hasattr(self, "names"): 
      self.names = {}
    for k, v in self.name.items():

       # seabird code, seabird name, optional units
       m = re.match('(.+): (.+?)(?: \[(.+)\]|(?:\s+)?$)', v)
       if m:
         self.names[k] = (m.group(1), m.group(2), m.group(3))
       else:
         raise ValueError("Unknown column descriptor format: %s" % (v,))

    # convert the dictionary into a list for easy zipping
    namelist = [self.names[i] for i in range(len(self.names))]
    self.namelist = namelist

  def Validate(self):
    ''' Validate the parsed cast.'''
    self.ValidateNames()
    self.ValidateSpans()

  def ValidateSpans(self):
    '''Checks the span lines for proper format.
Assumes a datcnv output has already been parsed.'''
  
    if not hasattr(self, "span"): 
      self.span = {}
 
    for k, v in self.span.items():
       m = re.match('(.+), (.+)', v)
       if m:
         # minval, maxval
         self.span[k] = (float(m.group(1)), float(m.group(2)))
       else:
         raise ValueError("Unknown span format: %s" % (v,))

  def ValidateBottles(self):
    '''Gets the bottle column names. Currently uses a hard-coded list.'''
    self.bottlename = ['bottle','firesequence','time','startscan','endscan'] 

  def IterateBottle(self,key):
    '''Returns a generator for bottles.  Each item is a dictionary'''
    self.ValidateBottles()
    self.ValidateNames()
    depthname = []
    index = None
    for i in range(len(self.namelist)):
      if re.match(key.lower(),self.namelist[i][1].lower()):
       index = int(i)
    for b in self.bottles:
      row = dict(zip(self.bottlename, b))
      if index is not None:
	s = int(row['startscan']) - 1
        e = int(row['endscan']) 
        n = e - s -1
        d = []
        for r in self.rows[s:e]:
          d.append(r[index])
        d = map(float,d)
        row['depth'] = sum(d)/n
      else: 
        row['depth'] = None
      yield row

    # convert the dictionary into a list for easy zipping
  def Iterate(self):
    '''Returns a generator.  Each item is a dictionary'''
    self.Validate()
    cols = ["%s: %s" % (code, name) for (code, name, units) in self.namelist]
    for r in self.rows:
      yield dict(zip(cols, r))


if __name__ == '__main__':    
  
  d = SeabirdCast('/home/cruise/ingest/ctd/cast006.cnv', '/home/cruise/ingest/ctd/cast006.CON')
  #d = SeabirdCast('/home/cruise/ingest/ctd/cast125.cnv', '/home/cruise/ingest/ctd/cast125.CON')
  d.GetConverted()
  d.Parse()
  #for r in d.Iterate():
   # print r
  
  d.GetBottle()
  d.ParseBottle()
  for r in d.IterateBottle('Pressure'):
    print r
  
