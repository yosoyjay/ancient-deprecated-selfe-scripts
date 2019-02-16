import sys, urllib, urlparse
import xml.xpath as xpath
import xml
import cgi
import StringIO 
import time, datetime
import xml.dom.minidom
from xml.dom.ext.reader import Sax2
from xml.dom.ext import PrettyPrint

namespaces = {
  "xsi":"http://www.w3.org/2001/XMLSchema-instance", 
  "swe":"http://www.opengis.net/swe/0", 
  "gml":"http://www.opengis.net/gml", 
  "sos":"http://www.opengis.net/sos/0", 
  "om":"http://www.opengis.net/om", 
  "ows":"http://www.opengeospatial.net/ows", 
  "xlink":"http://www.w3.org/1999/xlink" 
}

def Quote(s):
  return urllib.quote(s, "/:")

def UnQuote(s):
  return urllib.unquote(s)

class SOSException(Exception):
  pass

def ExtractValueFromNode(node):
  if node.nodeName == "#text":
    # text node
    value = node.nodeValue
  elif node.nodeValue == None:
    # element node
    value = node.localName
  else:
    # attribute node, or something we don't want
    value = node.nodeValue
  #v = value.strip()
  v = value.encode('utf-8')
  return v

class XMLRequest:
  def __init__(self, url):
    self.result = None
    self.url = url

  def Clear(self):
    self.result = None

  def GetValues(self, xpath):
    if not self.result: self.FetchXML()
    values = self.XPath(xpath)
    return [ExtractValueFromNode(o) for o in values]

  def FetchXML(self):
    response = urllib.urlopen(self.url)
    xmltext = response.read()

    reader = Sax2.Reader()
    d = reader.fromString(xmltext)
    context = xml.xpath.Context.Context(d)
    context.setNamespaces(namespaces)
    query = "/"
    e = xml.xpath.Compile(query)
    result = e.evaluate(context)
    self.result = result[0]

  def XPath(self, xpath, node=None):
    if not node: node = self.result
    if type(xpath) in (str, unicode):
      xpath = xml.xpath.Compile(xpath)

    c = xml.xpath.Context.Context(node)
    c.setNamespaces(namespaces)

    nodes = xpath.evaluate(c)
    return nodes

  def AsText(self):
    s = StringIO.StringIO()
    PrettyPrint(self.result,s)
    return s.getvalue()

class SOSRequest(XMLRequest):
  def IsException(self):
    return self.GetValues('//Exception | //exception')

  def ExceptionText(self):
    return self.GetValues('//ExceptionText')

  def FetchXML(self):
    XMLRequest.FetchXML(self)
    if self.IsException():
      raise SOSException(cgi.escape(self.AsText()))

class GetCapabilities(SOSRequest):
  def __init__(self, rooturl):
    url = rooturl 
    SOSRequest.__init__(self, url)

  def GetOfferings(self):
    return self.GetValues('//sos:ObservationOffering/@gml:id')

  def GetOfferingsContaining(self, pattern):
    return self.GetValues('//sos:ObservationOffering[contains(gml:description/text(),"%s")]/@gml:id' % (pattern,))

  def GetObservedProperties(self, offering):
    offpath = '//sos:ObservationOffering[@gml:id="%s"]' % (offering,)
    return self.GetValues(offpath + '/sos:observedProperty/@xlink:href')

  def GetTimePeriod(self, offering):
    offpath = '//sos:ObservationOffering[@gml:id="%s"]' % (offering,)
    begins = self.GetValues(offpath + '/sos:time/gml:TimePeriod/gml:beginPosition/text()')
    ends = self.GetValues(offpath + '/sos:time/gml:TimePeriod/gml:endPosition/text()')
    if not begins:
      begins = self.GetValues(offpath + '/sos:time/gml:TimePeriod/gml:beginPosition/@indeterminatePosition')

    if not ends:
      ends = self.GetValues(offpath + '/sos:time/gml:TimePeriod/gml:endPosition/@indeterminatePosition')

    return begins[0], ends[0]

  def GetEnvelope(self, offering):
    offpath = '//sos:ObservationOffering[@gml:id="%s"]' % (offering,)
    lowercorners = self.GetValues(offpath + '/gml:boundedBy/gml:Envelope/gml:lowerCorner/text()')
    uppercorners = self.GetValues(offpath + '/gml:boundedBy/gml:Envelope/gml:upperCorner/text()')
    return tuple(lowercorners[0].split()), tuple(uppercorners[0].split())

  def GetOperations(self):
    return self.GetValues('//ows:Operation/@name')

  def GetParameters(self, operation):
    oppath = '//ows:Operation[@name="%s"]' % (operation,)
    return self.GetValues(oppath + '/ows:Parameter/@name')

  def GetParameterValues(self, operation, parameter):
    oppath = '//ows:Operation[@name="%s"]' % (operation,)
    parampath = '/ows:Parameter[@name="%s"]' % (parameter,) 
    return self.GetValues(oppath + parampath + '/ows:Value/text()')

class DescribeSensor(SOSRequest):
  def __init__(self, sensorid, rooturl):
    self.sensorid = sensorid
    o = urlparse.urlparse(rooturl)
    url = "%s://%s/%s" % o[0:3] + "?request=DescribeSensor&SystemId=%s" % (Quote(sensorid),)
    SOSRequest.__init__(self, url)

  def OutputList(self):
    return self.GetValues("//OutputList/output/@name")

  def RecordTypeOf(self, name):
    return self.GetValues('//OutputList/output[@name="%s"]/swe:DataRecord/swe:field/@name' % (name,))
    

class GetObservation(SOSRequest):
  def GetResult(self):
    blocksep = self.GetValues('//swe:AsciiBlock/@blockSeparator')
    tokensep = self.GetValues('//swe:AsciiBlock/@tokenSeparator')
    rawresults = self.GetValues('//om:result/text()')

    results = []
    for r in rawresults: 
      b = blocksep[0]
      t = tokensep[0]

      for row in r.strip().split(b):
        results.append(tuple([f for f in row.split(t)]))

    return results

class GetObservation(GetObservation):
  def __init__(self, rooturl, offering, observedProp):
    self.offering = offering
    self.observedProp = observedProp
    self.conditions = []
    self.conditions.append(("request","GetObservation"))
    self.conditions.append(("offering",Quote(offering)))
    self.conditions.append(("observedProperty",Quote(observedProp)))
    self.rooturl = rooturl

  def FetchXML(self):
    conds = ["%s=%s" % c for c in self.conditions]
    o = urlparse.urlparse(self.rooturl)
    self.url = "%s://%s/%s" % o[0:3] + "?" + "&".join(conds)
    SOSRequest.__init__(self, self.url)
    SOSRequest.FetchXML(self)
     
  def AddTimeConstraint(self, timeinstant):
    self.conditions.append(('eventTime', timeinstant))

  def AddTimeRangeConstraint(self, timeranges):
    ranges = []
    for pair in timeranges: 
      ranges.append('/'.join(pair))
    timesyntax = ','.join(ranges)
    
    self.conditions.append(('eventTime', Quote(timesyntax)))

  def AddBoundingBoxConstraint(self, p1, p2):
    self.conditions.append(('bbox',"%s,%s,%s,%s" % (p1+p2)))



def GetStationData(getcapsurl, station, variable, hours="24"):
  # print the data for the given station and time period
  r = GetCapabilities(getcapsurl)
  r.FetchXML()

  endtime = datetime.datetime.now()
  delta = datetime.timedelta(hours=float(hours))
  starttime = endtime - delta

  ro = GetObservation(getcapsurl, station, variable)
  range = [d.strftime("%Y-%m-%d %H:%M:%S") for d in (starttime, endtime)]
  ro.AddTimeRangeConstraint([tuple(range)])
  ro.FetchXML()
  result = ro.GetResult()
  for r in result:
    d, t = r[0].split('T')
    year, month, day = d.split('-')
    t, z = t.split('+')
    hour, minute, second = t.split(':')
    print year, month, day, hour, minute, r[4]

def GetStationVariables(getcapsurl, station):
  # print the data for the given station and time period
  r = GetCapabilities(getcapsurl)
  r.FetchXML()

  props = r.GetObservedProperties(station)
  for p in props:
    print p

def ListStations(getcapsurl):
  # List all stations along with their lat/lon envelope
  r = GetCapabilities(getcapsurl)
  r.FetchXML()

  # use all offerings
  offs = r.GetOfferings()
  print "Offering_url, min_lon, min_lat, max_lon, max_lat "
  for o in offs:
    box = r.GetEnvelope(o)
    print "%s, %s, %s, %s, %s" % (o, box[0][0], box[0][1], box[1][0], box[1][1])


def Check(name, cond, env):
  sys.stdout.write(name)
  sys.stdout.flush()
  if eval(cond, env):
    sys.stdout.write('success')
  else:
    sys.stdout.write('FAILED')
  print ""
 
def ResultsReturned(rqst):
  try:
    rqst.FetchXML()
    rs = rqst.GetResult()
    success = len(rs)
  except SOSException:
    success = -1
  return success

def RunTest():
  if len(sys.argv) == 2:
    getcapsurl = sys.argv[1]
    detail = 0
  elif len(sys.argv) == 3:
    getcapsurl = sys.argv[1]
    detail = sys.argv[2]
  else:
    getcapsurl = 'http://data.stccmop.org/ws/sos.py'

  sys.stdout.write("Fetching GetCapabilites document ...")
  sys.stdout.flush()
  r = GetCapabilities(getcapsurl)
  r.FetchXML()
  sys.stdout.write("success\n")

  ops = r.GetOperations()
  Check("GetCapabilities advertised? ... ", "'GetCapabilities' in ops", locals())
  Check("DescribeSensor advertised? ... ", "'DescribeSensor' in ops", locals())
  Check("GetObservation advertised? ... ", "'GetObservation' in ops", locals())

  ps = r.GetParameters('GetCapabilities')
  Check("GetCapabilities supports 'service' argument? ... ", "'service' in ps", locals())
  ds = r.GetParameterValues('GetCapabilities', 'service')
  Check("Domain for service is 'SOS'? ...", "['SOS'] == ds", locals())

  ds = r.GetParameterValues('GetCapabilities', 'version')
  Check("Domain for version is '0.0.31'? ...", "['0.0.31'] == ds", locals())

  Check("GetCapabilities supports 'version' argument? ... ", "'version' in ps", locals())
  ps = r.GetParameters('DescribeSensor')
  Check("DescribeSensor supports 'SensorId' argument? ... ", "'SensorId' in ps", locals())
  sensors = r.GetParameterValues('DescribeSensor', 'SensorId')
  Check("Domain for 'SensorId' is populated? ... ", "len(sensors) > 0", locals())

  ps = r.GetParameters('GetObservation')
  Check("GetObservation supports 'offering' argument? ... ", "'offering' in ps", locals())
  Check("GetObservation supports 'observedProperty' argument? ... ", "'observedProperty' in ps", locals())
  Check("GetObservation supports 'eventTime' argument? ... ", "'eventTime' in ps", locals())
  Check("GetObservation supports 'bbox' argument? ... ", "'bbox' in ps", locals())

  ds = r.GetParameterValues('GetObservation', 'offering')
  Check("Domain for 'offering' is populated? ... ", "len(ds) > 0", locals())
  check = "len(r.GetParameterValues('GetObservation', 'observedProperty')) > 0"
  Check("Domain for 'observedProperty' is populated? ... ", check, locals())


  for s in sensors[-detail:]:
    
    sys.stdout.write("Test DescribeSensor for '%s' ... " % (s,))
    try:
      d = DescribeSensor(s, getcapsurl)
      d.FetchXML()
      sys.stdout.write("success")
    except SOSException:
      sys.stdout.write("FAILED")
    sys.stdout.write("\n")
    
    ol = d.OutputList()
    Check("Sensor '%s' provides OutputList? ..." % (s,), "len(ol) > 0", locals())
 
    if detail > 2:
      for o in ol:
        rt = d.RecordTypeOf(o)
        goodtype = ['time', 'latitude', 'longitude', 'depth']
        Check("OOSTethys record structrue for output '%s' at sensor '%s'? ... " % (o,s), "rt[:4] == goodtype", locals())

    
    d = DescribeSensor("foobar", getcapsurl)
    try:
      d.FetchXML()
      success = False
    except SOSException:
      success = True
    
    Check("Proper exception returned for bad sensor? ...", "success", locals())

  # use all offerings
  offs = r.GetOfferings()
  #offs = r.GetOfferingsContaining('ship')

  # just check the first one
  for o in offs[-detail:]:
    try:
      obs = r.GetObservedProperties(o)
      success = len(obs) > 0
    except:
      success = False
    Check("Offering '%s' advertises observedProperties? ..." % (o,), "success", locals())

    try:
      env = r.GetEnvelope(o)
      success = len(env) > 0
    except:
      success = False
    Check("Offering '%s' provides envelope? ..." % (o,), "success", locals())

    try:
      t = r.GetTimePeriod(o)
      success = len(t) > 0
    except:
      success = False
    Check("Offering '%s' provides time period? ..." % (o,), "success", locals())

    # just check the first one
    for ob in obs[-detail:]:
      orqst = GetObservation(getcapsurl, o, ob)
      vars = locals()
      vars.update(globals())
      Check("Get Latest for Offering '%s' and observedproperty '%s' ..." % (o,ob), "ResultsReturned(orqst)> 0", vars)

      hours = 24
      endtime = datetime.datetime.now()
      delta = datetime.timedelta(hours=float(hours))
      starttime = endtime - delta

      range = [d.strftime("%Y-%m-%d %H:%M:%S") for d in (starttime, endtime)]
      orqst.AddTimeRangeConstraint([tuple(range)])

      n = ResultsReturned(orqst)
      Check("Data Available in last 24 hours for Offering '%s' and observedproperty '%s'? ..." % (o,ob), "n > 0", locals())

      env = [(float(e), float(f)) for e,f in env]
      ll = (env[0][1]-0.001, env[0][0]-0.001)
      ur = (env[1][1]+0.001, env[1][0]+0.001)
      orqst.AddBoundingBoxConstraint(ll, ur)
      n = ResultsReturned(orqst)
      Check("Data Available in bounding box advertised by GetCapabilities in last 24 hours? ...", "n > 0", locals())

      r = GetObservation(getcapsurl, o, ob)
      r.AddTimeRangeConstraint([("foo", "bar")])

      n = ResultsReturned(r)
      Check("Exception returned for bogus time format?...", "n == -1", locals())


if __name__ == '__main__':
  RunTest()
