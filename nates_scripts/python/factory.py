import sys, os, os.path, traceback

import cmop.webproduct as wp
import cmop.xmlutil as xml
import cmop
import cmop.pylabutil as pylabutil
import traceback

matplotlib = wp.matplotlib
matplotlib.use('Agg')
import pylab
from pylab import *

import cmop.db as db
import tempfile

cdb = db.DB(user="pturner", password="whitehouse")

import matplotlib.numerix.ma as ma
import numpy

noproduct = '''
<i>No product specified</i>
'''

DELIM = ","

BUILDERHTML = '''
<html>
<head>%s</head>
<body>
<h2>Data Product Builder</h2>
<div class="menucontainer">%s</div>
<div class="formcontainer">
<form>
<input type="text" name="name"></input>
Configure Parameters...
<input type="text" name="name"></input>
<input type="text" name="default"></input>
<input type="text" name="type"></input>
Extract with SQL...
<textarea name="sql"></textarea>
Present with PyLab...
<textarea name="pylab"></textarea>
<input type="submit" name="Save" value="Save"></input>
</form>
</div>
</body>
</html>
'''

HTML = '''<html><head>%s</head><body>%s</body></html>'''

DEFAULTWIDTH = 4.5
DEFAULTHEIGHT = 3.7

JQUERY = "/scripts/javascripts/jquery.js"
CSS = "/scripts/css/factory.css"
JAVASCRIPT = "ws/product/factory.js.py"
GETPRODUCT = 'ws/product/factory.py/getproduct'
GETASCII = 'ws/product/factory.py/getascii'
GETASCII = 'ws/product/factory.py/describe'
GETFULLFORM = 'ws/product/factory.py/getfullform'
PROCESSFORM = 'ws/product/factory.py/getfullform'
GETPAGE = 'ws/product/factory.py/getpage'
GETFORM = 'ws/product/factory.py/getform'
GETSPEC = 'ws/product/factory.py/getspec'

WRAPPER = '<div style="width:%sin" class="product %s">%s</div>'

STYLE = '''
<LINK REL="stylesheet" HREF="http://%s/%s" TYPE="text/css">
''' % ("%s", CSS)

JQUERY = '''<script src="/scripts/javascripts/jquery.js" type="text/javascript"></script>'''
SCRIPTS = '''
<script type="text/javascript">JQ=jQuery.noConflict(); //rename $ function</script>
<script src="http://%s/%s" type="text/javascript"></script>
''' % ("%s", JAVASCRIPT)

def MakeParamDiv(classes, attrs, content):
  attrstr = " ".join(['%s="%s"' % (k,v) for k,v in attrs])
  jsclass = " ".join(classes)
  return '''
<div %s class="pad parameter %s">%s</div>
''' % (attrstr, jsclass, content)

def MakeInput(label, classes, product, name, content, right=False):
  jsclass = " ".join(classes)
  input = '''
%s <input class="%s" size="15" product="%s" name="%s" type="text" value="%s"/>
''' % (label, jsclass, product, name, content)
  return MakeParamDiv([right],[],input)

OPTION = '''<option value="%s" %s >%s</option>'''

def MakeSelect(label, attrs, classes, product, name, content, right=False):
  jsclass =  " ".join(classes)
  attrstr = " ".join(["%s=%s" % (k,v) for k,v in attrs])
  sel = '''
%s <select %s class="%s" product="%s" name="%s">
%s
</select>
''' % (label, attrstr, jsclass, product, name, content)
  return MakeParamDiv([right],[],sel)

DESCRIBE = '''<description>%s</description>'''

def MakeHidden(product, name, content):
  return '''
<input product="%s" type="hidden" name="%s" value="%s"/>
''' % (product, name, content)


LINK = '''<a href="%s">%s</a>'''

P = '''<p>'''
BR = '''<br/>'''

SPAN = '''<span>%s</span>'''

PLOTURL = '''http://%s/%s?%s''' % ("%s", GETPRODUCT, "%s")
ASCIIURL = '''http://%s/%s?%s''' % ("%s", GETASCII, "%s")

DW = '''<div class="container left" style="width:%sin;">%s</div>'''

FULLFORM = '''
<form action="processfullform" style="width:%sin;" class="factoryform" name="getform">
<input product="%s" type="hidden" name="product" value="%s"/>
%s
<div class="container left" style="width:%sin;">
<input type="submit" name="request" value="Product"/>
<input type="submit" name="request" value="ASCII"/>
</div>
</form>
'''

PRODUCTFORM = '''<form class="factoryform" style="width:%sin;" name="getform">%s</form>'''


IMG = '''
<a class="display" href="%s">
<img class="pad display" src="%s"/>
</a>
'''


SHOWVALUES = "(show values)"

RESET = '''<form><input type="submit" value="Reset Form"></input>%s</form>'''
UPDATE = '''<input type="submit" value="Update Form"></input><p>''' 
SUBMITPROD = '''<input type="submit" value="Update Form"></input><p>''' 
REFRESH = '''<a class="pad left" href="javascript:RefreshProduct('%s');">Refresh Product</a>'''
MAKELINK = '''<a class="pad right" href="javascript:MakeRequest('%s', '%s')">%s</a>'''

def emptyplot():
   # make the active figure an empty plot
    # with zero-length arrays,
    # scatter throws an error instead of 
    # generating an empty plot,
    # so use plot()
   title("No data to display")
   p = plot([], [])

class Environment:
  pass

def selected(d,default): 
  if d == default:
    return 'selected'
  else:
    return ''

def columnize(right):
  if right: 
    style = '''right'''
  else: 
    style = '''left'''
  return style

def MakeOptionsHTML(options, default):
  options = [OPTION % (d,selected(d,default),d) for d in options]
  return "\n".join(options)

class FactoryError(Exception):
  pass

class Parameter:
  sel = "SELECT * FROM product.factoryproductparameter WHERE product = '%s' AND parameter = '%s'"
  upd = '''UPDATE product.factoryproductparameter SET type='%s', ordinal=%s, domain=$P$%s$P$, "default"=$E$%s$E$ WHERE product = '%s' AND parameter = '%s';'''
  ins = '''
INSERT INTO product.factoryproductparameter VALUES 
('%s', '%s', '%s', $P$%s$P$, $E$%s$E$, %s);
'''

  selattr = '''
SELECT * 
  FROM product.parameterattribute 
 WHERE product = '%s' AND parameter = '%s'
   AND attribute = '%s';
'''
  updattr = '''
UPDATE product.parameterattribute SET value='%s' 
WHERE product = '%s' AND parameter = '%s' AND attribute='%s'; 
'''
  insattr = '''
INSERT INTO product.parameterattribute
VALUES ('%s', '%s', '%s', '%s');
'''

  def __init__(self, name, value):
    self.name = name
    self.value = value
    self.env = {}
    self.attributes = {}
 
  def Value(self):
    return self.value

  def DisplayOnly(self):
    return False

  def ComputeDomain(self, env):
    pass

  def JavascriptClasses(self):
    return [self.__class__.__name__]

  def ComputeDefault(self, env):
    return self.value

  def AddAttribute(self, name, val):
    self.attributes[name]=val

  def DisplayName(self):
    return self.attributes.get('display', self.name)

  def Configure(self, env):
    ''' Configure a parameter in a given environment. (Execute queries, etc.)'''
    self.env = env.copy()
    self.domain = self.ComputeDomain(env)
    self.default = self.ComputeDefault(env)
    #self.context = self.ComputeContext(env)
    #self.description = self.ComputeDescription(env)

  def Store(self, product, ordinal):
    results = cdb.execQuery(self.sel % (product, self.name))
    if results:
      cdb.execCommand(self.upd % (self.__class__.__name__, ordinal, "", self.Value(), product, self.name))
    else:
      cdb.execCommand(self.ins % (product, self.name, self.__class__.__name__, "", self.Value(), ordinal))

    self.StoreAttributes(product)

  def StoreAttributes(self, product):

    for k,v in self.attributes.items():
      results = cdb.execQuery(self.selattr % (product, self.name, k))
      if results:
        cdb.execCommand(self.updattr % (v, product, self.name, k))
      else:
        cdb.execCommand(self.insattr % (product, self.name, k, v))

  def Default(self):
    '''
Returns the default value of the parameter if the parameter has been previously configured.
Behavior is undefined otherwise.
'''
    return self.default

  def Domain(self):
    if not hasattr(self, 'domain'):
      raise ValueError("Parameter Domain not available; parameter has not been configured")
    else:
      return self.domain

  def AsHTML(self, product, right=False, default=None):
    if not default: default = self.Default()
    jsclass = self.JavascriptClasses()
    input = MakeInput(self.DisplayName(), jsclass, product, self.name, default, columnize(right))
    return input

  def __repr__(self):
    return "%s" % (self.name)

class DynamicDefaultParameter(Parameter):

  def ComputeDefault(self, env):
    '''
Compute the default value of the parameter given an environment
of string/value pairs
'''
    rs = cdb.execQuery(self.value % env)
    if rs:
      return rs[0][0]
    else:
      return ''

  def AsHTML(self, product, right=False, default=None):
    if not default: default = self.Default()
    jsclass = [columnize(right)] + self.JavascriptClasses()
    input = MakeInput(self.DisplayName(), jsclass, product, self.name, default, columnize(right))
    return input

  def JavascriptClasses(self):
    return [self.__class__.__name__, "refreshable"]

class HiddenParameter(DynamicDefaultParameter):
  '''
A HiddenParameter ibehaves like a DynamicDefaultParameter but is not rendered in HTML or other user interfaces.
'''

  def AsHTML(self, product, right=False, default=None):
    #return HIDDEN % (product, self.name, self.Default())
    return ''

class CalculatedParameter(DynamicDefaultParameter):
  '''
A CalculatedParameter behaves like a DynamicDefaultParameter but is not editable.
'''

  def DisplayOnly(self):
    '''Indicates that this parameter is not used as an input to any other parameters.
Used for long metadata values that are cumbersome to trnasmit as part of a url.
'''
    return True

  def DisplayName(self):
    return self.attributes.get('display', '')

  def AsHTML(self, product, right=False, default=None):
    val = self.Default()
    attrs = [('product', product), ('name',self.name)]
    jsclass = self.JavascriptClasses() 
    div = "%s: %s" % (self.DisplayName(), MakeParamDiv(jsclass, attrs, val))
    return MakeParamDiv([columnize(right)], [], div)

class SelectParameter(Parameter):
  def __init__(self, name, rawdomain):
    Parameter.__init__(self, name, "")
    self.rawdomain = rawdomain
    self.env = {}

  def ComputeDomain(self, env):
    return self.rawdomain.split(',')

  def ComputeDefault(self, env):
    '''
Returns the default value of the parameter if 
the parameter has been previously configured.
Behavior is undefined otherwise.
For a select parameter, the default value is the first 
value in the select list.
'''
    if self.domain:
      return self.domain[0]
    else:
      return SHOWVALUES

  def Value(self):
    return self.rawdomain

  def Store(self, product, ordinal):
    results = cdb.execQuery(self.sel % (product, self.name))
    if results:
      cdb.execCommand(self.upd % (self.__class__.__name__, ordinal, self.Value(), "", product, self.name))
    else:
      cdb.execCommand(self.ins % (product, self.name, self.__class__.__name__, self.Value(), "", ordinal))

    self.StoreAttributes(product)

  def AsHTML(self, product, right=False, default=None):
    if not default: default = self.Default()
    display = self.DisplayName()
    optionhtml = MakeOptionsHTML(self.domain, default)
    jsclass = self.JavascriptClasses()
    sel = MakeSelect(display, [], jsclass, product, self.name, optionhtml, columnize(right))
    return sel

class SQLSelectParameter(SelectParameter):
  def __init__(self, name, sql):
    self.name = name
    self.rawdomain = sql
    self.env = {}
    self.attributes = {}

  def JavascriptClasses(self):
    return [self.__class__.__name__, "refreshable"]

  def ComputeDomain(self, env):
    '''
Compute the domain of the parameter.
'''
    results = cdb.execQuery(self.rawdomain % self.env)
    return [r[0] for r in results]


class Product:
  def __init__(self, name, title=None):
    self.name = name
    if title: self.title = title 
    else: self.title = name
    self.params = []
    self.matlab = None
    self.describe = ''
    self.sql = None
    self.boundvars = {'product':self.name}
    self.arrays = None
    self.width = DEFAULTWIDTH
    self.height = DEFAULTHEIGHT

  def Store(self):
    sel = "SELECT * FROM product.factoryproduct WHERE product = '%s'"
    upd = '''UPDATE product.factoryproduct 
                SET extractwith=$SQL$%s$SQL$, 
                    plotwith=$PYLAB$%s$PYLAB$, 
                    describewith=$DESCR$%s$DESCR$ 
              WHERE product = '%s';'''
    ins = '''
INSERT INTO product.factoryproduct (product, extractwith, plotwith, describewith)
VALUES ('%s', $SQL$%s$SQL$, $PYLAB$%s$PYLAB$, $DESCR$%s$DESCR$);
'''

    results = cdb.execQuery(sel % (self.name,))
    if results:
      cdb.execCommand(upd % (self.sql, self.matlab, self.describe, self.name))
      ret = 2
    else:
      cdb.execCommand(ins % (self.name, self.sql, self.matlab, self.describe))
      ret = 1

    delete = '''DELETE FROM product.factoryproductparameter WHERE product = '%s' ''' % (self.name,)
    cdb.execCommand(delete)

    for i,p in enumerate(self.params):
      p.Store(self.name, i)

    return ret

  def Specification(self):
    ext = '''<extractwith><![CDATA[%s]]></extractwith>'''
    plot = '''<plotwith><![CDATA[%s]]></plotwith>'''
    describe = '''<describewith><![CDATA[%s]]></describewith>'''
    prod = '''<product name="%s">%s</product>'''
    param = '''<parameter name="%s" type="%s"><![CDATA[%s]]>%s</parameter>'''
    def tag(p):
      attrs = "\n".join(["<%s>%s</%s>" % (attr, val, attr) for attr, val in p.attributes.items()])
      return param % (p.name, p.__class__.__name__,p.Value(), attrs)
    params = "\n".join([tag(p) for p in self.params])
    extractwith = ext % (self.sql,)
    plotwith = plot % (self.matlab,)
    describewith = describe % (self.describe,)
    return prod % (self.name, params + extractwith + plotwith + describewith)

  def AsLink(self):
    form = "/%s?%s" % (GETFULLFORM,"product=%s" % (self.name,))
    spec = "/%s?%s" % (GETSPEC,"product=%s" % (self.name,))
    example = "/%s?%s" % (GETPAGE,"product=%s" % (self.name,))
    formlink = LINK % (form, "query form")
    speclink = LINK % (spec, "specification")
    exlink = LINK % (example, "product form")
    return "%s (%s) (%s) (%s)" % (self.name, formlink, speclink, exlink)

  def SetSize(width, height):  
    self.width = width
    self.height = height  

  def Extract(self):
    groundsql = self.GroundSQL()
    self.fields, tuples = cdb.execQuery(groundsql, True)
    self.arrays = [array(x) for x in zip(*tuples)]

  def Environment(self):
    return self.env

  def GroundSQL(self):
    if not self.sql: 
      raise FactoryError("Extraction Failed for product '%s': SQL not specified" % (self.name,))

    try:
      groundsql = self.sql % self.Environment()
    except KeyError, e:
      raise FactoryError("Missing a required parameter: %s" %(e,))

    return groundsql

  def ExtractWith(self, sql):
    self.sql = sql

  def PlotWith(self, matlab):
    self.matlab = matlab

  def DescribeWith(self, describe):
    self.describe = describe

  def Describe(self):
    if self.describe:
      return DESCRIBE % (self.describe % self.Environment(),)
    else:
      return DESCRIBE % ('',)

  def Plot(self):
    figure(figsize=(self.width,self.height))
    title(self.title, size='smaller')

    if not self.arrays:
      emptyplot()
      return

    d = self.PylabEnvironment()

    try:
      exec self.matlab in d
    except Exception, e:
      close() 
      e,v,t = sys.exc_info()
      raise e,v,t
 
  def PylabEnvironment(self): 
    d = dict(zip(self.fields, self.arrays))
    d.update(self.Environment())
    _ = Environment()
    _.__dict__ = d
    d.update(locals())
    d.update(pylab.__dict__)
    d['numpy'] = numpy
    d['cmop'] = cmop
    d['ma'] = ma
    d['wp'] = wp
    d['matplotlib'] = matplotlib
    d['pylabutil'] = pylabutil
    return d

  def AsASCII(self, **kwargs):
    p = self
    hdr = p.Describe()
    root = xml.ParseXML(hdr)
    result = ""
  
    results = xml.Xpath("//*", root)
    tab = ""
    for r in results:
      attr = r.nodeName
      valnode = xml.Xpath("text()", r)
      if valnode:
        val = xml.ExtractValueFromNode(valnode[0])
        result += "%s%s: %s\n" % (tab, attr, val)
      else:
        result += "%s:\n" % (attr,)
        tab += "  "

    result += "--\n"
    yield result
    for t in p.AsTuples(**kwargs):
      yield "%s\n" % (DELIM.join(t),)
  
  def AsImage(self, **kwargs):
    self.Configure(**kwargs)
    try:
      self.Make(**kwargs)
    except:
      wp.errorplot("Unknown Error: ")
    cmop.debug("Making image...", 5)
    img = wp.MakeImage()
    cmop.debug("Closing...", 5)
    close()
    cmop.debug("Returning image...", 5)
    return img

  def Configure(self, **actualparams):
    '''
Bind parameters to values supplied by the user.  
Adjusts three internal structures:
  boundvars : a dictionary of user-supplied string, value pairs
  params : a list of parameter objects not supplied by the user
  env : a dictionary of string, value pairs with a key for every paramter required by the product
        If the user did not supply a value, the parameter's default is used.
'''

    # boundvars is the environment supplied by the user
    self.boundvars = dict([(k,v) for k,v in actualparams.items() if v != SHOWVALUES])

    # add a bound variable representing the product itself
    self.boundvars['product'] = self.name

    # params is the list of parameter objects not supplied by the user
    self.unboundparams = [p for p in self.params if not p.name in self.boundvars]

    # env is the full environment needed to generate a product instance
    # including any defaults

    self.env = self.boundvars.copy()
    # Assumes document order agrees with parameter dependenciesa
    for p in self.params:
      # Resolve all the dynamic parts of the parameter in the curent environment
      p.Configure(self.env)

      # Add the parameter value to the environment
      if not self.env.has_key(p.name):
        self.env[p.name] = p.Default()

  def AddParam(self, p):
    self.params.append(p)
    return p

  def Make(self, **kwargs):
    try:
      self.Extract()
    except FactoryError, SQLError:
      wp.errorplot("Error Extracting Data:")

    try:
      self.Plot()
    except:
      wp.errorplot("Error Plotting Data:")

  def ResetForm(self):
    # create a reset form
    def hide(k,v):
      return MakeHidden(self.name, k, v)
    basearghtml = "".join([hide(k, v) for k,v in self.env.items()])
    reset = RESET % (basearghtml,)
    return reset

  def AsTuples(self):
    self.Extract()

    yield self.fields
    for tup in zip(*self.arrays):
      yield [str(x) for x in tup]
    

  def MakeHTMLParameters(self):
    '''Return html input tags for unbound parameters. Assumes defaults are computed where necessary; that self.env represents a complete closure.'''
    inputs = ""
    for i, p in enumerate(self.unboundparams):
      inputs += p.AsHTML(self.name, right=0, default=self.env[p.name]) 

    return inputs

  def GetParam(self, paramname):
    for p in self.params:
      if p.name == paramname:
        return p

  def AsFullForm(self, server="no.server.specified", **kwargs):
    self.Configure(**kwargs)
    inputs = self.MakeHTMLParameters()
    width = kwargs.get('width', '11')
    #meat = SCRIPTS % (server,) + inputs
    meat = inputs
    body = FULLFORM % (width, self.name, self.name, meat, width)
    html = WRAPPER % (self.width,self.name, body)
    head = SCRIPTS % (server,) + STYLE % (server,)
    html = HTML % (head, html)
    return html
 
  def AsForm(self,server="no.server.com", **kwargs):
    ''' 
Return an HTML form representation of a configured product.
Form elements supplied in kwargs are rendered as hidden elements.
'''
 
    # update the parameter bindings for this product
    self.Configure(**kwargs)

    # create hidden inputs for bound parameters (plus one for the name of the product)
    def hide(k,v):
      return MakeHidden(self.name, k, v)
  
    hiddens = "".join([hide(k,v) for k,v in self.boundvars.items()])

    # construct HTML for the unbound parameters
    inputs = self.MakeHTMLParameters()
 
    reset = self.ResetForm()
  
    # Build the url out of true parameters, skipping those used for display purposes 
    args = []
    for k,v in self.env.items():
      p = self.GetParam(k)
      if p:
       if p.DisplayOnly(): 
         continue
      args.append((k,v))

    # create the plot url
    qs = '&'.join(["%s=%s" % t for t in args])
    url = PLOTURL % (server,qs)
    asciilink = MAKELINK % (self.name, 'getascii', 'Download ASCII')
    describelink = MAKELINK % (self.name, 'describe', 'Metadata')

    img = IMG % (url, url)
    scripts = SCRIPTS % (server,) + STYLE % (server,)
    links = DW % (self.width, REFRESH % (self.name,) + asciilink + describelink)
    form = inputs + hiddens #UPDATE + inputs + hiddens
    body = scripts + PRODUCTFORM % (self.width,form) + BR + links + img
    #body = scripts + reset + PRODUCTFORM % (self.widthi, form) + refresh + P + img
    html = WRAPPER % (self.width,self.name, body)
    return html

def LoadParams(product):
  paramsql = '''
  SELECT type, parameter::varchar(30), 
         CASE WHEN type='Parameter' THEN "default"
              WHEN type='DynamicDefaultParameter' THEN "default"
              WHEN type='HiddenParameter' THEN "default"
              WHEN type='CalculatedParameter' THEN "default"
              ELSE domain
          END as value
    FROM product.factoryproductparameter 
   WHERE product = '%s' 
   ORDER BY ordinal'''

  ps = cdb.execQuery(paramsql % (product.name,))
  msg = ""
  for param in ps:
    msg += '''%s('%s',"""%s""")''' % param
    p = eval('''%s('%s',"""%s""")''' % param)
    product.AddParam(p)

    LoadParameterAttributes(product, p)


def LoadProduct(name):
  sql = '''
SELECT product::varchar(30), extractwith, plotwith, describewith
  FROM product.factoryproduct
 WHERE product = '%s'
''' % (name,)

  results = cdb.execQuery(sql)
  if len(results) > 1: 
     raise ValueError("Duplicate product '%s' in database? check key constraints." % (name,))
  if len(results) == 0:
     raise ValueError("Product '%s' not found." % (name,))

  r = results[0]
  p = Product(r[0])
  p.ExtractWith(r[1])
  p.PlotWith(r[2])
  p.DescribeWith(r[3])

  LoadParams(p)

  return p

def LoadParameterAttributes(prod, param):
  sql = '''
SELECT attribute, value 
  FROM product.parameterattribute 
 WHERE product = '%s' 
   AND parameter = '%s'
'''
  results = cdb.execQuery(sql % (prod.name,param.name))
  for name, val in results:
    param.AddAttribute(name.strip(), val.strip())

def LoadAllProducts():
  sql = '''
SELECT product::varchar(30), extractwith, plotwith, describewith
  FROM product.factoryproduct'''
  results = cdb.execQuery(sql)
  
  prods = {}
  for r in results:
    p = Product(r[0])
    p.ExtractWith(r[1])
    p.PlotWith(r[2])
    p.DescribeWith(r[3])
    LoadParams(p)

    prods[r[0]] = p

  return prods

