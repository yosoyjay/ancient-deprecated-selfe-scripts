# A library for generating data products on the web using pylab
import tempfile, os, traceback
import cmop

fullwritedir = "/tmp/"

try:
    from mod_python import apache, util, psp
    # running from a webserver
    os.environ["HOME"] = fullwritedir
    import matplotlib
    matplotlib.use('Agg')
except:
    import matplotlib

from pylab import *

class filestring:
  def __init__(self):
    self.data = ""

  def write(self, val): 
    self.data += val


def errorplot(prefix, width=10):
  msg = prefix
  s = filestring()
  traceback.print_exc(file=s)
  if s.data:
    msg += s.data
  
  fontsize = 10
  lines = ["      %s" % (c,) for c in msg.split('\n')]
  n = len(lines) + 1
  ht = n/4.0
  f = figure(figsize=(width,ht))
  text(0,0,"\n".join(lines),fontsize=fontsize, clip_on=False)

def wait(req):
  f = file("wait.gif")
  img = f.read()
  f.close()
  req.content_type = 'image/gif'
  return img

def MakeImage():
  # save the active figure to a tempfile
  # should be able to pass in a file object, but alas
  cmop.debug("Creating tempfile...", 5)
  (f, name) = tempfile.mkstemp(".png", dir=fullwritedir)
  os.close(f)
  cmop.debug("Saving image to tempfile %s" % (name,), 5)
  savefig(name, format='png')

  # close the active figure
  close()

  f = file(name)
  img = f.read()
  f.close()
  os.remove(name)

  cmop.debug("Returning image...", 5)
  return img

