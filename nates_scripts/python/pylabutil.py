import pylab 
from matplotlib.cbook import silent_list
import datetime
import numpy
import urllib, tempfile
from matplotlib.dates import WeekdayLocator,DateFormatter,SecondLocator,MinuteLocator,HourLocator,DayLocator,MonthLocator,YearLocator
def secondaxes(ax=None,position=True,sharedside = 'x'):
    """
    Make a new axes overlay ax (or the current axes if ax is None)
    sharing the specified axis.  The ticks for axn will be placed 
      on the specified side,
    and the axn instance is returned.
    
    ax = axes to overlay
    position = 'top' or 'bottom' for ticks
    sharedside = 'x' or 'y' for shared axis of graph 
       (forces agreement along this axis)
    """
    if ax is None:
        ax=pylab.gca()

    if sharedside == 'x':
        axn = pylab.gcf().add_axes(ax.get_position(), sharex=ax, frameon=False)
    elif sharedside == 'y':
        axn = pylab.gcf().add_axes(ax.get_position(), sharey=ax, frameon=False)
    pylab.draw_if_interactive()
    if sharedside == 'x':
      choose_axes(axn,position=position,side='y')
    if sharedside == 'y':
      choose_axes(axn,position=position,side='x')
    return axn

def adjust_axis(axn=None,c='r',step=0,side='x',position=True):
    '''
    adjusts position and color of x or y ticks

    axn = selected axis
    c = new color for ticks
    step = layer of offset
    side = x or y
    position = True or False True = top or right
    '''
    if axn is None:
      axn = pylab.gca()
    pylab.axes(axn)
    #locator = ticker.AutoLocator()
    #locator.autoscale()
    #eval('axn.' + side + 'axis.set_major_locator(locator)')
    #formatter = ticker.FixedFormatter()
    #reval('axn.' + side + 'axis.set_major_formatter(formatter)')
    choose_axes(axn,position=position,side=side)
    yt = pylab.getp(axn,side + 'ticklabels')
    ticklabels = pylab.getp(axn,side + 'ticklabels')
    pylab.setp(ticklabels,'color',c)
    pylab.setp(eval('axn.' + side + 'axis.label'),color=c)
    yt = pylab.getp(axn,side + 'ticklines')
    pylab.setp(yt,color=c)
    if step > 0:
        dpad = 0 
        dlen = 25 * step
	if side == 'x' and position:
	  drot = 90
	  ddir = 1
	elif side == 'x':
	  drot = -90
	  ddir = -1
  	elif side == 'y' and position:
	  drot = 0
          ddir = 1
        else:
          ddir = -1
          drot = 180
          dpad = -5
	eval('axn.'  + side + 'axis.get_major_ticks()')
	ticklabels = pylab.getp(axn,side + 'ticklabels')
        pylab.setp(ticklabels,dashdirection=ddir, dashlength=dlen, dashrotation=drot,dashpad=dpad)
    eval('pylab.setp(axn.get_' + side + 'gridlines(),c=\'' + c + '\')')
    pylab.draw_if_interactive()

def choose_axes(axn,position=True,side='x'):
    '''Set display of axis features to be one sided
       
       side specifies whether to act on x or y axis
       if position is true then use right or top of graph
    '''
    if axn is None:
        axn=pylab.gca()
    if side == 'y':
        if position:
           position = 'right'
        else:
           position = 'left'
    elif side == 'x':
        if position:
           position = 'top'
        else:
           position = 'bottom'
    eval('axn.' + side + 'axis.tick_' + position + '()')
    eval('axn.' + side + 'axis.set_label_position(position)')
    pylab.draw_if_interactive()

def adaptive_date_ticks(axx,dtrange=None,format = None, formatm = None, nticks = 0, fit = None,label=True,lformat=None,ltime=None,debug = False):
  '''Takes an axis and a time range in seconds and chooses the optimal time display

	ax is the xaxis or yaxis, 
        trange is the time range of the data in seconds, 
            defaults to the existing plot data range
        format is the preferred time format 
            e.g. '%Y' vs '%y' or '%b %d, %Y' vs '%m/%d'
        formatm is the preferred time format for the minor ticks
            e.g. '%Y' vs '%y' or '%b %d, %Y' vs '%m/%d'
        nticks forces the number of ticks, 
            defaults to number already in figure
        axname allows time on y axis (default is x axis)
	fit determines whether the image should:
            fit the data range exactly (fit='exact')
            fit the nearest tick mark outside the data range (fit='tick')
	    use the default fit (default or fit=None)
        label determines whether to add a label to the axis
        lformat allows the user to specify the format of the label
        ltime is the date number to use to generate the label,
            defaults to middle of data range'''
  if dtrange is None:
     bounds = pylab.getp(pylab.getp(axx,'data_interval'),'bounds')
     dtrange = pylab.diff(bounds)
  trange = dtrange * 86400
  if label  and ltime is None:
     ltime = pylab.mean(bounds)
  if nticks == 0:
     nticks = len(axx.get_ticklocs())
  if trange < 60: #image covers less than 1 minute
     bt = int(60/nticks)
     tloc = MinuteLocator()
     tlocm = SecondLocator(bysecond=range(bt,60,bt))
     tfor = ':%S'
     tform = '%M:%S'
     tlabel = '%b %d, %Y %H:%M'
     if debug: print 'seconds'
  elif trange/nticks < 60:
     tloc = MinuteLocator()
     tlocm = SecondLocator(bysecond=range(15,60,15))
     tform = ':%S'
     tfor = '%M:00'
     tlabel = '%b %d, %Y %H:' 
     if debug: print 'half minutes'
  elif trange/nticks <90:
     tloc = MinuteLocator(interval=1)
     tlocm = SecondLocator(bysecond=30)
     tform = ''
     tfor = '%M'
     tlabel = '%b %d, %Y %H:' 
     if debug: print 'minutes'
  elif trange/nticks < 120:
     tloc = MinuteLocator(byminute=range(0,60,2))
     tlocm = SecondLocator(bysecond=(0,30))
     tfor = '%M'
     tform = ''
     tlabel = '%b %d, %Y %H:'
     if debug: print '2 minutes'
  elif trange/nticks < 240:
     tloc = MinuteLocator(interval=5)
     tlocm = MinuteLocator(interval=1)
     tfor = '%H:%M'
     tform = ''
     tlabel = '%b %d, %Y'
     if debug: print '5 minutes'
  elif trange/nticks < 600:
     tloc = MinuteLocator(interval=10)
     tlocm = MinuteLocator(interval=1)
     tform = ''
     tfor = '%H:%M'
     tlabel = '%b %d, %Y'
     if debug: print '10 minutes'
  elif trange/nticks < 1200:
     tloc = HourLocator(interval=1)
     tlocm = MinuteLocator(byminute=range(15,60,15))
     tfor = '%H:%M'
     tform = ':%M'
     tlabel = '%b %d, %Y'
     if debug: print '30 minutes'
  elif trange/nticks < 2400:
     tloc = MinuteLocator(byminute=0)
     tlocm = MinuteLocator(byminute=30)
     tfor = '%H:00'
     tform = ''
     tlabel = '%b %d, %Y'
     if debug: print 'hour'
  elif trange/nticks < 4500:
     tloc = HourLocator(interval=2)
     tlocm = HourLocator(interval=1)
     tform = ''
     tfor = '%H:00'
     tlabel = '%b %d, %Y'
     if debug: print '2 hours'
  elif trange < 86400:
     tloc = HourLocator(byhour=range(0,24,3))
     tlocm = HourLocator()
     tform = ''
     tfor = '%H:00'
     tlabel = '%b %d, %Y'
     if debug: print '3 hours'
  elif trange < 86400*2:
     tloc = HourLocator(byhour=range(0,24,6))
     tlocm = HourLocator(byhour=range(0,24,2))
     tform = ''
     tfor = '%H:00\n%m/%d'
     tlabel = '%b %d, %Y'
     if debug: print '6 hours'
  elif dtrange < 2.5:
     tloc = DayLocator()
     tlocm = HourLocator(byhour=range(6,24,6))
     tfor = '%m/%d'
     tform = '%H:%M'
     tlabel = '%Y'
     if debug: print '1 day/ 6 hours'
  elif dtrange < 3:
     tloc = DayLocator()
     tlocm = HourLocator(byhour=12)
     tfor = '%m/%d'
     tform = '%H:%M'
     tlabel = '%Y'
     if debug: print '1 day/ 12 hours'
  elif dtrange/nticks < 0.5:
     tloc = DayLocator(interval=1)
     tlocm = HourLocator(byhour=range(0,24,12))
     tform = ''
     tfor = '%m/%d'
     tlabel = '%Y'
     if debug: print '1 day'
  elif dtrange < 31*6:
     tloc = MonthLocator()
     if dtrange/nticks < 1:
       interv = 1
       endv = 32
     elif dtrange/nticks < 2:
       interv = 2
       endv = 31
     elif dtrange/nticks < 3:
       interv = 3
       endv = 30
     elif dtrange/nticks < 4:
       interv = 5
       endv =30 
     elif dtrange/nticks < 6:
       interv = 7
       endv=28
     elif dtrange/nticks < 8:
       interv = 10
       endv =30
     else: 
       interv = 15
       endv=30
     tlocm = MonthLocator(bymonthday=range(interv,endv,interv))
     tfor = "%d\n%b '%y"
     tform = '%d'
     tlabel = 'date'
     if debug: print 'month'
  elif dtrange/nticks < 30:
     tloc = MonthLocator(bymonth=range(1,13,2))
     tlocm = MonthLocator(bymonth=range(2,13,2))
     tfor = "%b\n%y"
     tform = "%b"
     tlabel = 'date'
     if debug: print '2 month'
  elif dtrange < 366:
     tloc = MonthLocator(bymonth=range(1,13,6))
     tlocm = MonthLocator(interval=1)
     tfor = "%b\n'%y"
     tform = "%b"
     tlabel = 'date'
     if debug: print '3 month/month'
  elif dtrange < 730:
     tloc = MonthLocator(bymonth=range(1,13,6))
     tlocm = MonthLocator(bymonth=range(1,13,2))
     tfor = "%b\n'%y"
     tform = "%b"
     tlabel = 'date'
  elif dtrange < 365*3:
     tloc = MonthLocator(bymonth=range(1,13,6))
     tlocm = MonthLocator(bymonth=range(1,13,2))
     tfor = "%b\n'%y"
     tform = ""
     tlabel = 'date'
     if debug: print '6 month/2 month'
  elif dtrange < 365*4:
     tloc = MonthLocator(bymonth=1)
     tlocm = MonthLocator(bymonth=7)
     tfor = "%b\n'%y"
     tform = "%b"
     tlabel = 'date'
     if debug: print '6 month/2 month'
  elif dtrange/nticks < 240:
     tloc = YearLocator()
     tlocm = MonthLocator(bymonth=(1,7))
     tform = ""
     tfor = "%Y"
     tlabel = 'Year'
     if debug: print 'year/6 month'
  elif dtrange/nticks < 365:
     tloc = YearLocator(base=2)
     tlocm = MonthLocator(bymonth=(1,7))
     tfor = '%Y'
     tform = ''
     tlabel = 'year'
     if debug: print '12 month'
  elif dtrange/nticks < 580:
     tloc = YearLocator(base=2)
     tlocm = YearLocator()
     tfor = '%Y'
     tform = ''
     tlabel = 'year'
     if debug: print '2 year/year'
  elif dtrange < 20*365:
     tloc = YearLocator(base=5)
     tlocm = YearLocator()
     tfor = '%Y'
     tform = ''
     tlabel = 'year'
     if debug: print '5 year/year'
  elif dtrange < 100*365:
     tloc = YearLocator(base=10)
     tlocm = YearLocator(base=2)
     tfor = '%Y'
     tform = ''
     tlabel = 'year'
     if debug: print '10 year'
  else:
     interv = 2.7*dtrange/(nticks*365.25)
     print interv
     interv = int(numpy.round(pylab.matplotlib.numerix.power(10,int(pylab.matplotlib.numerix.log10(interv)*3)/3.)/(pylab.matplotlib.numerix.power(10,int(pylab.matplotlib.numerix.log10(interv)*3)/3)))*(pylab.matplotlib.numerix.power(10,int(pylab.matplotlib.numerix.log10(interv)*3)/3)*1.))
     print interv
     tloc = YearLocator(base=interv)
     tlocm = YearLocator(base=int(interv/5))
     tfor = '%Y'
     tform = ''
     tlabel = 'year'
     if debug: print 'variable long scale'
  if format is None:
     format = tfor
  if formatm is None:
     formatm = tform
  if lformat is None:
     lformat = tlabel
  axx.set_major_locator(tloc)
  axx.set_minor_locator(tlocm)
  axx.set_major_formatter(DateFormatter(format))
  axx.set_minor_formatter(DateFormatter(formatm))
  if label:
     if ltime is not None: #and lformat.find('%')<0:
        print 'reached here'
        dateform = pylab.DateFormatter(lformat)
        print lformat 
        print pylab.num2date(ltime)
        tlabel = dateform.strftime(pylab.num2date(ltime),lformat)
        print tlabel
     elif lformat.find('%')>0:
        tlabel = 'time'
     #if debug:
      # str = "range %8.4f range per tick %4.4f ticks %d\ntick format %s label %s" %  (dtrange, dtrange/nticks, nticks, tform, tlabel)
      # print str
      # ax.set_title(str) 
     pylab.setp(axx.get_label(),text=tlabel)
  enlarge_allticklines(axx,factor=1.5)
  
  pylab.draw_if_interactive()

def my_date_ticker_factory(span, tz=None, numticks=5, debug = False):
    """
    Create a date locator with numticks (approx) and a date formatter
    for span in days.  Return value is (locator, formatter)


    """

    if span==0: span = 1/24.

    minutes = span*24*60
    minutes30 = span*24*30
    minutes10 = span*24*10
    hours  = span*24
    hours2  = span*12
    hours3  = span*8
    hours6  = span*4 
    days   = span
    days2   = span/2
    days3   = span/3
    weeks  = span/7.
    halfmonths = span/14. # approx
    months = span/31. # approx
    m3 = span/(3*31.) # approx
    m6 = span/182. # approx
    years  = span/365.

    if years>numticks:
    
        locator = YearLocator(int(years/numticks), tz=tz)  # define
        fmt = '%Y'
        if debug: 
           print 'years %d > numticks %d, year interval = %d' % (years,numticks,int(years/numticks))
    elif months>numticks:
        locator = MonthLocator(interval=int(months/numticks),tz=tz)
        fmt = '%b %Y'
        if debug: 
           print 'months %d > numticks %d, month interval = %d' % (months,numticks,1)
    elif weeks>numticks:
        locator = MonthLocator(bymonthday=range(1,28,int(weeks/numticks)*7),tz=tz)
        fmt = '%b %d'
        if debug: 
           print 'weeks %d > numticks %d, week interval = %d' % (weeks,numticks,1)
    elif days>numticks:
        locator = DayLocator(interval=int(math.ceil(days/numticks)), tz=tz)
        fmt = '%b %d'
        if debug: 
           print 'days %d > numticks %d, day interval = %d' % (days,numticks,int(math.ceil(days/numticks)))
    elif hours>numticks:
        locator = HourLocator(interval=int(math.ceil(hours/numticks)), tz=tz)
        fmt = '%H:%M\n%b %d'
        if debug: 
           print 'hours %d > numticks %d, hour interval = %d' % (hours,numticks,int(math.ceil(hours/numticks)))
    elif minutes>numticks:
        locator = MinuteLocator(interval=int(math.ceil(minutes/numticks)), tz=tz)
        fmt = '%H:%M:%S'
        if debug: 
           print 'minutes %d > numticks %d, minute interval = %d' % (minutes,numticks,int(math.ceil(minutes/numticks)))
    else:
        locator = MinuteLocator(tz=tz)
        fmt = '%H:%M:%S'
        if debug: 
           print 'minutes %d < numticks %d, minute interval = %d' % (minutes,numticks,1)


    formatter = DateFormatter(fmt, tz=tz)
    return locator, formatter
def enlarge_minorticklines(axx,length=2):
        '''modify the tick lines in place'''
        lines = []
        for tick in axx.minorTicks:
            pylab.setp(tick.tick1line,ms=length)
            pylab.setp(tick.tick2line,ms=length)
def enlarge_allticklines(axx,factor=2):
        '''modify all tick lines in place'''
        lines = []
        for tick in axx.minorTicks:
            xt = pylab.getp(tick.tick1line,'ms')
            pylab.setp(tick.tick2line,ms=xt*factor)
            xt = pylab.getp(tick.tick2line,'ms')
            pylab.setp(tick.tick1line,ms=xt*factor)
        for tick in axx.majorTicks:
            xt = pylab.getp(tick.tick1line,'ms')
            pylab.setp(tick.tick2line,ms=xt*factor)
            xt = pylab.getp(tick.tick2line,'ms')
            pylab.setp(tick.tick1line,ms=xt*factor)


def drawcoastline(bmap, coastsegs, **kwargs):
  bmap.coastsegs = coastsegs
  bmap.drawcoastlines(**kwargs)

def plot_from_url(url, ax = None, suffix = '.png', plotflag = True):
  if ax is None:
    ax = pylab.gca()
  fname = tempfile.mktemp(suffix)
  im = urllib.urlretrieve(url,fname)
  d = pylab.imread(fname)
  if plotflag:
    ax.imshow(d)
    ax.set_axis_off()
    pylab.draw_if_interactive()
  return d

def figure_title(f,titlestr, **kwargs):
  f.text(0.5,0.9,titlestr,ha='center',**kwargs)
