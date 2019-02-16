import psycopg2
from numpy import *
import re
import pdb

class oe_db:
  '''a class for extracting basic data from the ourestuaries db'''

  def ___init___(self):
    pass

  def connect(self, host, dbname, user, pswd):
    self.host = host
    self.dbname = dbname
    self.user = user
    self.pswd = pswd

    self.srid = 32036   #this didn't seem to get initialized through __init__

#    self.conn = psycopg2.connect("dbname='"+dbname+"' user='"+user+"' host='"+host+"' password='"+pswd+"'")
    self.conn = psycopg2.connect("dbname="+dbname+" user="+user+" host="+host+" password="+pswd)
    self.curs = self.conn.cursor()


  def get_river_ith(self, river_names):
    ith = []
    for i in range(len(river_names)):
      execute_str = "SELECT fluxth_index FROM river_input WHERE id == '"+river_names[i]+"'"
      self.curs.execute(execute_str)
      rows = self.curs.fetchall()
      ith.append(rows[0])
    return ith


  def get_points_by_eid(self, eid):
    execute_str = "SELECT name, AsText(points), AsText(ac_norms), upstream_rivers FROM class_shapes WHERE type_name='point' and eid='"+eid+"'"
    self.curs.execute(execute_str)

    rows=self.curs.fetchall()
    npts = len(rows)
    names = []
    xy=zeros((npts, 2), float)
    norms=zeros((npts,2), float)
    rivers = []
    for i in range(npts):
      xy[i] = self.parse_pt(rows[i][1])
      names.append(rows[i][0])      
      if (rows[i][2]!=None):
        norms[i] = self.parse_pt(rows[i][2])
      else:
        norms[i] = None
      rivers.append(rows[0][3][0].split(','))
      
    return (names, xy, norms, rivers)


  def get_point(self, name, runid, username):
    execute_str = "SELECT name, AsText(points), AsText(ac_norms), upstream_rivers FROM class_shapes WHERE name='"+name+"' AND runid='"+runid+"' AND type_name='point' AND username='"+username+"'"
    self.curs.execute(execute_str)
    rows=self.curs.fetchall()
    pt = self.parse_pt(rows[0][1])
    name = rows[0][0]
    norms = self.parse_pt(rows[0][2])
    rivers = rows[0][3]
    return (name, pt, norms, rivers)


  def set_point(self, pt, norms, name, rivers, eid, runid, username):
    pt_text = "GeomFromText('POINT("+str(pt[0])+" "+str(pt[1])+")')"
    norms_text = "GeomFromText('POINT("+str(norms[0])+" "+str(norms[1])+")')"

    rivers_text = "{"
    for riv in rivers:
      rivers_text = rivers_text+riv+", "
    rivers_text = rivers_text.rstrip(", ")
    rivers_text = rivers_text+"}"
      
    execute_str = "SELECT name FROM class_shapes WHERE name='"+name+"' AND eid='"+eid+"' AND runid='"+runid+"' AND type_name='point' AND username='"+username+"'"
    self.curs.execute(execute_str)
    rows=self.curs.fetchall()
    if (len(rows)>0):
      execute_str = "UPDATE class_shapes SET points="+pt_text+", ac_norms="+norms_text+", upstream_rivers='"+rivers_text+"' WHERE name='"+name+"' AND eid='"+eid+"' AND runid='"+runid+"' AND type_name='point' AND username='"+username+"'"
    else:
      execute_str = "INSERT INTO class_shapes (runid, username, name, type_name, points, eid, ac_norms, upstream_rivers) VALUES ('"+runid+"', '"+username+"', '"+name+"', 'point', "+pt_text+", '"+eid+"', "+norms_text+", '"+rivers_text+"')"

#    print(execute_str)

    pdb.set_trace()
    self.curs2 = self.conn.cursor()

    err = self.curs2.execute(execute_str)

  

  def parse_pt(self, text):
    xy = re.findall("\d+\.*\d*", text)
    return array([float(xy[0]), float(xy[1])])


  def parse_linestring(self, text):
    tpts = text.split(',')
    pts = zeros((len(tpts), 2), float)
    for i in range(len(tpts)):
      pts[i,:] = self.parse_pt(tpts[i])
    return pts


  def set_transect(self, pts, norms, name, rivers, eid, runid, username):
    pts_text = "GeomFromText('LINESTRING("
    for i in range(pts.shape[0]):
      pts_text = pts_text+str(pts[i][0])+" "+str(pts[i][1])+","
    pts_text = pts_text.rstrip(",")      
    pts_text = pts_text + ")', "+str(self.srid)+")"


    norms_text = "GeomFromText('LINESTRING("
    for i in range(norms.shape[0]):
      norms_text = norms_text+str(norms[i][0])+" "+str(norms[i][1])+","
    norms_text = norms_text.rstrip(",")      
    norms_text = norms_text + ")', "+str(self.srid)+")"

    rivers_text = "{"
    for riv in rivers:
      rivers_text = rivers_text+riv+", "
    rivers_text = rivers_text.rstrip(", ")
    rivers_text = rivers_text+"}"
      
    execute_str = "SELECT name FROM class_shapes WHERE name='"+name+"' AND eid='"+eid+"' AND runid='"+runid+"' AND type_name='transect' AND username='"+username+"'"
    self.curs.execute(execute_str)
    rows=self.curs.fetchall()
    if (len(rows)>0):
      execute_str = "UPDATE class_shapes SET points="+pts_text+", ac_norms="+norms_text+", upstream_rivers='"+rivers_text+"' WHERE name='"+name+"' AND eid='"+eid+"' AND runid='"+runid+"' AND type_name='transect' AND username='"+username+"'"
    else:
      execute_str = "INSERT INTO class_shapes (runid, username, name, type_name, points, eid, ac_norms, upstream_rivers) VALUES ('"+runid+"', '"+username+"', '"+name+"', 'transect', "+pts_text+", '"+eid+"', "+norms_text+", '"+rivers_text+"')"

#    print(execute_str)

#    pdb.set_trace()
    self.curs2 = self.conn.cursor()

    err = self.curs2.execute(execute_str)


  def get_shape_info_by_eid(self, eid, type_name):
    execute_str = "SELECT name, runid, username FROM class_shapes WHERE type_name='"+type_name+"' and eid='"+eid+"'"
    self.curs.execute(execute_str)
    rows=self.curs.fetchall()
    npts = len(rows)
    names = []
    runids = []
    unames = []

    for i in range(npts):
      names.append(rows[i][0])
      runids.append(rows[i][1])
      unames.append(rows[i][2])

    return (names, runids, unames)

    
  def get_trans_or_cs(self, name, runid, username, type_name):
    execute_str = "SELECT AsText(points), AsText(ac_norms), upstream_rivers, eid FROM class_shapes WHERE type_name='"+type_name+"' and runid='"+runid+"' and name='"+name+"' and username='"+username+"'"
    self.curs.execute(execute_str)
    rows=self.curs.fetchall()
    rivers = []

    pts = self.parse_linestring(rows[0][0])
    if rows[0][1]==None:
      ac_norms = None
    else:
      ac_norms=self.parse_linestring(rows[0][1])
      
#    rvrs = re.findall("'.*'", rows[0][2][0])
#    for i in range(len(rvrs)):
#      arvr = rvrs[i].strip("'")
#      rivers.append(arvr)
    rivers = rows[0][2][0].split(',')
      
    return (pts, ac_norms, rivers, rows[0][3])
