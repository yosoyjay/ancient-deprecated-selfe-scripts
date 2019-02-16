import ourestuaries_db
import numpy


oe = ourestuaries_db.oe_db()
oe.connect("cdb01.stccmop.org", "ourestuaries", "hyde", "~good7poop")

(names, pts) = oe.get_points("yaqalsea")

npts = pts.shape[0]

print "npts: "+str(npts)+"\n"

for i in range(npts):
  print names[i]+": "+str(pts[i])
