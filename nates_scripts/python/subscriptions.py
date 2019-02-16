import sys
import cmop
import cmop.db as db

tables = [
    "cruise.ctdcast", "cruise.castobservation", 
    "cruise.tsg", "cruise.adcpconfiguration", 
    "cruise.adcpping", "cruise.adcpobservation"
]

class SubscriptionBroker:
  prefix = "unsent_to_"
  attrsql = '''
select attname 
  from pg_attribute a, 
       pg_class c,
       pg_namespace s
 where a.attrelid = c.oid
   and c.relnamespace = s.oid
   and s.nspname = '%s'
   and c.relname = '%s'
   and a.attnum > 0
   and a.atttypid != 0
'''

  def __init__(self, fromhost, tohost):
    self.fromhost = fromhost
    self.tohost = tohost
    self.dirtycolumn = self.column(self.tohost)
    self.fromdb = db.DB(user="writer", host=fromhost)
    self.todb = db.DB(user="writer", host=tohost)

  def column(self, host):
    return '%s%s' % (self.prefix,host.replace(".", "_"))

  def Clean(self, table):
    # Remove the subscription column
    # Don't do this as part of normal operation;
    # dropped columns are not permanently removed!
    msg = "SubscriptionBroker: CLEAN CALLED: %s dropping column from table %s on %s"
    cmop.info(msg % (self.tohost, table, self.fromhost))
    mod = '''ALTER TABLE %s DROP COLUMN "%s"'''
    modind = '''DROP INDEX %s_%s'''
    self.fromdb.begin()
    try:
      self.fromdb.execCommand(mod % (table, self.dirtycolumn))
      self.fromdb.execCommand(modind % (table.replace(".", "_"), self.dirtycolumn))
      self.fromdb.commit()
    except:
      self.fromdb.rollback()

  def UnSubscribe(self, table): 
    # new records will not be marked for broadcast,
    # but the subscription column will not be dropped
    # no effect if subscription column does not exist
    if self.SubscriptionColumnExists(table):
      msg = "SubscriptionBroker: %s UNsubscribing from table %s on %s"
      cmop.info(msg % (self.tohost, table, self.fromhost))
      mod = '''ALTER TABLE %s ALTER COLUMN "%s" DROP DEFAULT'''
      self.fromdb.execCommand(mod % (table, self.dirtycolumn))

  def Subscribed(self, table):
    checksql = self.attrsql + '''
   and attname = '%s'
   and atthasdef
'''
    sch, tab = table.split(".")
    check = checksql % (sch, tab, self.dirtycolumn)
    rs = self.fromdb.execQuery(check)
    return rs != []

  def Subscribe(self, table):
    # all new records will be marked for broadcast to current host
    if not self.Subscribed(table):
      msg = "SubscriptionBroker: %s subscribing to table %s on %s"
      cmop.info(msg % (self.tohost, table, self.fromhost))
      self.Configure(table)
      default = '''ALTER TABLE %s ALTER COLUMN "%s" SET DEFAULT True'''
      self.fromdb.execCommand(default % (table, self.dirtycolumn))

  def AttributeExists(self, db, table, attr):
    checksql = self.attrsql + '''
   and attname = '%s'
'''
    sch, tab = table.split(".")
    check = checksql % (sch, tab, attr)

    rs = db.execQuery(check)
    return rs != []

  def SubscriptionColumnExists(self, table):
    return self.AttributeExists(self.fromdb, table, self.dirtycolumn)

  def Configure(self, table):
    # Add a column on source table representing this subscription
    # also add an index
    if not self.SubscriptionColumnExists(table):
      if not self.tohost.Primarykey(table): 
        msg = "Table %s does not have a primary key defined; cannot transfer."
        raise TypeError(msg % (table,))

      msg = "SubscriptionBroker: Configuring %s subscribing to table %s on %s"
      cmop.info(msg % (self.tohost, table, self.fromhost))
      self.fromdb.begin()
      mod = '''ALTER TABLE %s ADD COLUMN "%s" bool; '''
      ind = '''CREATE INDEX "%s_%s" ON %s("%s")'''
      vac = '''VACUUM ANALYZE %s'''
      try:
        self.fromdb.execCommand(mod % (table, self.dirtycolumn))
        self.fromdb.execCommand(ind % (table.replace(".", "_"), 
                                       self.dirtycolumn, 
                                       table,
                                       self.dirtycolumn))
        self.fromdb.commit()
        self.fromdb.execCommand(vac % (table,))
      except:
        self.fromdb.rollback()
        raise

  def Transfer(self, table, filter="True", orderby=None, limit=None):
     if not self.SubscriptionColumnExists(table):
       return

     getattrs = self.fromdb.Attributes(table, self.prefix)
     setattrs = getattrs[:]

     if self.AttributeExists(self.todb, table, self.column(self.fromhost)):
       setattrs.append('"%s"' % (self.column(self.fromhost),))

     getattrsstr = ", ".join(getattrs)
     setattrsstr = ", ".join(setattrs)

     get = '''SELECT %s FROM %s WHERE "%s" AND %s''' 
     if orderby: get += " ORDER BY %s" % (orderby,)
     if limit: get += " LIMIT %s" % (limit,)
     get = get % (getattrsstr, table, self.dirtycolumn, filter)

     temptable = '''
CREATE TEMP TABLE %s_job AS (
     %s
)
'''
     temptable = temptable % (table, get)
     
     set = 'UPDATE %s SET "%s" = False FROM %s '
     def cond(t):
       return " AND ".join([db.equals(a,v) for a,v in zip(getattrs,t) if v != 'NULL'])
     set = set % (table, self.dirtycolumn, "%s")

  def Transfer2(self, table, filter="True", orderby=None, limit=None):
     if not self.SubscriptionColumnExists(table):
       return

     # get the data attributes
     getattrs = self.fromdb.Attributes(table, self.prefix)
     setattrs = getattrs[:]
    
     # if the reverse subscription is set, mark the new tuple as sent
     if self.AttributeExists(self.todb, table, self.column(self.fromhost)): 
       setattrs.append('"%s"' % (self.column(self.fromhost),))
       val = ['False']
     else:
       val = []
   
     getattrsstr = ", ".join(getattrs)


  def ValuesClause(self, xfertable, extravals=[]):
       '''
Read all tuples from the given table and prepare a VALUES clause from the results
Returns a tuple of the values clause as a string and the number of tuples involved.
'''
       select = '''SELECT * FROM %s''' % (xfertable,)
       cmop.debug("Reading tuples to transfer.", 8)
       rs = self.fromdb.execQuery(select)
       cnt = len(rs)
       if not rs:
         cmop.info("No tuples to transfer.")
         return "()", 0

       def preprow(r):
         vals = ["%s" % (db.quote(a),) for a in r]
         vals += extravals
         return "(%s)" % (", ".join(vals),)

       values = ", ".join([preprow(r) for r in rs])

       return values, cnt

  def Transfer(self, table, filter="True", orderby=None, limit=None):
     cmop.debug("Transferring %s" % (table,))

     if not self.SubscriptionColumnExists(table):
       cmop.debug("Subscription does not appear to be active; no subscription column found on %s" % (table,))
       return
   
     # get the data attributes
     keys = self.fromdb.PrimaryKey(table)
     if not keys: 
       raise TypeError("Table %s does not have a primary key defined; cannot transfer." % (table,))

     keyattrsstr = ", ".join(keys)
     getattrs = self.fromdb.Attributes(table, self.prefix)
     setattrs = getattrs[:]

     # if the reverse subscription is set, mark the new tuple as sent
     if self.AttributeExists(self.todb, table, self.column(self.fromhost)): 
       setattrs.append('"%s"' % (self.column(self.fromhost),))
       val = ['False']
     else:
       val = []

     getattrsstr = ", ".join(getattrs)
     setattrsstr = ", ".join(setattrs)

     get = '''SELECT %s FROM %s WHERE "%s" AND %s''' 
     if orderby: get += " ORDER BY %s" % (orderby,)
     if limit: get += " LIMIT %s" % (limit,)
     get = get % (getattrsstr, table, self.dirtycolumn, filter)

     xfertable = table.replace(".", "_") + "_xfer"

     create = '''CREATE TEMP TABLE %s AS  (%s)''' % (xfertable, get)

     insert = '''INSERT INTO %s (%s) VALUES %s''' % (xfertable, setattrsstr, "%s")

     drop = "DROP TABLE %s" % (xfertable,)

     def equal(col):
         return "%s.%s = %s.%s" % (table, col, xfertable, col)

     joincond = " AND ".join([equal(col) for col in keys])
     update = '''
UPDATE %s SET %s = False
  FROM %s
 WHERE %s
''' % (table, self.column(self.tohost), xfertable, joincond)

     createremote = '''
CREATE TEMP TABLE %s AS (
SELECT * FROM %s LIMIT 0
)''' % (xfertable, table)


     qualattrs = ", ".join(["%s.%s" % (xfertable, a) for a in setattrs])
     merge = '''
INSERT INTO %s (%s) ( 
  SELECT %s
    FROM %s LEFT JOIN %s
         ON (%s) 
   WHERE %s.%s IS NULL
)
''' % (table, setattrsstr, qualattrs, xfertable, table, joincond, table, keys[0])
     
     try:
       self.fromdb.begin()
       self.todb.begin()

       cmop.debug("Creating xfer table %s on %s" % (xfertable, self.fromhost))
       self.fromdb.execCommand(create)
     
       cmop.debug("Creating xfer table %s on %s" % (xfertable, self.tohost))
       self.todb.execCommand(createremote)

       cmop.debug("Extracting values from %s on %s" % (xfertable, self.fromhost))
       values, cnt = self.ValuesClause(xfertable, val)
 
       if cnt > 0: 
         cmop.debug("Inserting tuples on %s to %s" % (self.tohost,xfertable))
         self.todb.execCommand(insert % (values,))

         cmop.debug("Merging tuples on %s into %s" % (self.tohost,table))
         self.todb.execCommand(merge)
       
         cmop.debug("Marking tuples as sent on %s for %s" % (self.fromhost,table))
         self.fromdb.execCommand(update)
 
         cmop.debug("Dropping temp tables")
         self.fromdb.execCommand(drop)
         self.todb.execCommand(drop)
       
         cmop.debug("Committing on %s" % (self.tohost,))
         self.todb.commit()

         cmop.debug("Committing on %s" % (self.fromhost,))
         self.fromdb.commit()

     except:
       self.todb.rollback()
       self.fromdb.rollback()
       raise

     cmop.info("Transferred %s %s tuples from %s to %s (may include dupes)" % (cnt, table, self.fromhost, self.tohost))

def main():
  #cmop.DebugOn()
  broker = SubscriptionBroker('localhost', 'cdb02')
  
  for t in tables:
    broker.Transfer(t, limit=100)

if __name__ == '__main__':
  main()
