try: 
  import pg
except: 
  import _pg
  pg = _pg
import sys
import re
import popen2, time

import ConfigParser

import cmop
class loggerStub:
  '''A logging abstraction to make it easier to decouple this module from the cmop module.'''
  def debug(self, s, level=10):
    cmop.debug(s, level)

  def info(self, s):
    cmop.info(s)

  def profile(self, s):
    cmop.profile(s)

logger = loggerStub()

DELIM = ';'

class SQLError(Exception):
  pass

class DuplicateError(SQLError):
  pass

def commas(xs):
  return ', '.join([str(x) for x in xs])

defaults = cmop.config
class defaults:
    cfp = ConfigParser.ConfigParser()
    conffile = "/usr/local/cmop/etc/server.ini"
    try:
      cfp.read([conffile])
      hostname = cfp.get("dbserver","dbserver")
      dbname = cfp.get("dbserver","dbname")
      user = cfp.get("dbserver", "dbuser")
      password = cfp.get("dbserver", "dbpasswd")
    except:
      logger.debug("Error using config file '%s' ...using defaults" % (conffile,))
      dbname = 'cmop'
  #    hostname = 'cdb02.stccmop.org'
      hostname = 'localhost'
      user = 'reader'
      password = ''

def escape(v):
  return v.replace("'","''")

expr = re.compile("\s*('.*')|(\(.*\))\s*")
def quote(v):
  # format strings for postgres
  if type(v) == str:
    # add single quotes to strings, 
    # unless they are already there
    if not expr.match(v):
      return "'%s'" % (escape(v),)
    else:
      return v
  if v == None:
    return 'NULL'
  else:
    return v

def equals(a,v):
  # build an SQL equailty expression 
  # for attribute = value
  # floats are not necessarily equal
  if type(v) == float:
    return "trunc(%s::numeric,2)::numeric=trunc(%s,2)::numeric" % (a,v)
  else:
    return "%s=%s" % (a,v)

class DB:
    ''' 
    A Database abstraction class
    --works with postgres only, but exists to protect higher-level 
      apps from a db software change
    --Provides direct query (with results) and command (no results)
    --Provides abstraction of logging, configuration, and errorhandling for the database

    make sure user has permissions in pg_hba.conf on hostname
    '''   

    def __init__(self, dbname=defaults.dbname, host=defaults.hostname, \
                       user=defaults.user, password=defaults.password):
      self.dbname= dbname
      self.user= user
      self.host= host
      self.password= password
      self.dbconn = None
      self.connect(self.dbname,self.host,self.user,self.password)
      self.connect()
      self.suppressquery = False

    def execCommand(self, qry):
      '''Executes a SQL statement, ignoring the results'''
      #self.connect()
      logger.debug('Executing command (%s): \n%s' % (self.host,qry), 7)
      try:
        if not self.suppressquery:
          result = self.dbconn.query(qry)
      except pg.ProgrammingError, e:
        if "duplicate" in str(e):
          raise DuplicateError('%s on %s: %s \n SQL="%s"' % (self.dbname, self.host, e, qry))
        else:
          raise SQLError('%s on %s: %s \n SQL="%s"' % (self.dbname, self.host, e, qry))

    def close(self):
      if self.dbconn:
        self.dbconn.close()
        self.dbconn = None

    def reconnect(self): 
       self.close()
       self.connect(self.dbname,self.host,self.user,self.password)

    def execQuery(self, qry, withfields=False):
      '''
Executes a SQL statement, returning the results
      '''
      self.connect()
      logger.debug('Executing query (%s): \n%s' %(self.host,qry), 7)
      qry = qry.strip()
      if not len(qry):
        return []
      # any error handling we should do?
      #print qry
      #x = raw_input()
      #if x == '':
      #  x = qry
      #self.dbconn = None
      #self.dbconn.reset()
      #self.connect()
      try:
        self.response = self.dbconn.query(qry)
      except pg.ProgrammingError, e:
        raise SQLError('%s on %s: %s \n SQL=%s' % (self.dbname, self.host, e, qry.__repr__()))

      if self.response: 
        if withfields:
          return (self.response.listfields(), self.response.getresult())
        else: 
          return self.response.getresult()

    def connect(self, dbname=defaults.dbname, hostname=defaults.hostname, 
                      user=defaults.user, password=defaults.password):
      if not self.dbconn:
        self.dbconn = pg.DB (dbname, \
                                  hostname, \
                                  -1, \
                                  None, \
                                  None, \
                                  user, \
                                  password)
     
    def PrimaryKey(self, table):
      pair = table.split(".")
      if len(pair) == 1:
        sch, tab = 'public', pair[0]
      else:
        sch, tab = pair

      keysql = '''
select a.attname
  from pg_attribute a, pg_constraint c, pg_class t, pg_namespace s
 where contype = 'p'
   and conrelid = t.oid
   and t.relnamespace = s.oid 
   and a.attnum = ANY (c.conkey)
   and a.attrelid = t.oid
   and s.nspname = '%s'
   and t.relname = '%s'
'''
      rs = self.execQuery(keysql % (sch, tab))
      return [a[0] for a in rs]


    def Attributes(self, table, namefilter=None):
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
      pair = table.split(".")
      if len(pair) == 1:
        sch, tab = 'public', pair[0]
      else:
        sch, tab = pair

      if namefilter:
        attrs = (attrsql + ''' 
    and attname NOT LIKE '%s%%'
 ''') % (sch, tab, namefilter)

      rs = self.execQuery(attrs)
      return [a[0] for a in rs]

    def QuoteAsString(self, s):
      return "'%s'" % (s,)
    
    def appendTo(self, tablename, qry):
      '''
Append a query's results to a table, creating it if it doesn't exist.
     '''
      if self.checkTable(tablename): 
        qry = "INSERT INTO %s (%s)"
      else:
        qry = "CREATE TABLE %s AS (%s);" % (tablename, qry)

      self.execCommand(qry)

    def TupleExists(self, tablename, attrs, vals):
      conds = []
      # floating point numbers cannot be safely compared for equality
      for a, v in zip(attrs, vals):
        conds.append(equals(a, v))
      condition = " AND ".join(conds)
      sql = '''SELECT exists (SELECT %s FROM %s WHERE %s)''' % (attrs[0], tablename, condition)
      s = time.time()
      rs = self.execQuery(sql)
      logger.debug("Tuple Exists: %s (%s)" % (time.time() - s,rs[0][0]), 6)
      return rs[0][0] == 't'
   
    def IdempotentInsert(self, tablename, attrs, vals):
      exists = self.TupleExists(tablename, attrs, vals)
      if not exists:
        try:
          self.InsertTuple(tablename, attrs, vals)
        except DuplicateError: 
          pass 
        return True
      else:
        tup = zip(attrs,vals)
        logger.debug('IdempotentInsert: Tuple exists. (%s)' % (tup,), 6)
        return False

    def Upsert(self, tablename, keyattrs, keyvals, dataattrs, datavals):

      rs = self.TupleExists(keyattrs, keyvals)
      if not rs:
        self.InsertTuple(tablename, keyattrs+dataattrs, keyvals+datavals)
      else:
        upd = '''UPDATE %s SET %s WHERE %s'''
        setters = ", ".join(["%s=%s" % x for x in zip(dataattrs, datavals)])
        sql = upd % (tablename, setters, condition)
        self.execCommand(sql)

    def InsertTuple(self, tablename, attributes, values):
      assert(len(attributes) == len(values))
      ins = '''INSERT INTO %s (%s) VALUES (%s)'''
      sql = ins % (tablename, commas(attributes), commas(values))
      self.execCommand(sql)

    def InsertQry(self, tablename, qry):
      insert = '''INSERT INTO %s (%s)''' % (tablename, qry)
      self.execCommand(insert)

    def createTableAs(self,tablename, qry):
      '''
Drops and recreates a table based on a query's results.
     '''
      create = "CREATE TABLE %s AS (%s);" % (tablename, qry)

      if self.checkTable(tablename): 
        drop = "DROP TABLE %s; " % tablename
        create = drop + create

      self.execCommand(create)

    def dropTable(self, name):
      if self.checkTable(name):
        drop = "DROP TABLE %s CASCADE;" % name
        self.execCommand(drop)

    def begin(self):
      self.execCommand("BEGIN TRANSACTION")

    def rollback(self):
      self.execCommand("ROLLBACK")

    def commit(self):
      self.execCommand("COMMIT")

    def get_attnames(self, table):
      return self.dbconn.get_attnames(table)

    def checkTable(self, name):
      '''
returns None if table <name> does not exist
      '''
      check = "select relname from pg_class where relname = '%s'" % name
      result = self.execQuery(check)
      return result
    
    def Union(self, queries):
      return "\n UNION \n".join([q.replace(";", "") for q in queries])

    def createTable(self, name, fields, types, key=None, oids=False):
      '''
idempotent table creator.  Returns False if table already existed.
fields is a sequence of string column names and 
types is sequence of string type names. 
      '''
      exists = self.checkTable(name)
      if not exists:
        sql = "CREATE TABLE " + name + "("
        fieldstypes = ['%s %s' % ft for ft in zip(fields, types)]
        sql = sql + ','.join(fieldstypes)
        if key:
          sql = sql + ', PRIMARY KEY (%s)' % (','.join(key),)
        sql = sql + ')'
        if oids:
          sql = sql + ' WITH OIDS'
      
        self.execCommand(sql)
        return True
      else:
        return False

if __name__ == '__main__':
  print logger.debug("main")
