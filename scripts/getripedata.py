from warriorpy.shorthand import diriofile as df
from warriorpy.ripe import processing as p
from warriorpy.shorthand import linereader
from pymongo import MongoClient
import json
import logging
import logging.config

##################################################################
#                           LOGGING
##################################################################

logging.config.fileConfig('logging.conf', disable_existing_loggers=False)

# create logger
logger = logging.getLogger(__name__)
logger.debug(__name__+"logger loaded")

##################################################################
#                        GLOBAL AND SETUP
##################################################################

# paths
basedir = df.getdir(__file__)+'../'
rawdirlist = df.getlines(basedir+'state/datapaths.list')
datafiles = df.listfiles(basedir+rawdirlist[0], fullpath=True)

# database setup
mclient = MongoClient()
db = mclient.veracity
coll = db.ripe_meas30002_may2017

##################################################################
#                           CODE
##################################################################

def build_db(fnum):
    lr = linereader.linereader(datafiles[fnum])
    ind = 0
    while lr.gotmore:
        strline = lr.getnext()
        data = json.loads(strline)
        if p.isbadquery(data):
            continue
        try:
            parsed = p.parse_dns_json(data)
        except KeyError as e:
            logger.error('bad key: '+str(data))
            continue

        parsed['ind'] = ind
        ind += 1
        coll.insert_one(parsed)
        if ind % 50 == 0:
            print "entries so far: "+str(ind)
    print ind

def append_db(fnum):
    lr = linereader.linereader(datafiles[fnum])
    lastentry = coll.find().sort({'ind':-1}).limit(1)
    ind = lastentry['ind']
    print ind
    while lr.gotmore:
        strline = lr.getnext()
        data = json.loads(strline)
        if p.isbadquery(data):
            continue
        try:
            parsed = p.parse_dns_json(data)
        except KeyError as e:
            logger.error('bad key: '+str(data))
            continue

        parsed['ind'] = ind
        ind += 1
        coll.insert_one(parsed)
        if ind % 50 == 0:
            print "entries so far: "+str(ind)
    print ind


build_db(0)

