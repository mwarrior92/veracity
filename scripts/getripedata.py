from warriorpy.shorthand import diriofile as df
from warriorpy.ripe import processing as p
from warriorpy.ripe import ripe_api as ra
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
coll = db.m30002_may17_full
probe_cache = db.probe_data

##################################################################
#                           CODE
##################################################################


finished = False
# NOTE: this is not designed to be run in parrallel; multiple instances of this
# code will mess up the db entry pattern
def append_db(fnum):
    print "fnum is: "+str(fnum)
    lr = linereader.linereader(datafiles[fnum])
    lastentry = coll.find().sort([('ind',-1)]).limit(1)
    if lastentry.count() > 0:
        ind = lastentry[0]['ind'] + 1
        linenum = lastentry[0]['linenum'] + 1
    else:
        ind = 0
        linenum = 0
    print ind
    pos = 0

    # look for the linenum that corresponds to the first entry of this file
    if linenum > 0:
        while pos < 10 and lr.gotmore:
                strline = lr.getnext()
                pos += 1
                data = json.loads(strline)

                if p.isbadquery(data):
                    continue
                try:
                    logger.debug("parsing json...")
                    parsed = p.parse_dns_json(data)
                    # check if we've started imported this file yet
                    tmp = coll.find({
                        "time":parsed["time"],
                        "probe_id":parsed["probe_id"],
                        "domain": parsed["domain"]
                        }).limit(1)
                    # if we've seen it before, move pos to where file starts
                    if tmp.count() > 0:
                        print tmp[0]
                        print parsed
                        print "--------------------------------------------------------------------"
                        pos = tmp[0]['linenum']-pos+1
                    else:
                        # if we haven't seen it, begin appending from the end
                        pos = linenum
                except KeyError as e:
                    logger.error('bad key: '+str(data))
                # break once we've got the right '0' position for this file
                break
        if pos == 10:
            return

    lr = linereader.linereader(datafiles[fnum])

    print "pos is: " + str(pos)
    sincelastmiss = 0
    while lr.gotmore:
        strline = lr.getnext()
        pos += 1
        if linenum > pos:
            if pos % 1000 == 0:
                logger.debug("current index: "+str(pos))
            continue
        data = json.loads(strline)
        if p.isbadquery(data):
            continue
        try:
            logger.debug("parsing json...")
            parsed = p.parse_dns_json(data)
        except KeyError as e:
            logger.error('bad key: '+str(data))
            continue

        parsed['ind'] = ind
        parsed['linenum'] = linenum
        ind += 1
        linenum += 1
        # add additional info to probe
        # try cache first
        logger.debug("checking cache...")
        tmp = list(probe_cache.find({"probe_id": parsed["probe_id"]}).limit(1))
        if len(tmp) > 0:
            sincelastmiss += 1
            logger.debug("probe info cache hit")
            parsed['country'] = tmp[0]['country']
            if 'asn_v4' in tmp[0]:
                parsed['asn_v4'] = tmp[0]['asn_v4']
            if 'asn_v6' in tmp[0]:
                parsed['asn_v6'] = tmp[0]['asn_v6']
        else:
            print "missed on: "+str(parsed["probe_id"])+"; "\
            +str(sincelastmiss)+" since last incident"
            sincelastmiss = 0
            print "probe info db size: "+str(len(probe_cache.distinct("probe_id")))
            logger.debug("probe info cache miss")
            # if cache miss, do query to get probe info
            b, probe_info = ra.get_probe_info(parsed['probe_id'])
            if b:
                parsed['country'] = probe_info['country_code']
                if 'asn_v4' in probe_info:
                    parsed['asn_v4'] = probe_info['asn_v4']
                if 'asn_v6' in probe_info:
                    parsed['asn_v6'] = probe_info['asn_v6']
                if 'geometry' in probe_info:
                    parsed['geometry'] = probe_info['geometry']
                if 'prefix_v4' in probe_info:
                    parsed['prefix_v4'] = probe_info['prefix_v4']
                if 'prefix_v6' in probe_info:
                    parsed['prefix_v6'] = probe_info['prefix_v6']
            print "inserting..."
            probe_cache.insert_one(parsed)
        coll.insert_one(parsed)
        if ind % 500 == 0:
            print "entries so far: "+str(ind)
    global finished
    finished = True
    logger.info("finished with fnum "+str(fnum))
    return


if __name__ == "__main__":
    for index, datafile in enumerate(datafiles):
        print datafile
        logger.info("datafile: "+datafile)
        finished = False
        while not finished:
            try:
                append_db(index)
                print finished
            except Exception as e:
                logger.error(e)


def update_probe_info(fields):
    cache = list(probe_cache.find())
    for probe in cache:
        if len(set(probe).symmetric_difference(set(fields))) > 0:
            changed = dict()
            b, info = ra.get_probe_info(probe['probe_id'])
            if b:
                for field in fields:
                    if field in info:
                        changed[field] = info[field]
            if len(changed) > 0:
                probe_cache.update(
                        {'probe_id': probe['probe_id']},
                        {"$set": changed})
    return
