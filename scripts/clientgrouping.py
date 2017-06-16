from warriorpy.shorthand import diriofile as df
from pymongo import MongoClient
import json
from netaddr import IPNetwork as CIDR
from netaddr import IPAddress as IP

##################################################################
#                           LOGGING
##################################################################
import logging
import logging.config

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

def get_mask_group(mask):
    fmtmask = str(ipp.ip2int(mask))
    distinct_ip_list = coll.distinct("probe_ip",
            {"$where": "(this.probe_ip & "+fmtmask+")=="+fmtmask})
    print "num matched IPs in mask "+mask+": "+str(len(distinct_ip_list))
    return distinct_ip_list


def get_asn_group(asn):
    distinct_ip_list = coll.distinct("probe_ip", {"asn_v4": asn})
    return distinct_ip_list


def get_allocation_info(owner):
    # maybe move this to processing.py so it happens up front?
    pass
