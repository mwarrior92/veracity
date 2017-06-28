from warriorpy.shorthand import diriofile as df
from pymongo import MongoClient
import json
from netaddr import IPNetwork as CIDR
from netaddr import IPAddress as IP
import numpy as np
import veracity_vector as vv
from sklearn.cluster import DBSCAN
from sklearn import metrics

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

domains = df.getlines(basedir+'state/sites.csv')

# database setup
mclient = MongoClient()
db = mclient.veracity
coll = db.ripe_meas30002_may2017

ag_coll = db.ripe_meas30002_answergroups

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


def get_answer_groups(t, duration=30000, fmt=None):
    print "getting window"
    window = vv.get_window(t, duration, domains[:3])
    print "converting window to dict"
    dd = vv.window_to_dict(window)
    X = list()
    indl = list()
    # list of indices
    print "creating array"
    for probe in dd:
        X.append(np.array(vv.dict_to_vector(dd[probe])))
        indl.append(dd[probe]['ind'])
    X = np.array(X)
    logger.warning("performing db scan...")
    dbs = DBSCAN(eps=255).fit(X)
    labels = dbs.labels_
    print "labels: "+str(set(labels))
    return (labels, indl)

