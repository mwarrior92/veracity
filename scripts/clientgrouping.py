from warriorpy.shorthand import diriofile as df
from pymongo import MongoClient
import json
from netaddr import IPNetwork as CIDR
from netaddr import IPAddress as IP
import numpy as np
import veracity_vector as vv
from warriorpy.net_tools import ipparsing as ipp
from sklearn.cluster import DBSCAN
from sklearn import metrics
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist

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
plotsdir = df.rightdir(basedir+"plots/")

domains = df.getlines(basedir+'state/sites.csv')

# database setup
mclient = MongoClient()
db = mclient.veracity
coll = db.m30002_may17

ag = db.answergroups

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


def get_ipstructs(t, duration=30000, mask=32, fmt=None, country=None):
    print "getting window"
    if fmt is None:
        window = vv.get_window(t, duration, domains, country)
    elif type(fmt) is int:
        window = vv.get_window(t, duration, domains[:fmt], country)
    else:
        window = vv.get_window(t, duration, fmt, country)
    print "converting window to dict"
    dd = vv.window_to_dict(window)
    X = list()
    indl = list()
    # list of indices
    print "creating array"
    fmtmask = ipp.make_v4_prefix_mask(mask)
    wd = vv.get_weighting(t, duration, mask, fmt, country)
    for probe in dd:
        vec = vv.dict_to_ipstruct(dd[probe], fmtmask, weight=wd)
        X.append(vec)
        indl.append(dd[probe]['ind'])
    return np.array(X)


def get_hamming_distance(t, duration=30000, mask=32, fmt=None,
        method="average", country=None):
    X = get_ipstructs(t, duration, mask, fmt, country)
    dm = np.zeros((len(X) * (len(X) - 1)) // 2, dtype=np.double)
    k = 0
    for i in xrange(0, len(X)-1):
        for j in xrange(i + 1, len(X)):
            dm[k] = 1.0 - vv.distance_metric(X[i], X[j])
            k = k + 1
            if k % 1000 == 0:
                print k
    Z = linkage(dm, method)
    c, coph_dists = cophenet(Z, dm)
    print " Cophenetic Correlation Coefficient: "+str(c)
    plt.figure(figsize=(15, 10))
    plt.xlabel('sample index')
    plt.ylabel('distance')
    dendrogram(
        Z,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
    )
    plt.savefig(plotsdir+'hamming.pdf', bbox_inches='tight')

