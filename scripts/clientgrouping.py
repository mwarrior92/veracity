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

pd = db.probe_data

##################################################################
#                           CODE
##################################################################


def get_svl(t, duration=30000, mask=32, fmt=None, country_set=None,
        oddballs=False):
    window = vv.get_window(t, duration, country_set=country_set)
    dd = vv.window_to_dicts(window)
    anssets = vv.get_answer_space_dict(dd)
    sl = vv.sort_sites(anssets)
    fmt = vv.transform_fmt(fmt, sl)

    # remove any domains that only have 1 IP (since all nodes will see the
    # same thing)
    for dom in fmt:
        if len(anssets[dom]) < 2 or ('google' in dom and dom != 'google.com.'):
            del anssets[dom]
    fmt = list(set(anssets.keys()).intersection(set(fmt)))

    ps = vv.get_probe_space(dd, fmt)
    svl = vv.dicts_to_svl(dd, fmt, mask, oddballs)
    for dom in fmt:
        print "-----------"+dom+"-----------"
        tmp = sorted(anssets[dom])
        for val in tmp:
            if type(val) is int:
                print ipp.int2ip(val)
    return svl, fmt


def get_dendrogram(t, duration=30000, mask=32, fmt=None,
        method="average", country_set=None, p=0, oddballs=False,
        fname='dendrogram.pdf', X=None):

    if X is None:
        X, fmt = get_svl(t, duration, mask, fmt, country_set, oddballs)

    dm = np.zeros((len(X) * (len(X) - 1)) // 2, dtype=np.double)
    k = 0
    for i in xrange(0, len(X)-1):
        for j in xrange(i + 1, len(X)):
            dm[k] = 1.0 - vv.closeness(X[i], X[j])
            k = k + 1
            if k % 10000 == 0:
                print k
    Z = linkage(dm, method)
    c, coph_dists = cophenet(Z, dm)
    print " Cophenetic Correlation Coefficient: "+str(c)
    plt.figure(figsize=(15, 10))
    plt.xlabel('sample index')
    plt.ylabel('distance')
    d = dendrogram(
        Z,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
        truncate_mode="lastp",
        p=p,
    )
    plt.savefig(plotsdir+fname, bbox_inches='tight')
    return Z, d, X


def country_closeness(t, duration=30000, mask=32, fmt=None,
        method="average", country_set=None, oddballs=True):

    svl, fmt = get_svl(t, duration, mask, fmt, country_set, oddballs)

    # remove any countries that have less than 10 probes, as if they have near 1
    # probe, the analyses can be more skewed by outlier behavior
    csvl = vv.country_svl(svl)
    countries = csvl.keys()
    for c in countries:
        if len(csvl[c]) < 10:
            del csvl[c]

    for c in csvl:
        if c is None:
            continue
        dist_list = vv.get_dist_list(csvl[c])
        s = ""
        s += "--------"+c+"---------\n"
        s += "average:  "+str(vv.avg_dist(dist_list))+"\n"
        s += "max:      "+str(vv.max_dist(dist_list))+"\n"
        s += "min:      "+str(vv.min_dist(dist_list))+"\n"
        s += "median:      "+str(vv.median_dist(dist_list))+"\n"
        print s

        plt.figure(figsize=(15, 10))
        plt.xlabel('pairwise distance')
        plt.ylabel('# of pairs')
        plt.hist(dist_list, bins=50)
        plt.savefig(plotsdir+c+'_dist.pdf', bbox_inches='tight')


def asn_closeness(t, duration=30000, mask=32, fmt=None,
        method="average", country_set=None, oddballs=True):

    svl, fmt = get_svl(t, duration, mask, fmt, country_set, oddballs)

    # remove any countries that have less than 10 probes, as if they have near 1
    # probe, the analyses can be more skewed by outlier behavior
    asvl = vv.asn_svl(svl)
    asns = asvl.keys()
    for a in asns:
        if len(asvl[a]) < 10:
            del asvl[a]

    for a in asvl:
        if a is None:
            continue
        dist_list = vv.get_dist_list(asvl[a])
        s = ""
        s += "--------"+str(a)+"---------\n"
        s += "average:  "+str(vv.avg_dist(dist_list))+"\n"
        s += "max:      "+str(vv.max_dist(dist_list))+"\n"
        s += "min:      "+str(vv.min_dist(dist_list))+"\n"
        s += "median:      "+str(vv.median_dist(dist_list))+"\n"
        print s

        plt.figure(figsize=(15, 10))
        plt.xlabel('pairwise distance')
        plt.ylabel('# of pairs')
        plt.hist(dist_list, bins=50)
        plt.savefig(plotsdir+str(a)+'_dist.pdf', bbox_inches='tight')
