from warriorpy.shorthand import diriofile as df
import numpy as np
import veracity_vector as vv
from warriorpy.net_tools import ipparsing as ipp
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

##################################################################
#                           CODE
##################################################################


def get_dendrogram(t, duration=30000, mask=32, fmt=None,
        method="average", country_set=None, p=0, oddballs=False,
        fname='dendrogram.pdf', X=None):

    if X is None:
        X, fmt = vv.get_svl(t, duration, mask, fmt, country_set, oddballs)

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


def country_full_closeness_hist(t, duration=30000, mask=32, fmt=None,
        method="average", country_set=None, oddballs=True):

    svl, fmt = vv.get_svl(t, duration, mask, fmt, country_set, oddballs)

    # remove any countries that have less than 10 probes, as if they have near 1
    # probe, the analyses can be more skewed by outlier behavior
    csvl = vv.country_svl(svl)
    countries = csvl.keys()
    for c in countries:
        if len(csvl[c]) < 10:
            del csvl[c]

    for c in csvl:
        if c is None:
            print "found a 'None' in csvl..."
            continue
        dist_list = vv.get_dist_list(csvl[c])
        s = ""
        s += "--------"+c+"---------\n"
        s += "average:  "+str(np.mean(dist_list))+"\n"
        s += "max:      "+str(max(dist_list))+"\n"
        s += "min:      "+str(min(dist_list))+"\n"
        s += "median:      "+str(np.median(dist_list))+"\n"
        print s

        plt.figure(figsize=(15, 10))
        plt.xlabel('pairwise distance')
        plt.ylabel('# of pairs')
        plt.hist(dist_list, bins=50)
        plt.savefig(plotsdir+c+'_dist.pdf', bbox_inches='tight')


def asn_full_closeness_hist(t, duration=30000, mask=32, fmt=None,
        method="average", country_set=None, oddballs=True):

    svl, fmt = vv.get_svl(t, duration, mask, fmt, country_set, oddballs)

    # remove any countries that have less than 10 probes, as if they have near 1
    # probe, the analyses can be more skewed by outlier behavior
    asvl = vv.asn_svl(svl)
    asns = asvl.keys()
    for a in asns:
        if len(asvl[a]) < 10:
            del asvl[a]

    for a in asvl:
        if a is None:
            print "found a 'None' in asvl..."
            continue
        cl = vv.get_cl(asvl[a])
        s = ""
        s += "--------"+str(a)+"---------\n"
        s += "average:  "+str(vv.avg_dist(cl))+"\n"
        s += "max:      "+str(vv.max_dist(cl))+"\n"
        s += "min:      "+str(vv.min_dist(cl))+"\n"
        s += "median:      "+str(vv.median_dist(cl))+"\n"
        print s

        plt.figure(figsize=(15, 10))
        plt.xlabel('pairwise distance')
        plt.ylabel('# of pairs')
        plt.hist(dist_list, bins=50)
        plt.savefig(plotsdir+str(a)+'_dist.pdf', bbox_inches='tight')
