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

# TODO move this to make_plots.py
def get_dendrogram(t, duration=30000, mask=32, fmt=None,
        method="average", country_set=None, p=0, oddballs=False,
        fname='dendrogram.pdf', X=None):

    if X is None:
        X, fmt, _ = vv.get_svl(t, duration, mask, fmt, country_set, oddballs)

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
