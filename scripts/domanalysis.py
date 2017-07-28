from warriorpy.shorthand import diriofile as df
from pymongo import MongoClient
import cPickle as pickle
import json
import numpy as np
import veracity_vector as vv
from warriorpy.net_tools import ipparsing as ipp
from sklearn.cluster import DBSCAN
from sklearn import metrics
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist
import clientgrouping as cg
from collections import defaultdict
from statsmodels.distributions.empirical_distribution import ECDF
from warriorpy.shorthand import plotstuff as ps
from IPy import IP

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
coll = db.m30002_may17_full

ag = db.answergroups

##################################################################
#                           CODE
##################################################################


