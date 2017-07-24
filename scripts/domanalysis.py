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


def client_answer_space(t, duration=30000, mask=32, fmt=None, country=None):
    print "getting window"
    if fmt is None:
        fmt = domains
    elif type(fmt) is int:
        fmt = domains[-fmt:]

    window = vv.get_window(t, duration, fmt, country)
    print "converting window to dict"
    dd = vv.window_to_dict(window)
    infod = defaultdict(lambda: defaultdict(set))
    for probeip in dd:
        for dom in dd[probeip]:
            for ip in dd[probeip][dom]:
                if ip != 0 and IP(ipp.int2ip(ip)+"/32").iptype() == 'PUBLIC':
                    infod[probeip][ip].add(dom)

    print "plotting IPs/domains per client"
    # CDF of number of IPs seen by a client
    ippc = list() # IPs per client
    for probeip in infod:
        ips = float(len(infod[probeip]))
        doms = set.union(*[infod[probeip][z] for z in infod[probeip]])
        ippc.append(ips/float(len(doms)))
    ecdf = ECDF(ippc)
    x = list(ecdf.x)
    y = list(ecdf.y)
    fig, ax = plt.subplots(1, 1)
    ax.plot(x, y)
    ps.set_dim(fig, ax)
    plt.xlabel("IPs seen by probe")
    plt.ylabel("CDF of probes")
    plt.tight_layout()
    ps.saveme(ax, fig, plotsdir+"IPs_per_client")
    plt.close(fig)

    print "plotting shared IP"
    # bar of % of clients that saw shared (relative to other dom answers on same
    # client) IP for a dom
    total_clients = len(infod)
    domd = defaultdict(lambda: 0)
    for probeip in infod:
        for ans in infod[probeip]:
            if len(infod[probeip][ans]) > 1:
                for dom in infod[probeip][ans]:
                    domd[dom] += 1
    bars = list()
    for dom in fmt:
        bars.append(domd[dom])
    fig, ax = plt.subplots(1, 1)
    barlocs = range(len(bars)+1)[1:]
    ax.bar(barlocs, bars, tick_label=[str(z) for z in range(len(bars)+1)[1:]],
            align='center')
    # NOTE: tight layout will hang if a dimension is using log but doesn't cover
    # wide enough of a span of values (for example, if the only value is 1)
    ps.set_dim(fig, ax, ylog=True)
    plt.xlabel("domain")
    plt.ylabel("% multi-use clients")
    ps.saveme(ax, fig, plotsdir+"multi_use")
    plt.close(fig)

    print "plotting max domains per IP per client"
    # CDF of max number of domains per IP for a client
    dpip = list() # domains per ip
    for probeip in infod:
        dpip.append(max([len(infod[probeip][z]) for z in infod[probeip]]))
    ecdf = ECDF(dpip)
    x = list(ecdf.x)
    y = list(ecdf.y)
    fig, ax = plt.subplots(1, 1)
    ax.plot(x, y)
    ps.set_dim(fig, ax)
    plt.xlabel("max domains per IP for client")
    plt.ylabel("CDF of clients")
    plt.tight_layout()
    ps.saveme(ax, fig, plotsdir+"domains_per_IP")
    plt.close(fig)

    print "plotting clients per IP"
    # CDF of number of clients per IP for a domain
    ipd = defaultdict(lambda: defaultdict(set))
    for probeip in infod:
        for ans in infod[probeip]:
            for dom in infod[probeip][ans]:
                ipd[dom][ans].add(probeip)
    dd = defaultdict(list)
    for dom in ipd:
        for ans in ipd[dom]:
            dd[dom].append(len(ipd[dom][ans]))
    small = list()
    med = list()
    large = list()
    massive = list()
    print fmt
    for ind, dom in enumerate(fmt):
        if dom in dd:
            if ind > 39:
                massive.append((dd[dom], dom))
            elif ind > 25:
                large.append((dd[dom], dom))
            elif ind > 10:
                med.append((dd[dom], dom))
            else:
                small.append((dd[dom], dom))
        else:
            print "missing: "+dom+"!!!!!!!!!!!!!!!!!"

    fig, ax = plt.subplots(1, 1)
    for item in small:
        ecdf = ECDF(item[0])
        x = list(ecdf.x)
        y = list(ecdf.y)
        ax.plot(x, y, label=item[1])
    ps.set_dim(fig, ax, xlog=True)
    plt.xlabel("# of clients")
    plt.ylabel("CDF of IPs")
    lgd = ps.legend_setup(ax, 4, "top center", True)
    filename = plotsdir+"small_clients_per_IP"
    pickle.dump(ax, open(filename+'.axp', 'w'))
    pickle.dump(fig, open(filename+'.figp', 'w'))
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    fig, ax = plt.subplots(1, 1)
    for item in med:
        ecdf = ECDF(item[0])
        x = list(ecdf.x)
        y = list(ecdf.y)
        ax.plot(x, y, label=item[1])
    ps.set_dim(fig, ax, xlog=True)
    plt.xlabel("# of clients")
    plt.ylabel("CDF of IPs")
    lgd = ps.legend_setup(ax, 4, "top center", True)
    filename = plotsdir+"med_clients_per_IP"
    pickle.dump(ax, open(filename+'.axp', 'w'))
    pickle.dump(fig, open(filename+'.figp', 'w'))
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    fig, ax = plt.subplots(1, 1)
    for item in large:
        ecdf = ECDF(item[0])
        x = list(ecdf.x)
        y = list(ecdf.y)
        ax.plot(x, y, label=item[1])
    ps.set_dim(fig, ax, xlog=True)
    plt.xlabel("# of clients")
    plt.ylabel("CDF of IPs")
    lgd = ps.legend_setup(ax, 4, "top center", True)
    filename = plotsdir+"large_clients_per_IP"
    pickle.dump(ax, open(filename+'.axp', 'w'))
    pickle.dump(fig, open(filename+'.figp', 'w'))
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    fig, ax = plt.subplots(1, 1)
    for item in massive:
        ecdf = ECDF(item[0])
        x = list(ecdf.x)
        y = list(ecdf.y)
        ax.plot(x, y, label=item[1])
    ps.set_dim(fig, ax, xlog=True)
    plt.xlabel("# of clients")
    plt.ylabel("CDF of IPs")
    lgd = ps.legend_setup(ax, 4, "top center", True)
    filename = plotsdir+"massive_clients_per_IP"
    pickle.dump(ax, open(filename+'.axp', 'w'))
    pickle.dump(fig, open(filename+'.figp', 'w'))
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    # graph of domain pairings, where each edge represents the number of
    # pairings between two domains and each vertex is a domain

    # Apriori analysis of domain grouping


def plot_answer_space_histogram(t, duration=30000, mask=32, fmt=None,
        country=None):
    anssets = vv.get_answer_space_dict(t, duration, mask, fmt, country)
    spacelist = list()
    for dom in anssets:
        spacelist.append(len(anssets[dom]))
    plt.hist(spacelist, bins=len(set(spacelist)), log=True)
    plt.xlabel('number of IPs')
    plt.ylabel('number of domains')
    plt.savefig(plotsdir+'histogram.png', bbox_inches='tight')
    n = 0
    outlist = [z+": "+str(len(anssets[z]))+"\n" for z in anssets if len(anssets[z]) > n]

    print "number of sites with "+str(n)+" servers: "+str(len(outlist))

    n = 9
    outlist = [z+": "+str(len(anssets[z]))+"\n" for z in anssets if len(anssets[z]) > n]

    print "number of sites with "+str(n)+" servers: "+str(len(outlist))


    n = 50
    outlist = [z+": "+str(len(anssets[z]))+"\n" for z in anssets if len(anssets[z]) > n]

    print "number of sites with "+str(n)+" servers: "+str(len(outlist))

    n = 0
    outlist = [(z,len(anssets[z])) for z in anssets if len(anssets[z]) >
            n]

    sl = sorted(outlist, key=lambda z: z[1])
    for val in sl:
        print val

    overlap = list()
    overlap2 = list()
    if fmt is None:
        fmt = domains
    elif type(fmt) is int:
        fmt = domains[-fmt:]
    for x, a in enumerate(fmt):
        for y, b in enumerate(fmt):
            if y != x:
                inters = anssets[a].intersection(anssets[b])
                if 0 in inters:
                    inters.remove(0)
                if '0.0.0.0' in inters:
                    inters.remove('0.0.0.0')
                if len(inters) > 0:
                    overlap.append((a, b, len(inters)))
                    overlap2.append((a, b, inters))
    overlap = sorted(overlap, key=lambda z: z[2])
    overlap2 = sorted(overlap2, key=lambda z: len(z[2]))
    for val in overlap:
        print val

    outdict = defaultdict(list)
    for val in overlap2:
        outdict[val[0]].append(val[1])
    return outdict, anssets

