from warriorpy.shorthand import diriofile as df
from warriorpy.shorthand import plotstuff as ps
from matplotlib import pyplot as plt
import numpy as np
import veracity_vector as vv
import metric_validation as mv
import clientcompare as cc
import vgraphs as vg
import networkx as nx
from collections import defaultdict
from statsmodels.distributions.empirical_distribution import ECDF
from warriorpy.net_tools import ipparsing as ipp
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster
import sys
import math
import gc

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
statedir = df.rightdir(basedir+'state/')
rawdirlist = df.getlines(basedir+'state/datapaths.list')
datafiles = df.listfiles(basedir+rawdirlist[0], fullpath=True)
plotsdir = df.rightdir(basedir+"plots/")
ccachef = df.rightdir(statedir+"pickles/")+"ccache.pickle"

##################################################################
#                           CODE
##################################################################


def sausage_linkage(start_time, method="single", fname="", Zf=False, xlim=None, **kwas):
    '''
    :param start_time: int indicating the earliest query the window should include
    :param method: the linkage method to be used
    :param fname: string to be appended to end of plot file name
    :param Zf: boolean -> if True, will load from file (see code for file name).
    :param xlim: x axis limits for plot. Accepts formats: None, [a, b],
    :param **kwas: keyword arguments for vv.get_svl()

    plots the change in distribution of components with respect to maximum
    distance threshold

    '''

    Z, X = get_zx(start_time, method, fname, Zf, **kwas)
    y1 = defaultdict(list)
    y2 = list()
    x = np.arange(.0, 1.0, .01)
    t = len(X)
    for dist in x:
        labels = fcluster(Z, dist, criterion='distance')
        clusters = [[X[z] for z, c in enumerate(labels) if c == y] \
                for y in set(labels)]
        groups = [c for c in clusters if len(c) > 1]
        loners = [c for c in clusters if len(c) == 1]
        i = len(loners)
        g = len(groups)
        s = max([len(z) for z in clusters])
        y1['# groups'].append(g)
        y1['# individuals'].append(i)
        y1['max component size'].append(s)
        try:
            y2.append((t-i-s)/(g-1))
        except ZeroDivisionError:
            y2.append(0)

    fig, ax = plt.subplots(1, 1)
    for label in y1:
        ax.plot(x, y1[label], label=label)
    ax.set_xlabel('maximum distance')
    ax.set_ylabel('# components')

    ax2 = ax.twinx()
    ax2.plot(x, y2, 'k', label='# components')
    ax2.set_ylabel('avg component size')

    ps.set_dim(fig, ax, xdim=13, ydim=7.5, xlim=xlim)
    lgd = ps.legend_setup(ax, 5, "top center", True)
    filename = plotsdir+"linkage_"+method+fname
    plt.savefig(filename+'.pdf', bbox_inches='tight')
    plt.savefig(filename+'.png', bbox_inches='tight')


def plot_fmeasure(start_time, method="average", fname="", Zf=False, **kwas):
    '''
    :param start_time: int indicating the earliest query the window should include
    :param method: the linkage method to be used
    :param fname: string to be appended to end of plot file name
    :param **kwas: keyword arguments for vv.get_svl()
    '''
    Z, X = get_zx(start_time, method, fname, Zf, **kwas)    # thresholds = np.arange(0.


def get_zx(start_time, method="single", fname="", Zf=False, **kwas):
    '''
    :param start_time: int indicating the earliest query the window should include
    :param method: the linkage method to be used
    :param fname: string to be appended to end of plot file name
    :param Zf: boolean -> if True, will load from file (see code for file name).
        NOTE: it will save the most recent version you calculated. Make sure the
        right version of the file exists before setting Zf to true
    :param **kwas: keyword arguments for vv.get_svl()
    :return: linkage, dendrogram's output, svl

    computes and plots dendrogram with respect to distance between clients
    '''
    if Zf is False:
        kwas['start_time'] = start_time
        X, fmt, _, ccache = vv.get_svl(**kwas)
        logger.warning("svl len: "+str(len(X)))

        dm = np.zeros((len(X) * (len(X) - 1)) // 2, dtype=np.double)
        k = 0
        for i in xrange(0, len(X)-1):
            for j in xrange(i + 1, len(X)):
                dm[k] = 1.0 - ccache[X[i]][X[j]]
                k = k + 1
        ccache.dump()
        Z = linkage(dm, method)
        df.pickleout(statedir+'pickles/'+'Z_'+method+fname+'.pickle', (Z, dm, X))
        logger.warning('dumped Z to ' \
                +statedir+'pickles/'+'Z_'+method+fname+'.pickle')
    else:
        Z, dm, X = df.picklein(statedir+'pickles/'+'Z_'+method+fname+'.pickle')
        logger.warning('loaded Z from '+statedir+'pickles/'+'Z_'+method+fname+'.pickle')
    c, coph_dists = cophenet(Z, dm)

    return Z, X


def get_dendrogram(start_time, method="average", p=0, fname="", Zf=False, **kwas):
    '''
    :param start_time: int indicating the earliest query the window should include
    :param method: the linkage method to be used
    :param p: dendrogram will show last p merged clusters
    :param fname: string to be appended to end of plot file name
    :param Zf: boolean -> if True, will load from file (see code for file name).
        NOTE: it will save the most recent version you calculated. Make sure the
        right version of the file exists before setting Zf to true
    :param **kwas: keyword arguments for vv.get_svl()
    :return: linkage, dendrogram's output, svl

    computes and plots dendrogram with respect to distance between clients
    '''

    Z, X = get_zx(start_time, method, fname, Zf, **kwas)
    plt.figure(figsize=(15, 10))
    plt.xlabel('sample index')
    plt.ylabel('distance')
    sys.setrecursionlimit(10000)
    d = dendrogram(
        Z,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
        truncate_mode="lastp",
        p=p,
        no_labels=True,
    )
    filename = plotsdir+"dendrogram"+fname
    plt.savefig(filename+".png", bbox_inches='tight')
    plt.savefig(filename+".pdf", bbox_inches='tight')

    return Z, d, X


def plot_optimizing_window(start_time, duration, fname="", xlim=None,
        maxdur=90000*15, incr=30000, **kwas):
    '''
    :param start_time: int indicating the earliest query the window should include
    :param fname: string to be appended to end of plot file name
    :param xlim: x axis limits for plot. Accepts formats: None, [a, b],
    :param maxdur: the outer bound of the duration range to be covered
    :param incr: the number of seconds to increment the duration by in each loop
    :param **kwas: keyword arguments for vv.get_svl()

    makes line plot varying the duration (x axis) vs the closeness to one's self
    from a different point in time (e.g., for a 10 second duration, self A would
    be time 0-9, and self B would be time 10-19)
    '''

    allvals = list()
    allbars = list()
    allx = list()
    dur = duration
    kwas['return_ccache'] = False
    while dur < maxdur:
        print "getting svls..."
        kwas['duration'] = dur
        kwas['start_time'] = start_time
        svl, __, __ = vv.get_svl(**kwas)
        logger.warning("svl len: "+str(len(svl)))
        svl1 = dict()
        for sv in svl:
            svl1[sv.id] = sv
        kwas['start_time'] = start_time+dur
        svl, __, __ = vv.get_svl(**kwas)
        logger.warning("svl len: "+str(len(svl)))
        svl2 = dict()
        for sv in svl:
            svl2[sv.id] = sv

        print "calculating closeness for subnets...", dur
        vals = list()
        for pid in svl1:
            if pid in svl2:
                vals.append(vv.closeness(svl1[pid], svl2[pid]))

        allvals.append(np.mean(vals))
        allbars.append(np.std(vals))
        allx.append(float(dur)/(60.0*60.0*8.0))
        dur += incr



    fig, ax = plt.subplots(1, 1)
    ax.errorbar(allx, allvals, yerr=allbars)
    ps.set_dim(fig, ax, xdim=13, ydim=7.5, xlim=xlim)
    plt.xlabel("# 8 hour cycles in block duration")
    plt.ylabel("average self closeness")
    lgd = ps.legend_setup(ax, 4, "top center", True)
    filename = plotsdir+"avg_self_closeness"+fname
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    print "saving data..."
    outstr = df.overwrite(statedir+fname+'_avg_self_closeness.csv',
            df.list2col(allvals))


def plot_closeness(start_time, duration, fname="", xlim=[.6, 1.0], loops=15, **kwas):
    '''
    :param start_time: int indicating the earliest query the window should include
    :param fname: string to be appended to end of plot file name
    :param xlim: x axis limits for plot. Accepts formats: None, [a, b],
    :param loops: number of time blocks
    :param **kwas: keyword arguments for vv.get_svl()

    plots:
        1) CDF for pairwise closeness of each pair
        2) CDF for the average pairwise closeness experienced by each probe
        across all other probes
    '''
    means = defaultdict(list)
    vals = list()
    kwas['duration'] = duration
    for l in xrange(0, loops):
        print "getting svl..."
        kwas['start_time'] = start_time+duration*l
        svl, fmt, __, ccache = vv.get_svl(**kwas)
        logger.warning("svl len: "+str(len(svl)))
        print len(svl)

        print "calculating closeness for resolvers..."
        for i in xrange(0, len(svl)-1):
            for j in xrange(i + 1, len(svl)):
                vals.append(ccache[svl[i]][svl[j]])
                means[svl[i].get_id()].append(vals[-1])
                means[svl[j].get_id()].append(vals[-1])
        ccache.dump()
        del ccache, svl, fmt
        gc.collect()

    print "plotting..."
    fig, ax = plt.subplots(1, 1)

    ecdf = ECDF(vals)
    x = list(ecdf.x)
    y = list(ecdf.y)
    ax.plot(x, y, label="pairwise")

    ecdf = ECDF([np.mean(means[z]) for z in means])
    x = list(ecdf.x)
    y = list(ecdf.y)
    ax.plot(x, y, label="average (per client)")

    ps.set_dim(fig, ax, xdim=13, ydim=7.5, xlim=xlim)
    plt.xlabel("pairwise probe closeness")
    plt.ylabel("CDF of pairs")
    lgd = ps.legend_setup(ax, 4, "top center", True)
    filename = plotsdir+"closeness_diff_desc"+fname
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    print "saving data..."
    df.overwrite(statedir+'overall_closeness'+fname+'.csv',
        df.list2col(vals))
    df.overwrite(statedir+'overall_avg_closeness'+fname+'.csv',
        df.list2col([(z, np.mean(means[z])) for z in means]))


def plot_closeness_diff_desc(start_time, fname="", xlim=[.6, 1.0], rmask=16, **kwas):
    '''
    :param start_time: int indicating the earliest query the window should include
    :param fname: string to be appended to end of plot file name
    :param xlim: x axis limits for plot. Accepts formats: None, [a, b],
    :param rmask: mask for resolver IPs
    :param **kwas: keyword arguments for vv.get_svl()

    for each descriptor (ASN, country, registered prefix, /24 subnet), plot the
    CDF of the pairwise closeness of clients, such that the clients in a pair
    come from different groups in the descriptor (e.g., different countries
        for the country descriptor)
    '''
    print "getting svl..."
    kwas['start_time'] = start_time
    svl, fmt, __, ccache = vv.get_svl(**kwas)
    logger.warning("svl len: "+str(len(svl)))
    print len(svl)

    print "getting descriptor lists..."
    csvl = vv.country_svl(svl)
    asvl = vv.asn_svl(svl)
    ssvl = vv.subnet_svl(svl)
    #osvl = vv.owner_svl(svl)
    psvl = vv.prefix_svl(svl)
    lsvl = vv.ldns_svl(svl, rmask, False)
    fmtmask = ipp.make_v4_prefix_mask(rmask)
    to_remove = [
            '208.67.222.123',   # OpenDNS
            '208.67.220.123',
            '8.8.8.8',          # Google Public DNS
            '8.8.4.4',
            '64.6.64.6',        # Verisign
            '64.6.65.6']
    # remove massive public DNS providers
    for ip in to_remove:
        tmp = ipp.ip2int(ip) & fmtmask
        if tmp in lsvl:
            del lsvl[tmp]

    print "calculating closeness for resolvers..."
    lvals = list()
    resolvers = [c for c in lsvl if len(lsvl[c]) > 1]
    for i in xrange(0, len(resolvers)-1):
        for a in lsvl[resolvers[i]]:
            for j in xrange(i+1, len(resolvers)):
                for b in lsvl[resolvers[j]]:
                    lvals.append(ccache[a][b])

    print "calculating closeness for countries..."
    cvals = list()
    countries = [c for c in csvl if len(csvl[c]) > 1]
    for i in xrange(0, len(countries)-1):
        for a in csvl[countries[i]]:
            for j in xrange(i+1, len(countries)):
                for b in csvl[countries[j]]:
                    cvals.append(ccache[a][b])
    print "calculating closeness for ASNs..."
    avals = list()
    asns = [a for a in asvl if len(asvl[a]) > 1]
    for i in xrange(0, len(asns)-1):
        for a in asvl[asns[i]]:
            for j in xrange(i+1, len(asns)):
                for b in asvl[asns[j]]:
                    avals.append(ccache[a][b])
    print "calculating closeness for subnets..."
    svals = list()
    subnets = [s for s in ssvl if len(ssvl[s]) > 1]
    for i in xrange(0, len(subnets)-1):
        for a in ssvl[subnets[i]]:
            for j in xrange(i+1, len(subnets)):
                for b in ssvl[subnets[j]]:
                    svals.append(ccache[a][b])
    '''
    print "calculating closeness for owners..."
    ovals = list()
    owners = [o for o in osvl if len(osvl[o]) > 1]
    for i in xrange(0, len(owners)-1):
        for a in osvl[owners[i]]:
            for j in xrange(i+1, len(owners)):
                for b in osvl[owners[j]]:
                    ovals.append(ccache[a][b])
    '''
    print "calculating closeness for prefixes..."
    pvals = list()
    prefixes = [p for p in psvl if len(psvl[p]) > 1]
    for i in xrange(0, len(prefixes)-1):
        for a in psvl[prefixes[i]]:
            for j in xrange(i+1, len(prefixes)):
                for b in psvl[prefixes[j]]:
                    pvals.append(ccache[a][b])

    print "plotting..."
    #vals = [cvals, avals, svals, ovals, pvals]
    #labels = ['country', 'ASN', 'subnet', 'owner', 'prefix']
    vals = [cvals, avals, svals, pvals, lvals]
    labels = ['country', 'ASN', 'subnet', 'prefix', 'resolver']

    fig, ax = plt.subplots(1, 1)
    for i in xrange(0, len(vals)):
        ecdf = ECDF(vals[i])
        x = list(ecdf.x)
        y = list(ecdf.y)
        ax.plot(x, y, label=labels[i])
    ps.set_dim(fig, ax, xdim=13, ydim=7.5, xlim=xlim)
    plt.xlabel("pairwise probe closeness")
    plt.ylabel("CDF of pairs")
    lgd = ps.legend_setup(ax, 4, "top center", True)
    filename = plotsdir+"closeness_diff_desc"+fname
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    print "saving data..."
    for i in xrange(0, len(vals)):
        outstr = df.overwrite(statedir+labels[i]+'_diff.csv',
                df.list2col(vals[i]))
    ccache.dump()


def plot_closeness_same_desc(start_time, fname="", xlim=[.6, 1.0], rmask=16, **kwas):
    '''
    :param start_time: int indicating the earliest query the window should include
    :param fname: string to be appended to end of plot file name
    :param xlim: x axis limits for plot. Accepts formats: None, [a, b],
    :param rmask: mask for resolver IPs
    :param **kwas: keyword arguments for vv.get_svl()

    for each descriptor (ASN, country, registered prefix, /24 subnet), plot the
    CDF of the pairwise closeness of clients, such that the clients in a pair
    come from the same groups in the descriptor (e.g., same country for the
        country descriptor)
    '''
    print "getting svl..."
    kwas['start_time'] = start_time
    svl, fmt, __, ccache = vv.get_svl(**kwas)
    logger.warning("svl len: "+str(len(svl)))

    print "getting descriptor lists..."
    csvl = vv.country_svl(svl)
    asvl = vv.asn_svl(svl)
    ssvl = vv.subnet_svl(svl)
    #osvl = vv.owner_svl(svl)
    psvl = vv.prefix_svl(svl)
    lsvl = vv.ldns_svl(svl, rmask, False)
    fmtmask = ipp.make_v4_prefix_mask(rmask)
    to_remove = [
            '208.67.222.123',   # OpenDNS
            '208.67.220.123',
            '8.8.8.8',          # Google Public DNS
            '8.8.4.4',
            '64.6.64.6',        # Verisign
            '64.6.65.6']
    # remove massive public DNS providers
    for ip in to_remove:
        tmp = ipp.ip2int(ip) & fmtmask
        if tmp in lsvl:
            del lsvl[tmp]

    print "calculating closeness for resolvers..."
    lvals = list()
    resolvers = lsvl.keys()
    for k in resolvers:
        ksvl = lsvl[k]
        for a in xrange(0, len(ksvl)-1):
            for b in xrange(a+1, len(ksvl)):
                lvals.append(ccache[ksvl[a]][ksvl[b]])

    print "calculating closeness for countries..."
    cvals = list()
    countries = csvl.keys()
    for k in countries:
        ksvl = csvl[k]
        for a in xrange(0, len(ksvl)-1):
            for b in xrange(a+1, len(ksvl)):
                cvals.append(ccache[ksvl[a]][ksvl[b]])
    print "calculating closeness for ASNs..."
    avals = list()
    asns = asvl.keys()
    for k in asns:
        ksvl = asvl[k]
        for a in xrange(0, len(ksvl)-1):
            for b in xrange(a+1, len(ksvl)):
                avals.append(ccache[ksvl[a]][ksvl[b]])
    print "calculating closeness for subnets..."
    svals = list()
    subnets = ssvl.keys()
    for k in subnets:
        ksvl = ssvl[k]
        for a in xrange(0, len(ksvl)-1):
            for b in xrange(a+1, len(ksvl)):
                svals.append(ccache[ksvl[a]][ksvl[b]])
    '''
    print "calculating closeness for owners..."
    ovals = list()
    owners = osvl.keys()
    for k in owners:
        ksvl = osvl[k]
        for a in xrange(0, len(ksvl)-1):
            for b in xrange(a+1, len(ksvl)):
                ovals.append(ccache[ksvl[a]][ksvl[b]])
    '''
    print "calculating closeness for prefixes..."
    pvals = list()
    prefixes = psvl.keys()
    for k in prefixes:
        ksvl = psvl[k]
        for a in xrange(0, len(ksvl)-1):
            for b in xrange(a+1, len(ksvl)):
                pvals.append(ccache[ksvl[a]][ksvl[b]])

    print "plotting..."
    #vals = [cvals, avals, svals, ovals, pvals]
    #labels = ['country', 'ASN', 'subnet', 'owner', 'prefix']
    vals = [cvals, avals, svals, pvals, lvals]
    labels = ['country', 'ASN', 'subnet', 'prefix', 'resolver']

    fig, ax = plt.subplots(1, 1)
    for i in xrange(0, len(vals)):
        ecdf = ECDF(vals[i])
        x = list(ecdf.x)
        y = list(ecdf.y)
        ax.plot(x, y, label=labels[i])
    ps.set_dim(fig, ax, xdim=13, ydim=7.5, xlim=xlim)
    plt.xlabel("pairwise probe closeness")
    plt.ylabel("CDF of pairs")
    lgd = ps.legend_setup(ax, 4, "top center", True)
    filename = plotsdir+"closeness_same_desc"+fname
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    print "saving data..."
    for i in xrange(0, len(vals)):
        outstr = df.overwrite(statedir+labels[i]+'_same.csv',
                df.list2col(vals[i]))
    ccache.dump()


def plot_varying_mc(start_time, fname="", tmin=.75, tmax=1.01, tinc=.025, **kwas):
    '''
    :param start_time: int indicating the earliest query the window should include
    :param fname: string to be appended to end of plot file name
    :param tmin: start of threshold range (e.g., a of np.arange(a, b, c))
    :param tmax: end of threshold range (e.g., b of np.arange(a, b, c))
    :param tinc: the step size of the threshold range (e.g., c of np.arange(a,
        b, c)
    :param **kwas: keyword arguments for vv.get_svl()

    '''

    logger.info("getting svl")
    kwas['start_time'] = start_time
    svl, fmt, __, ccache = vv.get_svl(**kwas)
    mwl = np.arange(tmin, tmax, tinc)
    ccl = vg.get_cc_varying_mc(svl, mwl, ccache)

    fig, ax = plt.subplots(1, 1)
    labels = ['country', 'prefix', 'resolver', 'subnet', 'asn']
    x = mwl
    for label in labels:
        y = [z[label] for z in ccl]
        ax.plot(x,y, label=label)
    ax.set_xlabel('minimum closeness')
    ax.set_ylabel('% of component')

    ax2 = ax.twinx()
    y = [z['quantity'] for z in ccl]
    ax2.plot(x, y, 'k', label='# components')
    ax2.set_ylabel('# components')

    ps.set_dim(fig, ax, xdim=13, ydim=7.5)
    lgd = ps.legend_setup(ax, 5, "top center", True)
    plt.savefig(plotsdir+'components_'+fname+'.pdf', bbox_inches='tight')
    plt.savefig(plotsdir+'components_'+fname+'.png', bbox_inches='tight')
    ccache.dump()


def inv_hist(start_time, fname="", thresh=.35, **kwas):
    logger.info("getting svl...")
    kwas['start_time'] = start_time
    kwas['return_ccache'] = False
    svl, fmt, __ = vv.get_svl(**kwas)
    logger.info("getting ipsl...")
    ipsl, dompairs = get_ip_sets(svl)
    logger.info("getting pairing counts...")
    pc = vg.get_pairing_counts(ipsl)
    ipcount = len(pc)
    logger.info("building inv. graph...")
    G = vg.build_inv_graph(pc)
    vg.remove_far_edges(G, thresh)
    dd = vg.nodes_by_degree(G)
    vg.remove_degree_below(G, dd, 1)
    weights = [e[2] for e in G.edges_iter(data='weight')]
    cc = list(nx.connected_components(G))
    #print cc
    for c in cc:
        print "****************************"
        print str(len(c))
        print set([dom for ip in c for dom in dompairs[ip]])
        weights = [w for a,b,w in G.edges_iter(c, data='weight')]
        print "median weight: "+str(np.median(weights))
        print "average weight: "+str(np.mean(weights))
    print "num connected comps: "+str(len(cc))
    print "size of connected comps: "+str(np.median([len(z) for z in cc]))

    plt.figure(figsize=(15, 10))
    plt.xlabel('pairwise closeness')
    plt.ylabel('# of pairs (servers)')
    plt.hist(weights, bins=100)
    plt.savefig(plotsdir+fname+'inv_hist.pdf', bbox_inches='tight')


def plot_self_match(start_time, duration, fname="", loops=7, gap=0, thresh=10,
        **kwas):
    '''
    lines:  1) domain independent cdf of ALL matches
            2-n) cdf of matches for domain with answer space > thresh
            n-m) cdf of matches for ALL domains with answer space < thresh
    '''
    valsd = defaultdict(list) # {pid: [vals]}
    bigsvld = dict()
    kwas['duration'] = duration
    kwas['return_ccache'] = False
    for i in xrange(0, loops):
        print (t+2*duration*i, duration)
        kwas['start_time'] = start_time+2*(gap+duration)*i
        svld, allsvl, allfmt, anssets = mv.arrange_self_data(start_time, gap, 2, **kwas)

        pids, vals = mv.self_match(svld)
        for pid, val in zip(pids, vals):
            valsd[pid].append(val)
            bigsvld[pid] = svld[pid]

    results = list()
    for pid in valsd:
        results.append((pid, np.mean(valsd[pid])))

    results = sorted(results, key=lambda z: z[1], reverse=True)

    fig, ax = plt.subplots(1, 1)
    ecdf = ECDF([z[1] for z in results])
    x = list(ecdf.x)
    y = list(ecdf.y)
    ax.plot(x, y, color='k')
    ax.axvline(np.median(x), color='r', linestyle='--')
    ps.set_dim(fig, ax, xdim=13, ydim=7.5)
    plt.xlabel("closeness to self")
    plt.ylabel("CDF of clients")
    filename = plotsdir+"self_closeness"+fname
    fig.savefig(filename+'.png', bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_inches='tight')
    plt.close(fig)

    labels = defaultdict(list)
    for pid, val in results:
        labels['countries'].append((val, bigsvld[pid][0].get_country()))
        labels['subnets'].append((val, bigsvld[pid][0].get_subnet()))
        labels['prefixes'].append((val, bigsvld[pid][0].get_prefix()))
        labels['asns'].append((val, bigsvld[pid][0].get_asn()))
        labels['resolvers'].append((val, bigsvld[pid][0].get_ldns()))

    for k in labels:
        data = sorted([(y, np.mean([z[0] for z in labels[k] if z[1] == y]),
                len([z[0] for z in labels[k] if z[1] == y])) \
                for y in set([v[1] for v in labels[k]])], key=lambda x: x[1])
        df.overwrite(statedir+'self_closeness_'+k+fname+'.csv',
                df.list2col(data))

    print "saving data..."
    outstr = df.overwrite(statedir+'self_closeness'+fname+'.csv',
            df.list2col(results))


def plot_examine_self_diff(start_time, fname="", loops=2, gap=0, thresh=10,
        **kwas):
    '''
    lines:  1) domain independent cdf of ALL matches
            2-n) cdf of matches for domain with answer space > thresh
            n-m) cdf of matches for ALL domains with answer space < thresh
    '''
    kwas['start_time'] = start_time
    kwas['return_ccache'] = False
    svld, allsvl, allfmt, anssets = mv.arrange_self_data(start_time, gap, loops, **kwas)

    sm = mv.examine_self_diff(svld)
    vals = [[], []]
    labels = ['all', 'small']
    for dom in sm:
        vals[0] = vals[0] + sm[dom]
        if len(anssets[dom]) < thresh:
            vals[1] += sm[dom]
        else:
            vals.append(sm[dom])
            labels.append(dom)

    fig, ax = plt.subplots(1, 1)
    for i in xrange(0, len(vals)):
        ecdf = ECDF(vals[i])
        x = [z-32+mask for z in ecdf.x]
        y = list(ecdf.y)
        ax.plot(x, y, label=labels[i])
    ps.set_dim(fig, ax, xdim=13, ydim=7.5)
    plt.xlabel("self mask match by domain")
    plt.ylabel("CDF of clients")
    lgd = ps.legend_setup(ax, 4, "top center", True)
    filename = plotsdir+"self_mask"+fname
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    print "saving data..."
    for i in xrange(0, len(vals)):
        outstr = df.overwrite(statedir+labels[i]+'_self_jaccard.csv',
                df.list2col(vals[i]))


def plot_examine_diff_diff(start_time, fname="", loops=2, gap=0,
        thresh=10, **kwas):
    '''
    lines:  1) domain independent cdf of ALL matches
            2-n) cdf of matches for domain with answer space > thresh
            n-m) cdf of matches for ALL domains with answer space < thresh
    '''
    kwas['start_time'] = start_time
    kwas['return_ccache'] = False
    svld, allsvl, allfmt, anssets = mv.arrange_self_data(start_time, gap, loops, **kwas)

    sm = mv.examine_diff_diff(svld)
    vals = [[], []]
    labels = ['all', 'small']
    for dom in sm:
        vals[0] = vals[0] + sm[dom]
        if len(anssets[dom]) < thresh:
            vals[1] += sm[dom]
        else:
            vals.append(sm[dom])
            labels.append(dom)

    fig, ax = plt.subplots(1, 1)
    for i in xrange(0, len(vals)):
        print "*****************"+labels[i]+"*********************"
        print vals[i]
        ecdf = ECDF(vals[i])
        x = list(ecdf.x)
        y = list(ecdf.y)
        ax.plot(x, y, label=labels[i])
    ps.set_dim(fig, ax, xdim=13, ydim=7.5)
    plt.xlabel("diff mask match by domain")
    plt.ylabel("CDF of clients")
    lgd = ps.legend_setup(ax, 4, "top center", True)
    filename = plotsdir+"diff_mask"+fname
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    print "saving data..."
    for i in xrange(0, len(vals)):
        outstr = df.overwrite(statedir+labels[i]+'_diff_jaccard.csv',
                df.list2col(vals[i]))


def plot_measure_expansion(start_time, fname="", loops=10, gap=0, thresh=10,
        **kwas):

    kwas['start_time'] = start_time
    kwas['return_ccache'] = False
    svld, allsvl, allfmt, anssets = mv.arrange_self_data(start_time, gap, loops, **kwas)
    keys = svld.keys()

    counts = mv.measure_expansion(svld)
    domvals = defaultdict(lambda: defaultdict(list))
    allvals = defaultdict(list)
    smallvals = defaultdict(list)

    for c in counts:
        for dom in c:
            for i, val in enumerate(c[dom]['ratio']):
                allvals[i].append(val)
                if len(anssets[dom]) < thresh:
                    smallvals[i].append(val)
                else:
                    domvals[dom][i].append(val)

    labels = ['all', 'small'] + domvals.keys()
    vals = list()
    for i in labels:
        vals.append([])
    for i in sorted(allvals.keys()):
        vals[0].append(np.mean(allvals[i]))
        vals[1].append(np.mean(smallvals[i]))
        for j, dom in enumerate(labels[2:]):
            vals[j+2].append(np.mean(domvals[dom][i]))

    fig, ax = plt.subplots(1, 1)
    for i in xrange(0, len(vals)):
        ax.plot(vals[i], label=labels[i])
    ps.set_dim(fig, ax, xdim=13, ydim=7.5)
    plt.xlabel("cycle #")
    plt.ylabel("# new IPs / ans. size")
    lgd = ps.legend_setup(ax, 4, "top center", True)
    filename = plotsdir+"newvssize"+fname
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    print "saving data..."
    for i in xrange(0, len(vals)):
        outstr = df.overwrite(statedir+labels[i]+'newvssize.csv',
                df.list2col(vals[i]))


def plot_resolver_comparison(start_time, fname="", xlim=[.6, 1.0], rmask=16, **kwas):
    '''
    :param start_time: int indicating the earliest query the window should include
    :param fname: string to be appended to end of plot file name
    :param xlim: x axis limits for plot. Accepts formats: None, [a, b],
    :param rmask: mask for resolver IPs
    :param **kwas: keyword arguments for vv.get_svl()
    :returns: [country, ASN, subnet, prefix] pair dictionaries of closeness lists

    gets pairwise closeness of probes with different descriptors to find odd
    behavior (probes in difference descriptors with high closeness scores)

    NOTE: writes data to files for conveniece
    '''


    print("getting svl...")
    kwas['start_time'] = start_time
    svl, fmt, __, ccache = vv.get_svl(**kwas)
    logger.warning("svl len: "+str(len(svl)))

    nearbies = cc.nearby_probes_diff_ldns(svl, rmask)

    vals = defaultdict(list)
    fmtmask = ipp.make_v4_prefix_mask(rmask)
    for group in nearbies:
        for i in xrange(0, len(group)-1):
            for j in xrange(i+1, len(group)):
                a = group[i]
                b = group[j]
                closeness = ccache[a][b]
                if a.get_ldns() & fmtmask == b.get_ldns() & fmtmask:
                    vals['same LDNS'].append(closeness)
                else:
                    vals['diff LDNS'].append(closeness)
    ccache.dump()

    fig, ax = plt.subplots(1, 1)
    for l in vals:
        ecdf = ECDF(vals[l])
        x = list(ecdf.x)
        y = list(ecdf.y)
        ax.plot(x, y, label=l)
    ps.set_dim(fig, ax, xdim=13, ydim=7.5, xlim=xlim)
    plt.xlabel("pairwise probe closeness")
    plt.ylabel("CDF of pairs")
    lgd = ps.legend_setup(ax, 4, "top center", True)
    filename = plotsdir+"closeness_ldns"+fname
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)



