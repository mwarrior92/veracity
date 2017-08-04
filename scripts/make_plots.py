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


def get_dendrogram(t, duration=30000, mask=32, fmt=None,
        method="average", country_set=None, p=0, oddballs=True,
        fname='.pdf', maxmissing=0):
    '''
    :param t: int indicating the earliest query the window should include
    :param duration: int indication the span of time covered by the window,
        in seconds
    :param mask: int, prefix mask to use over domain IPs
    :param fmt: see transform fmt
    :param country_set: the set of countries the window should include queries from.
        If None, then all countries will be inluded
    :param method: the linkage method to be used
    :param p: dendrogram will show last p merged clusters
    :param oddballs: if True, will include non-public IPs (10.x.x.x, etc); if
        False, will only include public IPs
    :param fname: string to be appended to end of plot file name
    :param maxmissing: the maximum number of domains that can be missing
        measurements for a client to still be included in this measurement
    :return: linkage, dendrogram's output, svl

    computes and plots dendrogram
    '''

    X, fmt, _ = vv.get_svl(t, duration, mask, fmt, country_set, oddballs, maxmissing)
    ccache = vv.init_ccache(None, ccachef, t, duration, mask, fmt, oddballs, maxmissing)
    logger.warning("svl len: "+str(len(X)))

    dm = np.zeros((len(X) * (len(X) - 1)) // 2, dtype=np.double)
    k = 0
    for i in xrange(0, len(X)-1):
        for j in xrange(i + 1, len(X)):
            dm[k] = 1.0 - ccache[X[i]][X[j]]
            k = k + 1
            if k % 10000 == 0:
                ccache.dump()
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
    filename = plotsdir+"dendrogram"+fname
    plt.savefig(filename+".png", bbox_inches='tight')
    plt.savefig(filename+".pdf", bbox_inches='tight')

    return Z, d, X


def plot_optimizing_window(t, duration=30000, mask=32, fmt=None,
        country_set=None, oddballs=True, fname="", xlim=None,
        maxdur=90000*15, incr=30000, maxmissing=0):
    '''
    :param t: int indicating the earliest query the window should include
    :param duration: int indication the span of time covered by the window,
        in seconds
    :param mask: int, prefix mask to use over domain IPs
    :param fmt: see transform fmt
    :param country_set: the set of countries the window should include queries from.
        If None, then all countries will be inluded
    :param oddballs: if True, will include non-public IPs (10.x.x.x, etc); if
        False, will only include public IPs
    :param fname: string to be appended to end of plot file name
    :param xlim: x axis limits for plot. Accepts formats: None, [a, b],
    :param maxmissing: the maximum number of domains that can be missing
        measurements for a client to still be included in this measurement
        [a, None], [None, b], where None values result in the default/automatic
        limits
    :param maxdur: the outer bound of the duration range to be covered
    :param incr: the number of seconds to increment the duration by in each loop
    :param maxmissing: the maximum number of domains that can be missing
        measurements for a client to still be included in this measurement

    makes line plot varying the duration (x axis) vs the closeness to one's self
    from a different point in time (e.g., for a 10 second duration, self A would
    be time 0-9, and self B would be time 10-19)
    '''

    allvals = list()
    allbars = list()
    allx = list()
    dur = duration
    while dur < maxdur:
        print "getting svls..."
        svl, _, _ = vv.get_svl(t, dur, mask, fmt, country_set, oddballs, maxmissing)
        logger.warning("svl len: "+str(len(svl)))
        svl1 = dict()
        for sv in svl:
            svl1[sv.id] = sv
        svl, _, _ = vv.get_svl(t+dur, dur, mask, fmt, country_set, oddballs, maxmissing)
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


def plot_closeness_diff_desc(t, duration=30000, mask=32, fmt=None,
        country_set=None, oddballs=True, fname="", xlim=[.6, 1.0], maxmissing=0):
    '''
    :param t: int indicating the earliest query the window should include
    :param duration: int indication the span of time covered by the window,
        in seconds
    :param mask: int, prefix mask to use over domain IPs
    :param fmt: see transform fmt
    :param country_set: the set of countries the window should include queries from.
        If None, then all countries will be inluded
    :param method: the linkage method to be used
    :param p: dendrogram will show last p merged clusters
    :param oddballs: if True, will include non-public IPs (10.x.x.x, etc); if
        False, will only include public IPs
    :param fname: string to be appended to end of plot file name
    :param maxmissing: the maximum number of domains that can be missing
        measurements for a client to still be included in this measurement

    for each descriptor (ASN, country, registered prefix, /24 subnet), plot the
    CDF of the pairwise closeness of clients, such that the clients in a pair
    come from different groups in the descriptor (e.g., different countries
        for the country descriptor)
    '''
    print "getting svl..."
    svl, fmt, _ = vv.get_svl(t, duration, mask, fmt, country_set, oddballs, maxmissing)
    logger.warning("svl len: "+str(len(svl)))
    print len(svl)

    print "getting descriptor lists..."
    csvl = vv.country_svl(svl)
    asvl = vv.asn_svl(svl)
    ssvl = vv.subnet_svl(svl)
    #osvl = vv.owner_svl(svl)
    psvl = vv.prefix_svl(svl)

    ccache = vv.init_ccache(None, ccachef, t, duration, mask, fmt, oddballs, maxmissing)

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
    vals = [cvals, avals, svals, pvals]
    labels = ['country', 'ASN', 'subnet', 'prefix']

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


def plot_closeness_same_desc(t, duration=30000, mask=32, fmt=None,
        country_set=None, oddballs=True, fname="", xlim=[.6, 1.0], maxmissing=0):
    '''
    :param t: int indicating the earliest query the window should include
    :param duration: int indication the span of time covered by the window,
        in seconds
    :param mask: int, prefix mask to use over domain IPs
    :param fmt: see transform fmt
    :param country_set: the set of countries the window should include queries from.
        If None, then all countries will be inluded
    :param method: the linkage method to be used
    :param p: dendrogram will show last p merged clusters
    :param oddballs: if True, will include non-public IPs (10.x.x.x, etc); if
        False, will only include public IPs
    :param fname: string to be appended to end of plot file name
    :param maxmissing: the maximum number of domains that can be missing
        measurements for a client to still be included in this measurement

    for each descriptor (ASN, country, registered prefix, /24 subnet), plot the
    CDF of the pairwise closeness of clients, such that the clients in a pair
    come from the same groups in the descriptor (e.g., same country for the
        country descriptor)
    '''
    print "getting svl..."
    svl, fmt, _ = vv.get_svl(t, duration, mask, fmt, country_set, oddballs, maxmissing)
    logger.warning("svl len: "+str(len(svl)))

    print "getting descriptor lists..."
    csvl = vv.country_svl(svl)
    asvl = vv.asn_svl(svl)
    ssvl = vv.subnet_svl(svl)
    #osvl = vv.owner_svl(svl)
    psvl = vv.prefix_svl(svl)

    ccache = vv.init_ccache(None, ccachef, t, duration, mask, fmt, oddballs, maxmissing)

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
    vals = [cvals, avals, svals, pvals]
    labels = ['country', 'ASN', 'subnet', 'prefix']

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


def plot_csize_vs_mc(t, duration=30000, mask=32, fmt=None, country_set=None,
        oddballs=True, fname="", maxmissing=0, tmin=.5, tmax=1.01,
        tinc=.01):
    '''
    :param t: int indicating the earliest query the window should include
    :param duration: int indication the span of time covered by the window,
        in seconds
    :param mask: int, prefix mask to use over domain IPs
    :param fmt: see transform fmt
    :param country_set: the set of countries the window should include queries from.
        If None, then all countries will be inluded
    :param oddballs: if True, will include non-public IPs (10.x.x.x, etc); if
        False, will only include public IPs
    :param fname: string to be appended to end of plot file name
    :param maxmissing: the maximum number of domains that can be missing
        measurements for a client to still be included in this measurement
    :param tmin: start of threshold range (e.g., a of np.arange(a, b, c))
    :param tmax: end of threshold range (e.g., b of np.arange(a, b, c))
    :param tinc: the step size of the threshold range (e.g., c of np.arange(a,
        b, c)

    plots the average cluster
    '''

    print "getting svl"
    svl, fmt, _ = vv.get_svl(t, duration, mask, fmt, country_set, oddballs, maxmissing)
    ccache = vv.init_ccache(None, ccachef, t, duration, mask, fmt, oddballs, maxmissing)
    print "calling csize vs mc"
    x, y = vg.csize_vs_mc(svl, np.arange(tmin, tmax, tinc), ccache)
    plt.figure(figsize=(15, 10))
    plt.xlabel('minimum closeness')
    plt.ylabel('size of biggest component')
    plt.plot(x,y)
    plt.savefig(plotsdir+'components_'+fname+'.pdf', bbox_inches='tight')
    plt.savefig(plotsdir+'components_'+fname+'.png', bbox_inches='tight')
    ccache.dump()


def plot_ccount_vs_mc(t, duration=30000, mask=32, fmt=None, country_set=None,
        oddballs=True, fname="", maxmissing=0, ccache=None, tmin=.5, tmax=1.01,
        tinc=.01):

    svl, fmt, _ = vv.get_svl(t, duration, mask, fmt, country_set, oddballs, maxmissing)
    ccache = vv.init_ccache(None, ccachef, t, duration, mask, fmt, oddballs, maxmissing)
    x, y = vg.ccount_vs_mc(svl, np.arange(tmin, tmax, tinc), ccache)
    plt.figure(figsize=(15, 10))
    plt.xlabel('minimum closeness')
    plt.ylabel('# of components')
    plt.plot(x,y)
    plt.savefig(plotsdir+'components_'+fname+'.pdf', bbox_inches='tight')
    ccache.dump()


def inv_hist(t, duration=30000, mask=32, fmt=None, country_set=None,
        oddballs=True, fname="", thresh=.35, maxmissing=0):

    logger.info("getting svl...")
    svl, fmt, _ = vv.get_svl(t, duration, mask, fmt, country_set, oddballs, maxmissing)
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


def plot_self_match(t, duration=6*30000, mask=32, fmt=None,
        country_set=None, oddballs=True, fname="", loops=2,
        gap=1, thresh=10, maxmissing=0):
    '''
    lines:  1) domain independent cdf of ALL matches
            2-n) cdf of matches for domain with answer space > thresh
            n-m) cdf of matches for ALL domains with answer space < thresh
    '''
    svld, allsvl, allfmt, anssets = mv.arrange_self_data(t, duration, gap, loops, mask,
            fmt, country_set, oddballs, maxmissing)
    keys = svld.keys()

    sm = mv.self_match(svld)
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
        x = list(ecdf.x)
        y = list(ecdf.y)
        ax.plot(x, y, label=labels[i])
    ps.set_dim(fig, ax, xdim=13, ydim=7.5)
    plt.xlabel("self jaccard diff by domain")
    plt.ylabel("CDF of clients")
    lgd = ps.legend_setup(ax, 4, "top center", True)
    filename = plotsdir+"self_jaccard"+fname
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    print "saving data..."
    for i in xrange(0, len(vals)):
        outstr = df.overwrite(statedir+labels[i]+'_self_jaccard.csv',
                df.list2col(vals[i]))


def plot_examine_self_diff(t, duration=30000, mask=32, fmt=None,
        country_set=None, oddballs=True, fname="", loops=2,
        gap=0, thresh=10, maxmissing=0):
    '''
    lines:  1) domain independent cdf of ALL matches
            2-n) cdf of matches for domain with answer space > thresh
            n-m) cdf of matches for ALL domains with answer space < thresh
    '''
    svld, allsvl, allfmt, anssets = mv.arrange_self_data(t, duration, gap, loops, mask,
            fmt, country_set, oddballs, maxmissing)
    keys = svld.keys()

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


def plot_examine_diff_diff(t, duration=30000, mask=32, fmt=None,
        country_set=None, oddballs=True, fname="", loops=2,
        gap=0, thresh=10, maxmissing=0):
    '''
    lines:  1) domain independent cdf of ALL matches
            2-n) cdf of matches for domain with answer space > thresh
            n-m) cdf of matches for ALL domains with answer space < thresh
    '''
    svld, allsvl, allfmt, anssets = mv.arrange_self_data(t, duration, gap, loops, mask,
            fmt, country_set, oddballs, maxmissing)
    keys = svld.keys()

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


def plot_measure_expansion(t, duration=30000, mask=32, fmt=None,
        country_set=None, oddballs=True, fname="", loops=10,
        gap=0, thresh=10, maxmissing=0):
    svld, allsvl, allfmt, anssets = mv.arrange_self_data(t, duration, gap, loops, mask,
            fmt, country_set, oddballs, maxmissing)
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


def plot_resolver_comparison(t, duration=30000, mask=32, fmt=None,
        country_set=None, oddballs=True, fname="", xlim=[.6, 1.0],
        maxmissing=0, rmask=16):
    '''
    :param t: int indicating the earliest query the window should include
    :param duration: int indication the span of time covered by the window,
        in seconds
    :param mask: int, prefix mask to use over domain IPs
    :param fmt: see transform fmt
    :param country_set: the set of countries the window should include queries from.
        If None, then all countries will be inluded
    :param method: the linkage method to be used
    :param p: dendrogram will show last p merged clusters
    :param oddballs: if True, will include non-public IPs (10.x.x.x, etc); if
        False, will only include public IPs
    :param fname: string to be appended to end of plot file name
    :param maxmissing: the maximum number of domains that can be missing
        measurements for a client to still be included in this measurement
    :returns: [country, ASN, subnet, prefix] pair dictionaries of closeness lists

    gets pairwise closeness of probes with different descriptors to find odd
    behavior (probes in difference descriptors with high closeness scores)

    NOTE: writes data to files for conveniece
    '''


    print("getting svl...")
    svl, fmt, _ = vv.get_svl(t, duration, mask, fmt, country_set, oddballs, maxmissing)
    logger.warning("svl len: "+str(len(svl)))
    ccache = vv.init_ccache(None, ccachef, t, duration, mask, fmt, oddballs, maxmissing)

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



