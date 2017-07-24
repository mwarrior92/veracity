from warriorpy.shorthand import diriofile as df
from warriorpy.shorthand import plotstuff as ps
from matplotlib import pyplot as plt
import numpy as np
import veracity_vector as vv
import metric_validation as mv
import vgraphs as vg
import networkx as nx
from collections import defaultdict
from statsmodels.distributions.empirical_distribution import ECDF

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

##################################################################
#                           CODE
##################################################################


def plot_closeness_diff_desc(t, duration=30000, mask=32, fmt=None,
        country_set=None, oddballs=True, fname="", ccache=None):
    print "getting svl..."
    svl, fmt = vv.get_svl(t, duration, mask, fmt, country_set, oddballs)
    print len(svl)

    print "getting descriptor lists..."
    csvl = vv.country_svl(svl)
    asvl = vv.asn_svl(svl)
    ssvl = vv.subnet_svl(svl)
    #osvl = vv.owner_svl(svl)
    psvl = vv.prefix_svl(svl)

    ccache = vv.init_ccache(ccache)

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
    ps.set_dim(fig, ax, xdim=13, ydim=7.5, xlim=[.6, 1.0])
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

    return ccache


def plot_closeness_same_desc(t, duration=30000, mask=32, fmt=None,
        country_set=None, oddballs=True, fname="", ccache=None):
    print "getting svl..."
    svl, fmt = vv.get_svl(t, duration, mask, fmt, country_set, oddballs)

    print "getting descriptor lists..."
    csvl = vv.country_svl(svl)
    asvl = vv.asn_svl(svl)
    ssvl = vv.subnet_svl(svl)
    #osvl = vv.owner_svl(svl)
    psvl = vv.prefix_svl(svl)

    ccache = vv.init_ccache(ccache)

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
    ps.set_dim(fig, ax, xdim=13, ydim=7.5, xlim=[.6, 1.0])
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

    return ccache


def plot_csize_vs_mc(t, duration=30000, mask=32, fmt=None, country_set=None,
        oddballs=True, fname=""):

    print "getting svl"
    svl, fmt = vv.get_svl(t, duration, mask, fmt, country_set, oddballs)
    print "calling csize vs mc"
    x, y = vg.csize_vs_mc(svl, np.arange(.45, .65, .01))
    plt.figure(figsize=(15, 10))
    plt.xlabel('minimum closeness')
    plt.ylabel('size of biggest component')
    plt.plot(x,y)
    plt.savefig(plotsdir+'components_'+fname+'.pdf', bbox_inches='tight')


def plot_ccount_vs_mc(t, duration=30000, mask=32, fmt=None, country_set=None,
        oddballs=True, fname=""):

    svl, fmt = vv.get_svl(t, duration, mask, fmt, country_set, oddballs)
    x, y = vg.ccount_vs_mc(svl, np.arange(.45, .48, .01))
    plt.figure(figsize=(15, 10))
    plt.xlabel('minimum closeness')
    plt.ylabel('# of components')
    plt.plot(x,y)
    plt.savefig(plotsdir+'components_'+fname+'.pdf', bbox_inches='tight')


def inv_hist(t, duration=30000, mask=32, fmt=None, country_set=None,
        oddballs=True, fname="", thresh=.35):

    logger.info("getting svl...")
    svl, fmt = vv.get_svl(t, duration, mask, fmt, country_set, oddballs)
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
        country_set=None, oddballs=True, fname="", ccache=None, loops=2,
        gap=1, thresh=10):
    '''
    lines:  1) domain independent cdf of ALL matches
            2-n) cdf of matches for domain with answer space > thresh
            n-m) cdf of matches for ALL domains with answer space < thresh
    '''
    svld, allsvl, allfmt = mv.arrange_self_data(t, duration, gap, loops, mask,
            fmt, country_set, oddballs)
    keys = svld.keys()
    anssets = vv.get_answer_space_dict(allsvl, allfmt)

    #anssets = vv.get_answer_space_dict(allsvl, allfmt)
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
        print "*****************"+labels[i]+"*********************"
        print vals[i]
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
        country_set=None, oddballs=True, fname="", ccache=None, loops=2,
        gap=0, thresh=10):
    '''
    lines:  1) domain independent cdf of ALL matches
            2-n) cdf of matches for domain with answer space > thresh
            n-m) cdf of matches for ALL domains with answer space < thresh
    '''
    svld, allsvl, allfmt = mv.arrange_self_data(t, duration, gap, loops, mask,
            fmt, country_set, oddballs)
    keys = svld.keys()
    anssets = vv.get_answer_space_dict(allsvl, allfmt)

    #anssets = vv.get_answer_space_dict(allsvl, allfmt)
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
        print "*****************"+labels[i]+"*********************"
        print vals[i]
        ecdf = ECDF(vals[i])
        x = list(ecdf.x)
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
        country_set=None, oddballs=True, fname="", ccache=None, loops=2,
        gap=0, thresh=10):
    '''
    lines:  1) domain independent cdf of ALL matches
            2-n) cdf of matches for domain with answer space > thresh
            n-m) cdf of matches for ALL domains with answer space < thresh
    '''
    svld, allsvl, allfmt = mv.arrange_self_data(t, duration, gap, loops, mask,
            fmt, country_set, oddballs)
    keys = svld.keys()
    anssets = vv.get_answer_space_dict(allsvl, allfmt)

    #anssets = vv.get_answer_space_dict(allsvl, allfmt)
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
        country_set=None, oddballs=True, fname="", ccache=None, loops=10,
        gap=0, thresh=10):
    svld, allsvl, allfmt = mv.arrange_self_data(t, duration, gap, loops, mask,
            fmt, country_set, oddballs)
    keys = svld.keys()
    anssets = vv.get_answer_space_dict(allsvl, allfmt)

    vals = [[], []]
    labels = ['all', 'small']
    for dom in sm:
        vals[0] = vals[0] + sm[dom]
        if len(anssets[dom]) < thresh:
            vals[1] += sm[dom]
        else:
            vals.append(sm[dom])
            labels.append(dom)
    counts = mv.measure_expansion(svld):
    for c in counts:
        for dom in c:
            pass

    return

