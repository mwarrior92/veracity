from warriorpy.shorthand import diriofile as df
from warriorpy.shorthand import plotstuff as ps
from matplotlib import pyplot as plt
import numpy as np
import veracity_vector as vv
import metric_validation as mv
import clientcompare as cc
import matplotlib.cm as cm
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
import copy

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


def plot_fmeasure(start_time, method="complete", fname="", Zf=False, **kwas):
    '''
    :param start_time: int indicating the earliest query the window should include
    :param method: the linkage method to be used
    :param fname: string to be appended to end of plot file name
    :param **kwas: keyword arguments for vv.get_svl()

    scatter plot:
        x -> max distance threshold
        y -> f-measure
        lineplot -> # components (y)

    other output:
        list of desc. that shared optimal components
    '''

    Z, svl = get_zx(start_time, method, fname, Zf, **kwas)
    dsvl = dict()
    dsvl['country'] = vv.country_svl(svl)
    dsvl['asn'] = vv.asn_svl(svl)
    dsvl['subnet'] = vv.subnet_svl(svl)
    dsvl['prefix'] = vv.prefix_svl(svl)
    dsvl['ldns'] = vv.ldns_svl(svl, rmask, False)
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
    vals = defaultdict(list)
    grouping = dict()
    count = list()
    for max_dist in np.arange(0, 1.01, .01):
        data = defaultdict(lambda: defaultdict(list))
        labels = fcluster(Z, max_dist, criterion='distance')
        clusters = [[(c, svl[z]) for z, c in enumerate(labels) if c == y] \
                for y in set(labels)]
        if len(clusters) > 1 and len(clusters) < len(svl):
            count.append(max_dist, len(clusters))
            for c, blob in clusters:
                for desc in dsvl:
                    cluster = [getattr(sv,"get_"+desc)() for sv in blob]
                    for d in set(cluster):
                        localcount = float(len([z for z in cluster if z == d]))
                        localsize = float(len(cluster))
                        globalsize = float(len(dsvl[desc][d]))
                        precision = localcount / localsize
                        recall = localcount / globalsize
                        fmeasure = (2*precision*recall)/(precision+recall)
                        data[desc][d].append((fmeasure, c, max_dist))
    for desc in data:
        for d in data[desc]:
            maxf, maxc, maxd = max(data[desc][d], key=lambda z: z[0])
            vals[desc].append((maxf, maxd))
            grouping[(desc, maxc, maxd)].append((d, maxf))

    print "plotting..."
    fig, ax = plt.subplots(1, 1)

    vals['resolver'] = vals.pop('ldns')
    colors = iter(cm.rainbow(np.linspace(0, 1, len(vals))))
    for desc in vals:
        y, x = zip(*vals[desc])
        heatmap, xedges, yedges = np.histogram2d(x, y, bins=50)
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        plt.clf()
        plt.imshow(heatmap.T, extent=extent, origin='lower')
        cb = PLT.colorbar()
        cb.set_label('# of '+make_plural(desc))
        plt.xlabel("max distance threshold")
        plt.ylabel("F-measure")
        ax.grid(b=True, which='both', color='b', linestyle='-')
        filename = plotsdir+"fmeasure_"+desc+fname
        fig.savefig(filename+'.png', bbox_inches='tight')
        fig.savefig(filename+'.pdf', bbox_inches='tight')
        plt.close(fig)

    fig, ax = plt.subplots(1, 1)
    x, y = zip(*count)
    ax.plot(x, y)
    plt.xlabel("max distance threshold")
    ax.set_ylabel('# components')
    ps.set_dim(fig, ax, xdim=13, ydim=7.5, xlim=xlim)
    lgd = ps.legend_setup(ax, 4, "top center", True)
    filename = plotsdir+"component_count"+fname
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    print "saving data..."
    df.overwrite(plotsdir+'fmeasure_groups'+fname+'.csv',
        df.list2col(vals))

    groups = [(z[0], groups[z]) for z in groups if len(groups[z]) > 1]
    if len(groups) > 0:
        groups = sorted(groups, key=lambda z: z[0]+str(z[1])+str(z[2]))
        df.overwrite(plotsdir+"fmeasure_groups"+fname+".csv", df.list2col(groups))


def make_plural(desc):
        if desc == 'country':
            return 'countries'
        if desc == 'subnet':
            return 'subnets'
        if desc == 'asn':
            return 'ASNs'
        if desc == 'prefix':
            return 'prefixes'
        if desc == 'resolver' or desc == 'ldns':
            return 'resolvers'


def change_in_components(start_time, duration, comp_count=100.0,
        method="complete", fname="", loops=15, Zf=False, **kwas):
    '''
    :param start_time: int indicating the earliest query the window should include
    :param method: the linkage method to be used
    :param fname: string to be appended to end of plot file name
    :param **kwas: keyword arguments for vv.get_svl()
    '''
    component_tuples = dict()
    kwas['duration'] = duration
    for i in xrange(0, loops):
        Z, svl = get_zx(start_time, method, fname, Zf, **kwas)
        best = -1
        best_set = None
        prev = best
        for max_dist in np.arange(0, 1.01, .01):
            data = defaultdict(lambda: defaultdict(list))
            labels = fcluster(Z, max_dist, criterion='distance')
            clusters = [(str(y)+"_"+str(max_dist), set([svl[z].get_id() for z,
                c in enumerate(labels) if c == y])) for y in set(labels)]
            count = float(len([z for z in clusters if len(z)>1]))
            if abs(1-(count/comp_count)) > best:
                best = abs(1-(count/comp_count))
                best_set = (clusters, max_dist, i)
            # stop it when it starts moving away from target
            if count > prev and prev > comp_count:
                break
            prev = count
        component_tuples[best_set[2]] = best_set[:2]
        start_time += duration

    df.pickleout(plotsdir+"component_tuples"+fname+".pickle",
            component_tuples)

    for i in xrange(0, len(component_tuples)-1):
        listA = component_tuples[i][0]
        compA = dict()
        labelA = dict()
        for j in xrange(0, len(listA)):
            labelA[j], compA[j] = zip(*listA[j])
            compA[j] = set(compA[j])
        for k in xrange(i+1, len(component_tuples)):
            listB = component_tuples[k][0]
            for l in xrange(0, len(listB)):
                labelB[l], compB[l] = zip(*listB[l])
            for m in xrange(0, len(listA)):
                for n in xrange(0, len(listB)):
                    anb = float(len(compA[m].intersection(compB[n])))
                    aub = float(len(compA[m].union(compB[n])))
                    jac = anb/aub
                    if jac > M[(i,k, m)][0]: # [(runA ind, runB ind, compA ind)]
                        M[(i,k,m)] = (jac, n) # (jaccard, compB ind)
                                              # use component_tuples.pickle to
                                              # map indices to data
                    if jac == 1.0:
                        break
                M2[k-i].append(m[(i,k,m)][0])
    df.pickleout(plotsdir+"component_matches"+fname+".pickle", M)
    df.pickleout(plotsdir+"component_changes"+fname+".pickle", M2)

    x = list()
    y = list()
    std = list()
    for dist in M2:
        x.append(dist)
        y.append(np.mean(M2[dist]))
        std.append(np.std(M2[dist]))

    fig, ax = plt.subplots(1, 1)
    ax.errorbar(x, y, yerr=std)
    ps.set_dim(fig, ax, xdim=13, ydim=7.5, xlim=xlim)
    plt.xlabel("distance in days")
    plt.ylabel("mean component match")
    filename = plotsdir+"change_in_components"+fname
    fig.savefig(filename+'.png', bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_inches='tight')
    plt.close(fig)


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
        df.pickleout(plotsdir+'pickles/'+'Z_'+method+fname+'.pickle', (Z, dm, X))
        logger.warning('dumped Z to ' \
                +plotsdir+'pickles/'+'Z_'+method+fname+'.pickle')
    else:
        Z, dm, X = df.picklein(plotsdir+'pickles/'+'Z_'+method+fname+'.pickle')
        logger.warning('loaded Z from '+plotsdir+'pickles/'+'Z_'+method+fname+'.pickle')
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
    outstr = df.overwrite(plotsdir+fname+'_avg_self_closeness.csv',
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

    NOTE: plot 3.1
    '''
    means = defaultdict(list)
    vals = list()
    kwas['duration'] = duration
    for l in xrange(0, loops):
        print "getting svl..."
        kwas['start_time'] = start_time+duration*l
        svl, __, __, ccache = vv.get_svl(**kwas)
        logger.warning("svl len: "+str(len(svl)))
        print len(svl)

        print "calculating closeness for resolvers..."
        for i in xrange(0, len(svl)-1):
            for j in xrange(i + 1, len(svl)):
                vals.append(ccache[svl[i]][svl[j]])
                means[svl[i].get_id()].append(vals[-1])
                means[svl[j].get_id()].append(vals[-1])
        ccache.dump()
        del ccache, svl, __
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
    filename = plotsdir+"overall_closeness"+fname
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    print "saving data..."
    df.overwrite(plotsdir+'overall_closeness'+fname+'.csv',
        df.list2col(vals))
    df.overwrite(plotsdir+'overall_avg_closeness'+fname+'.csv',
        df.list2col([(z, np.mean(means[z])) for z in means]))


def plot_closeness_diff_desc(start_time, duration, fname="", xlim=[.6, 1.0],
        rmask=16, loops=31, **kwas):
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

    NOTE: plot 4.2
    '''
    lvals = list()
    cvals = list()
    avals = list()
    svals = list()
    pvals = list()
    kwas['duration'] = duration
    for l in xrange(0, loops):
        print "getting svl..."
        kwas['start_time'] = start_time+duration*l
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
        resolvers = [c for c in lsvl if len(lsvl[c]) > 1]
        for i in xrange(0, len(resolvers)-1):
            for a in lsvl[resolvers[i]]:
                for j in xrange(i+1, len(resolvers)):
                    for b in lsvl[resolvers[j]]:
                        lvals.append(ccache[a][b])

        print "calculating closeness for countries..."
        countries = [c for c in csvl if len(csvl[c]) > 1]
        for i in xrange(0, len(countries)-1):
            for a in csvl[countries[i]]:
                for j in xrange(i+1, len(countries)):
                    for b in csvl[countries[j]]:
                        cvals.append(ccache[a][b])
        print "calculating closeness for ASNs..."
        asns = [a for a in asvl if len(asvl[a]) > 1]
        for i in xrange(0, len(asns)-1):
            for a in asvl[asns[i]]:
                for j in xrange(i+1, len(asns)):
                    for b in asvl[asns[j]]:
                        avals.append(ccache[a][b])
        print "calculating closeness for subnets..."
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
        outstr = df.overwrite(plotsdir+labels[i]+'_diff.csv',
                df.list2col(vals[i]))
    ccache.dump()


def plot_closeness_same_desc(start_time, duration, fname="", xlim=[.6, 1.0], rmask=16,
        loops=31, **kwas):
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

    NOTE: plot 4.1
    '''
    lvals = list()
    cvals = list()
    avals = list()
    svals = list()
    pvals = list()
    kwas['duration'] = duration
    for l in xrange(0, loops):
        print "getting svl..."
        kwas['start_time'] = start_time+duration*l
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
        resolvers = lsvl.keys()
        for k in resolvers:
            ksvl = lsvl[k]
            for a in xrange(0, len(ksvl)-1):
                for b in xrange(a+1, len(ksvl)):
                    lvals.append(ccache[ksvl[a]][ksvl[b]])

        print "calculating closeness for countries..."
        countries = csvl.keys()
        for k in countries:
            ksvl = csvl[k]
            for a in xrange(0, len(ksvl)-1):
                for b in xrange(a+1, len(ksvl)):
                    cvals.append(ccache[ksvl[a]][ksvl[b]])
        print "calculating closeness for ASNs..."
        asns = asvl.keys()
        for k in asns:
            ksvl = asvl[k]
            for a in xrange(0, len(ksvl)-1):
                for b in xrange(a+1, len(ksvl)):
                    avals.append(ccache[ksvl[a]][ksvl[b]])
        print "calculating closeness for subnets..."
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
        print type(vals[i][0])
        print labels[i], "\n"
        print len(vals[i])
        ecdf = ECDF(np.array(vals[i]))
        x = list(ecdf.x)
        y = list(ecdf.y)
        ax.plot(x, y, label=labels[i])
    ps.set_dim(fig, ax, xdim=13, ydim=7.5)
    plt.xlabel("pairwise probe closeness")
    plt.ylabel("CDF of pairs")
    lgd = ps.legend_setup(ax, 4, "top center", True)
    filename = plotsdir+"closeness_same_desc"+fname
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    print "saving data..."
    for i in xrange(0, len(vals)):
        outstr = df.overwrite(plotsdir+labels[i]+'_same.csv',
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
    :param start_time: int indicating the earliest query the window should include
    :param gap: the gap (in seconds) between each iteration's dataset
    :param loops: the number of iterations (datasets) / 2 (since we need 2 per
        comparison in this case)
    :param **kwas: keyword arguments for vv.get_svl()

    plots CDF:
        x -> closeness to self for back-to-back iteration windows (e.g., days 1-2
        vs days 3-4, days 5-6 vs days 7-8, ...)
        y -> CDF of clients

    NOTE: plot 3.2
    '''
    valsd = defaultdict(list) # {pid: [vals]}
    bigsvld = dict()
    kwas['duration'] = duration
    kwas['return_ccache'] = False
    for i in xrange(0, loops):
        print (start_time+2*duration*i, duration)
        tmp_start = start_time+2*(gap+duration)*i
        svld, allsvl, allfmt, anssets = mv.arrange_self_data(tmp_start, gap, 2, **kwas)

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
        df.overwrite(plotsdir+'self_closeness_'+k+fname+'.csv',
                df.list2col(data))

    print "saving data..."
    outstr = df.overwrite(plotsdir+'self_closeness'+fname+'.csv',
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
        outstr = df.overwrite(plotsdir+labels[i]+'_self_jaccard.csv',
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
        outstr = df.overwrite(plotsdir+labels[i]+'_diff_jaccard.csv',
                df.list2col(vals[i]))


def plot_measure_expansion(start_time, fname="", loops=31, gap=0, thresh=10,
        **kwas):
    '''
    :param start_time: int indicating the earliest query the window should include
    :param gap: the gap (in seconds) between each iteration's dataset
    :param loops: the number of iterations (datasets)
    :param **kwas: keyword arguments for vv.get_svl()

    line plot:
        x -> n (such that it represents the nth iteration)
        y -> # of new IPs observed by a client on nth iteration
        line -> each line corresponds to one domain
    '''

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
    marker = ps.get_markers()
    style = ps.get_styles()
    for i in xrange(0, len(vals)):
        ax.plot(vals[i], label=labels[i], fillstyle='full', marker=next(marker),
                markerfacecolor='white', markevery=6, linestyle=next(style))
    ps.set_dim(fig, ax)
    plt.xlabel("cycle #")
    plt.ylabel("# new IPs / ans. size")
    ax.grid(b=True, which='major', color='b', linestyle='-')
    lgd = ps.legend_setup(ax, 3, "top right", True)
    filename = plotsdir+"expansion"+fname
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    print "saving data..."
    for i in xrange(0, len(vals)):
        outstr = df.overwrite(plotsdir+labels[i]+'newvssize.csv',
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



