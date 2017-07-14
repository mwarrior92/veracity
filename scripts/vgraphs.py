from warriorpy.shorthand import diriofile as df
from pymongo import MongoClient
import json
import numpy as np
import veracity_vector as vv
from warriorpy.net_tools import ipparsing as ipp
from matplotlib import pyplot as plt
import networkx as nx
from collections import defaultdict

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

def add_closeness_edges(svl, ccache, G=None):
    '''
    svl     -> list of smartvecs
    ccache  -> closeness_cache
    G       -> preexisting graph
    '''
    if G is None:
        G = nx.Graph()
    for i in xrange(0, len(svl)-1):
        for j in xrange(i+1, len(svl)):
            G.add_edge((svl[i].ip, i), (svl[j].ip, j),
                    weight=ccache[svl[i]][svl[j]])
    return G


def remove_far_edges(G, min_closeness):
    edges = [(a, b) for (a, b, d) in G.edges_iter(data=True) \
                                    if d['weight'] < min_closeness]

    G.remove_edges_from(edges)


def get_connected_comps(svl, min_closeness, ccache):
    '''
    return set of connected components
    '''
    G = add_closeness_edges(svl, ccache)
    remove_far_edges(G, min_closeness)
    return nx.connected_components(G)


def get_cc_varying_mc(svl, mcl, ccache=None):
    ccl = list()
    ccache = vv.init_ccache(ccache)
    for min_closeness in mcl:
        print min_closeness
        ccl.append(get_connected_comps(svl, min_closeness, ccache))
    return ccl


def csize_vs_mc(svl, mcl, ccache=None):
    X = list()
    Y = list()
    print "calling get_cc_varying_mc..."
    ccl = get_cc_varying_mc(svl, mcl, ccache)
    for i in xrange(0, len(mcl)):
        Y.append(max([len(z) for z in ccl[i]]))
        X.append(mcl[i])
    return X, Y


def ccount_vs_mc(svl, mcl, ccache=None):
    X = list()
    Y = list()
    ccl = get_cc_varying_mc(svl, mcl, ccache)
    for i in xrange(0, len(mcl)):
        tmp = list(ccl[i])
        Y.append(len(tmp))
        if len(tmp) > 1:
            for g in xrange(0, len(tmp)):
                print "************* group "+str(g)+" ****************"
                countries = [svl[z[1]].get_country() for z in tmp[g]]
                for c in set(countries):
                    count = float(len([z for z in countries if z==c]))
                    print c+": "+str(100.0*count/float(len(tmp[g])))+"%"
            print "---------------------------------------------\n"
        X.append(mcl[i])
    return X, Y


def plot_csize_vs_mc(t, duration=30000, mask=32, fmt=None, country_set=None,
        oddballs=True, fname=""):

    print "getting svl"
    svl, fmt = vv.get_svl(t, duration, mask, fmt, country_set, oddballs)
    print "calling csize vs mc"
    x, y = csize_vs_mc(svl, np.arange(.45, .65, .01))
    plt.figure(figsize=(15, 10))
    plt.xlabel('minimum closeness')
    plt.ylabel('size of biggest component')
    plt.plot(x,y)
    plt.savefig(plotsdir+'components_'+fname+'.pdf', bbox_inches='tight')


def plot_ccount_vs_mc(t, duration=30000, mask=32, fmt=None, country_set=None,
        oddballs=True, fname=""):

    svl, fmt = vv.get_svl(t, duration, mask, fmt, country_set, oddballs)
    x, y = ccount_vs_mc(svl, np.arange(.45, .48, .01))
    plt.figure(figsize=(15, 10))
    plt.xlabel('minimum closeness')
    plt.ylabel('# of components')
    plt.plot(x,y)
    plt.savefig(plotsdir+'components_'+fname+'.pdf', bbox_inches='tight')


def get_ip_sets(svl):
    ipsl = list() # ip set list
    dompairs = defaultdict(set)
    for sv in svl:
        ipsl.append(set())
        for dom in sv:
            for ip in sv[dom]:
                ipsl[-1].add(ip)
                dompairs[ip].add(dom)
    return ipsl, dompairs


def get_pairing_counts(ipsl):
    pc = defaultdict(lambda: defaultdict(int))
    for ips in ipsl:
        ipl = list(ips)
        for a in xrange(0, len(ipl)-1):
            for b in xrange(a+1, len(ipl)):
                pc[ipl[a]][ipl[b]] += 1
                pc[ipl[b]][ipl[a]] += 1
            pc[ipl[a]]['total'] += 1
        pc[ipl[-1]]['total'] += 1
    return pc


def build_inv_graph(pc):
    G = nx.Graph()
    for a in pc:
        for b in [z for z in pc[a] if z != 'total']:
            # 2 * number of probes that saw both a and b
            n = float(2*pc[a][b])
            # number of probes that saw a + number of probes that saw b
            d = float(pc[a]['total']+pc[b]['total'])
            G.add_edge(a, b, weight=n/d)
    return G


def nodes_by_degree(G):
    nd = G.degree()
    dd = defaultdict(list)
    for node in nd:
        dd[nd[node]].append(node)
    return dd


def remove_degree_above(G, dd, max_degree):
    for d in [z for z in dd if z > max_degree]:
        G.remove_nodes_from(dd[d])


def remove_degree_below(G, dd, min_degree):
    for d in [z for z in dd if z < min_degree]:
        G.remove_nodes_from(dd[d])


def inv_hist(t, duration=30000, mask=32, fmt=None, country_set=None,
        oddballs=True, fname="", thresh=.35):

    logger.info("getting svl...")
    svl, fmt = vv.get_svl(t, duration, mask, fmt, country_set, oddballs)
    logger.info("getting ipsl...")
    ipsl, dompairs = get_ip_sets(svl)
    logger.info("getting pairing counts...")
    pc = get_pairing_counts(ipsl)
    ipcount = len(pc)
    logger.info("building inv. graph...")
    G = build_inv_graph(pc)
    #dd = nodes_by_degree(G)
    #degree = int(.39*ipcount)
    #remove_degree_above(G, dd, degree)
    remove_far_edges(G, thresh)
    dd = nodes_by_degree(G)
    remove_degree_below(G, dd, 1)
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
