from warriorpy.shorthand import diriofile as df
import numpy as np
import veracity_vector as vv
from warriorpy.net_tools import ipparsing as ipp
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

def add_closeness_edges(svl, ccache, G=None, minweight=0):
    '''
    :param svl: list of smartvecs
    :param ccache: closeness_cache
    :param G: preexisting graph
    :result: graph

    builds graph using pairwise closeness as edge weights
    '''
    if G is None:
        G = nx.Graph()
    for i in xrange(0, len(svl)-1):
        for j in xrange(i+1, len(svl)):
            w = ccache[svl[i]][svl[j]]
            if w >= minweight:
                G.add_edge(svl[i], svl[j], weight=w)
    return G


def remove_far_edges(G, min_weight):
    '''
    :param G: graph
    :param min_weight: float, all edges below this value will be removed
    :return: updated graph
    '''
    edges = [(a, b) for (a, b, d) in G.edges_iter(data=True) \
                                    if d['weight'] < min_weight]

    G.remove_edges_from(edges)


def component_analysis(c):
    countries = list()
    prefixes = list()
    resolvers = list()
    subnets = list()
    asns = list()
    for sv in c:
        countries.append(sv.get_country())
        prefixes.append(sv.get_prefix())
        resolvers.append(sv.get_ldns())
        subnets.append(sv.get_subnet(24))
        asns.append(sv.get_asn())

    distro = [(z, countries.count(z)) for z in set(countries)]
    country = max(distro, key=lambda z: z[1])
    country = (country[0], float(country[1])/float(len(countries)))

    distro = [(z, prefixes.count(z)) for z in set(prefixes)]
    prefix = max(distro, key=lambda z: z[1])
    prefix = (prefix[0], float(prefix[1])/float(len(countries)))

    distro = [(z, resolvers.count(z)) for z in set(resolvers)]
    resolver = max(distro, key=lambda z: z[1])
    resolver = (resolver[0], float(resolver[1])/float(len(countries)))

    distro = [(z, subnets.count(z)) for z in set(subnets)]
    subnet = max(distro, key=lambda z: z[1])
    subnet = (subnet[0], float(subnet[1])/float(len(countries)))

    distro = [(z, asns.count(z)) for z in set(asns)]
    asn = max(distro, key=lambda z: z[1])
    asn = (asn[0], float(asn[1])/float(len(countries)))

    return country, prefix, resolver, subnet, asn


def get_cc_varying_mc(svl, mwl, ccache):
    '''
    :param svl: list of smartvecs
    :param mwl: list of min_weight thresholds for which edges are kept
    :param ccache: closeness_cache
    :return: list of list of connected components, ordered such that each item
        in the outer list corresponds to the min_weight in the same position
        in mwl
    '''
    ccl = list()
    G = add_closeness_edges(svl, ccache, minweight=mwl[0])
    cc = nx.connected_components(G)

    countries = list()
    prefixes = list()
    resolvers = list()
    subnets = list()
    asns = list()
    lens = list()
    count = 0
    tinycs = list()
    for c in cc:
        if len(c) > 1:
            country, prefix, resolver, subnet, asn = component_analysis(c)
            countries.append(country[1])
            prefixes.append(prefix[1])
            resolvers.append(resolver[1])
            subnets.append(subnet[1])
            asns.append(asn[1])
        lens.append(len(c))
        count += 1
        # remove tiny components to avoid watering down averages
        if len(c) < 3:
            tinycs += list(c)
    G.remove_nodes_from(tinycs)
    vals = {
        'country': np.mean(countries),
        'prefix': np.mean(prefixes),
        'resolver': np.mean(resolvers),
        'subnet': np.mean(subnets),
        'asn': np.mean(asns),
        'len': np.median(lens, overwrite_input=True),
        'min_closeness': mwl[0],
        'quantity': count}
    print vals
    ccl.append(vals)


    for min_closeness in mwl[1:]:
        count = 0
        remove_far_edges(G, min_closeness)
        cc = nx.connected_components(G)

        countries = list()
        prefixes = list()
        resolvers = list()
        subnets = list()
        asns = list()
        lens = list()
        tinycs = list()

        for c in cc:
            if len(c) > 1:
                country, prefix, resolver, subnet, asn = component_analysis(c)
                countries.append(country[1])
                prefixes.append(prefix[1])
                resolvers.append(resolver[1])
                subnets.append(subnet[1])
                asns.append(asn[1])
            lens.append(len(c))
            count += 1
            if len(c) < 3:
                tinycs += list(c)
        G.remove_nodes_from(tinycs)
        vals = {
            'country': np.mean(countries),
            'prefix': np.mean(prefixes),
            'resolver': np.mean(resolvers),
            'subnet': np.mean(subnets),
            'asn': np.mean(asns),
            'len': np.median(lens, overwrite_input=True),
            'min_closeness': min_closeness,
            'quantity': count}
        print vals
        ccl.append(vals)
    return ccl


def get_ip_sets(svl):
    '''
    :param svl: list of smartvecs
    :return: (list of sets of IPs, dict {ip: set(domains)})

    each item in ipsl is the set of IPs seen by one probe
    '''
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
    '''
    :param ipsl: output from get_ip_sets()
    :return: dict {ip_A: {ip_B: # times ip_A and ip_B were seen by same probe}}
    '''
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
    '''
    :param pc: output from get_pairing_counts()
    :return: graph of ip associations, based on pc
    '''
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
    '''
    :param G: graph
    :return: dict {degree: [nodes with degree==degree]}

    NOTE: this should always be run immediately before dd will be used
    '''
    nd = G.degree()
    dd = defaultdict(list)
    for node in nd:
        dd[nd[node]].append(node)
    return dd


def remove_degree_above(G, dd, max_degree):
    '''
    :param G: graph
    :param dd: output from nodes_by_degree()
    :param max_degree: will remove nodes with degree higher than this
    '''
    for d in [z for z in dd if z > max_degree]:
        G.remove_nodes_from(dd[d])


def remove_degree_below(G, dd, min_degree):
    '''
    :param G: graph
    :param dd: output from nodes_by_degree()
    :param min_degree: will remove nodes with degree lower than this
    '''
    for d in [z for z in dd if z < min_degree]:
        G.remove_nodes_from(dd[d])


