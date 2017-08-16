from warriorpy.shorthand import diriofile as df
from warriorpy.net_tools import dns_tools as dt
from warriorpy.net_tools import ipparsing as ipp
from warriorpy.shorthand import plotstuff as ps
from matplotlib import pyplot as plt
from math import log
from math import factorial as fact
import numpy as np
import veracity_vector as vv
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


def get_domain_matrix(start_time, **kwas):
    '''
    :param start_time: int indicating the earliest query the window should include
    :param **kwas: keyword arguments for vv.get_svl()
    :return: (m) matrix of client pairs vs domains,
             (fmt) list of domains

    other outputs:
        -> csv with pairs vs domains matrix (m)
        -> csv with list of domain pair correlations (corrs)
        -> csv with list of mean Jaccard for each domain (means)
    '''
    kwas['start_time'] = start_time
    kwas['return_ccache'] = False
    svl, fmt, anssets = vv.get_svl(**kwas)
    combs = fact(len(svl))/(fact(2)*fact(len(svl)-2))
    m = np.zeros(combs, len(fmt)))
    p = 0
    for i in xrange(0, len(svl)-1):
        a = svl[i]
        for j in xrange(i+1, len(svl)):
            b = svl[j]
            for k in xrange(0, len(fmt)):
                domtotal = sum([a[dom][z] for z in a[dom]]+[b[dom][z] for z in b[dom]])
                overlap = set(a[dom]).intersection(set(b[dom]))
                aweight = [a[dom][z] for z in a[dom] if z in overlap]
                bweight = [b[dom][z] for z in b[dom] if z in overlap]
                m[p,k] = sum(aweight+bweight)/domtotal
            p += 1

    df.overwrite(plotsdir+"dommatrix"+fname+".csv", df.list2col(fmt))
    df.append(plotsdir+"dommatrix"+fname+".csv", df.list2col(m))

    C = np.corrcoef(m, rowvar=False)
    corrs = list()
    for i in xrange(0, len(fmt)-1):
        for j in xrange(i+1, len(fmt)):
            corrs.append((fmt[i]+"_"+fmt[j], C[i, j]))
    corrs = sorted(corrs, key=lambda z: z[1])
    means = sorted(fmt, zip(np.mean(m, axis=0)), key=lambda z: z[1])

    df.overwrite(plotsdir+"domcorr"+fname+".csv",df.list2col(corrs))
    df.overwrite(plotsdir+"dommean"+fname+".csv",df.list2col(means))

    meand = dict(means)
    # get mean jaccard vs # IPs seen
    mj_ni = [(meand[dom], len(anssets[dom])) for dom in meand]

    fig, ax = plt.subplots(1, 1)

    colors = iter(cm.rainbow(np.linspace(0, 1, len(mj_ni))))
    for x, y in mj_ni:
        ax.scatter(x, y, color=next(colors))
    plt.xlabel("mean jaccard")
    plt.ylabel("# IPs observed")
    plt.grid()
    ax2 = ax.twinx()
    x, y = zip(*count)
    ax2.plot(x, y)
    ax2.set_ylabel('# components')
    ps.set_dim(fig, ax, xdim=13, ydim=7.5, ylog=True)
    filename = plotsdir+"jaccard_vs_ipspace"+fname
    fig.savefig(filename+'.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(filename+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)

    return m, fmt


def arrange_self_data(start_time, gap=0, loops=2, **kwas):
    '''
    :param start_time: int indicating the earliest query the window should include
    :param gap: the gap (in seconds) between each iteration's dataset
    :param loops: the number of iterations (datasets)
    :param **kwas: keyword arguments for vv.get_svl()
    '''
    svld = defaultdict(list) # dict {id: [svl]}
    allsvl = list()
    allfmt = set()
    anssets = defaultdict(set)

    kwas['return_ccache'] = False
    for l in xrange(0, loops):
        kwas['start_time'] = start_time+l*(gap+duration)
        svl, fmt2, anssets2 = vv.get_svl(**kwas)
        logger.warning("svl len: "+str(len(svl)))
        allfmt |= set(fmt2)
        for dom in anssets2:
            anssets[dom] |= set(anssets2[dom])
        for i in xrange(0, len(svl)):
            svld[svl[i].get_id()].append(svl[i])
            allsvl.append(svl[i])

    return svld, allsvl, list(allfmt), anssets


def dom_traits(svld, anssets):
    '''
    :param svld: output from arrange_self_data()
    :return: dict {<trait>: {<trait val>: {dom: [clients]}}}, where 'clients'
    refers to the list of clients that witnessed said trait value from said
    domain

    traits include:
        ansm -> order of magnitude (bins) of the number of answers
        spacem -> order of magnitude (log bins) of the number of IPs seen
            across all probes
        prefix -> number of bits in registered prefix (TODO)
        shared -> bool: does dom'full IP set intersect with that of another dom
            NOTE: this is technically difficult to know; for simplicity, I'll
            start by only comparing to what I see from my data; manual
            inspection or something else is necessary to really say this or not
        lpm -> longest prefix match between IPs in dom's full answer space
        spm -> shorted prefix match between IPs in dom's full answer space
        pilpm -> the (average) longest prefix match for EACH IP. For example,
            the 4 bit IPs [1100, 1101, 0110, 0111] would have a pilpm of 3: the
            some IP's match '110', and the rest match '011', which both have
            a length of 3 bits and in turn average to (3+3+3+3/4) = 3 bits
    '''
    # {<trait>: {<trait val>: {dom: [clients]}}}
    dtd = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    checked = set()
    for a in anssets:
        if len(anssets[a]) == 0:
            continue
        checked.add(a)
        dtd['spacem'][int(round(log(len(anssets[a]), 2)))][a]
        for b in [z for z in anssets if z not in checked]:
            overlap = anssets[a].symmetric_difference(anssets[b])
            overlap = [ip for ip in overlap if ipp.is_public(ip)]
            if len(overlap) > 0:
                dtd['shared'][True][a] = True
                dtd['shared'][True][b] = True
        if a not in dtd['shared'][True]:
            dtd['shared'][False][a] = True
        mpms = list()
        for ip in anssets[a]:
            matches = [ipp.prefix_match(ip, z) for z in anssets[a] if \
                                z != ip]
            mpms.append(max(matches))
        dtd['pilpm'][int(round(log(np.mean(mpms),2)))][a] = True

    for pid in svld:
        svl = svld[pid]
        for sv in svl:
            for dom in sv:
                dtd['ansm'][len(sv.vec[dom])][dom].append(sv)
                ipstrs = [ipp.int2ip(z) for z in sv.vec[dom]]
                cidrs = [dt.get_cidr(ip) for ip in ipstrs if ipp.is_public(ip)]
                masks = [int(cidr.split("/")[1]) for cidr in cidrs if cidr is \
                            not None]
                for mask in masks:
                    dtd['prefix'][mask][dom].append(sv)
                pms = list()
                for ip in sv.vec[dom]:
                    matches = [ipp.prefix_match(ip, z) for z in anssets[dom] if \
                                z != ip]
                    pms += matches
                dtd['lpm'][max(pms)][dom].append(sv)
                dtd['spm'][min(pms)][dom].append(sv)


    return dtd


def self_match(svld):
    results = list()
    pids = svld.keys()
    for pid in pids:
        tmpsm = defaultdict(list)
        for i in xrange(0, len(svld[pid])-1):
            for j in xrange(i+1, len(svld[pid])):
                A = svld[pid][i]
                B = svld[pid][j]
                results.append(vv.closeness(A, B))
    return pids, results


def measure_expansion(svld):
    allcounts = list()
    for pid in svld:
        count = defaultdict(dict)
        for sv in svld[pid]:
            for dom in sv:
                ips = set(sv.vec[dom])
                if dom not in count:
                    count[dom]['ips'] = set()
                    count[dom]['new_ips'] = list()
                    count[dom]['old_ips'] = list()
                    count[dom]['size'] = list()
                    count[dom]['ratio'] = list()
                intsn = count[dom]['ips'].intersection(ips)
                count[dom]['new_ips'].append(len(ips) - len(intsn))
                count[dom]['old_ips'].append(len(intsn))
                count[dom]['ips'] |= ips
                count[dom]['size'].append(len(count[dom]['ips']))
                count[dom]['ratio'].append(float(count[dom]['new_ips'][-1])/float(len(ips)))
        allcounts.append(count)

    return allcounts


def examine_self_diff(svld):
    sm = defaultdict(list) # {dom: [self match value]}
                            # one self match value per client per domain
    for pid in svld:
        tmpsm = defaultdict(list)
        for i in xrange(0, len(svld[id])-1):
            for j in xrange(i+1, len(svld[pid])):
                A = svld[pid][i]
                B = svld[pid][j]
                for dom in [z for z in A if z in B]:
                    a = set(A.vec[dom])
                    b = set(B.vec[dom])
                    n = a.symmetric_difference(b)
                    if len(n) > 0:
                        matches = list()
                        for ip in n:
                            if ip not in a:
                                for ip2 in [z for z in a]:
                                    matches.append(ipp.prefix_match(ip, ip2))
                            elif ip not in b:
                                for ip2 in [z for z in b]:
                                    matches.append(ipp.prefix_match(ip, ip2))
                            tmpsm[dom].append(max(matches))
        for dom in tmpsm:
            sm[dom].append(np.mean(tmpsm[dom]))
    return sm


def examine_diff_diff(svld):
    sm = defaultdict(list) # {dom: [self match value]}
                            # one self match value per client per domain
    ids = svld.keys()
    for p in xrange(0, len(ids)-1):
        for q in xrange(p+1, len(ids)):
            tmpsm = defaultdict(list)
            for i in xrange(0, len(svld[ids[p]])):
                for j in xrange(0, len(svld[ids[q]])):
                    A = svld[ids[p]][i]
                    B = svld[ids[q]][j]
                    for dom in [z for z in A if z in B]:
                        a = set(A.vec[dom])
                        b = set(B.vec[dom])
                        n = a.symmetric_difference(b)
                        if len(n) > 0:
                            matches = list()
                            for ip in n:
                                if ip not in a:
                                    for ip2 in [z for z in a]:
                                        matches.append(ipp.prefix_match(ip, ip2))
                                elif ip not in b:
                                    for ip2 in [z for z in b]:
                                        matches.append(ipp.prefix_match(ip, ip2))
                                tmpsm[dom].append(max(matches))
        for dom in tmpsm:
            sm[dom].append(np.mean(tmpsm[dom]))
    return sm


def examine_self_diff(svld):
    sm = defaultdict(list) # {dom: [self match value]}
                            # one self match value per client per domain
    for pid in svld:
        tmpsm = defaultdict(list)
        for i in xrange(0, len(svld[pid])-1):
            for j in xrange(i+1, len(svld[pid])):
                A = svld[pid][i]
                B = svld[pid][j]
                for dom in [z for z in A if z in B]:
                    a = set(A.vec[dom])
                    b = set(B.vec[dom])
                    n = a.symmetric_difference(b)
                    if len(n) > 0:
                        matches = list()
                        for ip in n:
                            if ip not in a:
                                for ip2 in [z for z in a]:
                                    matches.append(ipp.prefix_match(ip, ip2))
                            elif ip not in b:
                                for ip2 in [z for z in b]:
                                    matches.append(ipp.prefix_match(ip, ip2))
                            tmpsm[dom].append(max(matches))
        for dom in tmpsm:
            sm[dom].append(np.mean(tmpsm[dom]))
    return sm


def examine_diff_diff(svld):
    sm = defaultdict(list) # {dom: [self match value]}
                            # one self match value per client per domain
    ids = svld.keys()
    for p in xrange(0, len(ids)-1):
        for q in xrange(p+1, len(ids)):
            tmpsm = defaultdict(list)
            for i in xrange(0, len(svld[ids[p]])):
                for j in xrange(0, len(svld[ids[q]])):
                    A = svld[ids[p]][i]
                    B = svld[ids[q]][j]
                    for dom in [z for z in A if z in B]:
                        a = set(A.vec[dom])
                        b = set(B.vec[dom])
                        n = a.symmetric_difference(b)
                        if len(n) > 0:
                            matches = list()
                            for ip in n:
                                if ip not in a:
                                    for ip2 in [z for z in a]:
                                        matches.append(ipp.prefix_match(ip, ip2))
                                elif ip not in b:
                                    for ip2 in [z for z in b]:
                                        matches.append(ipp.prefix_match(ip, ip2))
                                tmpsm[dom].append(max(matches))
        for dom in tmpsm:
            sm[dom].append(np.mean(tmpsm[dom]))
    return sm

