from warriorpy.shorthand import diriofile as df
from warriorpy.net_tools import dns_tools as dt
from warriorpy.net_tools import ipparsing as ipp
from warriorpy.shorthand import plotstuff as ps
from matplotlib import pyplot as plt
from math import log
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


def arrange_self_data(t, duration=30000, gap=1000, loops=2, mask=32,
        fmt=None, country_set=None, oddballs=True, fname="", ccache=None):

    svld = defaultdict(list) # dict {ip: [svl]}
    allsvl = list()
    allfmt = set()

    for l in xrange(0, loops):
        svl, fmt2 = vv.get_svl(t, duration, mask, fmt, country_set, oddballs)
        for i in xrange(0, len(svl)):
            svld[svl[i].ip].append(svl[i])
            allsvl.append(svl[i])
            allfmt |= set(fmt2)

    return svld, allsvl, list(allfmt)


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

    for ip in svld:
        svl = svld[ip]
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
    sm = defaultdict(list) # {dom: [self match value]}
                            # one self match value per client per domain
    for ip in svld:
        tmpsm = dict()
        for i in xrange(0, len(svld[ip])-1):
            for j in xrange(i+1, len(svld[ip])):
                a = svld[ip][i]
                b = svld[ip][j]
                for dom in [z for z in a if z in b]:
                    n = 2.0*float(a.vec[dom].symmetric_difference(b.vec[dom]))
                    d = float(len(a.vec[dom])+len(b.vec[dom]))
                    tmpsm[dom].append(n/d)
        for dom in tmpsm:
            sm[dom].append(np.mean(tmpsm[dom]))
    return sm

