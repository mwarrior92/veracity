from __future__ import print_function
from warriorpy.shorthand import diriofile as df
from warriorpy.shorthand import easymath as em
from warriorpy.shorthand import plotstuff as ps
from matplotlib import pyplot as plt
import numpy as np
import veracity_vector as vv
import metric_validation as mv
import vgraphs as vg
import networkx as nx
from collections import defaultdict
from statsmodels.distributions.empirical_distribution import ECDF
from warriorpy.net_tools import ipparsing as ipp
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
import sys
import math

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


def distance(c, d): # c: closeness, d: physical distance (km)
    return math.sqrt(c**2 + math.log(max([d,1.0]), 40075)**2)


def closest_diff_desc(t, duration=30000, mask=32, fmt=None,
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
    :returns: [country, ASN, subnet, prefix] pair dictionaries of closeness lists

    gets pairwise closeness of probes with different descriptors to find odd
    behavior (probes in difference descriptors with high closeness scores)

    NOTE: writes data to files for conveniece
    '''
    print("getting svl...")
    svl, fmt, _ = vv.get_svl(t, duration, mask, fmt, country_set, oddballs, maxmissing)
    logger.warning("svl len: "+str(len(svl)))

    print("getting descriptor lists...")
    #csvl = vv.country_svl(svl)
    asvl = vv.asn_svl(svl)
    #ssvl = vv.subnet_svl(svl)
    psvl = vv.prefix_svl(svl)

    ccache = vv.init_ccache(None, ccachef, t, duration, mask, fmt, oddballs, maxmissing)

    idc = defaultdict(list)
    # {idA_idB:closeness}
    iic = dict()
    # {asnA_asnB: [closeness]}
    ddc = defaultdict(list)
    print("\n\ncalculating closeness for ASNs...")
    asns = [c for c in asvl if len(asvl[c]) > 1]
    for i in xrange(0, len(asns)-1):
        print(asns[i],end=", ")
        sys.stdout.flush()
        for a in asvl[asns[i]]:
            for j in xrange(i+1, len(asns)):
                for b in asvl[asns[j]]:
                    closeness = ccache[a][b]
                    ad = str(a.get_asn())
                    bd = str(b.get_asn())
                    aid = str(a.get_id())
                    bid = str(b.get_id())
                    dist = em.latlong_distance_km(a.get_coordinates(),
                            b.get_coordinates())
                    dist = distance(closeness, dist)
                    idc[aid+"_"+bd].append((closeness, dist))
                    idc[bid+"_"+ad].append((closeness, dist))
                    iic["_".join(sorted([aid, bid]))] = (closeness, dist)
                    ddc["_".join(sorted([ad, bd]))].append((closeness,
                        dist))
    ccache.dump()


    idac = sorted([(k, np.mean([q[0] for q in idc[k]]), np.mean([q[1] for q in \
            idc[k]])) for k in idc], key=lambda z: z[2], reverse=True)
    idac = [(z[0], z[1]) for z in idac]
    filename = plotsdir+"asn_idac"+fname+".csv"
    df.overwrite(filename, df.list2col(idac))

    ddac = sorted([(k, np.mean([q[0] for q in ddc[k]]), np.mean([q[1] for q in \
            ddc[k]])) for k in ddc], key=lambda z: z[2], reverse=True)
    ddac = [(z[0], z[1]) for z in ddac]
    filename = plotsdir+"asn_ddac"+fname+".csv"
    df.overwrite(filename, df.list2col(ddac))

    iic = sorted([(k, iic[k][0], iic[k][1]) for k in iic], reverse=True,
            key=lambda z: z[2])
    iic = [(z[0], z[1]) for z in iic]
    filename = plotsdir+"asn_iic"+fname+".csv"
    df.overwrite(filename, df.list2col(iic))

    # {idA_prefixB: [closeness]}
    idc = defaultdict(list)
    # {idA_idB:closeness}
    iic = dict()
    # {prefixA_prefixB: [closeness]}
    ddc = defaultdict(list)
    print("\n\ncalculating closeness for prefixes...")
    prefixes = [c for c in psvl if len(psvl[c]) > 1]
    for i in xrange(0, len(prefixes)-1):
        print(prefixes[i],end=", ")
        sys.stdout.flush()
        for a in psvl[prefixes[i]]:
            for j in xrange(i+1, len(prefixes)):
                for b in psvl[prefixes[j]]:
                    closeness = ccache[a][b]
                    ad = str(a.get_prefix())
                    bd = str(b.get_prefix())
                    aid = str(a.get_id())
                    bid = str(b.get_id())
                    dist = em.latlong_distance_km(a.get_coordinates(),
                            b.get_coordinates())
                    dist = distance(closeness, dist)
                    idc[aid+"_"+bd].append((closeness, dist))
                    idc[bid+"_"+ad].append((closeness, dist))
                    iic["_".join(sorted([aid, bid]))] = (closeness, dist)
                    ddc["_".join(sorted([ad, bd]))].append((closeness,
                        dist))
    ccache.dump()


    idac = sorted([(k, np.mean([q[0] for q in idc[k]]), np.mean([q[1] for q in \
            idc[k]])) for k in idc], key=lambda z: z[2], reverse=True)
    idac = [(z[0], z[1]) for z in idac]
    filename = plotsdir+"prefix_idac"+fname+".csv"
    df.overwrite(filename, df.list2col(idac))

    ddac = sorted([(k, np.mean([q[0] for q in ddc[k]]), np.mean([q[1] for q in \
            ddc[k]])) for k in ddc], key=lambda z: z[2], reverse=True)
    ddac = [(z[0], z[1]) for z in ddac]
    filename = plotsdir+"prefix_ddac"+fname+".csv"
    df.overwrite(filename, df.list2col(ddac))

    iic = sorted([(k, iic[k][0], iic[k][1]) for k in iic], reverse=True,
            key=lambda z: z[2])
    iic = [(z[0], z[1]) for z in iic]
    filename = plotsdir+"prefix_iic"+fname+".csv"
    df.overwrite(filename, df.list2col(iic))

    svd = dict()
    for sv in svl:
        svd[sv.get_id()] = sv

    return svd


def farthest_same_desc(t, duration=30000, mask=32, fmt=None,
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
    :returns: [country, ASN, subnet, prefix] pair dictionaries of closeness lists

    gets pairwise closeness of probes with different descriptors to find odd
    behavior (probes in difference descriptors with high closeness scores)

    NOTE: writes data to files for conveniece
    '''
    print("getting svl...")
    svl, fmt, _ = vv.get_svl(t, duration, mask, fmt, country_set, oddballs, maxmissing)
    logger.warning("svl len: "+str(len(svl)))

    print("getting descriptor lists...")
    csvl = vv.country_svl(svl)
    asvl = vv.asn_svl(svl)
    ssvl = vv.subnet_svl(svl)
    psvl = vv.prefix_svl(svl)

    ccache = vv.init_ccache(None, ccachef, t, duration, mask, fmt, oddballs, maxmissing)

    # {idA_countryB: [closeness]}
    idc = defaultdict(list)
    # {idA_idB:closeness}
    iic = dict()
    # {countryA_countryB: [closeness]}
    ddc = defaultdict(list)
    print("calculating closeness for countries...")
    countries = [c for c in csvl if len(csvl[c]) > 1]
    for i in xrange(0, len(countries)):
        print(countries[i],end=", ")
        sys.stdout.flush()
        for j in xrange(0, len(csvl[countries[i]])-1):
            a = csvl[countries[i]][j]
            for l in xrange(j+1, len(csvl[countries[i]])):
                b = csvl[countries[i]][l]
                closeness = ccache[a][b]
                ad = a.get_country()
                bd = b.get_country()
                aid = str(a.get_id())
                bid = str(b.get_id())
                dist = em.latlong_distance_km(a.get_coordinates(),
                        b.get_coordinates())
                dist = distance(closeness, dist)
                idc[aid+"_"+bd].append((closeness, dist))
                idc[bid+"_"+ad].append((closeness, dist))
                iic["_".join(sorted([aid, bid]))] = (closeness, dist)
                ddc["_".join(sorted([ad, bd]))].append((closeness,
                    dist))
    ccache.dump()


    idac = sorted([(k, np.mean([q[0] for q in idc[k]]), np.mean([q[1] for q in \
            idc[k]])) for k in idc], key=lambda z: z[2])
    idac = [(z[0], z[1]) for z in idac]
    filename = plotsdir+"country_sdidac"+fname+".csv"
    df.overwrite(filename, df.list2col(idac))

    ddac = sorted([(k, np.mean([q[0] for q in ddc[k]]), np.mean([q[1] for q in \
            ddc[k]])) for k in ddc], key=lambda z: z[2])
    ddac = [(z[0], z[1]) for z in ddac]
    filename = plotsdir+"country_sdac"+fname+".csv"
    df.overwrite(filename, df.list2col(ddac))

    iic = sorted([(k, iic[k][0], iic[k][1]) for k in iic],
            key=lambda z: z[2])
    iic = [(z[0], z[1]) for z in iic]
    filename = plotsdir+"country_sdiic"+fname+".csv"
    df.overwrite(filename, df.list2col(iic))


    # {idA_subnetB: [closeness]}
    idc = defaultdict(list)
    # {idA_idB:closeness}
    iic = dict()
    # {subnetA_subnetB: [closeness]}
    ddc = defaultdict(list)
    print("\n\ncalculating closeness for subnets...")
    subnets = [c for c in ssvl if len(ssvl[c]) > 1]
    for i in xrange(0, len(subnets)):
        print(subnets[i],end=", ")
        sys.stdout.flush()
        for j in xrange(0, len(ssvl[subnets[i]])-1):
            a = ssvl[subnets[i]][j]
            for l in xrange(j+1, len(ssvl[subnets[i]])):
                b = ssvl[subnets[i]][l]
                closeness = ccache[a][b]
                ad = str(a.get_subnet())
                bd = str(b.get_subnet())
                aid = str(a.get_id())
                bid = str(b.get_id())
                dist = em.latlong_distance_km(a.get_coordinates(),
                        b.get_coordinates())
                dist = distance(closeness, dist)
                idc[aid+"_"+bd].append((closeness, dist))
                idc[bid+"_"+ad].append((closeness, dist))
                iic["_".join(sorted([aid, bid]))] = (closeness, dist)
                ddc["_".join(sorted([ad, bd]))].append((closeness,
                    dist))
    ccache.dump()


    idac = sorted([(k, np.mean([q[0] for q in idc[k]]), np.mean([q[1] for q in \
            idc[k]])) for k in idc], key=lambda z: z[2])
    idac = [(z[0], z[1]) for z in idac]
    filename = plotsdir+"subnet_sdidac"+fname+".csv"
    df.overwrite(filename, df.list2col(idac))

    ddac = sorted([(k, np.mean([q[0] for q in ddc[k]]), np.mean([q[1] for q in \
            ddc[k]])) for k in ddc], key=lambda z: z[2])
    ddac = [(z[0], z[1]) for z in ddac]
    filename = plotsdir+"subnet_sdac"+fname+".csv"
    df.overwrite(filename, df.list2col(ddac))

    iic = sorted([(k, iic[k][0], iic[k][1]) for k in iic],
            key=lambda z: z[2])
    iic = [(z[0], z[1]) for z in iic]
    filename = plotsdir+"subnet_sdiic"+fname+".csv"
    df.overwrite(filename, df.list2col(iic))

    # {idA_asnB: [closeness]}
    idc = defaultdict(list)
    # {idA_idB:closeness}
    iic = dict()
    # {asnA_asnB: [closeness]}
    ddc = defaultdict(list)
    print("\n\ncalculating closeness for ASNs...")
    asns = [c for c in asvl if len(asvl[c]) > 1]
    for i in xrange(0, len(asns)):
        print(asns[i],end=", ")
        sys.stdout.flush()
        for j in xrange(0, len(asvl[asns[i]])-1):
            a = asvl[asns[i]][j]
            for l in xrange(j+1, len(asvl[asns[i]])):
                b = asvl[asns[i]][l]
                closeness = ccache[a][b]
                ad = str(a.get_asn())
                bd = str(b.get_asn())
                aid = str(a.get_id())
                bid = str(b.get_id())
                dist = em.latlong_distance_km(a.get_coordinates(),
                        b.get_coordinates())
                dist = distance(closeness, dist)
                idc[aid+"_"+bd].append((closeness, dist))
                idc[bid+"_"+ad].append((closeness, dist))
                iic["_".join(sorted([aid, bid]))] = (closeness, dist)
                ddc["_".join(sorted([ad, bd]))].append((closeness,
                    dist))
    ccache.dump()


    idac = sorted([(k, np.mean([q[0] for q in idc[k]]), np.mean([q[1] for q in \
            idc[k]])) for k in idc], key=lambda z: z[2])
    idac = [(z[0], z[1]) for z in idac]
    filename = plotsdir+"asn_sdidac"+fname+".csv"
    df.overwrite(filename, df.list2col(idac))

    ddac = sorted([(k, np.mean([q[0] for q in ddc[k]]), np.mean([q[1] for q in \
            ddc[k]])) for k in ddc], key=lambda z: z[2])
    ddac = [(z[0], z[1]) for z in ddac]
    filename = plotsdir+"asn_sdac"+fname+".csv"
    df.overwrite(filename, df.list2col(ddac))

    iic = sorted([(k, iic[k][0], iic[k][1]) for k in iic],
            key=lambda z: z[2])
    iic = [(z[0], z[1]) for z in iic]
    filename = plotsdir+"asn_sdiic"+fname+".csv"
    df.overwrite(filename, df.list2col(iic))

    # {idA_prefixB: [closeness]}
    idc = defaultdict(list)
    # {idA_idB:closeness}
    iic = dict()
    # {prefixA_prefixB: [closeness]}
    ddc = defaultdict(list)
    print("\n\ncalculating closeness for prefixes...")
    prefixes = [c for c in psvl if len(psvl[c]) > 1]
    for i in xrange(0, len(countries)):
        print(prefixes[i],end=", ")
        sys.stdout.flush()
        for j in xrange(0, len(psvl[prefixes[i]])-1):
            a = psvl[prefixes[i]][j]
            for l in xrange(j+1, len(psvl[prefixes[i]])):
                b = psvl[prefixes[i]][l]
                closeness = ccache[a][b]
                ad = a.get_prefix()
                bd = b.get_prefix()
                aid = str(a.get_id())
                bid = str(b.get_id())
                dist = em.latlong_distance_km(a.get_coordinates(),
                        b.get_coordinates())
                dist = distance(closeness, dist)
                idc[aid+"_"+bd].append((closeness, dist))
                idc[bid+"_"+ad].append((closeness, dist))
                iic["_".join(sorted([aid, bid]))] = (closeness, dist)
                ddc["_".join(sorted([ad, bd]))].append((closeness,
                    dist))
    ccache.dump()


    idac = sorted([(k, np.mean([q[0] for q in idc[k]]), np.mean([q[1] for q in \
            idc[k]])) for k in idc], key=lambda z: z[2])
    idac = [(z[0], z[1]) for z in idac]
    filename = plotsdir+"prefix_sdidac"+fname+".csv"
    df.overwrite(filename, df.list2col(idac))

    ddac = sorted([(k, np.mean([q[0] for q in ddc[k]]), np.mean([q[1] for q in \
            ddc[k]])) for k in ddc], key=lambda z: z[2])
    ddac = [(z[0], z[1]) for z in ddac]
    filename = plotsdir+"prefix_sdac"+fname+".csv"
    df.overwrite(filename, df.list2col(ddac))

    iic = sorted([(k, iic[k][0], iic[k][1]) for k in iic],
            key=lambda z: z[2])
    iic = [(z[0], z[1]) for z in iic]
    filename = plotsdir+"prefix_sdiic"+fname+".csv"
    df.overwrite(filename, df.list2col(iic))

    svd = dict()
    for sv in svl:
        svd[sv.get_id()] = sv

    return svd


# get probes with different local DNS resolvers that have most other features
# (IP/subnet, country, ASN, owner) in common (such that local DNS is the only
# distinguishing factor) NOTE: will need to do test queries from these probes to
# our own server to distinguish resolvers that are actually distinct from those
# that are actually just forwarding requests to the same upstream resolvers
def nearby_probes_diff_ldns(svl, rmask=16):

    print("reducing svl to only probes with public ldns")
    svl = [sv for sv in svl if ipp.is_public(sv.get_ldns())]
    print("getting asn descriptor list...")
    asvl = vv.asn_svl(svl)

    nearbies = list()
    strnearbies = list()

    for asn in asvl:
        if asn is None:
            continue
        ssvl = vv.subnet_svl(asvl[asn], 16)
        for subnet in ssvl:
            if subnet is None:
                continue
            if len(ssvl[subnet]) > 1:
                csvl = vv.country_svl(ssvl[subnet])
                for country in csvl:
                    if country is None:
                        continue
                    if len(csvl[country]) > 1:
                        osvl = vv.owner_svl(csvl[country])
                        for owner in osvl:
                            if owner is None:
                                continue
                            if len(osvl[owner]) > 1:
                                resolvers = [z.get_ldns() for z in osvl[owner]]
                                # collapse redundant resolvers in the same
                                # subnet (e.g., 8.8.8.8 and 8.8.4.4 -> 8.8.x.x)
                                fmtmask = ipp.make_v4_prefix_mask(rmask)
                                r2 = defaultdict(list)
                                for r in resolvers:
                                    r2[r & fmtmask].append(r)
                                # keep resolvers with at least 1 probe //2 probes
                                keep = list()
                                k2 = list()
                                for z in r2:
                                    if len(r2[z]) > 1:
                                        # // keep += r2[z]
                                        k2.append(z)
                                    keep += r2[z]
                                #if len(k2) > 0:
                                if len(r2) > 1:
                                    print("has different resolvers!!!")
                                if len(keep) > 1: # if there's stuff to compare
                                    print("found some!!")
                                    print("asn: "+str(asn))
                                    svl = [sv for sv in osvl[owner] if sv.get_ldns() \
                                            in keep]
                                    nearbies.append(svl)
                                    strnearbies.append(["|~"+ipp.int2ip(z.get_ldns())+\
                                            "__"+ipp.int2ip(z.get_ip())+"~|" \
                                            for z in svl])
    df.overwrite(plotsdir+"nearbies.csv", df.list2col(strnearbies))
    logger.warning("nearbies: "+str(len(nearbies)))
    return nearbies



