from warriorpy.shorthand import diriofile as df
from warriorpy.net_tools import ipparsing as ipp
import numpy as np
from pymongo import MongoClient
from IPy import IP
from collections import defaultdict

"""
class for answer vectors

purpose:
    to create an n-dimensional "space" for determining the "closeness" of clients

assumptions:

    FIRST ANSWER: I only work with the first answer from each answer set (i.e. the answer that is
    actually being used

    TIME BLOCKS: It looks like the the probes complete a cycle through the set of domains for
    meas. 30002 every ~30000 seconds, which is between 8-9 hours. I therefor assume that any
    incidents within the same 30000 second block are essentially 'simultaneous'.

"""

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

domains = df.getlines(basedir+'state/sites.csv')

# database setup
mclient = MongoClient()
db = mclient.veracity
data = db.m30002_may17
probe_cache = db.probe_data

##################################################################
#                           CODE
##################################################################


# Probes happen to query certain domains multiple times back-to-back for
# some reason; in these cases, I will take only their first query
def get_window(start_time, duration, domain_set, country=None):
    cmd = [
            {"$match": {
                "time": {
                        "$gte": start_time,
                        "$lt": start_time+duration
                        },
                "domain": {"$in": domain_set},
                "ipv4.answer_ip_list": {"$exists": True}
                }
            },
            {"$sort": {"ind": 1}},
            {"$group": {
                "_id": {
                    "probe_ip": "$probe_ip",
                    "probe_id": "$probe_id",
                    "domain": "$domain"
                    },
                "answers": {"$first": "$ipv4.answer_ip_list"},
                "ind": {"$first": "$ind"}
                },
            }
        ]
    if country is not None:
        cmd[0]["$match"]["country"] = country
    res = data.aggregate(cmd, allowDiskUse=True)
    return res


# convert from set of raw documents from single time block to dictionary of
# domains vs answer-sets. If there is more than one answer for a domain, I take
# the first answer
def window_to_dict(window):
    d = defaultdict(lambda: defaultdict(list))
    for doc in window:
        if doc['_id']['probe_ip'] not in d:
            d[doc['_id']['probe_ip']]['ind'].append(doc['ind'])
        if len(doc['answers']) > 0:
            d[doc['_id']['probe_ip']][doc['_id']['domain']] += doc['answers']
    return d


def dict_to_ipstruct(d, fmt=None, oddballs=False, weight=None, mask=32):
    '''
    NOTE: each time an ip is seen, its respective domain's weight is added to
    its entry in the vector; in other words, the more often an IP is seen, the
    more important it is
    '''
    fmtmask = ipp.make_v4_prefix_mask(mask)
    vector = defaultdict(int)
    if fmt is None:
        fmt = domains
    elif type(fmt) is int:
        fmt = domains[:fmt]
    for domain in fmt:
        if domain in d:
            for ip in d[domain]:
                ipm = ip & fmtmask
                ipstr = ipp.int2ip(ipm)
                if (ipm != 0 and IP(ipstr+"/32").iptype() == "PUBLIC") or oddballs:
                    vector[ipm] += weight[domain]/float(len(d[domain]))
    return vector


# generate vector from dict of {domain: answer} entries
def dict_to_vector(d, fmt=None):
    vector = list()
    if fmt is None:
        fmt = domains
    elif type(fmt) is int:
        fmt = domains[:fmt]
    for domain in fmt:
        if domain in d:
            vector.append(d[domain][0])
        else:
            vector.append(0)
    return vector


def weight_vector(vec, wvec):
    pass


def distance_metric(a, b):
    '''
    (a n b) / (b u a), where each entry is weighted using a tuple (ip, weight)

    NOTE: since every value contributes to sums twice (once from a and once from
    b), the weights are effectively half-weights. This is to account for the
    fact that different domains - or different numbers of domains - may contribute
    the same IP for a and b.
    '''
    totalval = sum([a[z] for z in a]+[b[z] for z in b])
    overlap = set(a.keys()).intersection(set(b.keys()))
    aweight = [a[z] for z in a if z in overlap]
    bweight = [b[z] for z in b if z in overlap]
    return float(sum(aweight+bweight))/float(totalval)


def get_vectors(t, duration=30000, mask=32, fmt=None, country=None):
    print "getting window"
    if fmt is None:
        window = get_window(t, duration, domains, country)
    elif type(fmt) is int:
        window = get_window(t, duration, domains[:fmt], country)
    else:
        window = get_window(t, duration, fmt, country)
    print "converting window to dict"
    dd = window_to_dict(window)
    X = list()
    indl = list()
    # list of indices
    print "creating array"
    fmtmask = ipp.make_v4_prefix_mask(mask)
    for probe in dd:
        vec = np.array(dict_to_vector(dd[probe])) & fmtmask
        X.append(vec)
        indl.append(dd[probe]['ind'])
    return np.array(X)




def get_answer_space_dict(t, duration=30000, mask=32, fmt=None, country=None):
    X = get_vectors(t, duration, mask, fmt, country)
    anssets = defaultdict(set)
    if fmt is None:
        fmt = domains
    elif type(fmt) is int:
        fmt = domains[:fmt]
    for vec in X:
        for ind, dom in enumerate(fmt):
            anssets[dom].add(ipp.int2ip(vec[ind]))
    return anssets


def sort_sites(t, duration=30000, mask=32, fmt=None, country=None):
    anssets = get_answer_space_dict(t, duration, mask, fmt, country)
    spacelist = list()
    for dom in anssets:
        spacelist.append(len(anssets[dom]))

    outlist = [(z,len(anssets[z])) for z in anssets.keys()]

    sl = sorted(outlist, key=lambda z: z[1])
    sl = [z[0] for z in sl]

    slcol = df.list2col(sl)
    df.overwrite(basedir+'state/sites.csv', slcol)
    df.overwrite(basedir+'state/site_sizes.csv', df.list2col(outlist))
    for val in sl:
        print val
    return outlist


def get_weighting(t, duration=30000, mask=32, fmt=None, country=None):
    anssets = get_answer_space_dict(t, duration, mask, fmt, country)
    spacelist = list()
    for dom in anssets:
        spacelist.append(len(anssets[dom]))

    sl = [(z,len(anssets[z])) for z in anssets.keys()]

    wd = dict()
    for item in sl:
        wd[item[0]] = float(item[1])
    return wd

