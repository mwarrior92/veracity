from warriorpy.shorthand import diriofile as df
from warriorpy.net_tools import ipparsing as ipp
from warriorpy.net_tools import dns_tools as dt
import numpy as np
from pymongo import MongoClient
from IPy import IP
from collections import defaultdict
import math

"""
class for answer vectors

purpose:
    to create an n-dimensional "space" for determining the "closeness" of clients

assumptions:

    ENTIRE ANSWER SETS: I user the entire answer set from each query (as opposed
    to only the first answer, for example)

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
pdata = db.probe_data

##################################################################
#                           CODE
##################################################################


# Probes happen to query certain domains multiple times back-to-back for
# some reason; in these cases, I will take only their first query
def get_window(start_time, duration, domain_set=None, country_set=None):
    cmd = [
            {"$match": {
                "time": {
                        "$gte": start_time,
                        "$lt": start_time+duration
                        },
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
    if country_set is not None:
        cmd[0]["$match"]["country"] = {"$in": country_set}
    if domain_set is not None:
        cmd[0]["$match"]["domain"] = {"$in": domain_set}
    res = data.aggregate(cmd, allowDiskUse=True)
    return res


# convert from set of raw documents from single time block to dictionary of
# domains vs answer-sets. If there is more than one answer for a domain
def window_to_dicts(window):
    d = defaultdict(lambda: defaultdict(list))
    for doc in list(window):
        if len(doc['answers']) > 0:
            d[doc['_id']['probe_ip']]['inds_'+doc['_id']['domain']].append(doc['ind'])
            d[doc['_id']['probe_ip']][doc['_id']['domain']] += doc['answers']
    return d


def transform_fmt(fmt, doms=None):
    if doms is None:
        doms = domains

    if fmt is None:
        fmt = doms
    elif type(fmt) is int:
        if fmt > 0:
            fmt = doms[:fmt]
        elif fmt < 0:
            fmt = doms[fmt:]
    return fmt


def dicts_to_svl(dd, fmt=None, mask=32, oddballs=False):
    '''
    convert dicts to list of smartvecs
    '''
    fmt = transform_fmt(fmt)
    svl = list()
    for probe_ip in dd:
        if probe_ip == "":
            continue
        tmp = smartvec(dd[probe_ip], probe_ip, fmt,
            mask, oddballs)
        if len(set(tmp.vec).symmetric_difference(set(fmt))) < min([3, len(fmt)]):
            # if vec isn't missing anything, keep it
            svl.append(tmp)
    return svl


def country_svl(svl):
    csvl = defaultdict(list) # {country: svl}
    for sv in svl:
        csvl[sv.get_country()].append(sv)
    return csvl


def asn_svl(svl):
    asvl = defaultdict(list) # {asn: svl}
    for sv in svl:
        asvl[sv.get_asn()].append(sv)
    return asvl


def subnet_svl(svl, mask=24):
    ssvl = defaultdict(list) # {subnet: svl}
    for sv in svl:
        ssvl[sv.get_subnet(mask)].append(sv)
    return ssvl


def owner_svl(svl):
    osvl = defaultdict(list) # {owner: svl}
    for sv in svl:
        osvl[sv.get_owner()].append(sv)
    return osvl


def prefix_svl(svl):
    psvl = defaultdict(list) # {prefix_v4: svl}
    for sv in svl:
        psvl[sv.get_prefix()].append(sv)
    return psvl


class smartvec:
    def __init__(self, d, probe_ip, fmt=None, mask=32, oddballs=False):
        '''
        d is a dictionary with indices stored in key 'inds_'+domain and answer lists stored
        in keys of their respective domain (e.g. 'google.com.')
        '''
        self.vec = defaultdict(lambda: defaultdict(float)) # {domain: {ans: cum. weight}}
        fmt = transform_fmt(fmt)
        fmtmask = ipp.make_v4_prefix_mask(mask)
        for dom in fmt:
            if dom in d:
                query_count = float(len(d['inds_'+dom]))
                answer_count = float(len(d[dom]))
                for ip in d[dom]:
                    ipm = ip & fmtmask
                    ipstr = ipp.int2ip(ipm)
                    if IP(ipstr+"/32").iptype() == "PUBLIC" or oddballs:
                        self.vec[dom][ipm] += query_count / answer_count
        self.ip = probe_ip


    def __repr__(self):
        return "smartvec(ip="+ipp.int2ip(self.ip)+")"


    def __str__(self):
        s = ""
        for dom in self.vec:
            s += dom+": "+", ".join([ipp.int2ip(ip) for ip in self.vec[dom]])+"\n\n"
        return s


    def get_ip(self):
        return ipp.int2ip(self.ip)


    def get_answers(self, dom):
        if dom in self.vec:
            return [ipp.int2ip(ip) for ip in self.vec[dom]]
        return None


    def get_subnet(self, mask):
        fmtmask = ipp.make_v4_prefix_mask(mask)
        return ipp.int2ip(self.ip & fmtmask)


    def get_probe_info(self):
        tmp = list(pdata.find({'probe_ip': self.ip}).limit(1))
        if len(tmp) > 0:
            return tmp[0]
        else:
            return None


    def get_country(self):
        tmp = self.get_probe_info()
        if tmp is not None:
            return tmp['country']
        else:
            return None


    def get_asn(self):
        tmp = self.get_probe_info()
        if tmp is not None:
            return tmp['asn_v4']
        else:
            return None


    def get_id(self):
        tmp = self.get_probe_info()
        if tmp is not None:
            return tmp['probe_id']
        else:
            return None


    def get_prefix(self):
        tmp = self.get_probe_info()
        if tmp is not None:
            return tmp['prefix_v4']
        else:
            return None


    def get_owner(self):
        return dt.get_owner_name(ipp.int2ip(self.ip))


    def sees_private_self(self):
        tmp = self.get_probe_info()
        if tmp is not None:
            ip = tmp['ipv4']['perceived_ip']
            if IP(ipp.int2ip(ip)+"/32").iptype() == "PUBLIC":
                return False
        return True


    def sees_private_ldns(self):
        tmp = self.get_probe_info()
        if tmp is not None:
            ip = tmp['ipv4']['perceived_ldns']
            if IP(ipp.int2ip(ip)+"/32").iptype() == "PUBLIC":
                return False
        return True


    def __len__(self):
        return len(self.vec)


    def __iter__(self):
        for item in dict.__iter__(self.vec):
            yield item


    def __contains__(self, item):
        return item in self.vec


    def __getitem__(self, item):
        return self.vec[item]


def closeness(a, b):
        '''
        (a n b) / (b u a), where each entry is weighted using a tuple (ip, weight)

        NOTE: since every value contributes to sums twice (once from a and once from
        b), the weights are effectively half-weights. This is to account for the
        fact that different domains - or different numbers of domains - may contribute
        the same IP for a and b.

        NOTE: each domain is normalized, so domains with a lot of IPs per query
        response won't skew the results
        '''
        n = 0.0 # numerator
        d = float(len(a)) # denominator
        for dom in [j for j in a if j in b]:
            domtotal = sum([a[dom][z] for z in a[dom]]+[b[dom][z] for z in b[dom]])
            overlap = set(a[dom]).intersection(set(b[dom]))
            aweight = [a[dom][z] for z in a[dom] if z in overlap]
            bweight = [b[dom][z] for z in b[dom] if z in overlap]
            n += sum(aweight+bweight)/domtotal
        return n/d


def get_dist_list(svl):
    # NOTE: omits self comparisons and associative (redundant) comparisons;
    # can't use get_dist_list output with linkage
    dist_list = list()
    checked = set()
    for ii, i in enumerate(svl):
        checked.add(ii)
        for jj, j in enumerate(svl):
            if jj not in checked:
                dist_list.append(1.0-closeness(i, j))
    return dist_list


def avg_dist(dist_list):
    return np.mean(dist_list)


def max_dist(dist_list):
    return max(dist_list)


def min_dist(dist_list):
    return min(dist_list)


def median_dist(dist_list):
    return np.median(dist_list)


def get_answer_space_dict(dd, fmt=None):
        anssets = defaultdict(set)
        fmt = transform_fmt(fmt)
        for probe in dd:
            for dom in fmt:
                if dom in dd[probe]:
                    anssets[dom] |= set(dd[probe][dom])
        return anssets


def get_probe_space(dd, fmt=None):
    ps = defaultdict(int)
    fmt = transform_fmt(fmt)
    for probe in dd:
        for dom in fmt:
            if dom in dd[probe]:
                ps[dom] += 1
    print ps
    return ps


def get_weighting(anssets):
        '''
        use output from get_answer_space_dict() as input for this
        '''
        sl = [(z,len(anssets[z])) for z in anssets]

        wd = dict()
        for item in sl:
            wd[item[0]] = 1.0 /float(item[1]) #math.log(float(item[1])+2, 2)
        return wd


def sort_sites(anssets):
        '''
        use output from get_answer_space_dict() as input for this
        '''
        spacelist = [(z,len(anssets[z])) for z in anssets]

        sl = sorted(spacelist, key=lambda z: z[1])
        for val in sl:
            print val
        sl = [z[0] for z in sl]

        return sl
