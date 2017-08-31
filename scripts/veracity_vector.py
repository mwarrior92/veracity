from warriorpy.shorthand import diriofile as df
from warriorpy.shorthand import vinegar as vngr
from warriorpy.net_tools import ipparsing as ipp
from warriorpy.net_tools import dns_tools as dt
import numpy as np
from pymongo import MongoClient
from IPy import IP
from collections import defaultdict
import math
import cPickle as pickle
import pprint
import numbers
import time
import copy

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
statedir = df.rightdir(basedir+"state/")
rawdirlist = df.getlines(basedir+'state/datapaths.list')
datafiles = df.listfiles(basedir+rawdirlist[0], fullpath=True)
vngr.set_cache_dir(df.rightdir(statedir+"pickles"))

# database setup
mclient = MongoClient()
db = mclient.veracity
data = db.m30002_may17_full
pdata = db.probe_data

#printing
pp = pprint.PrettyPrinter(indent=4)

# vector processor
vproc = df.picklein(statedir+'nmf_rank20_pkr01.pickle')

##################################################################
#                           CODE
##################################################################


# Probes happen to query certain domains multiple times back-to-back for
# some reason; in these cases, I will take only their first query

def get_window(start_time, duration, domain_set=None, country_set=None):
    '''
    :param start_time: int indicating the earliest query the window should include
    :param duration: int indication the span of time covered by the window,
        in seconds
    :param domain_set: the list of domains the window should include queries for.
        If None, then all domains will be included
    :param country_set: the set of countries the window should include queries from.
        If None, then all countries will be inluded
    :return: dictionary, where each entry is a mongodb doc containing a DNS query
    '''

    match = {"$match": {
                "time": {
                        "$gte": start_time,
                        "$lt": start_time+duration
                        },
                "ipv4.answer_ip_list": {"$exists": True}
                }
            }
    group1 = {"$group": {
                "_id": {
                    "probe_id": "$probe_id",
                    "domain": "$domain"
                    },
                "answers": { "$push": "$ipv4.answer_ip_list" },
                "ind": { "$addToSet": "$ind"},
                "probe_ip": {"$addToSet": "$probe_ip"},
                "ldns": {"$addToSet": "$ipv4.perceived_ldns"}
                }
             }
    unwind = {"$unwind": "$answers"}
    group2 = {"$group": {
                "_id": "$_id",
                "answers": {"$push": "$answers"},
                "ind": {"$first": "$ind"},
                "probe_ip": {"$first": "$probe_ip"},
                "ldns": {"$first": "$ldns"}
                }
             }
    group3 = {"$group": {
                "_id": "$_id.probe_id",
                "domains": {"$push": "$_id.domain"},
                "answers": {"$push": "$answers"},
                "ind": {"$push": "$ind"},
                "probe_ip": {"$addToSet": "$probe_ip"},
                "ldns": {"$addToSet": "$ldns"}
                }
             }

    cmd = [match, group1, unwind, unwind, group2, group3]
    if type(country_set) is list:
        cmd[0]["$match"]["country"] = {"$in": country_set}
    if type(domain_set) is list:
        cmd[0]["$match"]["domain"] = {"$in": domain_set}
    res = data.aggregate(cmd, allowDiskUse=True)
    return res


def transform_fmt(fmt, doms=None):
    '''
    :param fmt: None, int, or list of domains. If None, will assign list of
        domains read from file; if int n > 0, will assign first n domains from list
        read from file; if int n < 0, will assign last n domains read from file; if
        list of domains, will use directly (ignoring list of domains read from file)
    :param doms: list of all domains; if None, reads from file
    :return: list of domains
    '''

    if fmt is None:
        fmt = doms
    elif type(fmt) is int:
        if fmt > 0:
            fmt = doms[:fmt]
        elif fmt < 0:
            fmt = doms[fmt:]
    return fmt


def dicts_to_svl(dd, mask=32, oddballs=False):
    '''
    :param dd:  output from window_to_dicts()
    :param fmt: see transform fmt
    :param mask: int, prefix mask to use over domain IPs
    :param oddballs: if True, will include non-public IPs (10.x.x.x, etc); if
        False, will only include public IPs
    :return: list of smartvecs

    convert dicts to list of smartvecs
    '''
    svl = list()
    anssets = defaultdict(set)
    for probe in dd:
        if len(probe['probe_ip']) == 0:
            continue
        if len(probe['probe_ip'][0]) == 0:
            continue
        noip = True
        for i in xrange(0, len(probe['probe_ip'])):
            if noip:
                for j in xrange(0, len(probe['probe_ip'][i])):
                    if isinstance(probe['probe_ip'][i][j], numbers.Number):
                        noip = False
                        break
        if noip:
            continue
        svl.append(smartvec(probe, mask, oddballs))
        for i, dom in enumerate(probe['domains']):
            anssets[dom] |= set(probe['answers'][i])
    doms = sorted(anssets.keys(), key=lambda z: len(anssets[z]))
    return svl, doms, anssets


def country_svl(svl):
    '''
    :param svl: list of smartvecs
    :return: dict: {country: list of smartvecs}
    '''
    csvl = defaultdict(list) # {country: svl}
    for sv in svl:
        csvl[sv.get_country()].append(sv)
    return csvl


def asn_svl(svl):
    '''
    :param svl: list of smartvecs
    :return: dict: {asn: list of smartvecs}
    '''
    asvl = defaultdict(list) # {asn: svl}
    for sv in svl:
        asvl[sv.get_asn()].append(sv)
    return asvl


def subnet_svl(svl, mask=24):
    '''
    :param svl: list of smartvecs
    :param mask: subnets used will be of format /mask
    :return: dict: {subnet: list of smartvecs}, where 'subnet' refers to a /mask
        prefix
    '''
    ssvl = defaultdict(list) # {subnet: svl}
    for sv in svl:
        ssvl[sv.get_subnet(mask)].append(sv)
    return ssvl


def owner_svl(svl):
    '''
    :param svl: list of smartvecs
    :return: dict: {owner: list of smartvecs}
    '''
    osvl = defaultdict(list) # {owner: svl}
    for sv in svl:
        try:
            osvl[sv.get_owner()].append(sv)
        except:
            # I now generic except is bad code, but sending queries over the
            # Internet bugs out a lot so I need a catch all
            pass
    return osvl


def prefix_svl(svl):
    '''
    :param svl: list of smartvecs
    :return: dict: {prefix: list of smartvecs}, where prefix is the size of the
        block known to registries
    '''
    psvl = defaultdict(list) # {prefix_v4: svl}
    for sv in svl:
        psvl[sv.get_prefix()].append(sv)
    return psvl


def ldns_svl(svl, rmask=32, allowPrivate=True):
    '''
    :param svl: list of smartvecs
    :return: dict: {country: list of smartvecs}
    '''
    fmtmask = ipp.make_v4_prefix_mask(rmask)
    csvl = defaultdict(list) # {country: svl}
    for sv in svl:
        ip = sv.get_ldns()
        if isinstance(ip, numbers.Number):
            if ipp.is_public(ip) or allowPrivate:
                csvl[ip & fmtmask].append(sv)
    return csvl


def ddfloat():
    return defaultdict(float)


class smartvec:
    def __init__(self, d, mask=32, oddballs=False):
        '''
        :param d: a single subdict d from the dict dd output from
            window_to_dicts()
        :param mask: prefix mask t use on domain IPs
        :param oddballs: if True, will include non-public IPs (10.x.x.x, etc); if
            False, will only include public IPs

        container for a set of query results from a single probe

        NOTE: mask is only applied to domain IPs (not the probe's IP)
        NOTE: after mask is applied, original IPs are not preserved; for
            comparisons, separate sv's should be constructed
        '''
        self.vec = defaultdict(ddfloat) # {domain: {ans: cum. weight}}
        fmtmask = ipp.make_v4_prefix_mask(mask)
        self.query_count = defaultdict(float)
        self.answer_count = defaultdict(float)
        self.inds = defaultdict(list)
        for i, dom in enumerate(d['domains']):
            query_count = float(len(d['ind'][i]))
            answer_count = float(len(d['answers'][i]))
            qaratio = query_count / answer_count
            self.query_count[dom] = query_count
            self.answer_count[dom] = answer_count
            self.inds[dom] = d['ind'][i]
            for ip in d['answers'][i]:
                ipm = ip & fmtmask
                ipstr = ipp.int2ip(ipm)
                if IP(ipstr+"/32").iptype() == "PUBLIC" or oddballs:
                    self.vec[dom][ipm] += qaratio
            if len(self.vec[dom]) < 1:
                del self.vec[dom]
            else:
                self.vec[dom] = dict(self.vec[dom])
        self.vec = dict(self.vec)
        self.mask = mask
        self.ip = set()
        for ipset in d['probe_ip']:
            for ip in ipset:
                if not isinstance(ip, numbers.Number):
                    continue
                self.ip.add(ip)
        self.ldns = set()
        for ipset in d['ldns']:
            for ip in ipset:
                if not isinstance(ip, numbers.Number):
                    continue
                self.ldns.add(ip)
        self.ldns = list(self.ldns)
        self.id = d['_id']
        self.probe_info = None
        self.owner = None


    def absorb(self, sv):
        # merge two smartvectors
        for dom in sv:
            if dom not in self.vec:
                self.vec[dom] = dict()
            if dom not in self.inds:
                self.inds[dom] = list()
            self.inds[dom] += sv.inds[dom]
            if dom in self.vec:
                R1 = self.query_count[dom] / self.answer_count[dom]
            else:
                R1 = 1
            R2 = sv.query_count[dom] / sv.answer_count[dom]
            if dom not in self.query_count:
                self.query_count[dom] = 0
            Q1 = self.query_count[dom]
            if dom not in sv.query_count:
                sv.query_count[dom] = 0
            Q2 = sv.query_count[dom]
            if dom not in self.answer_count:
                self.answer_count[dom] = 0
            if dom not in sv.answer_count:
                sv.answer_count[dom] = 0
            A = self.answer_count[dom]+sv.answer_count[dom]
            self.answer_count[dom] = A
            self.query_count[dom] = Q1 + Q2
            for ip in set(sv.vec[dom].keys()+self.vec[dom].keys()):
                if ip not in self.vec[dom]:
                    self.vec[dom][ip] = 0
                self.vec[dom][ip] *= (R1*Q1/A)
                if ip not in sv.vec[dom]:
                    sv.vec[dom][ip] = 0
                self.vec[dom][ip] += (sv.vec[dom][ip]*(R2*Q2/A))
        self.ip = self.ip.union(sv.ip)
        self.ldns = sorted(list(set(self.ldns+sv.ldns)))


    def __repr__(self):
        return "smartvec(ip="+ipp.int2ip(self.get_ip())+")"


    def __str__(self):
        s = ""
        for dom in self.vec:
            s += dom+": "+", ".join([ipp.int2ip(ip) for ip in self.vec[dom]])+"\n\n"
        return s


    def get_ip(self):
        return list(self.ip)[0]


    def get_ip_str(self):
        return ipp.int2ip(self.get_ip())


    def get_answers(self, dom):
        if dom in self.vec:
            return [ipp.int2ip(ip) for ip in self.vec[dom]]
        return None


    def get_subnet(self, mask=24):
        fmtmask = ipp.make_v4_prefix_mask(mask)
        return ipp.int2ip(self.get_ip() & fmtmask)


    def get_probe_info(self):
        if self.probe_info is None:
            tmp = list(pdata.find({'probe_id': self.get_id()}).limit(1))
            if len(tmp) > 0:
                self.probe_info = tmp[0]
                return tmp[0]
            else:
                return None
        else:
            return self.probe_info


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


    def get_ldns(self):
        for ip in self.ldns:
            if ipp.is_public(ip):
                return ip
        return self.ldns[0]


    def get_id(self):
        return self.id


    def get_prefix(self):
        tmp = self.get_probe_info()
        if tmp is not None:
            return tmp['prefix_v4']
        else:
            return None


    def get_owner(self):
        if hasattr(self, "owner"):
            if self.owner is not None:
                return self.owner
        self.owner = dt.get_owner_name(ipp.int2ip(self.get_ip()))

        return self.owner

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


    def get_real_ldns(self):
        if sees_private_ldns:
            pass
        else:
            pass

    def get_coordinates(self):
        tmp = self.get_probe_info()
        if tmp is not None:
            try:
                return tmp['geometry']['coordinates']
            except KeyError:
                return None
        else:
            return None


    def uses_ecs(self):
        pass


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
        :param a: smartvec
        :param b: smartvec
        :return: closeness between a and b

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
            domtotal = sum([a[dom][z] for z in a[dom]] + \
                           [b[dom][z] for z in b[dom]])
            overlap = set(a[dom]).intersection(set(b[dom]))
            aweight = [a[dom][z] for z in a[dom] if z in overlap]
            bweight = [b[dom][z] for z in b[dom] if z in overlap]
            n += sum(aweight+bweight)/domtotal
        return n/d


def get_cl(svl, cc=None):
    '''
    :param svl: list of smartvecs
    :return: list of pairwise closeness measures between input smartvecs

    NOTE: omits self comparisons and associative (redundant) comparisons;
    '''
    cc = init_ccache(cc)
    cl = list()
    for i in xrange(0, len(svl)-1):
        for j in xrange(i+1, len(svl)):
            cl.append(cc[svl[i]][svl[j]])
    return cl


def get_probe_space(dd, fmt=None):
    '''
    :param dd: output from window_to_dicts()
    :param fmt: see transform fmt
    :return: dict {domain: number of probes to query domain}
    '''
    ps = defaultdict(int)
    fmt = transform_fmt(fmt)
    for probe in dd:
        for dom in fmt:
            if dom in dd[probe]:
                ps[dom] += 1
    return ps


def get_weighting(anssets):
    '''
    :param anssets: output from get_answer_space_dict()
    :return: dict {domain: weight}
    '''
    sl = [(z,len(anssets[z])) for z in anssets]

    wd = dict()
    for item in sl:
        wd[item[0]] = 1.0 /float(item[1]) #math.log(float(item[1])+2, 2)
    return wd


@vngr.cache_me_outside
def get_svl(start_time, duration=30000, mask=32, fmt=None, country_set=None,
        oddballs=True, maxmissing=0, return_ccache=True,
        ccachef=df.rightdir(statedir+"pickles/")+"ccache.pickle",
        mindomsize=2):
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
    :return: list of smartvecs
    '''
    # if the duration is too long, mongo times out the cursor before we can
    # finish processing, so we should break it into smaller jobs with
    # get_big_svl()
    if duration > 24*60*60:
        logger.warning("using get_big_svl(), which assumes that 'duration' is a multple of 24 hours")
        return get_big_svl(start_time, duration, mask, fmt, country_set,
            oddballs, maxmissing, return_ccache, ccachef, mindomsize)
    logger.info("->window...")
    window = get_window(start_time, duration, country_set=country_set, domain_set=fmt)

    logger.info("->svl...")
    svl, doms, anssets = dicts_to_svl(window, mask, oddballs)

    logger.debug(str(doms))

    fmt = transform_fmt(fmt, doms)
    # remove any domains that only have 1 IP (since all nodes will see the
    # same thing)
    for dom in fmt:
        if len(anssets[dom]) < mindomsize \
          or len([sv for sv in svl if dom not in sv]) > 0.5*len(svl):
            del anssets[dom]
    fmt = sorted(list(set(anssets.keys()).intersection(set(fmt))))
    svl = reduce_svl(svl, fmt, maxmissing)

    for dom in fmt:
        logger.debug("-----------"+dom+"-----------")
        tmp = sorted(anssets[dom])
        for val in tmp:
            if type(val) is int:
                logger.debug(ipp.int2ip(val))

    if return_ccache:
        ccache = init_ccache(None, ccachef, start_time, duration, mask, fmt, oddballs, maxmissing)
        return svl, fmt, dict(anssets), ccache
    else:
        return svl, fmt, dict(anssets)


def get_big_svl(start_time, duration=30000, mask=32, fmt=None, country_set=None,
        oddballs=True, maxmissing=0, return_ccache=True,
        ccachef=df.rightdir(statedir+"pickles/")+"ccache.pickle",
        mindomsize=1):
    '''
    see get_svl()
    '''
    dur = 60*60*24
    end_time = start_time + duration

    anssets = defaultdict(set)
    myfmt = set()

    svl, tmp_fmt, tmp_anssets = get_svl(start_time, duration=dur,
            mask=mask, fmt=None, country_set=country_set,
            oddballs=oddballs, maxmissing=1000, return_ccache=False,
            ccachef=ccachef, mindomsize=1)
    svld = dict()
    for sv in svl:
        svld[sv.get_id()] = sv
    for dom in tmp_anssets.keys():
        anssets[dom] |= tmp_anssets[dom]
    myfmt |= set(tmp_fmt)
    start_time += dur

    while start_time+dur <= end_time:
        tmp_svl, tmp_fmt, tmp_anssets = get_svl(start_time, duration=dur,
            mask=mask, fmt=None, country_set=country_set,
            oddballs=oddballs, maxmissing=1000, return_ccache=False,
            ccachef=ccachef, mindomsize=1)

        for dom in tmp_anssets.keys():
            anssets[dom] |= tmp_anssets[dom]
        myfmt |= set(tmp_fmt)

        for sv in tmp_svl:
            svid = sv.get_id()
            if svid in svld:
                svld[svid].absorb(sv)
            else:
                svld[svid] = sv

        start_time += dur

    svl = svld.values()
    if fmt is None:
        fmt = list(myfmt)
    for dom in fmt:
        if len(anssets[dom]) < mindomsize \
          or len([sv for sv in svl if dom not in sv]) > 0.5*len(svl):
            del anssets[dom]
    svl = reduce_svl(svl, fmt, maxmissing)

    if return_ccache:
        ccache = init_ccache(None, ccachef, start_time, duration, mask, fmt, oddballs, maxmissing)
        return svl, fmt, dict(anssets), ccache
    else:
        return svl, fmt, dict(anssets)


def reduce_svl(svl, fmt, maxmissing=0):
    '''
    :return: new svl that only contains clients that are missing upto mxli
    '''
    out_svl = list()
    sfmt = set(fmt)
    for sv in svl:
        sdoms = set(sv.vec.keys())
        if len(sfmt.difference(sdoms)) <= maxmissing:
            out_svl.append(sv)
    return out_svl


class closeness_cache:
    '''
    cache so that closeness calculations don't have to be repeated
    '''
    def __init__(self, f):
        self.cache = dict()
        self.tmp_item = None
        self.f = f
        self.hash = "random"
        self.indirect_cache = dict()
        self.changed = False
        self.hits = 0
        self.misses = 0

    def get_closeness(self, a, b):
        key = tuple(sorted([a.id, b.id]))
        if key in self.cache:
            self.hits += 1
            return self.cache[key]
        else:
            self.misses += 1
            tmp = closeness(a, b)
            self.cache[key] = tmp
            self.changed = True
            return tmp

    def __getitem__(self, b):
        # if this works, it should let you do cc[a][b] to get closeness(a,b)
        if self.tmp_item is None:
            self.tmp_item = b
            return self
        else:
            a = self.tmp_item
            self.tmp_item = None
            return self.get_closeness(a, b)

    def load(self, *args, **kwargs):
        cache = df.picklein(self.f)
        invals = df.make_hashable(list(args)+zip(kwargs.keys(), kwargs.values()))
        self.hash = invals
        if cache is not None:
            self.indirect_cache = cache
            if invals in cache:
                tmp = df.picklein(cache[invals])
                if tmp is not None:
                    self.cache = tmp
                    self.hits = 0
                    self.misses = 0
                    logger.warning("loaded "+str(cache[invals]))
                    return
        logger.warning("fresh cache")

    def dump(self):
        if self.changed:
            invals = self.hash
            cache = self.indirect_cache
            logger.warning("dumping, DO NOT ctrl-C!")
            if invals not in cache:
                cache[invals] = df.rightdir(statedir+"pickles")+"closeness"+str(time.time())+".pickle"
                df.pickleout(self.f, cache)
            df.pickleout(cache[invals], self.cache)
            logger.warning("safe to ctrl-C")
            self.changed = False
            time.sleep(2)
        if self.hits + self.misses > 0:
            logger.debug("hit rate: "+ \
                    str(float(self.hits)/float(self.hits+self.misses)) \
                    +", total hits: "+str(self.hits))


def init_ccache(ccache, f, *args, **kwargs):
    if ccache is None:
        c = closeness_cache(f)
        c.load(*args, **kwargs)
        return c
    else:
        return ccache


