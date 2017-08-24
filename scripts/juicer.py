import ripe.atlas.cousteau as rac
from warriorpy.shorthand import diriofile as df
from warriorpy.ripe import ripe_api as ra
import filter_builder as fb
from pymongo import MongoClient

##################################################################
#                           LOGGING
##################################################################
import logging
import logging.config


logging.config.fileConfig('logging.conf', disable_existing_loggers=False)

# create logger
logger = logging.getLogger(__name__)
logger.debug(__name__+"logger loaded")

# database setup
mclient = MongoClient()
db = mclient.veracity
pdata = db.probe_data

# fields that are specific to probes
pfields = ['country', 'asn_v4', 'asn_v6', 'geometry', 'probe_ip']

##################################################################
#                        GLOBAL AND SETUP
##################################################################


# paths
basedir = df.getdir(__file__)+'../'
statedir = df.rightdir(basedir+"state/")

# external parameters
params = df.getlines(statedir+'params')
creation_key = params[0]
destruction_key = params[1]


##################################################################
#                            CODE
##################################################################


def get_measurements(filters):
    '''
    see https://atlas.ripe.net/docs/api/v2/reference/#!/measurements/ for
    details other optional query params
    '''
    measurements = rac.MeasurementRequest(**filters)
    return measurements


def prime_measurements(measurements, *fields):
    cache = dict()
    out = list()
    def get_probe_info(pid):
        if pid in cache:
            return cache[pid]
        tmp = list(pdata.find({'probe_id': pid}).limit(1))
        if len(tmp) > 0:
            probe_info = tmp[0]
            cache[pid] = tmp[0]
            return tmp[0]
        else:
            b, probe_info = ra.get_probe_info(pid)
            if b:
                pdata.insert_one(probe_info)
                cache[pid] = probe_info
                return probe_info
            else:
                return {}

    for meas in measurements:
        for f in fields:
            if f == 'domain':
                if meas['type'] == 'dns':
                    if 'query_argument' in meas:
                        if meas['query_argument'] is not None and\
                        meas['query_argument'] != 'null':
                            vd = meas['query_argument'].split('.')
                            vd = [z for z in vd if len(z)>20]
                            if len(vd) == 0:
                                meas[f] = meas['query_argument']
                            else:
                                meas[f] = 'RANDOM.'+'.'.join(vd[1:])
                            if meas['query_argument'] in ['.', 'ip.server',
                                    'hostname.bind']:
                                meas[f] = 'OPERATOR_QUERY'
                        else:
                            meas[f] = 'NONE'
                    else:
                        meas[f] = 'NONE'
                else:
                    vd = meas['target'].split('.')
                    vd = [z for z in vd if not z.isdigit()]
                    if len(vd) > 0:
                        meas[f] = meas['target']
                    else:
                        meas[f] = 'direct_to_IP'
            elif f in pfields:
                probes = [z['id'] for z in meas['probes']]
                for probe in probes:
                    pinfo = get_probe_info(probe)
                    if f in pinfo:
                        if probe not in meas:
                            meas[probe] = dict()
                        meas[probe][f] = pinfo[f]
        out.append(meas)

    return out


def filter_measurements(measurements, filters):
    out = list()
    for meas in measurements:
        good = True
        for key in filters:
            if key == 'start_time':
                if key in meas:
                    if filters[key] <= meas[key]:
                        continue
                    else:
                        break
            elif key == 'participant_count' or key == 'stop_time':
                if key in meas:
                    if meas[key] is not None and meas[key] != 'null':
                        if filters[key] >= meas[key]:
                            continue
                        else:
                            break
            elif 'time' in key or 'interval' in key:
                # TODO handling times here is complicated;
                # better to use list comprehensions on whole list later
                pass
            elif key not in meas:
                good = False
                break
            if type(meas[key]) is list:
                if len(filters[key].intersection(meas[key]))<len(filters[key]):
                    good = False
                    break
            elif meas[key] not in filters[key]:
                good = False
                break
        if good:
            out.append(meas)
    return out


def dmaximize_sample(measurements, dvals, dfield, mfield=None, mprobe=None,
        minsize=10):
    '''
    :param dvals: list of discrete values
    :param valfield: field to check
    :param mfield: meas field trait to maximize
    :param mprobe: probe field trait to maximize

    this aims to maximize the quantity of measurements (or probes performing
    said measurements) that are distinct by some field (mfield or mprobe) by
    choosing the optimal subset of dvals

    NOTE: you may need to prime the measurements data to include the field
    you're looking for
    '''
    mem = dict()

    def get_intersection(s):
        if s in mem:
            return mem[s]
        elif len(s) > 2:
            h = len(s)/2
            a = s[0:h]
            b = s[h:]
            A = get_intersection(a)
            B = get_intersection(b)
            tmp = A.intersection(B)
            mem[s] = tmp
            return tmp
        elif len(s) == 2:
            if dfield in pfields:
                A = set([z for z in measurements if z[dfield] == dvals[s[0]]])
                B = set([z for z in measurements if z[dfield] == dvals[s[1]]])
            else:
                A = set()
                B = set()
                for meas in measurements:
                    for probe in meas['probes']:
                        if meas[probe][dfield] == dvals[s[0]]:
                            A.add(probe)
                        elif meas[probe][dfield] == dvals[s[1]]:
                            B.add(probe)
            tmp = A.intersection(B)
            mem[s] = tmp
            return tmp
        else:
            tmp = set([z for z in measurements if z[dfield] == dvals[s[0]]])
            mem[s] = tmp
            return tmp


    mask = [False]*len(dvals)
    sets = defaultdict(list)
    for __ in xrange(0, 2**len(mask)):
        # increment
        t = True
        i = 0
        while t:
            t = mask[i]
            mask[i] = not mask[i]
            i = (i+1)%len(mask)

        inds = tuple([i for i, z in enumerate(mask) if z])
        subset = get_intersection(inds)
        sets[len(subset)].append(inds)

    sizes = sorted(sets.keys(), reverse=True)
    for size in sizes:
        if len(max(inds, key=lambda i: len(i))) >= minsize:
            tmp = max(inds, key=lambda i: len(i))
            vals = [dvals[z] for z in tmp]
            return size, vals
    return None, None
