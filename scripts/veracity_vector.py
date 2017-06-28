from warriorpy.shorthand import diriofile as df
from warriorpy.net_tools import ipparsing as ipp
from pymongo import MongoClient

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
def get_window(start_time, duration, domain_set):
    res = data.aggregate([
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
                    ], allowDiskUse=True)
    return res


# convert from set of raw documents from single time block to dictionary of
# domains vs answer-sets. If there is more than one answer for a domain, I take
# the first answer
def window_to_dict(window):
    d = dict()
    for doc in window:
        if doc['_id']['probe_ip'] not in d:
            d[doc['_id']['probe_ip']] = dict()
            d[doc['_id']['probe_ip']]['ind'] = doc['ind']
        if len(doc['answers']) > 0:
            d[doc['_id']['probe_ip']][doc['_id']['domain']] = doc['answers'][0]
        else:
            d[doc['_id']['probe_ip']][doc['_id']['domain']] = 0
    return d


# generate vector from dict of {domain: answer} entries
def dict_to_vector(d, fmt=None):
    vector = list()
    if fmt is None:
        fmt = domains[:3]
    for domain in fmt:
        if domain in d:
            vector.append(d[domain])
        else:
            vector.append(0)
    return vector


def weight_vector(vec, wvec):
    pass
