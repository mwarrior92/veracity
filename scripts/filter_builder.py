import datetime
from warriorpy.shorthand import diriofile as df

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
#                            CODE
##################################################################

def datetime_to_unixsecs(*datetime_args):
    t = datetime.datetime(*datetime_args)
    et = datetime.datetime(1970,1,1)
    return (t-et).total_seconds()


class query_filter:
    def __init__(self, meas_type, af=4, public_only=True, status=[2, 4, 5, 8]):
        # IPv4

        self.filter = dict()
        self.filter['af'] = af
        # ongoing or stopped measurements
        self.filter['status'] = status
        self.filter['type'] = meas_type
        self.filter['is_public'] = public_only
        self.filter['optional_fields'] = 'probes'


    def set_time(self, lenience=300, *datetime_args):
        t = datetime_to_unixsecs(*datetime_args)
        self.filter['start_time__lte'] = t+lenience
        self.filter['stop_time__gte'] = t-lenience


    def set_start_time(self, *datetime_args):
        self.filter['start_time__lte'] = datetime_to_unixsecs(*datetime_args)


    def set_stop_time(self, *datetime_args):
        self.filter['stop_time__gte'] = datetime_to_unixsecs(*datetime_args)


    def set_time_window(self, duration, *datetime_args):
        '''
        duration -> length of time window
        '''
        t = datetime_to_unixsecs(*datetime_args)
        self.filter['start_time__lte'] = t+duration
        #self.filter['start_time__lte'] = t
        #self.filter['stop_time__gte'] = t-duration


    def remove_filter(self, key):
        if key in self.filter:
            del self.filter[key]


    def allow_private(self):
        self.remove_filter('is_public')


    def set_af(self, af):
        self.filter['af'] = af


    def set_interval(self, target_interval=240, lenience=60):
        '''
        lenience -> how much can results deviate from target interval?
        '''
        self.filter['interval__gte'] = target_interval - lenience
        self.filter['interval__lte'] = target_interval + lenience


    def set_status(self, status):
        self.filter['status'] = status


    def set_dst(self, target=None, prefix=None, suffix=None, contains=None):
        if target is not None:
            self.filter['target'] = target
        if prefix is not None:
            self.filter['target__startswith'] = target
        if suffix is not None:
            self.filter['target__endswith'] = target
        if contains is not None:
            self.filter['target__contains'] = target


    def set_nameserver(self, **kwas):
        '''
        see set_dst; set_nameserver() is only for clarity / disambiguation. See
        the "query_argument" response field for comparison
        '''
        self.set_dst(**kwas)


    def manual_set(self, **kwas):
        '''
        use this to set any field(s)
        '''
        for key in kwas:
            self.filter[key] = kwas[key]


    def get_filter(self):
        return self.filter


class result_filter:
    # for most of the values, use a set; this allows for use of "in"
    # if you have multple values
    def __init__(self, meas_type):
        self.filter = dict()
        self.type = meas_type


    def set_min_probes(self, size):
        self.filter['participant_count'] = size


    def set_domain(d):
        # d should be a set of domains
        if self.type == 'dns':
            self.filter['query_argument'] = d
        else:
            self.filter['target'] = d


    def set_time(self, lenience=300, *datetime_args):
        t = datetime_to_unixsecs(*datetime_args)
        self.filter['start_time'] = t-lenience
        self.filter['stop_time'] = t+lenience


    def set_start_time(self, *datetime_args):
        self.filter['start_time'] = datetime_to_unixsecs(*datetime_args)


    def set_stop_time(self, *datetime_args):
        self.filter['stop_time'] = datetime_to_unixsecs(*datetime_args)


    def set_time_window(self, duration, *datetime_args):
        '''
        duration -> length of time window
        '''
        t = datetime_to_unixsecs(*datetime_args)
        self.filter['start_time'] = t+duration
        self.filter['stop_time'] = t


    def allow_private(self):
        del self.filter['is_public']


    def set_af(self, af):
        self.filter['af'] = af


    def set_interval(self, target_interval=240, lenience=60):
        '''
        lenience -> how much can results deviate from target interval?
        '''
        self.filter['interval__gte'] = target_interval - lenience
        self.filter['interval__lte'] = target_interval + lenience


    def set_status(self, status):
        self.filter['status'] = status


    def set_dst(self, target=None, prefix=None, suffix=None, contains=None):
        if target is not None:
            self.filter['target'] = target
        if prefix is not None:
            self.filter['target__startswith'] = target
        if suffix is not None:
            self.filter['target__endswith'] = target
        if contains is not None:
            self.filter['target__contains'] = target


    def set_nameserver(self, **kwas):
        '''
        see set_dst; set_nameserver() is only for clarity / disambiguation. See
        the "query_argument" response field for comparison
        '''
        self.set_dst(**kwas)


    def manual_set(self, **kwas):
        '''
        use this to set any field(s)
        '''
        for key in kwas:
            self.filter[key] = kwas[key]


    def get_filter(self):
        return self.filter


