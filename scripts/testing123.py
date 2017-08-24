import juicer as jcr
import filter_builder as fb
from warriorpy.shorthand import diriofile as df

qfb = fb.query_filter('dns')
qfb.set_time_window(60*60*24, 2017, 7, 21)
qfb.manual_set(resolve_on_probe=True)

qfs = qfb.get_filter()

print qfs,"\n"

msms = jcr.get_measurements(qfs)

print "priming..."
msms = jcr.prime_measurements(msms, 'domain')


rfb = fb.result_filter('dns')
rfb.set_min_probes(1000)
rfb.set_time_window(60*60*24, 2017, 7, 21)
rfs = rfb.get_filter()

print rfs
print "filtering..."
msms = jcr.filter_measurements(msms, rfs)


g = None
doms = set()
for m in msms:
    if m['domain'] is not None:
        if 'ripe' not in m['domain']:
            doms.add((m['domain'], m['id']))

df.overwrite('doms.csv', df.list2col(sorted(doms, key=lambda z: z[0])))
df.pickleout('msms.pickle', msms)
