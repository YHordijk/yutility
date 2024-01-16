from tcutility import results, log
import os

data = []
calc_dir = 'tmp/functional_test'
for d in os.listdir(calc_dir):
    res = results.read(os.path.join(calc_dir, d))
    status = res.status.code
    emoji = {'R': log.Emojis.wait, 'PD': log.Emojis.sleep, 'S': log.Emojis.good, 'F': log.Emojis.fail, 'W': log.Emojis.warning, 'U': log.Emojis.question}.get(status, '') + ' '
    data.append((d.split('.')[0], emoji + res.status.name, '{: >6.2f}'.format(res.timing.total or 0), d.split('.')[1]))

data = sorted(data, key=lambda x: int(x[0]))
log.table(data, ['ID', 'STATUS', 'TIMING (s)', 'FUNCTIONAL'])