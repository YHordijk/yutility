import sys
from datetime import datetime
from time import perf_counter
import numpy as np
import itertools
import json
from math import floor, ceil

logfile = sys.stdout
tab_level = 0
max_width = 200
print_date = True


emojis = {
    'wait':     '🕒', 
    'good':     '✅', 
    'cancel':   '🛑', 
    'sleep':    '💤',
    'fail':     '❌',
    'send':     '📤',
    'receive':  '📥',
    'empty':    '⠀⠀',
    'finish':   '🏁',
    'warning':  '⚠️',
    'question': '❔ ',
}

rarrow = '─>'
larrow = '<─'
lrarrow = rlarrow = '<─>'
angstrom = 'Å' 


class NoPrint:
    def __init__(self, stdout=None, stderr=None):
        self.stdout = stdout
        self.stderr = stderr
        if self.stdout is None:
            self.stdout = sys.stdout
        if self.stderr is None:
            self.stderr = sys.stderr

    def __enter__(self):
        sys.stdout = None
        sys.stderr = None

    def __exit__(self, *args, **kwargs):
        sys.stdout = self.stdout
        sys.stderr = self.stderr


def time_stamp():
    now = datetime.now()
    return f'[{now.year}/{str(now.month).zfill(2)}/{str(now.day).zfill(2)} {str(now.hour).zfill(2)}:{str(now.minute).zfill(2)}:{str(now.second).zfill(2)}] '


def log(message='', /, end='\n'):
    if type(message) is dict:
        message = json.dumps(message, indent=4, sort_keys=True)
    message = str(message)
    message = message.split('\n')
    for m in message:
        # print(m)
        # m = m.encode('utf-8')
        # print(m)
        if max_width > 0 and len(m) > max_width:
            m = m[:max_width - 4] + ' ...' 
        if print_date:
            print(time_stamp() + '\t'*tab_level + m, file=logfile, end=end, flush=True)
        else:
            print('\t'*tab_level + m, file=logfile, end=end, flush=True)
    # logfile.flush()


arrow_prefix = []
def arrow(txt, tags=['split']):
    s = ''
    for tag in arrow_prefix + tags:
        if tag == 'start':
            s += '┯ '
        elif tag == 'startinv':
            s += '┷ '
        elif tag == 'end':
            s += '╰─> '
        elif tag == 'straight':
            s += '│   '
        elif tag == 'split':
            s += '├─> '
        elif tag == 'skip':
            s += '    '
        elif tag == 'vert':
            s += '────'
        elif tag == 'endvert':
            s += '╰──> '
        elif tag == 'splitvert':
            s += '├──> '
    log(s + txt)


def print_list(items, sep='   ', header=None, return_str=False, hline=[]):
    if return_str:
        returns = ''
    # if two-dimensional
    if type(items[0]) in [list, set, tuple]:
        column_lens = [max([len(str(r[i])) for r in items]) for i in range(len(items[0]))]
        if type(sep) in [list, tuple]:
            sep_lens = [len(s) for s in sep]
        else:
            sep_lens = [len(sep) for _ in range(len(items[0]))]

        hline = [h % len(items[0]) - 1 for h in hline]

        if header is not None:
            header_lens = [len(head) - head.count('⠀')//2 for head in header]
            column_lens = [max(column_lens[i], header_lens[i]) for i, head in enumerate(header)]
            line = ''
            for i in range(len(header)):
                line += str(header[i]).ljust(column_lens[i])
                if i < len(header) - 1:
                    line += ' '*sep_lens[i]
            if return_str:
                returns += line + '\n'
                returns += '─' * (sum(column_lens) + sum(sep_lens)) + '\n'
            else:
                log(line)
                log('─' * (sum(column_lens) + sum(sep_lens)))

        for col, item in enumerate(items):
            line = ''
            for i in range(len(item)):
                line += str(item[i]).ljust(column_lens[i])
                if i < len(item) - 1:
                    if type(sep) in [list, tuple]:
                        s = sep[i]
                    else:
                        s = sep
                    line += s

            if return_str:
                returns += line + '\n'
                if col in hline:
                    returns += '─' * (sum(column_lens) + sum(sep_lens))
            else:
                log(line)
                if col in hline:
                    log('─' * (sum(column_lens) + sum(sep_lens)))
    else:
        for i, item in enumerate(items):
            log(str(item))

    if return_str:
        return returns


loading_bar_start_time = 0
def loading_bar(i, N, Nsegments=50, Nsteps=10, start_char='├', end_char='│', fill_char='─', empty_char=' ', center_char='>', comment=''):
    N = max(1, N)

    global loading_bar_start_time
    if loading_bar_start_time == 0:
        loading_bar_start_time = perf_counter()

    if i % (N//min(N, Nsteps)) == 0 or i == N:
        segment = int(i/N*Nsegments)
        fill_seg = fill_char*segment
        if segment == Nsegments:
            center_seg = fill_char
        elif i == 0:
            center_seg = empty_char
            
        else:
            center_seg = center_char
        empty_seg = empty_char*(Nsegments-segment)

        if i == 0:
            eta = '???'
        else:
            eta = f'{(perf_counter() - loading_bar_start_time)/max(i,1) * (N-i):.1f}'
        log(f'{i:{int(np.log10(N))+1}}/{N} {start_char}{fill_seg + center_seg + empty_seg}{end_char} {i/N:6.1%} ETA: {eta}s {comment}', end='\r')

        if i == N:
            log()
            loading_bar_start_time = 0


def print_image(a, draw_edge=True, round_edge=True):
    dither_chars = [' ', '░', '▒', '▓', '█']
    a = (a - a.min())/(a.max() - a.min()) * (len(dither_chars)-1)
    a = a.astype(int)
    width = a.shape[1]
    if draw_edge:
        s = '╭' if round_edge else '┌'
        s += '─' * width
        s += '╮' if round_edge else '┐'
        log(s)

    for row in a:
        s = '│' if draw_edge else ''
        for y in row:
            d = dither_chars[y]
            s += d
        s += '│' if draw_edge else ''
        log(s)

    if draw_edge:
        s = '╰' if round_edge else '└'
        s += '─' * width
        s += '╯' if round_edge else '┘'
        log(s)


def print_matrix(a, xlabels=[], ylabels=[], round_edge=True, dither=False, empty_where=None):
    a = np.array(a)
    if empty_where is None:
        empty_where = np.zeros_like(a)
    if round_edge:
        corners = ['╭', '╮', '╯', '╰']
    else:
        corners = ['┌', '┐', '┘', '└']

    dither_chars = [' ', '░', '▒', '▓', '█']

    if xlabels == []:
        ylabel_max_len = 0
    else:
        ylabel_max_len = max(len(yl) for yl in ylabels) + 1
    rows = []
    for i, x in enumerate(a):
        if i == 0:
            row1 = ' '*ylabel_max_len + corners[0]
        if ylabels == []:
            row2 = ' '*ylabel_max_len + '│'
        else:
            row2 = ylabels[i].rjust(ylabel_max_len) + '│'
        if i == a.shape[0]-1:
            row3 = ' '*ylabel_max_len + corners[3]
        else:
            row3 = ' '*ylabel_max_len + '├'

        for j, y in enumerate(x):
            if j == (a.shape[1]-1) and i == 0:
                row1 += '───' + corners[1]
            elif i == 0:
                row1 += '───┬'

            if not empty_where[i, j]:
                if dither:
                    d = dither_chars[int(y/a.max() * (len(dither_chars))-1)]
                    row2 += f' {d} │'
                else:
                    row2 += f' {y} │'
            else:
                row2 += '   │'

            if i == a.shape[0]-1:
                if j == (a.shape[1]-1):
                    row3 += '───' + corners[2]
                else:
                    row3 += '───┴'
            else:
                if j == (a.shape[1]-1):
                    row3 += '───┤'
                else:
                    row3 += '───┼'

        if i == 0:
            rows.append(row1)

        rows.append(row2)
        rows.append(row3)

    for r in rows:
        log(r)

    if xlabels is not None:
        for yl in itertools.zip_longest(*xlabels, fillvalue=' '):
            log(' '*ylabel_max_len + ''.join([f'  {l} ' for l in yl]))


def boxed_text(txt, round_edge=True, align='left', double_edge=True, title=None):
    straights = ['│', '─']
    if round_edge and not double_edge:
        corners = ['╭', '╮', '╯', '╰']
    elif double_edge:
        corners = ['╔', '╗', '╝', '╚']
        straights = ['║', '═']
    else:
        corners = ['┌', '┐', '┘', '└']

    txts = txt.split('\n')
    maxlen = max(len(txt) for txt in txts)

    if title is not None:
        s = corners[0] + straights[1]*floor((maxlen+2 - len(title))/2-1) + ' ' + title + ' ' + straights[1]*ceil((maxlen+2 - len(title))/2-1) + corners[1] + '\n'
    else:
        s = corners[0] + straights[1]*(maxlen+2) + corners[1] + '\n'
    for txt in txts:
        if align == 'left':
            s += f'{straights[0]} ' + txt.ljust(maxlen) + f' {straights[0]}\n'
        if align == 'right':
            s += f'{straights[0]} ' + txt.rjust(maxlen) + f' {straights[0]}\n'
        if align == 'center':
            s += f'{straights[0]} ' + txt.center(maxlen) + f' {straights[0]}\n'
    s += corners[3] + straights[1]*(maxlen+2) + corners[2] + '\n'

    log(s)
