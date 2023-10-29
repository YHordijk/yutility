import sys
from datetime import datetime
from time import perf_counter
import numpy as np
import itertools
import json
from yutility import dictfunc
from strip_ansi import strip_ansi

logfile = sys.stdout
tab_level = 0
max_width = 0
print_date = True
use_colors = False


class Emojis:
    wait = '🕒'
    good = '✅'
    cancel = '🛑'
    sleep = '💤'
    fail = '❌'
    send = '📤'
    receive = '📥'
    empty = '⠀⠀'
    finish = '🏁'
    warning = '⚠️'
    question = '❔'
    info = 'ℹ️'
    rarrow = '─>'
    larrow = '<─'
    lrarrow = rlarrow = '<─>'
    angstrom = 'Å' 

    def __getitem__(cls, key):
        return getattr(cls, key)


emojis = Emojis()

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


def log(message='', end='\n', print_time_stamp=True):
    if isinstance(message, dict):
        lst = dictfunc.dict_to_list(message)
        lst_ = []
        for line in lst:
            lst_.append([])
            for item in line:
                if isinstance(item, (str, int, float, bool)):
                    lst_[-1].append(item)
                elif isinstance(item, list):
                    lst_[-1].append([element if isinstance(item, (str, int, float, bool)) else repr(element) for element in item])
                else:
                    lst_[-1].append(repr(item))
        message = dictfunc.list_to_dict(lst_)
        message = json.dumps(message, indent=4, sort_keys=True)
    message = str(message)
    if not use_colors:
        message = strip_ansi(message)
    message = message.split('\n')
    for m in message:
        if max_width > 0 and len(m) > max_width:
            m = m[:max_width - 4] + ' ...' 

        m = '\t'*tab_level + m
        if print_date and print_time_stamp:
            m = time_stamp() + m

        print(m, file=logfile, end=end, flush=True)


arrow_prefix = []
def arrow(txt='', tags=['straight']):
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

        hline = [h+len(items) for h in hline if h < 0]

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
    Nsteps = N if logfile.isatty() else Nsteps

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


def loading_bar2(sequence, comment='', Nsegments=50, Nsteps=10, start_char='├', end_char='│', fill_char='─', empty_char=' ', center_char='>'):
    N = len(sequence)
    Nsteps = N if logfile.isatty() else Nsteps
    Ndigits = int(np.log10(N))+1
    max_length = 0

    loading_bar_start_time = perf_counter()
    for i, val in enumerate(sequence):
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

            s = f'{i:{Ndigits}}/{N} {start_char if i > 0 else "|"}{fill_seg + center_seg + empty_seg}{end_char} {i/N:6.1%} ETA: {eta}s {comment}'.ljust(max_length)
            max_length = max(len(s), max_length)
            log(s, end='\r')
        
        yield val
    s = bcolors.OKGREEN + f'{i+1:{Ndigits}}/{N} {start_char}{fill_char*(Nsegments+1)}┤ 100.0% {comment}'.ljust(max_length) + bcolors.ENDC
    log(s, end='\r')
    log()


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

    cell_widths = [max(len(str(x)) for x in a[:, i]) for i in range(a.shape[1])]

    rows = []
    for i, x in enumerate(a):
        if i == 0:
            row1 = corners[0]
        if ylabels == []:
            row2 = '│'
        else:
            row2 = '│'
        if i == a.shape[0]-1:
            row3 = corners[3]
        else:
            row3 = '├'

        for j, y in enumerate(x):
            cw = cell_widths[j]
            if j == (a.shape[1]-1) and i == 0:
                row1 += '─'*cw + '──' + corners[1]
            elif i == 0:
                row1 += '─'*cw + '──┬'

            if not empty_where[i, j]:
                if dither:
                    d = dither_chars[int(y/a.max() * (len(dither_chars))-1)]
                    row2 += f' {d.center(cw)} │'
                else:
                    row2 += f' {str(y).center(cw)} │'
            else:
                row2 += ' '*cw + '  │'

            if j == (a.shape[1]-1) and ylabels != []:
                row2 += ' '*cw + ' ' + ylabels[i]

            if i == a.shape[0]-1:
                if j == (a.shape[1]-1):
                    row3 += '─'*cw + '──' + corners[2]
                else:
                    row3 += '─'*cw + '──┴'
            else:
                if j == (a.shape[1]-1):
                    row3 += '─'*cw + '──┤'
                else:
                    row3 += '─'*cw + '──┼'

        if i == 0:
            rows.append(row1)

        rows.append(row2)
        rows.append(row3)

    for r in rows:
        log(r)

    if xlabels is not None:
        for yl in itertools.zip_longest(*xlabels, fillvalue=' '):
            log(' ' + ''.join([f' {l} '.center(cell_widths[j]+3) for j, l in enumerate(yl)]))


def boxed_text(txt, round_edge=True, align='left', double_edge=False, title=None, title_align='left', color=''):
    endc = bcolors.ENDC if color else ''

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

    # build first row
    if title is not None:
        if title_align == 'left':
            s = corners[0] + (' ' + title + ' ').ljust(maxlen+2, straights[1]) + corners[1] + '\n'
        elif title_align == 'right':
            s = corners[0] + (' ' + title + ' ').rjust(maxlen+2, straights[1]) + corners[1] + '\n'
        else:
            s = corners[0] + (' ' + title + ' ').center(maxlen+2, straights[1]) + corners[1] + '\n'
    else:
        s = corners[0] + straights[1]*(maxlen+2) + corners[1] + '\n'
    # build main body of box
    for txt in txts:
        if align == 'left':
            s += f'{straights[0]} ' + color + txt.ljust(maxlen) + endc + f' {straights[0]}\n'
        if align == 'right':
            s += f'{straights[0]} ' + color + txt.rjust(maxlen) + endc + f' {straights[0]}\n'
        if align == 'center':
            s += f'{straights[0]} ' + color + txt.center(maxlen) + endc + f' {straights[0]}\n'
    # build final row
    s += corners[3] + straights[1]*(maxlen+2) + corners[2] + '\n'

    log(s.removesuffix('\n'))


class BoxedText:
    def __init__(self, title=None):
        self.lines = []
        self.title = title

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        ...


def info(txt):
    boxed_text(txt, round_edge=True, double_edge=False, title='Info', color=bcolors.OKBLUE)


def warn(txt):
    boxed_text(txt, double_edge=True, title='Warning', title_align='center', color=bcolors.WARNING)


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


if __name__ == '__main__':
    boxed_text('testing, 1, 2, 3, 4, 5, 6, 7, 8\nsecond row, 1, 2, 3, 4, 5, 6\n...\n...\n\n...\n...\nLast row here', title='ReactionRunner')
    mat = np.arange(3 * 3).reshape(3, 3)
    print_matrix(mat, xlabels=['a', 'ab', 'abc'], ylabels=['a', 'ab', 'abc'])
    mat = np.arange(3 * 4).reshape(3, 4)
    mat[0, 0] = 41203
    print_matrix(mat, xlabels=['a', 'ab', 'abc', 'abcd'], ylabels=['a', 'ab', 'abc'], round_edge=False)
    X, Y = np.meshgrid(np.linspace(-10, 10, 40), np.linspace(-10, 10, 20))
    mat = np.exp(-(X**2 + Y**2)/100)
    print_image(mat)
    mat = np.cos(X/2) * np.sin(Y/2)
    print_image(mat, draw_edge=False)
    log()
    arrow('Start', ['start'])
    arrow()
    arrow('First step', ['split'])
    arrow('First substep of first step', ['straight', 'split'])
    arrow('Second substep of first step', ['straight', 'split'])
    arrow('Third and final substep of first step', ['straight', 'end'])
    arrow()
    arrow('Second step', ['split'])
    arrow('Substep of second step', ['straight', 'split'])
    arrow('Substep of substep of second step', ['straight', 'straight', 'end'])
    arrow('Substep of substep of substep of second step', ['straight', 'straight', 'skip', 'end'])
    arrow('Substep of second step', ['straight', 'end'])
    arrow()
    arrow('Final step', ['split'])
    arrow()
    arrow(f'{emojis.good} The end', ['startinv'])

    info('This is important info')
    warn('This is an important warning!')

    import time
    for x in loading_bar2(range(100), 'Sleeping test'):
        time.sleep(.01)
