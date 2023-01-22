import numpy as np
from PIL import Image

# [' ', '░', '▒', '▓', '█']


def dither_text(array, chars=[' ', '▒', '█'], size=(200, 200)):
    chars = chars[::-1]
    block_size = int(array.shape[1] / size[1]), int(array.shape[0] / size[0])
    array = np.array(Image.fromarray(array).resize(
        (block_size[0] * size[0], block_size[1] * size[1])))

    cond = np.zeros(size)
    for x in range(size[1]):
        for y in range(size[0]):
            cond[y,
                 x] = np.sum(array[x * block_size[1]:(x + 1) * block_size[1],
                                   y * block_size[0]:(y + 1) * block_size[0]])

    cond = cond / cond.max()
    cond = (cond * (len(chars) - 1)).astype(int)

    file = open('dithered.txt', 'w+')
    for y in range(size[1]):
        line = ''.join([chars[i] for i in cond[:, y]])
        print(line, file=file)
        print(line)


def dither_braille(array, size=(10, 10), thresh=0.5):
    chars = []
    for c in [0, 1, 2, 3, 8, 9, 10, 11]:
        s = 8 * (c % 2)
        for s in [0, 8]:
            chars.append([chr(int('28' + hex(c)[2:] + hex(i)[2:],
                                  16)) for i in range(s,
                                                      8 + s)] + [chr(int('28' + hex(c + 4)[2:] + hex(i)[2:],
                                                                         16)) for i in range(s,
                                                                                             8 + s)])

    chars[0][0] = '⢀'

    block_size = int(array.shape[1] /
                     size[1]), int(array.shape[0] /
                                   size[0])  # width, height
    # height should be 4n, width should be 2m
    block_size = max(4, block_size[1] -
                     (block_size[1] %
                      4)), max(2, block_size[0] -
                               (block_size[0] %
                                2))

    array = np.array(Image.fromarray(array).resize(
        (block_size[1] * size[1], block_size[0] * size[0])))
    array = 1 - array / array.max()
    cs = []
    for x in range(size[0]):
        cs_ = []
        for y in range(size[1]):
            sub = array[x *
                        block_size[0]:(x +
                                       1) *
                        block_size[0], y *
                        block_size[1]:(y +
                                       1) *
                        block_size[1]]
            ls = sub[:, :sub.shape[1] // 2]
            rs = sub[:, sub.shape[1] // 2:]
            s = ls.shape[0] // 4
            lss = [np.sum(ls[i * s:(i + 1) * s, :]) //
                   (ls.shape[1] * thresh) for i in range(4)]
            l = int(sum(min(v, 1) * 2**j for j, v in enumerate(lss)))
            rss = [np.sum(rs[i * s:(i + 1) * s, :]) //
                   (rs.shape[1] * thresh) for i in range(4)]
            r = int(sum(min(v, 1) * 2**j for j, v in enumerate(rss)))
            cs_.append(chars[r][l])
        cs.append(''.join(cs_))
    file = open('dithered.txt', 'w+')
    [print(c, file=file) for c in cs]
    [print(c) for c in cs]


# inp = r"Z:\Downloads\IMG-20220813-WA0007.jpg"
inp = '/Users/yumanhordijk/PhD/ychem/utility/Sxazg5f.png'
img = np.array(Image.open(inp).convert('L'))

dither_text(img, size=(120, 20))
dither_braille(img, size=(50, 200))
