
def lev(a, b):
    if len(b) == 0:
        return len(a)
    if len(a) == 0:
        return len(b)

    def tail(s):
        return s[1:]

    if a[0] == b[0]:
        return lev(tail(a), tail(b))

    return 1 + min(lev(tail(a), b), lev(a, tail(b)), lev(tail(a), tail(b)))


def get_closest(a, others, compare_func=lev):
	lowest_dist = None
	lowest_strs = []
	for other in others:
		dist = compare_func(a, other)
		if lowest_dist is None:
			lowest_dist = dist

		if dist == lowest_dist:
			lowest_strs.append(other)

		if dist < lowest_dist:
			lowest_dist = dist
			lowest_strs = [other]

	return lowest_strs