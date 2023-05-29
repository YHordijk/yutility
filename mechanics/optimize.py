from yutility import numdiff


def minimize(pos, function, h=1e-5, eps=1e-3, update_strength=1e-1, max_iter=1000):
    last_E = 0
    for i in range(max_iter):
        E = function(pos)
        F = -numdiff.gradient(function, pos, h=h)

        if abs(E - last_E) < eps:
            break

        pos += F * update_strength
        last_E = E

    print(pos)
    return pos
