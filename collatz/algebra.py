import re

M = MatrixSpace(QQ, 2)

A, B = M([[3/2, 1/2], [0, 1]]), M([[1/2, 0], [0, 1]])
a, b = A^-1, B^-1

def sign(x):
    return (x > 0) - (x < 0)

def mat(ex):
    return eval("*".join(ex))

def reduce(ex):
    reduced = ex.replace("Aa", "").replace("aA", "").replace("Bb", "").replace("bB", "")
    return reduce(reduced) if reduced != ex else ex

def inv(ex):
    return ex.swapcase()[::-1]

def rot(ex, i=1):
    return ex[i:] + ex[:i]

def power(ex, x):
    return ex * int(x) if x >= 0 else inv(ex) * -int(x)

def fixed(m):
    return list(m)[0][1] / (1 - list(m)[0][0])

def fixer(n, i, c="bABa"):
    p, q, r = euler_phi((n/1).denom()), log(list(mat(c))[0][1].denom(), 2), sign(list(mat(c))[0][1])
    a, b = list([B, A][i]^p)[0]
    s = (n*(1 - a) - b) * 2^p

    return reduce("".join(map(power, *zip(*[("B", p - q), (c, r * s), ("B", q - p), ("A", p)] if i else [("B", p - q), (c, r * s), ("B", q)]))))

def comm(n):
    x, y = fixer(n, 1), fixer(n, 0)
    *map(print, [x, y]),
    return reduce(x + y + inv(x) + inv(y))
