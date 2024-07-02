from collections import defaultdict
from itertools import permutations, product

N = 41
Q.<q> = PowerSeriesRing(QQ, default_prec=N)


# Modular forms
def G(k):
    return eisenstein_series_qexp(k, prec=N)

delta = delta_qexp(prec=N)


# Master class for sums of MacMahonesque functions w/ products
# +: polynomial addition
# *: (noncommutative) polynomial multiplication
# @: harmonic product (only defined on sums of letters)
# %: quasi-shuffle product
class _U(dict):
    def __init__(self, dct=None):
        if isinstance(dct, (tuple, list)):
            super().__init__({tuple(dct): 1})
            
        else:
            try:
                int(dct)
                super().__init__({(): dct})
            except:
                super().__init__({k: v for k, v in dict(dct or {}).items() if v})
    
    def __add__(self, other):
        other = _U(other)
        return _U({k: self.get(k, 0) + other.get(k, 0) for k in [*self, *other]})
    
    def __radd__(self, other):
        return _U.__add__(_U(other), self)
    
    def __matmul__(self, other):
        if any(len(k) > 1 for k in [*self, *other]):
            raise NotImplementedError
            
        def c(i, j, m):
            return factorial(i) * factorial(j) / factorial(i + j + 1) if m == i + j else ((-1)^i * binomial(i, m + 1) + (-1)^j * binomial(j, m + 1)) * bernoulli(i + j - m) / (i + j - m)
        
        return sum([self[(i,)] * other[(j,)] * c(i, j, m) * U(m + 1) for i, *_ in self for j, *_ in other for m in range(i + j + 1)], start=_U())
    
    def __rmatmul__(self, other):
        return _U.__matmul__(_U(other), self)
    
    def __mod__(self, other):
        other = _U(other)
        
        if [*self.keys()] == [()]:
            return other * self[()]
        elif [*other.keys()] == [()]:
            return self * other[()]
            
        return sum([self[(x, *w)] * other[(y, *v)] * ((x,) * (_U(w) % (y, *v)) + (y,) * ((x, *w) % _U(v)) + (U(x) @ U(y)) * (_U(w) % _U(v)))
                    for x, *w in self for y, *v in other], start=_U())
    
    def __rmod__(self, other):
        return _U.__mod__(_U(other), self)
    
    def __mul__(self, other):
        other = _U(other)
        return sum([{k + l: v * w} for k, v in self.items() for l, w in other.items()], start=_U())
    
    def __rmul__(self, other):
        return _U.__mul__(_U(other), self)
    
    def __neg__(self):
        return _U({k: -v for k, v in self.items()})
    
    def __sub__(self, other):
        return self + -other
    
    def __rsub__(self, other):
        return _U.__sub__(_U(other), self)
    
    def qexp(self):
        return sum(v * prod(m^i for m, i in zip(reversed(p.to_exp_dict().values()), k)) * q^n
                   for n in range(1, N) for p in Partitions(n) for k, v in self.items() if len(p.to_exp_dict()) == len(k))


# Use these to actually initialize
def U(*a):
    return _U(a)

def Usym(*a):
    return sum(U(*perm) for perm in permutations(a))
        

# For generating all quasimodular forms of weight <= k
# MUCH faster than collecting all the Usym_a's
def Mtwiddly(k, subspace=0):
    bound = k // 2 + 1
    series = defaultdict(Q)
    
    if not subspace:
        series[()] = Q(1)
    
    for n in range(1, N):
        buckets = [[] for _ in range(bound)]
        
        for p in Partitions(n):
            dct = p.to_exp_dict()
            
            if len(dct) < bound:
                buckets[len(dct)].append(dct)
                
        def s(a):
            if len(b := [*{*a}]) == 1:
                return sum(prod(m^b[0] for m in dct.values()) for dct in buckets[len(a)])
            else:
                return sum(prod(m^i for m, i in zip(reversed(dct.values()), perm)) for dct in buckets[len(a)] for perm in permutations(a))
                
        for l in range(bound):
            for v in range(1, k - 2*l, 2):
                a = v, *[1] * l
                if not subspace or len(a) == subspace:
                    series[a[::-1]] += s(a) * q^n
        
        for l in range(bound):
            for v in range(3, k - 2*l - 4, 2):
                for w in range(v, k - v - 2*l - 1, 2):
                    a = w, v, *[1] * l
                    if not subspace or len(a) == subspace:
                        series[a[::-1]] += s(a) * q^n
                 
    return dict(series)

# Matrix corresponding to a space of forms
def mat(space, base=QQ):
    return matrix(base, [series.padded_list(N) for series in space.values()])


# Returns the vector, but also prints it all pretty
def find_linear_combo(space, value):
    A = mat(space).transpose()
    b = vector(QQ, value.padded_list(N))
    
    sol = A.solve_right(b)
    denom = lcm(coeff.denom() for coeff in sol if coeff)
    
    zipped = [tup for tup in zip(space.keys(), sol * denom, space.values()) if tup[1]]
    
    padding = 1 + max(len(str(index)) for index, *_ in zipped)
    for index, coeff, series in zipped:
        print(f"{f'{index}':{padding}}: {int(coeff):+10} * ({series})")
        
    print(padding * " " + f"= {denom:+10} * ({value})\n\n")
    
    for _, coeff, _ in zipped:
        print(f"{abs(int(coeff)):{12 + padding}} =", factor(abs(coeff)))
    
    print()
    return denom, sol


# Finds a "good basis" without doing completely reduced echelon form
def find_good_basis(space, rescale=False):
    m = mat(space)
    keys = [*space.keys()]
    
    offset = 0
    while not m.transpose()[offset]:
        offset += 1
    
    row = 0
    while row < m.nrows():
        for col in range(row):
            if m[row][col + offset]:
                m.add_multiple_of_row(row, col, factor := -m[row][col + offset] / m[col][col + offset])
                
        if m[row][row + offset] < 0:
            m.rescale_row(row, -1)
            
        if (d := (m[row][row + offset] / 1).denom()) > 1:
            m.rescale_row(row, d)
            
        if not any(m[row]):
            keys.pop(row)
            m = m.delete_rows([row])
                
            row -= 1
            
        elif rescale and m[row][row + offset] != 1:
            m.rescale_row(row, 1/gcd(m[row][offset:]))
            
        row += 1

    for key, row in zip(keys, m):
        sol = mat(space).transpose().solve_right(row)
        print(key)
                
    return m


def stagger(k, rescale=False):
    M = matrix(QQ, [[1], *[[0]] * (N - 1)])
    for subspace in range(1, k//2 + 1):
        M = M.augment(find_good_basis(Mtwiddly(k, subspace), rescale=rescale).transpose())
    
    return M.transpose()
	
	
find_good_basis(Mtwiddly(16))
