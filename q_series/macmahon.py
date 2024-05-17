from collections import defaultdict
from itertools import permutations, product

N = 40
Q.<q> = PowerSeriesRing(QQ, default_prec=N)


# Modular forms
def E(k):
    return eisenstein_series_qexp(k, normalization="constant", prec=N)

delta = delta_qexp(prec=N)
    

# For playing with individual M_a's and U_a's
def M(a, n):
    return sum(prod(m^i for m, i in zip(reversed(p.to_exp_dict().values()), a)) for p in Partitions(n) if len(p.to_exp_dict()) == len(a))

def test(n):
    return (n^2 - 3*n + 2) * M([1], n) - 8 * M([1, 1], n)

def U(*a):
    return sum(M(a, n) * q^n for n in range(1, N))

def Usym(*a):
    return sum(U(*perm) for perm in permutations(a))


# For generating all quasimodular forms of weight <= k
# MUCH faster than collecting all the Usym_a's
def Mtwiddly(k):
    bound = k // 2 + 1
    series = defaultdict(Q)
    series[()] = Q(1)
    
    for n in range(1, N):
        buckets = [[] for _ in range(bound)]
        
        for p in Partitions(n):
            dct = p.to_exp_dict()
            
            if len(dct) < bound:
                buckets[len(dct)].append(dct)
        
        for l in range(1, bound):
            for a in product(range(1, k - l + 1, 2), repeat=l):
                if l + sum(a) <= k:
                    series[tuple(sorted(a))] += sum(prod(m^i for m, i in zip(reversed(dct.values()), a)) for dct in buckets[l]) * q^n
                 
    return series


# Returns the vector, but also prints it all pretty
def find_linear_combo(space, value):
    A = matrix(QQ, [series.padded_list(N) for series in space.values()]).transpose()
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


# Bachmann
def B(*a):
    return prod(1 / factorial(s - 1) for s in a) * U(*[s - 1 for s in reversed(a)])
