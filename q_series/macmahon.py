from collections import defaultdict
from itertools import permutations, product

N = 36
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

# Bachmann's notation
def B(*a):
    return prod(1 / factorial(s - 1) for s in a) * U(*[s - 1 for s in reversed(a)])


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
                
        for l in range(bound):
            for v in range(1, k - 2*l, 2):
                a = v, *[1] * l
                series[a[::-1]] += sum(prod(m^i for m, i in zip(reversed(dct.values()), perm)) for dct in buckets[len(a)] for perm in permutations(a)) * q^n
        
        for l in range(bound):
            for v in range(3, k - 2*l - 4, 2):
                for w in range(v, k - v - 2*l - 1, 2):
                    a = w, v, *[1] * l
                    series[a[::-1]] += sum(prod(m^i for m, i in zip(reversed(dct.values()), perm)) for dct in buckets[len(a)] for perm in permutations(a)) * q^n
                 
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
def find_good_basis(space):
    m = mat(space)
    keys = [*space.keys()]
    
    row = 0
    dim = m.nrows()
    while row < dim:
        for col in range(row):
            if m[row][col]:
                m.add_multiple_of_row(row, col, factor := -m[row][col] / m[col][col])
                
        if m[row][row] < 0:
            m.rescale_row(row, -1)
            
        if (d := (m[row][row] / 1).denom()) > 1:
            m.rescale_row(row, d)
            
        if not any(m[row]):
            keys.pop(row)
                  
            for k in range(row, m.nrows() - 1):
                m.swap_rows(k, k + 1)
                
            row -= 1
            dim -= 1
            
        row += 1

    for key, row in zip(keys, m[:dim]):
        sol = mat(space).transpose().solve_right(row)
        print(key, sol)
                
    return m[:dim]
