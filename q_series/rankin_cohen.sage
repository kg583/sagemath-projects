# Ring of forms

from sage.modular.etaproducts import qexp_eta

N = 50

Q.<q> = PowerSeriesRing(QQ)

M.<E2, E4, E6, eta> = LaurentPolynomialRing(QQ)
Delta = (E4^3 - E6^2) / 1728

def derivative(ex):
    derivatives = [
        (E2^2 - E4) / 12,
        (E2 * E4 - E6) / 3,
        (E2 * E6 - E4^2) / 2,
        E2 * eta / 24
    ]
    
    gens = M.gens()
    ex = M(ex)
    return sum(coeff * prod(h^i for h, i in zip(gens, exponents) if g != h) * exponents[j] * g^(exponents[j] - 1) * derivatives[j] for coeff, exponents in zip(ex.coefficients(), ex.exponents()) for j, g in enumerate(gens))

def D(k, ex):
    for _ in range(k):
        ex = derivative(ex)
        
    return ex

def Eqexp(k):
    return eisenstein_series_qexp(k, prec=N, normalization="constant")

def E(k):
    if k == 0:
        return 1
    
    if k == 2:
        return E2
    
    gens = [E4^a * E6^((k - 4*a) // 6) for a in range(int(k / 4) + 1) if (k - 4*a) % 6 == 0]
    mat = matrix(QQ, [qexp(series).padded_list(N) for series in gens]).transpose()
    vec = vector(QQ, Eqexp(k).padded_list(N))
    
    return sum(coeff * gen for coeff, gen in zip(gens, mat.solve_right(vec)))

def nat_basis(ex):
    def name_cusp_form(k, i):
        name = r"Delta"
        if i > 1:
            name += f"^{i}"
        if k - 12 * i > 0:
            name += f" * E{k - 12*i}"
            
        return name
    
    k = weight(ex)
    gens = {f"E{k}": E(k)} | {name_cusp_form(k, i): delta^i * E(k - 12*i) for i in range(1, int(k / 12) + (k % 12 != 2))}
    mat = matrix(QQ, [qexp(series).padded_list(N) for series in gens.values()]).transpose()
    vec = vector(QQ, qexp(ex).padded_list(N))
    
    return {name: value for name, value in zip(gens, mat.solve_right(vec))}
    

def weight(ex):
    return max(sum(a * b for a, b in zip(exponents, (2, 4, 6, 1/2))) for exponents in M(ex).exponents())

def g(k):
    return k*(3*k - 1) // 2

def R(n):
    b = int(sqrt(24*n + 1))
    return range((1 - b) // 6 - 2, (1 + b) // 6 + 2)

# The brackets

def rcalg(nu, f=eta, g=eta^-1):
    k, l = weight(f), weight(g)
    return sum((-1)^r * Rational(gamma(k + nu) * gamma(l + nu) / factorial(s := nu - r) / gamma(k + r) / factorial(r) / gamma(l + s)) * D(r, f) * D(s, g) for r in range(nu + 1))

def rcqexp(nu):
    return 1/24^nu * sum(q^n * sum((-1)^r * Rational(gamma(nu + 1/2) * gamma(nu - 1/2) / factorial(s := nu - r) / gamma(r + 1/2) / factorial(r) / gamma(s - 1/2)) * \
                                   sum((-1)^k * (6*k - 1)^(2*r) * (24*n - (6*k - 1)^2)^s * Partitions(n - g(k)).cardinality()
                                       for k in R(n) if g(k) <= n)
                                   for r in range(nu + 1))
                         for n in range(N))

# The formulas

def P(nu, m, n, conjectured=False):
    if conjectured:
        return (-1)^nu * binomial(2 * nu - 2, nu - 2) * sum(m^(nu - j) * (-6 * n)^j * 2 * nu * (2 * nu - 1) / (2 * nu - j) / (2 * nu - j - 1) * binomial(2 * nu - j, j) for j in range(nu + 1))
    else:
        return sum((-1)^r * Rational(gamma(nu + 1/2) * gamma(nu - 1/2) / factorial(s := nu - r) / gamma(r + 1/2) / factorial(r) / gamma(s - 1/2)) * m^r * (24*n - m)^s for r in range(nu + 1))

def recur(nu, n):
    return (sum((-1)^(k-1) * P(nu, (6*k - 1)^2, n) * Partitions(n - g(k)).cardinality() for k in R(n) if k != 0) - (-1)^nu * (4*nu / bernoulli(2*nu)) * binomial(2*nu - 2, nu - 2) * sigma(n, 2*nu - 1)) / P(nu, 1, n)


# woah they're the same!!!!!
nu = 5

print(qexp(rcalg(nu)), rcqexp(nu), sep="\n")
recur(nu, 17), Partitions(17).cardinality()
