# Ring of forms

from sage.modular.etaproducts import qexp_eta

N = 50

Q.<q> = PowerSeriesRing(QQ)

M.<E2, E4, E6, eta> = LaurentPolynomialRing(QQ)

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


def qexp(ex):
    return ex.subs(E2=eisenstein_series_qexp(2, prec=N, normalization="constant"),
                   E4=eisenstein_series_qexp(4, prec=N, normalization="constant"),
                   E6=eisenstein_series_qexp(6, prec=N, normalization="constant"),
                   eta=qexp_eta(Q, prec=N))
    

def weight(ex):
    return max(sum(a * b for a, b in zip(exponents, (2, 4, 6, 1/2))) for exponents in M(ex).exponents())

def g(k):
    return k*(3*k - 1) // 2

def R(n):
    b = int(sqrt(24*n + 1))
    return range((1 - b) // 6 - 2, (1 + b) // 6 + 2)

# The brackets

def rcalg(f, g, nu):
    k, l = weight(f), weight(g)
    return sum((-1)^r * Rational(gamma(k + nu) * gamma(l + nu) / factorial(s := nu - r) / gamma(k + r) / factorial(r) / gamma(l + s)) * D(r, f) * D(s, g) for r in range(nu + 1))

def rcqexp(nu):
    return 1/24^nu * sum(q^n * sum((-1)^r * Rational(gamma(nu + 1/2) * gamma(nu - 1/2) / factorial(s := nu - r) / gamma(r + 1/2) / factorial(r) / gamma(s - 1/2)) * \
                                   sum((-1)^k * (6*k - 1)^(2*r) * (24*n - (6*k - 1)^2)^s * Partitions(n - g(k)).cardinality()
                                       for k in R(n) if g(k) <= n)
                                   for r in range(nu + 1))
                         for n in range(N))

def P(nu, m, n):
    return sum((-1)^r * Rational(gamma(nu + 1/2) * gamma(nu - 1/2) / factorial(s := nu - r) / gamma(r + 1/2) / factorial(r) / gamma(s - 1/2)) * m^r * (24*n - m)^s for r in range(nu + 1))

# poly == P(nu, *var("m n"))
def poly(nu):
    m, n = var("m n")
    return (-1)^nu * binomial(2 * nu - 2, nu - 2) * sum(m^(nu - j) * (-6 * n)^j * 2 * nu * (2 * nu - 1) / (2 * nu - j) / (2 * nu - j - 1) * binomial(2 * nu - j, j) for j in range(nu + 1)).expand()

def recur(nu, n):
    return (sum((-1)^(k-1) * poly(nu).subs(m=(6*k - 1)^2, n=n) * Partitions(n - g(k)).cardinality() for k in R(n) if k != 0) - (-1)^nu * (4*nu / bernoulli(2*nu)) * binomial(2*nu - 2, nu - 2) * sigma(n, 2*nu - 1)) / poly(nu).subs(m=1,n=n)
