from sage.modular.etaproducts import qexp_eta

N = 30

Q.<q> = PowerSeriesRing(QQ)

P.<E2, E4, E6, eta> = LaurentPolynomialRing(QQ)


def derivative(ex):
    derivatives = [
        (E2^2 - E4) / 12,
        (E2 * E4 - E6) / 3,
        (E2 * E6 - E4^2) / 2,
        -E2 * eta / 24
    ]
    
    gens = P.gens()
    ex = P(ex)
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
    return max(sum(a * b for a, b in zip(exponents, (2, 4, 6, 1/2))) for exponents in P(ex).exponents())


def rc(f, g, nu):
    k, l = weight(f), weight(g)
    return sum((-1)^r * Rational(gamma(k + nu) * gamma(l + nu) / factorial(s := nu - r) / gamma(k + r) / factorial(r) / gamma(l + s)) * D(r, f) * D(s, g) for r in range(nu + 1))
