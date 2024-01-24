from sage.modular.etaproducts import qexp_eta

PREC = 50

LL.<z> = LaurentSeriesRing(QQ, default_prec=PREC)
MM.<q> = PowerSeriesRing(LL, default_prec=PREC)


def coeffs(f):
    return dict(zip(f.exponents(), f.coefficients()))


def E(k):
    return eisenstein_series_qexp(k, prec=PREC, K=LL, normalization="constant")

def T(n, k, eps=None):
    def inner(f, *args, **kwargs):
        return hecke_operator_on_qexp(f, n, k, eps, prec=PREC, *args, **kwargs)
    
    return inner

eta = qexp_eta(MM, PREC)
delta = delta_qexp(prec=PREC, K=LL)
j = j_invariant_qexp(prec=PREC, K=LL)

J = E(4)^2 * E(6) / delta / (j - z)
