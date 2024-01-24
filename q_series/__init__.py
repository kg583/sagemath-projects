from sage.modular.etaproducts import qexp_eta

PREC = 50

LL.<z> = LaurentSeriesRing(ZZ, default_prec=PREC)
MM.<q> = PowerSeriesRing(LL, default_prec=PREC)


def coeffs(f):
    return dict(zip(f.exponents(), f.coefficients()))


def E(k, **kwargs):
    return eisenstein_series_qexp(k, prec=PREC, K=LL, normalization="constant")(**kwargs)

def eta(**kwargs):
    return qexp_eta(MM, PREC)(**kwargs)

def delta(**kwargs):
    return delta_qexp(prec=PREC, K=LL)(**kwargs)

def j(**kwargs):
    return j_invariant_qexp(prec=PREC, K=LL)(**kwargs)

def T(n, k, eps=None):
    def inner(f, *args, **kwargs):
        return hecke_operator_on_qexp(f, n, k, eps, prec=PREC, *args, **kwargs)
    
    return inner

def J(**kwargs):
    return (E(4)^2 * E(6) / delta / (j() - z))(**kwargs)
