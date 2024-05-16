from sage.modular.etaproducts import qexp_eta

PREC = 50

Z = ZZ
# or
Z.<z> = LaurentSeriesRing(ZZ, default_prec=PREC)

Q.<q> = PowerSeriesRing(Z, default_prec=PREC)


# Utilities
def coeffs(f):
    return dict(zip(f.exponents(), f.coefficients()))


# Modular forms
def E(k):
    return eisenstein_series_qexp(k, prec=PREC, K=Z, normalization="constant")

def delta():
    return delta_qexp(prec=PREC, K=Z)

def j():
    return j_invariant_qexp(prec=PREC, K=Z)

def eta():
    return qexp_eta(Q, PREC)

def T(n, k, eps=None):
    def inner(f, *args, **kwargs):
        return hecke_operator_on_qexp(f, n, k, eps, prec=PREC, *args, **kwargs)
    
    return inner
