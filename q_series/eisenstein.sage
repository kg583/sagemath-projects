# Quasimodular form ring
E.<G2, G4, G6> = QQ[]

# Differential operator ring
P.<D> = QQ[]

# Eisenstein series
def G(k):
    if k <= 6:
        return [G2, G4, G6][k//2 - 1]
    else:
        def d(j):
            return G(2*j + 4) * (2*j + 3) * factorial(j)
        
        n = (k - 8) // 2
        return sum(binomial(n, j) * d(j) * d(n - j) for j in range(n + 1)) * (3*n + 6) / (2*n + 9)

def weight(ex):
    return max(2*i + 4*j + 6*k for i, j, k in ex.exponents())
        
def qexp(ex):
    return ex.subs(G2=eisenstein_series_qexp(2, prec=20), G4=eisenstein_series_qexp(4, prec=20), G6=eisenstein_series_qexp(6, prec=20))

def diff_op(poly, ex):
    def deriv(e2, e4, e6):
        if e4 == e6 == 0:
            return (-2 * G2^2 + 5/6 * G4) * e2 * G2^abs(e2 - 1)
        elif e2 == e6 == 0:
            return (-8 * G2 * G4 + 7/10 * G6) * e4 * G4^abs(e4 - 1)
        elif e2 == e4 == 0:
            return (-12 * G2 * G6 + 400/7 * G4^2) * e6 * G6^abs(e6 - 1)
        else:
            return deriv(e2, 0, 0) * G4^e4 * G6^e6 + G2^e2 * deriv(0, e4, 0) * G6^e6 + G2^e2 * G4^e4 * deriv(0, 0, e6)
    
    derivatives = [ex]
    res = 0
    for i, d in enumerate(list(poly)):
        res += d * derivatives[i]
        derivatives.append(sum(c * deriv(*e) for e, c in zip(derivatives[-1].exponents(), derivatives[-1].coefficients())))
    
    return res
