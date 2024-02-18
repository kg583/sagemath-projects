# Cache to accelerate calculation of traces
from functools import cache

# Default power series ring and precision
# May need to increase for large weight traces
R.<q> = PowerSeriesRing(QQ)
PREC = 200

# Gets the coefficient of q^n in f
def coeff(f, n):
    return dict(zip(f.exponents(), f.coefficients())).get(n, 0)

# Common sqrt bound
def bound(n):
    return int(sqrt(n)) + 1

# The Eisenstein series normalized with constant Fourier coefficient 1
def Eis(k):
    return eisenstein_series_qexp(k, normalization='constant', prec=PREC)

# Delta and the j-function as q-expansions
Delta = (Eis(4)^3 - Eis(6)^2) / 1728
j = j_invariant_qexp(prec=PREC)

# Extraneous Sp factors as defined by Kaneko-Zagier
def delta(k):
    return k % 3

def epsilon(k):
    return (k % 4) / 2

# Space dimensions
def dimM(k):
    return dimS(k) + 1

def dimS(k):
    return max(k // 12 - (k % 12 == 2), 0)

# The divisor polynomial of f of weight k, modulo p if given
def divisor_polynomial(f, k, p=None):
    return poly_in_j(f / (Delta ^ dimS(k) * Eis(4) ^ delta(k) * Eis(6) ^ epsilon(k)), p).factor()

# Expresses f as a polynomial in j, modulo p if given
def poly_in_j(P, p=None):
    e = 1 - P.exponents()[0]
    coeffs = [0] * e
    while e:
        e -= 1
        coeffs[e] = coeff(P, -e)
        P -= coeffs[e] * j^e

    if p is not None:
        T.<x> = PolynomialRing(GF(p))
    else:
        T.<x> = PolynomialRing(QQ)

    return T(coeffs)

# Hurwitz-Kroncecker class numbers stored in a dictionary
# Requires the b-table of A259825 on the OEIS downloaded locally
H = {}
with open('b259825.txt') as file:
    for line in file:
        index, value = map(int, line.split())
        H.update({index: value / 12})

# Second-order recurrence component for the Eichler-Selberg trace formula
# Returns a functor which is called with the desired weight
@cache
def P(t, n):
    def adjust_index(k):
        return BinaryRecurrenceSequence(t, -n)(k - 1)
    return adjust_index

# Divisor sum in the Eichler-Selberg trace formula
def l(k, n):
    return 1/2 * sum(min(d, n // d)^(k - 1) for d in range(1, n + 1) if n % d == 0)

# Left-hand sum in the Eichler-Selberg trace formula
def s(k, n):
    bound = int(sqrt(n)) + 1
    return sum(P(t, n)(k) * H.get(4 * n - t^2, 0) for t in range(-2 * bound, 2 * bound))

# The trace of T(n) for weight k eigenforms using the Eichler-Selberg trace formula
def trace(k, n):
    return -1/2 * s(k, n) - l(k, n)

# T_k as a q-expansion
def trace_form(k):
    return R([trace(k, n) for n in range(1, PREC)]) * q

# T_k-hat as a q-expansion
def modified_trace_form(k, p):
    return R([trace(k, n) % p if n % p else 0 for n in range(1, PREC)]) * q


# Example 5.2 (other examples are computed similarly)
f = 3 * trace_form(28)
g = trace_form(76)
n = (76 - 28) / (5 - 1)
S5 = 1

f % 5
g % 5

divisor_polynomial(g, 76, 5) / divisor_polynomial(g, 28, 5)
x^(n // 3) * S5
