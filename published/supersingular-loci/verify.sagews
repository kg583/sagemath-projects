# Cache to accelerate calculation of traces
from functools import cache

# Space dimensions
def dimM(k):
    return dimS(k) + 1

def dimS(k):
    return max(k // 12 - (k % 12 == 2), 0)

# Maximum power of q required to verify the equivalence of two modular forms of weight k
# Given by Sturm's theorem
def max_coeff(k):
    return k // 12 + 1

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


# Verification of Theorem 2.1 (ii) and (iii) for p = 5
def prover5():
    max_ = 136

    for k in range(12, 136, 2):
        tk = [trace(k, n) % 5 for n in range(2, max_coeff(max_) + 1)]
        for c in range(1, 4):
            # (iii)
            if dimS(k) == dimS(l := k + 4 * c):
                assert(all(tk[n - 2] == trace(l, n) % 5 for n in range(2, max_coeff(l) + 1)))

            # (ii)
            if k - 24 in [0, 4, 6, 8, 10, 14]:
                l = k + 24 * c
                assert(all((c + 1) * tk[n - 2] % 5 == trace(l, n) % 5 for n in range(2, max_coeff(l) + 1)))

        # k % 20 or print(k)

    return True


# Verification of Theorem 2.1 (ii) and (iii) for p = 7
def prover7():
    max_ = 352

    for k in (i for i in range(12, max_, 2)):
        tk = [trace(k, n) % 7 for n in range(2, max_coeff(max_) + 1)]
        for c in range(1, 6):
            # (iii)
            if dimS(k) == dimS(l := k + 6 * c):
                assert(all(tk[n - 2] == trace(l, n) % 7 for n in range(2, max_coeff(l) + 1)))

            # (ii)
            if k - 48 in [0, 4, 6, 8, 10, 14]:
                l = k + 48 * c
                assert(all((c + 1) * tk[n - 2] % 7 == trace(l, n) % 7 for n in range(2, max_coeff(l) + 1)))

        # k % 20 or print(k)

    return True


# Verification of Theorem 2.1 (ii) and (iii) for p = 11
def prover11():
    max_ = 1336

    for k in (i for i in range(12, max_, 2)):
        tk = [trace(k, n) % 11 for n in range(2, max_coeff(max_) + 1)]
        for c in range(1, 10):
            # (iii)
            if dimS(k) == dimS(l := k + 10 * c):
                assert(all(tk[n - 2] == trace(l, n) % 11 for n in range(2, max_coeff(l) + 1)))

            # (ii)
            if k - 120 in [0, 4, 6, 8, 10, 14]:
                l = k + 120 * c
                assert(all((c + 1) * tk[n - 2] % 11 == trace(l, n) % 11 for n in range(2, max_coeff(l) + 1)))

        # k % 20 or print(k)

    return True
