# Master file for all crank calculations and searches

import itertools
import datetime

MAX_N = 75
Z.<z> = LaurentSeriesRing(ZZ)
Q.<q> = PowerSeriesRing(Z, default_prec=MAX_N)


# C(z, tau) and C_k(a_1, a_2, ...; z, tau)
# Note the inner function! Use like C(0, z, q)() or C(k, z, q)(a_1, a_2, ...), where the number of a_i's must match (k + (k % 2)) // 2
def C(k=0, var_zeta=z, var_q=q):
	def series(*args):
		if k == 0:
			return prod((1 - var_q^n) / (1 - var_zeta * var_q^n) / (1 - var_zeta^-1 * var_q^n) for n in range(1, MAX_N + 1))
		else:
			return C(0, 1, var_q)()^((k - (k % 2)) / 2) * prod(C(0, var_zeta^arg, var_q)() for arg in args)
			
	return series
			
		
# Infinite families for Conjecture 4.1		
def A(k, var_zeta=z, var_q=q):
	return C(k, var_zeta, var_q)(*range(2, (k + (k % 2)) / 2 + 2))
	

def B(k, var_zeta=z, var_q=q):
	return C(k, var_zeta, var_q)(*range(3, (k + (k % 2)) / 2 + 3))
	

# Creates a map from the power of z to its coefficient in the input polynomial
# Sage does not have a built-in for this for Laurent polynomials!
def laurent_dict(poly):
	coeffs = str(poly).replace(" - ", " + -").split(" + ")
	dictionary = {}
	for coeff in coeffs:
		dictionary.update({int((coeff + ("^1" if "z" in coeff else "*z^0")).split("^")[1]): int(("1*" + coeff + ("*z^0" if "z" not in coeff else "")).split("*")[-2])})
		
	return dictionary
	

# Determines if a Laurent polynomial is unimodal, assuming it is symmetric
def is_unimodal(poly):
	coeffs = laurent_dict(poly)
	
	for i in range(min(coeffs.keys()), 0):
		c, d = coeffs.get(i, 0), coeffs.get(i + 1, 0)
		if c > d:
			return False
			
	return True
	

# Finds the minimal value such that the given crank is unimodal in all larger [q^n]
def min_unimodal(crank):
	min_n = 1
	for n in range(1, MAX_N):
		if not is_unimodal(crank.coefficients()[n].laurent_polynomial()):
			min_n = n
				
	return min_n


# Formatting of min_unimodal return for printing
def format_min(min_n, file=None):
	if min_n == MAX_N - 1:
		print("Not unimodal", file=file)
	else:
		print(f"UNIMODAL for n>{min_n}", file=file)


# Searches all possible cranks of a given k
def find_cranks(k, file=None):
	start = datetime.datetime.now()
	for alist in itertools.combinations(range(1, k + 1), (k + 1) // 2):
		print(alist, end=": ", file=file)
		format_min(min_unimodal(C(k)(*alist)), file=file)
			
	print(f"Runtime: {datetime.datetime.now() - start}")


# Verifies Conjecture 4.1 for k up to max_k
def verifyAB(max_k, file=None):
	start = datetime.datetime.now()
	for k in range(3, max_k):
		print(f"A{k}", end=": ", file=file)
		format_min(min_unimodal(A(k)), file=file)
		print(f"B{k}", end=": ", file=file)
		format_min(min_unimodal(B(k)), file=file)
		
	print(f"Runtime: {datetime.datetime.now() - start}")
