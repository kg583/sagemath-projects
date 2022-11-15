from sage.rings.integer import Integer

from labeled_digraphs import LabeledDiGraph


class CollatzMapping:
	"""
	Implementation of generalized Collatz mappings
	
	Provides a number of algebraic values and objects useful for analyzing Collatz mappings
	Also supports calling as a function to iterate a mapping
	"""
	def __init__(self, *eqs):
		"""
		Initialize a mapping using a list of Expressions, provided in any order
		
		Each Expression must
			1) be linear in the provided variable, which may vary across each input,
			2) occupy a unique equivalence class modulo the degree of the mapping, and
			3) have linear coefficient that is nonzero module the degree.	
		"""
		self.d = Integer(len(eqs))
		self.m, self.r = [None] * self.d, [None] * self.d
		
		for eq in eqs:
			coeffs = {power: coeff for coeff, power in eq.coefficients()}
			
			try:
				m, r = Integer(coeffs[1]), Integer(-coeffs.get(0, 0))
			except KeyError:
				raise ValueError(f"Component '{eq}' is constant") from None
				
			try:
				i = r / m % self.d
			except ZeroDivisionError:
				raise ValueError(f"Linear coefficient of '{eq}' is divisible by the degree ({self.d})") from None
				
			if self.m[i] is not None:
				raise ValueError(f"Component '{eq}' clashes with earlier component on {i} mod {self.d}") from None
				
			self.m[i], self.r[i] = m, r
			
	def __call__(self, x):
		"""
		Return the mapping evaluated at x
		"""
		i = x % self.d
		return (self.m[i] * x - self.r[i]) / self.d
		
	def __getitem__(self, i):
		"""
		Return the ith component of the mapping
		"""
		return (self.m[i], self.r[i])
		
	def degree(self):
		"""
		Return the degree of the mapping
		"""
		return self.d
				
	def is_expansive(self):
		"""
		Check whether every linear coefficient exceeds the degree
		"""
		return all(self.m[i] > self.d for i in range(self.d))
		
	def is_non_negative(self):
		"""
		Check whether every constant coefficient is non-negative
		"""
		return all(self.r[i] >= 0 for i in range(self.d))
		
	def is_non_positive(self):
		"""
		Check whether every constant coefficient is non-positive
		"""
		return all(self.r[i] <= 0 for i in range(self.d))
		
	def is_reductive(self):
		"""
		Check whether every linear coefficient is bounded by the degree
		"""
		return all(self.m[i] < self.d for i in range(self.d))
		
	def is_regular(self):
		"""
		Check whether the 0th component is non-affine
		"""
		return self.r[0] == 0
		
	def is_simple(self):
		"""
		Check whether the 0th component is the identity
		"""
		return self.is_regular() and self.m[0] == 1
		
	def mu(self, v):
		"""
		Given an iteration vector, return its μ-vector
		"""
		return [self.m[int(vi)] for vi in v]
		
	def rho(self, v):
		"""
		Given an iteration vector, return its ρ-vector
		"""
		return [self.r[int(vi)] for vi in v]
		
	def pi(self, v):
		"""
		Return the product value of an iteration vector
		"""
		return prod(self.m[int(vi)] for vi in v)
		
	def N(self, v):
		"""
		Return the cycle numerator corresponding to an iteration vector
		"""
		return sum(self.r[int(v[j])] * self.pi(v[j+1:]) * self.d^j for j in range(len(v)))
		
	def D(self, v):
		"""
		Return the cycle denominator corresponding to an iteration vector
		"""
		return self.pi(v) - self.d ^ len(v)
		
	def x(self, v):
		"""
		Return the cycle value corresponding to an iteration vector
		"""
		return self.N(v) / self.D(v)
		
	def G(self, v):
		"""
		Return the cancellation of the cycle corresponding to an iteration vector
		"""
		return gcd(self.N(v), self.D(v))
		
	def delta(self, v):
		"""
		Return the reduced cycle denominator corresponding to an iteration vector
		"""
		return self.D(v) / self.G(v)
		
	def graph(self, p, abelian=False):
		"""
		Return the unary iteration graph of the mapping modulo p
		
		The vertex set is Z/pZ
		An edge labeled i connects vertices u and v if the ith component evaluated at u equals v
		If the graph is specified to be abelian, all constant coefficients are ignored
		"""
		G = LabeledDiGraph(loops=True, multiedges=True)
		for x in range(p):
			for i in range(self.d):
				G.add_edge(x, (self.m[i] * x - self.r[i] * (not abelian)) / self.d % p, i)
				
		return G
		
	def full_graph(self, p):
		"""
		Return the binary iteration graph of a mapping modulo p
		
		This graph is the tensor product of the unary iteration graph with itself
		This graph is isomorphic to the tensor product of the unary iteration graph with its abelianization
		"""
		G = LabeledDiGraph(loops=True, multiedges=True)
		for x in range(p):
			for y in range(p):
				for i in range(self.d):
					G.add_edge((x, y), ((self.m[i] * x - self.r[i]) / self.d % p, (self.m[i] * y - self.r[i]) / self.d % p), i)
					
		return G
