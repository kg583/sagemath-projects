from itertools import permutations

from sage.rings.integer import Integer

from labelled_digraphs import LabelledDiGraph


class CollatzMapping:
	"""
	Implementation of generalized Collatz mappings
	
	Provides a number of algebraic values and objects useful for analyzing Collatz mappings
	Also supports calling as a function to iterate a mapping
	"""
	def __init__(self, *eqs):
		"""
		Initialize a mapping using a list of Expressions or coefficient tuples, provided in any order
		
		Each Expression must
			1) be linear in the provided variable, which may vary across each input,
			2) occupy a unique equivalence class modulo the degree of the mapping, and
			3) have linear coefficient that is nonzero module the degree.	
		"""
		self.d = Integer(len(eqs))
		self.m, self.r = [None] * self.d, [None] * self.d
		
		for eq in eqs:
			try:
				coeffs = {power: coeff for coeff, power in eq.coefficients()}
			except AttributeError:
				coeffs = {1: eq[0], 0: eq[1]}
			
			try:
				m, r = Integer(coeffs[1]), Integer(-coeffs.get(0, 0))
			except KeyError:
				raise ValueError(f"Component '{eq}' is constant")
				
			try:
				i = r / m % self.d
			except ZeroDivisionError:
				raise ValueError(f"Linear coefficient of '{eq}' is divisible by the degree ({self.d})")
				
			if self.m[i] is not None:
				raise ValueError(f"Component '{eq}' clashes with earlier component on {i} mod {self.d}")
				
			self.m[i], self.r[i] = m, r
			
	def __call__(self, x, i=None):
		"""
		Return the mapping evaluated at x
		
		If an explicit index is given, it will be used instead of x modulo the degree
		"""
		i = i or x % self
		return (self[i][0] * x - self[i][1]) / self.d
		
	def __getitem__(self, i):
		"""
		Return the ith component of the mapping
		"""
		return (self.m[int(i)], self.r[int(i)])
	
	def __pow__(self, n):
		"""
		Return a mapping which is equivalent to iterating the mapping n times
		"""
		if n == 1:
			return self
		
		eqs = []
		for vec in permutations(range(self.d), r=int(n)):
			m, r = 1, 0
			for i in vec:
				m = self[i][0] * m / self.d
				r = (self[i][0] * r + self[i][1]) / self.d
				
			eqs.append((m * self.d ^ n, -r * self.d ^ n))
			
		return CollatzMapping(*eqs)
	
	def __rmod__(self, x):
		"""
		Return x modulo the degree
		"""
		return x % self.d
	
	def is_periodic(self, x, max_iterations=None):
		"""
		Return whether a given input is a periodic point under the mapping
		"""
		return not self.iterate(x, max_iterations)[0]
	
	is_cyclic = is_periodic
	
	def iterate(self, x, max_iterations=None):
		"""
		Iterate the mapping on a given input until a cycle is reached or max_iterations have been computed
		
		Separates output into lists of pre-periodic and periodic iterates
		If max_iterations is reached and a cycle is not detected, all iterates will be assumed pre-periodic
		"""
		iterates = [x]
		while max_iterations is None or len(iterates) < max_iterations:		
			x = self(x)
			if x in iterates:
				break
				
			iterates.append(x)
		
		else:
			return iterates, []
		
		loop = iterates.index(x)
		return iterates[:loop], iterates[loop:]
	
	def vector(self, x, length):
		"""
		Compute the iteration vector of a given input out to some length
		"""
		vector = []
		while len(vector) < length:
			vector.append(x % self.d)
			x = self(x)
			
		return vector
		
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
		Check whether the 0th component is affine
		"""
		return self.r[0] == 0
		
	def is_simple(self):
		"""
		Check whether the 0th component is the identity
		"""
		return self.is_regular() and self.m[0] == 1
	
	def expected_behavior(self):
		"""
		Return the conjectured behavior of this mapping
		"""
		return "convergent" if prod(self.m) < self.d ^ self.d else "divergent"
		
	def mu(self, v):
		"""
		Given an iteration vector, return its μ-vector
		"""
		return [self[vi][0] for vi in v]
		
	def rho(self, v):
		"""
		Given an iteration vector, return its ρ-vector
		"""
		return [self[vi][1] for vi in v]
		
	def pi(self, v):
		"""
		Return the product value of an iteration vector
		"""
		return prod(self[vi][0] for vi in v)
		
	def N(self, v):
		"""
		Return the cycle numerator corresponding to an iteration vector
		"""
		return sum(self[v[j]][1] * self.pi(v[j+1:]) * self.d^j for j in range(len(v)))
		
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
		G = LabelledDiGraph(loops=True, multiedges=True)
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
		G = LabelledDiGraph(loops=True, multiedges=True)
		for x in range(p):
			for y in range(p):
				for i in range(self.d):
					G.add_edge((x, y), self(x, i) % p, self(y, i) % p), i)
					
		return G
