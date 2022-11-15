from sage.rings.integer import Integer

from labeled_digraphs import LabeledDiGraph


class CollatzMapping:
	def __init__(self, *eqs):
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
		i = x % self.d
		return (self.m[i] * x - self.r[i]) / self.d
		
	def __getitem__(self, i):
		return (self.m[i], self.r[i])
		
	def degree(self):
		return self.d
				
	def is_expansive(self):
		return all(self.m[i] > self.d for i in range(self.d))
		
	def is_non_negative(self):
		return all(self.r[i] >= 0 for i in range(self.d))
		
	def is_non_positive(self):
		return all(self.r[i] <= 0 for i in range(self.d))
		
	def is_reductive(self):
		return all(self.m[i] < self.d for i in range(self.d))
		
	def is_regular(self):
		return self.r[0] == 0
		
	def is_simple(self):
		return self.is_regular and self.m[0] == 1
		
	def mu(self, v):
		return [self.m[int(vi)] for vi in v]
		
	def rho(self, v):
		return [self.r[int(vi)] for vi in v]
		
	def pi(self, v):
		return prod(self.m[int(vi)] for vi in v)
		
	def N(self, v):
		return sum(self.r[int(v[j])] * self.pi(v[j+1:]) * self.d^j for j in range(len(v)))
		
	def D(self, v):
		return self.pi(v) - self.d ^ len(v)
		
	def x(self, v):
		return self.N(v) / self.D(v)
		
	def G(self, v):
		return gcd(self.N(v), self.D(v))
		
	def delta(self, v):
		return self.D(v) / self.G(v)
		
	def graph(self, p, abelian=False):
		G = LabeledDiGraph(loops=True, multiedges=True)
		for x in range(p):
			for i in range(self.d):
				G.add_edge(x, (self.m[i] * x - self.r[i] * (not abelian)) / self.d % p, i)
				
		return G
		
	def full_graph(self, p):
		G = LabeledDiGraph(loops=True, multiedges=True)
		for x in range(p):
			for y in range(p):
				for i in range(self.d):
					G.add_edge((x, y), ((self.m[i] * x - self.r[i]) / self.d % p, self.m[i] * y / self.d % p), i)
					
		return G
