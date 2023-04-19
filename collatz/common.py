from copy import copy
from itertools import product

from sage.graphs.digraph import DiGraph


class LabelledDiGraph(DiGraph):
	"""
	A LabelledDiGraph is an extension of a standard DiGraph which assumes that all edges have some accompanying label.
	Multiple methods which are undefined or not specified for graphs with labels are implemented here.
	
	Initialization can be done with any of the data formats used for a standard DiGraph.
	Label existence is not checked during initialization to reduce overhead and permit piecewise construction of graphs.
	A dedicated backend to watch edge construction would be preferable.
	"""
	def is_cycle(self, path, start, simple=False):
		"""
		Check if a given path connects a vertex to itself.
		
		A path is a sequence of labels of incident edges.
		If the cycle is simple, the starting vertex is visited only at the start and end of the cycle. 
		"""
		return self.is_walk(path, start, start, simple)
	
	def is_regular(self, k=None):
		"""
		Check whether the graph is regular with respect to its labels.
		
		A labelled digraph is regular if and only if every vertex has in-degree and out-degree k with respect to each label.
		"""
		edges_in = {v: {l: 0 for l in self.edge_labels()} for v in self}
		edges_out = copy(edges_in)
		for o, t, l in self.edge_iterator(labels=True):
			edges_in[t][l] += 1
			edges_out[o][l] += 1
		
		if k is None:
			k = edges_in.items()[0][1].items()[0][1]
		
		return all(set(edges_in[v].values()) == {k} for v in self) and all(set(edges_out[v].values()) == {k} for v in self)
	
	def is_vertex_transitive(self, partition=None, verbosity=0, order=False, return_group=True, orbits=False):
		"""
		Check whether the graph is vertex transitive with respect to its labels.
		"""
		return super().is_vertex_transitive(partition=partition, verbosity=verbosity, edge_labels=True, return_group=return_group, orbits=False)
	
	def is_walk(self, path, start, end, simple=False):
		"""
		Check if a given path connects the starting and ending vertices.
		
		A path is a sequence of labels of incident edges.
		If the walk is simple, the ending vertex is not visited before the walk is complete.
		"""
		walk = self.walk(path, start)
		return walk[-1] == end and (end not in walk[1:-1] or not simple)
		
	def label_subgraph(self, label, vertices=None, edges=None, inplace=False, vertex_property=None, edge_property=None, algorithm=None, immutable=None):
		"""
		Return a subgraph of the graph where all edges have a specified label.
		
		Any of the usual arguments can be passed and will be respected, including additional edge properties.
		"""
		return self.subgraph(vertices=vertices, edges=edges, inplace=inplace, vertex_property=vertex_property,
							 edge_property=lambda e: e[2] == label and (edge_property or (lambda e: True))(e),
							 algorithm=algorithm, immutable=immutable)
		
	def labelled_incoming_edge_iterator(self, vertex, label):
		"""
		Return an iterator over all arriving edges from a vertex with a specified label.
		"""
		for edge in self.incoming_edge_iterator(vertex, labels=True):
			if edge[2] == label:
				yield edge
				
	def labelled_incoming_edges(self, vertex, label):
		"""
		Return a list of all arriving edges from a vertex with a specified label.
		"""
		return list(self.labelled_incoming_edge_iterator(vertex, label))
				
	def labelled_outgoing_edge_iterator(self, vertex, label):
		"""
		Return an iterator over all departing edges from a vertex with a specified label.
		"""
		for edge in self.outgoing_edge_iterator(vertex, labels=True):
			if edge[2] == label:
				yield edge
				
	def labelled_outgoing_edges(self, vertex, label):
		"""
		Return a list of all departing edges from a vertex with a specified label.
		"""
		return list(self.labelled_outgoing_edge_iterator(vertex, label))
		
	def labels_of(self, edges):
		"""
		Return a set of labels from a collection of edges, provided as three-tuples of origin, terminus, and label.
		"""
		return set([*zip(*edges)][2])
		
	def neighbor_in_by_label_iterator(self, vertex, label):
		"""
		Return an iterator over the in-neighbors of a vertex with a specified label.
		"""
		for edge in self.labelled_incoming_edge_iterator(vertex, label):
			yield edge[0]
			
	def neighbor_out_by_label_iterator(self, vertex, label):
		"""
		Return an iterator over the out-neighbors of a vertex with a specified label.
		"""
		for edge in self.labelled_outgoing_edge_iterator(vertex, label):
			yield edge[1]
			
	def neighbors_in_by_label(self, vertex, label):
		"""
		Return a list of the in-neighbors of a vertex with a specified label.
		"""
		return list(self.neighbor_in_by_label_iterator(vertex, label))
		
	def neighbors_out_by_label(self, vertex, label):
		"""
		Return a list of the out-neighbors of a vertex with a specified label.
		"""
		return list(self.neighbor_out_by_label_iterator(vertex, label))	
		
	def path_from_edges(self, edges):
		"""
		Return a list of labels given a sequence of edges.
		
		The edges must define a path through the graph, in that each edge's terminus is the next edge's origin.
		Edges are not required to have labels.
		"""
		path = []
		while edges[1:]:
			if edges[0][1] != edges[1][0]:
				raise ValueError(f"Edges {edges[0]} and {edges[1]} are not incident")
				
			edge = edges.pop(0)
			path.append(self.edge_label(*edge) if len(edge) < 3 else edge[2])
			
		return path + [edges[-1][2]]
		
	def path_from_vertices(self, vertices):
		"""
		Return a list of labels given a sequence of vertices.
		
		The vertices must define a path through the graph, in that each vertex has an oriented edge connecting it to the next vertex.
		"""
		path = []
		while vertices[1:]:
			try:
				path.append({t: l for o, t, l in self.outgoing_edge_iterator(vertices[0], labels=True)}[vertices[1]])
			
			except KeyError:
				raise ValueError(f"Vertices {vertices[0]} and {vertices[1]} are not adjacent")
				
			vertices.pop(0)
			
		return path
	
	def permutation_group(self):
		"""
		Compute the permutation group corresponding to the graph interpreted as a Cayley graph with labels as generators.
		
		The graph must be regular of degree one to be a Cayley graph.
		"""
		from sage.groups.perm_gps.permgroup_named import SymmetricGroup
		from sage.groups.perm_gps.permgroup import PermutationGroup

		S = SymmetricGroup(range(self.order()))
		vertex_ordering = {vertex: i for i, vertex in enumerate(self.vertices())}

		labels = self.edge_labels()
		images = {label: [None] * self.order() for label in labels}
		
		for label in labels:
			for vertex in self:
				neighbors_out = self.neighbors_out_by_label(vertex, label)
				if len(neighbors_out) != 1:
					raise ValueError("Graph is not regular of degree one and thus is not a Cayley graph")
					
				images[label][vertex_ordering[vertex]] = vertex_ordering[neighbors_out[0]]
				
		try:
			return PermutationGroup([*map(S, images.values())])
		except ValueError:
			raise ValueError("Graph is not a Cayley graph")
		
	def spanning_tree(self, start):
		"""
		Find a spanning tree of the graph from a given start node as a list of edges.
		"""
		seen = {start}
		stack = [start]
		tree = []
		
		while stack:
			for edge in self.outgoing_edge_iterator(stack.pop(0), labels=True):
				if edge[1] not in seen:
					seen.add(edge[1])
					stack.append(edge[1])
					tree.append(edge)
					
		return tree
		
	def tensor_product(self, other):
		"""
		Return the labeled tensor product (also called the Kronecker or categorical product) of the graph with another labeled digraph.
		
		The tensor product of two labeled digraphs is a graph with vertex set equal to the Cartesian product of the component vertex sets.
		An edge connecting (u, w) and (v, x) exists in the product graph if and only if u ~ w and v ~ x via edges with the same label.
		"""
		G = LabelledDiGraph(loops=(self.allows_loops() or other.allows_loops()), multiedges=(self.allows_multiple_edges() or other.allows_multiple_edges()))
		for u, w, l in self.edge_iterator(labels=True):
			for v, x, m in other.edge_iterator(labels=True):
				if l == m:
					G.add_edge((u, v), (w, x), l)
						
		return G

	categorical_product = tensor_product
	kronecker_product = tensor_product
				
	def walk(self, path, start):
		"""
		Return a list of vertices given a path in the graph and a starting vertex.
		
		The path must be traversible from the given starting vertex.
		If the graph is not regular, the choice of outgoing edge is based on the internal ordering.
		"""
		seq = [start]
		while path:
			label = path.pop(0)
			
			try:
				start = next(self.labelled_outgoing_edge_iterator(start, label))
				seq.append(start)
				
			except StopIteration:
				raise ValueError(f"There is no outgoing edge labeled {label} from vertex {start}")
				
		return seq


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
			3) have linear coefficient that is nonzero modulo the degree.	
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
		i = x % self if i is None else i
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
		for vec in product(range(self.d), repeat=int(n)):
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
					G.add_edge((x, y), (self(x, i) % p, self(y, i) % p), i)
					
		return G
	

def sigma(v, k=1):
	"""
	Rotate a vector k times to the left
	"""
	return v[k:] + v[:k]


C = CollatzMapping(3*x+1, x)
