from copy import copy

from sage.graphs.digraph import DiGraph


class LabeledDiGraph(DiGraph):
	"""
	A LabeledDiGraph is an extension of a standard DiGraph which assumes that all edges have some accompanying label.
	Multiple methods which are undefined or not specified for graphs with labels are implemented here.
	
	Initialization can be done with any of the data formats used for a standard DiGraph.
	Label existence is not checked during initialization to reduce overhead and permit piecewise construction of graphs.
	"""
	def is_regular(self):
		"""
		Check whether the graph is regular with respect to its labels.
		
		A labeled digraph is regular if and only if every vertex has the same in-degree and out-degree for each label.
		"""
		edges_in, edges_out = set(), set()
		for o, t, l in self.edge_iterator(labels=True):
			if (t, l) in edges_in:
				return False
			else:
				edges_in.add((t, l))
				
			if (o, l) in edges_out:
				return False
			else:
				edges_out.add((o, l))
					
		return True
		
	def label_subgraph(self, label, vertices=None, edges=None, inplace=False, vertex_property=None, edge_property=None, algorithm=None, immutable=None):
		"""
		Return a subgraph of the graph where all edges have a specified label.
		
		Any of the usual arguments can be passed and will be respected, including additional edge properties.
		"""
		return self.subgraph(vertices=vertices, edges=edges, inplace=inplace, vertex_property=vertex_property,
							 edge_property=lambda e: e[2] == label and edge_property(e),
							 algorithm=algorithm, immutable=immutable)
		
	def labeled_incoming_edge_iterator(self, vertex, label):
		"""
		Return an iterator over all arriving edges from a vertex with a specified label.
		"""
		for edge in self.incoming_edge_iterator(vertex, labels=True):
			if edge[2] == label:
				yield edge
				
	def labeled_incoming_edges(self, vertex, label):
		"""
		Return a list of all arriving edges from a vertex with a specified label.
		"""
		return list(self.labeled_incoming_edge_iterator(vertex, label))
				
	def labeled_outgoing_edge_iterator(self, vertex, label):
		"""
		Return an iterator over all departing edges from a vertex with a specified label.
		"""
		for edge in self.outgoing_edge_iterator(vertex, labels=True):
			if edge[2] == label:
				yield edge
				
	def labeled_outgoing_edges(self, vertex, label):
		"""
		Return a list of all departing edges from a vertex with a specified label.
		"""
		return list(self.labeled_outgoing_edge_iterator(vertex, label))
		
	def labels_of(self, edges):
		"""
		Return a set of labels from a collection of edges, provided as three-tuples of origin, terminus, and label.
		"""
		return set([*zip(*edges)][2])
		
	def neighbor_in_by_label_iterator(self, vertex, label):
		"""
		Return an iterator over the in-neighbors of a vertex with a specified label.
		"""
		for edge in self.labeled_incoming_edge_iterator(vertex, label):
			yield edge[0]
			
	def neighbor_out_by_label_iterator(self, vertex, label):
		"""
		Return an iterator over the out-neighbors of a vertex with a specified label.
		"""
		for edge in self.labeled_outgoing_edge_iterator(vertex, label):
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
		"""
		path = []
		while edges[1:]:
			if edges[0][1] != edges[1][0]:
				raise ValueError(f"Edges {edges[0]} and {edges[1]} are not incident") from None
				
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
				raise ValueError(f"Vertices {vertices[0]} and {vertices[1]} are not adjacent") from None
				
			vertices.pop(0)
			
		return path
	
	def permutation_group(self):
		"""
		Compute the permutation group corresponding to the graph interpreted as a Cayley graph with labels as generators.
		
		The graph must be regular of degree one to be a Cayley graph.
		"""
		S = SymmetricGroup(self.vertices())
		labels = self.edge_labels()
		images = {label: {vertex: None for vertex n self.vertices()} for label in labels}
		
		for label in labels:
			for vertex in self:
				neighbors_out = self.neighbors_out_by_label(vertex, label)
				if len(neighbors_out) != 1:
					raise ValueError("Graph is not regular of degree one and thus is not a Cayley graph") from None
					
				images[label][vertex] = neighbors_out[0]
				
		try:
			return PermutationGroup(*[S(images[label].values()) for label in images])
		except ValueError:
			raise ValueError("Graph is not a Cayley graph") from None
		
	def spanning_paths(self, start):
		"""
		Return a dictionary of paths connecting a given starting vertex to any other.
		
		Each entry is index by the terminal vertex, with paths stored as lists of labels.
		"""
		paths = {start: []}
		for edge in self.spanning_tree(start):
			paths[edge[1]] = paths[edge[0]] + [edge[2]]
			
		return path
		
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
	
	def stallings_fold(self, inplace=False):
		"""
		Fold the graph in the manner of Stallings, either inplace or as a copy.
		
		A graph is folded if no two incoming or outgoing edges from a vertex share the same label.
		The language of a graph composed labeled cycles starting at a central vertex is equivalent to the language of its folding.
		"""
		G = self if inplace else copy(self)
		
		folded = set()
		while len(folded) < G.order():
			vertices = iter(G)
			while (vertex := next(vertices)) in folded:
				pass
			
			incoming, outgoing = {}, {}
			for edge in G.incoming_edge_iterator(vertex):
				incoming[edge[2]] = incoming.get(edge[2], []) + [edge[0]]
				
			for edge in G.outgoing_edge_iterator(vertex):
				outgoing[edge[2]] = outgoing.get(edge[2], []) + [edge[1]]
				
			G.merge_vertices(incoming)
			G.merge_vertices(outgoing)
			folded.add(vertex)
			
		if not inplace:
			return G
		
	def tensor_product(self, other):
		"""
		Return the labeled tensor product (also called the Kronecker or categorical product) of the graph with another labeled digraph.
		
		The tensor product of two labeled digraphs is a graph with vertex set equal to the Cartesian product of the component vertex sets.
		An edge connecting (u, w) and (v, x) exists in the product graph if and only if u ~ w and v ~ x via edges with the same label.
		"""
		G = LabeledDiGraph(loops=(self.allows_loops() or other.allows_loops()), multiedges=(self.allows_multiple_edges() or other.allows_multiple_edges()))
		for u, w, l in self.edge_iterator(labels=True):
			for v, x, m in other.edge_iterator(labels=True):
				if l == m:
					P.add_edge((u, v), (w, x), l)
						
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
				start = next(self.labeled_outgoing_edge_iterator(start, label))
				seq.append(start)
				
			except StopIteration:
				raise ValueError(f"There is no outgoing edge labeled {label} from vertex {start}") from None
				
		return seq
		
	def is_walk(self, path, start, end, simple=False):
		"""
		Check if a given path connects the starting and ending vertices.
		
		If the walk is simple, the ending vertex is not visited before the walk is complete.
		"""
		walk = self.walk(path, start)
		return walk[-1] == end and (end not in walk[1:-1] or not simple)
		
	def is_cycle(self, path, start, simple=False):
		"""
		Check if a given path connects a vertex to itself.
		
		If the cycle is simple, the starting vertex is visited only at the start and end of the cycle. 
		"""
		return self.is_walk(path, start, start, simple)
