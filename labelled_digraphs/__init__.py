from copy import copy

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
		return super().is_vertex_transitive(partition=partition, verbosity=verbosity, edge_labels=True, return_group=return_group, orbits=orbits)
	
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
			
	def relabel(self, label_map):
		"""
		Relabel the graph using a given map between old and new labels
		"""
		for edge in self.edges():
			try:
				self.set_edge_label(edge[0], edge[1], label_map.get(edge[2], edge[2]))
			except AttributeError:
				self.set_edge_label(edge[0], edge[1], label_map[edge[2]])
		
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
