class LabeledDiGraph(DiGraph):
	def is_regular(self):
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
		
	def label_subgraph(self, label, vertices=None, edges=None, inplace=False, vertex_property=None, algorithm=None, immutable=None):
		return self.subgraph(vertices=vertices, edges=edges, inplace=inplace, vertex_property=vertex_property, edge_property=lambda e: e[2] == label,
							 algorithm=algorithm, immutable=immutable)
		
	def labeled_incoming_edge_iterator(self, vertex, label):
		for edge in self.incoming_edge_iterator(vertex, labels=True):
			if edge[2] == label:
				yield edge
				
	def labeled_incoming_edges(self, vertex, label):
		return list(self.labeled_incoming_edge_iterator(vertex, label))
				
	def labeled_outgoing_edge_iterator(self, vertex, label):
		for edge in self.outgoing_edge_iterator(vertex, labels=True):
			if edge[2] == label:
				yield edge
				
	def labeled_outgoing_edges(self, vertex, label):
		return list(self.labeled_outgoing_edge_iterator(vertex, label))
		
	def labels_of(self, edges):
		return set([*zip(*edges)][2])
		
	def neighbor_in_by_label_iterator(self, vertex, label):
		for edge in self.labeled_incoming_edge_iterator(vertex, label):
			yield edge[0]
			
	def neighbor_out_by_label_iterator(self, vertex, label):
		for edge in sefl.labeled_outgoing_edge_iterator(vertex, label):
			yield edge[1]
			
	def neighbors_in_by_label(self, vertex, label):
		return list(self.neighbor_in_by_label_iterator(vertex, label))
		
	def neighbors_out_by_label(self, vertex, label):
		return list(self.neighbor_out_by_label_iterator(vertex, label))	
		
	def path_from_edges(self, edges):
		path = []
		while edges[1:]:
			if edges[0][1] != edges[1][0]:
				raise ValueError(f"Edges {edges[0]} and {edges[1]} are not incident") from None
				
			edge = edges.pop(0)
			path.append(self.edge_label(*edge) if len(edge) < 3 else edge[2])
			
		return path + [edges[-1][2]]
		
	def path_from_vertices(self, vertices):
		path = []
		while vertices[1:]:
			try:
				path.append({t: l for o, t, l in self.outgoing_edge_iterator(vertices[0], labels=True)}[vertices[1]])
			
			except KeyError:
				raise ValueError(f"Vertices {vertices[0]} and {vertices[1]} are not adjacent") from None
				
			vertices.pop(0)
			
		return path
		
	def spanning_paths(self, start):
		paths = {start: []}
		tree = self.spanning_tree(start)
	
		for edge in tree:
			paths[edge[1]] = paths[edge[0]] + [edge[2]]
			
		return path
		
	def spanning_tree(self, start):
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
		G = LabeledDiGraph(loops=(self.allows_loops() or other.allows_loops()), multiedges=(self.allows_multiple_edges() or other.allows_multiple_edges()))
		for u, w, l in self.edge_iterator(labels=True):
			for v, x, m in other.edge_iterator(labels=True):
				if l == m:
					P.add_edge((u, v), (w, x), l)
						
		return G

	categorical_product = tensor_product
	kronecker_product = tensor_product
				
	def walk(self, path, start):
		seq = [start]
		while path:
			for edge in self.outgoing_edge_iterator(start, labels=True):
				if edge[2] == path.pop(0):
					seq.append(edge[1])
					start = edge[1]
					break
					
		return seq
		
	def is_walk(self, path, start, end):
		return self.walk(path, start)[-1] == end
		
	def is_cycle(self, path, start):
		return self.is_walk(path, start, start)
