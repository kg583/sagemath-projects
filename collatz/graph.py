from labelled_digraphs import LabelledDigraph


def collatz_graph(limit=5):
	"""
	Draws a graph of the Collatz function and its inverse iterated from the points 0 and -1
	The vertices are arranged in concentric circles based on their total stopping time
	"""
	f = [lambda x: x/2, lambda x: (3*x+1)/2]
	g = [lambda x: 2*x, lambda x: (2*x-1)/3]

	graph = LabelledDiGraph([(0, 0, "+0"), (-1, -1, "+1")], loops=True, multiedges=True)
	levels = {0: 0, -1: 0}
	
	
	for n in range(limit):
		for v in graph.vertices():
			for i in (0, 1):
				ed = (v, f[i](v), f"+{i}")
				if ed not in graph.edges():
					graph.add_edge(ed)
					
				if ed[1] not in levels:
					levels[ed[1]] = n + 1
					
				ed = (g[i](v), v, f"+{i}")
				if ed not in graph.edges():
					graph.add_edge(ed)
					
				if ed[0] not in levels:
					levels[ed[0]] = n + 1
			
	levels = [[k for k in levels if levels[k] == n] for n in range(limit + 1)]
	pos = {0: (-1, 0), -1: (1, 0)}
	for n in range(1, limit + 1):
		m = len(levels[n])
		i = 0
		for v in levels[n]:		
			ex = (n + 1) * e ^ (2 * pi * I * i / m)
			pos[v] = (real(ex), imag(ex))
			i += 1
				
	graph.show(pos=pos, color_by_label=True, figsize=[70, 70])
	return graph
