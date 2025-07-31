### tools to paste together (G, a)s and (H, b)s without gluing A \times B edges

from itertools import combinations, combinations_with_replacement, product

from sage.all import Graph

def read_graphs_from_file(filename):
    """
    Reads a list of graphs from a file. Each graph is expected to be in
    graph6 format, one graph per line.
    Returns a list of graphs.
    """
    graphs = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                graph = Graph(line)
                graphs.append(graph)
    return graphs

def write_graphs_to_file(graphs, filename):
    """
    Writes a list of graphs to a file in graph6 format.
    Each graph is written on a new line.
    """
    with open(filename, 'w') as f:
        for graph in graphs:
            f.write(graph.graph6_string() + '\n')

def unique_graphs(graph_list):
    seen = set()
    unique = []
    for G in graph_list:
        canon = G.canonical_label()
        # convert to string or another hashable type
        canon_str = canon.graph6_string()
        if canon_str not in seen:
            seen.add(canon_str)
            unique.append(G)
    return unique

def get_all_pastes(graphs, min_degree=0):
    """
    Given a list of graphs, returns a list of all possible pastes
    (G, a) and (H, b) where a and b are vertices in G and H respectively.
    Returns a dict where keys are d = |nbhd(a)| = |nbhd(b)|, and values are lists of pasted graphs.
    """
    results = {}
    graphs_to_point_candidates = {}
    for i, G in enumerate(graphs):
        a_candidates = [orbit[0] for orbit in G.automorphism_group().orbits()]
        graphs_to_point_candidates[i] = a_candidates
    
    # G1 and G2 are different graphs
    for i, j in combinations(range(len(graphs)), 2):
        G1 = graphs[i]
        G2 = graphs[j]
        a_candidates = graphs_to_point_candidates[i]
        b_candidates = graphs_to_point_candidates[j]
        for a, b in product(a_candidates, b_candidates):
            attempt = try_paste_together(G1, a, G2, b, min_degree)
            if attempt is None:
                continue
            pastes, d = attempt
            if d not in results:
                results[d] = []
            results[d].extend(pastes)
    # G1 and G2 are the same graph
    for i, G in enumerate(graphs):
        a_candidates = graphs_to_point_candidates[i]
        for a, b in combinations_with_replacement(a_candidates, 2):
            attempt = try_paste_together(G, a, G, b, min_degree)
            if attempt is None:
                continue
            pastes, d = attempt
            if d not in results:
                results[d] = []
            results[d].extend(pastes)

    return results

def try_paste_together(G, a, H, b, min_degree=0):
    """
    Try to paste together (G, a) and (H, b) without gluing A x B edges.
    Returns a 2-tuple:
    1. list of pasted graphs if successful, otherwise returns None.
    2. d = |nbhd(a)| = |nbhd(b)|, the size of the neighborhood of a and b.
    """
    a_nbhd = G.subgraph(G.neighbors(a))
    b_nbhd = H.subgraph(H.neighbors(b))

    if len(a_nbhd) != len(b_nbhd):
        return None
    if len(a_nbhd) < min_degree:
        return None

    is_iso, iso_map = b_nbhd.is_isomorphic(a_nbhd, certificate=True)
    if not is_iso:
        return None

    results = []
    
    # Create a new graph that combines G and H
    d = len(a_nbhd)
    n = len(G) + len(H) - d
    A = set(G.vertices()) - set(G.neighbors(a)) - {a}
    B = set(H.vertices()) - set(H.neighbors(b)) - {b}
    ## relabel so that the order is a, b, A, B, nbhd(a) = nbhd(b) = K
    G_relabel_map = {}
    H_relabel_map = {}
    G_relabel_map[a] = 0
    H_relabel_map[b] = 1
    for i, v in enumerate(A, start=2):
        G_relabel_map[v] = i
    for i, v in enumerate(B, start=2 + len(A)):
        H_relabel_map[v] = i
    for i, v in enumerate(a_nbhd.vertices(), start=2 + len(A) + len(B)):
        G_relabel_map[v] = i

    a_nbhd_automorphisms = a_nbhd.automorphism_group()
    for aut in a_nbhd_automorphisms:
        ## start F as a copy of  G
        F = G.copy()
        F.relabel(G_relabel_map)

        ## add b
        F.add_vertex(H_relabel_map[b])
        for u in [G_relabel_map[u] for u in G]:
            F.add_edge(u, H_relabel_map[b])

        ## add all vertices in B = H - nbhd(b)
        for v in B:
            label = H_relabel_map[v]
            F.add_vertex(label)
            F.add_edge(G_relabel_map[a], label)
            for u in H.neighbors(v):
                if u in iso_map:
                    K_nbhr = G_relabel_map[aut(iso_map[u])]
                    F.add_edge(label, K_nbhr)
                elif H_relabel_map[u] in F:
                    F.add_edge(label, H_relabel_map[u])

        results.append(F)

    return (results, d)