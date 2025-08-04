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

def write_graphs_to_file(graphs, filename, append=False):
    """
    Writes a list of graphs to a file in graph6 format.
    Each graph is written on a new line.
    """
    mode = 'a' if append else 'w'
    with open(filename, mode) as f:
        for graph in graphs:
            f.write(graph.graph6_string() + '\n')

def get_all_pastes(graphs, min_degree=0):
    """
    Given a list of graphs, returns a list of all possible pastes
    (G, a) and (H, b) where a and b are vertices in G and H respectively.
    Returns a dict where keys are d = |K|, and values are lists of pasted graphs.
    """
    Ks_to_graphs = {}
    for G in graphs:
        for a in [orbit[0] for orbit in G.automorphism_group().orbits()]:
            K = G.subgraph(G.neighbors(a))
            K_can = K.canonical_label()
            K_can_str = K_can.graph6_string()
            K_automorphisms = K.automorphism_group()
            if K_can_str not in Ks_to_graphs:
                Ks_to_graphs[K_can_str] = []
            Ks_to_graphs[K_can_str].append((a, G, K_automorphisms))

    results = {}
    print(f"Found {len(Ks_to_graphs)} unique K graphs.")

    for K_can_str, tuples in Ks_to_graphs.items():
        for (a, G, a_nbhd_autos), (b, H, b_nbhd_autos) in product(tuples, repeat=2):
            pastes, d = try_paste_together(G, a, a_nbhd_autos, H, b, min_degree)
            if pastes is not None:
                if d not in results:
                    results[d] = []
                for paste in pastes:
                    results[d].append(paste)
            num_pastes = sum(len(pastes) for pastes in results.values())
            if num_pastes >=  10000:
                for d in results:
                    write_graphs_to_file(results[d], f'pasted_graphs_d{d}.txt', append=True)
                print(f"Processed {num_pastes} pastes for K = \"{K_can_str}\".")
                results.clear()

    return results

def try_paste_together(G, a, a_nbhd_autos, H, b, min_degree=0):
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

    for aut in a_nbhd_autos:
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