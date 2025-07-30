### tools to paste together (G, a)s and (H, b)s without gluing A \times B edges

from itertools import product

def get_all_pastes(graphs, min_degree=0):
    """
    Given a list of graphs, returns a list of all possible pastes
    (G, a) and (H, b) where a and b are vertices in G and H respectively.
    Each paste is represented as a 2-tuple: (graph, d), where d is the size
    of the neighborhood of a and b.
    """
    results = []
    for G1, G2 in product(graphs, repeat=2):
        a_candidates = [orbit[0] for orbit in G1.automorphism_group().orbits()]
        b_candidates = [orbit[0] for orbit in G2.automorphism_group().orbits()]
        for a, b in product(a_candidates, b_candidates):
            attempt = try_paste_together(G1, a, G2, b, min_degree)
            if attempt is None:
                continue
            pastes, d = attempt
            for paste in pastes:
                results.append((paste, d))
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
    ## relabel so that the order is a, b, nbhd(a) = nbhd(b) = K, A, B
    G_relabel_map = {}
    H_relabel_map = {}
    G_relabel_map[a] = 0
    H_relabel_map[b] = 1
    for i, v in enumerate(a_nbhd.vertices(), start=2):
        G_relabel_map[v] = i
    for i, v in enumerate(A, start=2 + d):
        G_relabel_map[v] = i
    for i, v in enumerate(B, start=2 + len(A) + d):
        H_relabel_map[v] = i

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