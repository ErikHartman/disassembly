import networkx as nx
import math


def get_disassembly(P: dict, disassembly_indexes: dict):
    """
    P is a dict of {object:copy_number}
    disassembly_indexes is a dict of {object:disassembly index}
    """
    disassembly = 0
    n_t = sum(P.values())  # total number of copies in ensemble

    for sequence in P.keys():
        n_i_minus_one = P[sequence] - 1
        if n_i_minus_one > 0:
            disassembly += (math.e ** disassembly_indexes[sequence]) * (
                n_i_minus_one / n_t
            )
    return disassembly


def get_disassembly_indexes(G: nx.DiGraph, min_weight : float = 0.01):
    """
    The input is a weighted graph and a min_weight threshold.

    Edges with weight below min weight are removed.

    sequence -w-> subsequence

    When we caluclate disassembly index we want to know:

    Given that we have a sequence:
        What is the estimated length to get there

    """
    # Normalize node outputs to make them into probabilites
    G.remove_edges_from(nx.selfloop_edges(G))
    for node in G.nodes():
        in_edges = G.in_edges(node, data=True)
        total_in = sum([data["weight"] for s, t, data in in_edges])
        for key, target, data in in_edges:
            nx.set_edge_attributes(
                G, {(key, target): {"weight": data["weight"] / total_in}}
            )

    # Remove very improbable edges
    no_0_G = nx.DiGraph()
    for source, target, data in G.edges(data=True):
        if data["weight"] > min_weight:
            no_0_G.add_edge(source, target, **data)
    print(f"Removing low prob. (<{min_weight}) edges. Before: {G.number_of_edges()} | after: {no_0_G.number_of_edges()}")
    G = no_0_G
    
    disassembly_indexes = {sequence: 0 for sequence in G.nodes()}
    for i, (sequence, _) in enumerate(disassembly_indexes.items()):
        longest_object = sorted(list(disassembly_indexes.keys()), key=len)[-1]
        if sequence == longest_object:
            continue
        paths_from_sequence_to_longest_object = list(
            nx.all_simple_paths(G.reverse(), sequence, longest_object)
        )  # This is slow as hell for many paths

        print(
            f"\r {i} / {len(disassembly_indexes.keys())} | checking {len(paths_from_sequence_to_longest_object)} paths",
            end="",
        )

        weighted_length = 0
        for path in paths_from_sequence_to_longest_object:
            prob = 1
            length = len(path) - 1
            for i in range(len(path) - 1):
                if prob < 0.0001:
                    prob = 0
                    break
                prob *= G.get_edge_data(path[i + 1], path[i])["weight"] / sum(
                    [data["weight"] for _, _, data in G.in_edges(path[i], data=True)]
                )

            weighted_length += length * prob
        disassembly_indexes[sequence] = weighted_length

    return disassembly_indexes
