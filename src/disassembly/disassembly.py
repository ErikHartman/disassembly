import networkx as nx
import math

"""

TODO: Updaate MC approximation to start in terminal nodes

"""


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
                if prob < 0.001:
                    prob = 0
                    break
                prob *= G.get_edge_data(path[i + 1], path[i])["weight"] / sum(
                    [data["weight"] for _, _, data in G.in_edges(path[i], data=True)]
                )

            weighted_length += length * prob
        disassembly_indexes[sequence] = weighted_length
    print(f"\n Averaged DI: {sum(disassembly_indexes.values()) / len(disassembly_indexes.keys()):.2f}")
    return disassembly_indexes

import numpy as np
def get_disassembly_indexes_mc(G: nx.DiGraph, N_particles : int):
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
    terminal_nodes = [node for node in G.nodes() if G.out_degree(node) == 0]
    print(terminal_nodes)
    G = G.reverse()
    
    disassembly_indexes = {sequence: [] for sequence in G.nodes()}
    
    # release particle
    
    for _ in range(N_particles):
        starting_sequence = np.random.choice(terminal_nodes)
        sequence = starting_sequence
        path = []
        steps = 0
        while True:
            path.append(sequence)
            out_edges = G.out_edges(sequence, data=True)
            if len(out_edges) == 0:
                for i, node in enumerate(path):
                    disassembly_indexes[node].append(steps-i)
                break
            weights = np.array([weight["weight"] for _, _, weight in out_edges])
            sum_weights = sum(weights)
            weights = [w / sum_weights for w in weights]
            targets = [target for _, target, _ in out_edges]
            next_node = np.random.choice(targets, p=weights)
            sequence = next_node
            steps += 1

        
    return disassembly_indexes

