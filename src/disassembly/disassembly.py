"""
Functions for estimating disassembly.

"""


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


def get_disassembly_indexes(G: nx.DiGraph):
    # Normalize node outputs to make them into probabilites
    G.remove_edges_from(nx.selfloop_edges(G))
    for node in G.nodes():
        in_edges = G.in_edges(node, data=True)
        total_in = sum([data["weight"] for s, t, data in in_edges])
        for key, target, data in in_edges:
            nx.set_edge_attributes(
                G, {(key, target): {"weight": data["weight"] / total_in}}
            )

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
    print(
        f"\n Averaged DI: {sum(disassembly_indexes.values()) / len(disassembly_indexes.keys()):.2f}"
    )
    return disassembly_indexes


import numpy as np


def get_disassembly_indexes_mc(G: nx.DiGraph, N_particles: int):

    G.remove_edges_from(nx.selfloop_edges(G))
    terminal_nodes = [node for node in G.nodes() if G.out_degree(node) == 0]
    G = G.reverse()

    disassembly_indexes = {sequence: [] for sequence in G.nodes()}

    # Release particle

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
                    disassembly_indexes[node].append(steps - i)
                break
            weights = np.array([weight["weight"] for _, _, weight in out_edges])
            sum_weights = sum(weights)
            weights = [w / sum_weights for w in weights]
            targets = [target for _, target, _ in out_edges]
            next_node = np.random.choice(targets, p=weights)
            sequence = next_node
            steps += 1

    # Make sure all sequences are passed
    non_passed_sequences = [
        seq for seq in disassembly_indexes.keys() if len(disassembly_indexes[seq]) == 0
    ]
    if len(non_passed_sequences) > 0:
        for seq in non_passed_sequences:
            starting_sequence = seq
            sequence = starting_sequence
            path = []
            steps = 0
            while True:
                path.append(sequence)
                out_edges = G.out_edges(sequence, data=True)
                if len(out_edges) == 0:
                    for i, node in enumerate(path):
                        disassembly_indexes[node].append(steps - i)
                    break
                weights = np.array([weight["weight"] for _, _, weight in out_edges])
                sum_weights = sum(weights)
                weights = [w / sum_weights for w in weights]
                targets = [target for _, target, _ in out_edges]
                next_node = np.random.choice(targets, p=weights)
                sequence = next_node
                steps += 1

    mean_disassembly_indexes = {seq: np.mean(disassembly_indexes[seq]) for seq in disassembly_indexes.keys()}
    print(
        f"\n Averaged DI: {sum(mean_disassembly_indexes.values()) / len(mean_disassembly_indexes.keys()):.2f}"
    )

    return mean_disassembly_indexes
