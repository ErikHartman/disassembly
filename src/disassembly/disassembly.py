"""
Functions for estimating disassembly.

"""


import networkx as nx
import math
import numpy as np
import random



def get_disassembly(P: dict, disassembly_indexes: dict):
    """
    P is a dict of {object:copy_number}
    disassembly_indexes is a dict of {object:disassembly index}
    """
    disassembly = 0
    n_t = sum(P.values())  # total number of copies in ensemble

    for sequence in P.keys():
        if P[sequence]  > 0:
            disassembly += math.e ** (disassembly_indexes[sequence]) * (
                P[sequence]  / n_t
            )
    return disassembly


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
                if sum(weights) == 0:
                    print("weights == 0")
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



def get_disassembly_indexes_mc_joint(G: nx.DiGraph, N_particles: int):
    """
    For joint disassembly spaces we need to consider the contribution of a particle
    to the assembly path.
    
    """
    
    G.remove_edges_from(nx.selfloop_edges(G))
    G_orig = G.copy()
    terminal_nodes = [node for node in G.nodes() if G.out_degree(node) == 0]
    G = G.reverse()

    disassembly_indexes = {sequence: [] for sequence in G.nodes()}

    
    starting_sequences = np.random.choice(terminal_nodes, size=N_particles)
    # Release particle
    for particle in range(N_particles):
        starting_sequence = starting_sequences[particle]
        sequence: str = starting_sequence
        path = [] # List of nodes passed: n_4, n_3, n_2, n_1, n_0
        steps = []  # List of weights in path, w_43, w_32, w_21, w_10
        # The di for n_3 (i=1) is w_10 + w_21 + w_32, i.e., sum(steps[i,:])
        while True:
            path.append(sequence)
            out_edges = G.out_edges(sequence, data=True)
            if len(out_edges) == 0:
                for i, node in enumerate(path):
                    disassembly_indexes[node].append(sum(steps[i:]))
                break
            weights = np.array([weight["weight"] for _, _, weight in out_edges])
            sum_weights = sum(weights)
            weights = [w / sum_weights for w in weights]
            targets = [target for _, target, _ in out_edges]
            next_node = random.choices(targets, weights=weights)[0]
            weight_from_next_node_to_sequence = G_orig[next_node][sequence]["weight"] / sum(
                [
                    data["weight"]
                    for _,_, data in G_orig.out_edges(next_node, data=True)
                ]
            )
            if sequence.startswith(next_node) or sequence.endswith(next_node):
                steps.append(
                    2 * weight_from_next_node_to_sequence
                )  # if its not a terminal peptide, two cuts were technically needed
            else:
                steps.append(1 * weight_from_next_node_to_sequence)
            sequence = next_node

    # Make sure all sequences are passed
    non_passed_sequences = [
        seq for seq in disassembly_indexes.keys() if len(disassembly_indexes[seq]) == 0
    ]
    if len(non_passed_sequences) > 0:
        for seq in non_passed_sequences:
            starting_sequence = seq
            sequence = starting_sequence
            path = []
            steps = []
            while True:
                path.append(sequence)
                out_edges = G.out_edges(sequence, data=True)
                if len(out_edges) == 0:
                    for i, node in enumerate(path):
                        disassembly_indexes[node].append(sum(steps[i:]))
                    break
                weights = np.array([weight["weight"] for _, _, weight in out_edges])
                sum_weights = sum(weights)
                weights = [w / sum_weights for w in weights]
                targets = [target for _, target, _ in out_edges]
                next_node = random.choices(targets, weights=weights)[0]
                weight_from_next_node_to_sequence = G_orig[next_node][sequence]["weight"] / sum(
                    [
                        data["weight"]
                        for _,_, data in G_orig.out_edges(next_node, data=True)
                    ]
                )
                if sequence.startswith(next_node) or sequence.endswith(next_node):
                    steps.append(
                        2 * weight_from_next_node_to_sequence
                    )  # if its not a terminal peptide, two cuts were technically needed
                else:
                    steps.append(1 * weight_from_next_node_to_sequence)
                sequence = next_node

    mean_disassembly_indexes = {
        seq: np.mean(disassembly_indexes[seq]) for seq in disassembly_indexes.keys()
    }
    print(
        f"\n Average DI: {sum(mean_disassembly_indexes.values()) / len(mean_disassembly_indexes.keys()):.2f}"
    )

    return mean_disassembly_indexes


"""
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
"""