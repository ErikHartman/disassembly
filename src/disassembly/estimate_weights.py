"""
Functions for estimating weights in disassembly graph.

Main function is esimate_weights().

"""

import networkx as nx
import numpy as np
from disassembly.simulate_proteolysis import amino_acids
from disassembly.util import KL, normalize_dict, get_trend


def estimate_weights(
    P: dict,
    lr: float = 0.1,
    n_iterations: int = 100,
    early_stop: bool = True,
):
    keys = list(P.keys())
    P = normalize_dict(P)
    G = nx.DiGraph()
    G.add_nodes_from([(k, {"layer": len(k)}) for k in keys])
    for key1 in keys:
        for key2 in keys:
            if key1 == key2:
                G.add_edge(key2, key1, weight=np.random.uniform(0.5, 1))
            elif key1 in key2:  # key 1 = ABC, key 2 = ABCD
                G.add_edge(key2, key1, weight=np.random.uniform(0, 1))

    for node in G.nodes():
        out_edges = G.out_edges(node, data=True)

        total_out = sum([data["weight"] for s, t, data in out_edges])
        for key, target, data in out_edges:
            nx.set_edge_attributes(
                G, {(key, target): {"weight": data["weight"] / total_out}}
            )

    generated = {}
    kls = []
    weights = np.zeros((len(G.edges()), n_iterations), dtype=float)
    lr_cooldown = 100
    for i in range(n_iterations):
        lr_cooldown -= 1
        p_generated = generate_guess(G, keys)
        generated[i] = p_generated
        kl = KL(P.values(), p_generated.values())
        trend = get_trend(kls[-50:])

        print(
            f"\r {i} / {n_iterations} | {kl:.2f}, mean: {np.mean(kls[-25:]):.2f} | {trend} | nz: {len(np.nonzero(weights[:, i-1])[0])}",
            end="",
            flush=True,
        )

        if (
            lr_cooldown <= 0
            and (trend == "Increasing" or trend == "Plateau")
            and lr > 0.0001
        ):
            lr_cooldown = 75
            lr = lr / 2
            print(f"\nLearning rate decreased to {lr}")

        G = update_weights(G, P, p_generated, lr)

        kls.append(kl)
        weights[:, i] = [data["weight"] for _, _, data in G.edges(data=True)]

        if i > 200 and get_trend(kls[-100:]) == "Plateau" and early_stop:
            break

        if np.mean(kls[-100:]) < 0.01 and early_stop:
            break

    return G, kls, generated, weights


def generate_guess(G: nx.DiGraph, keys):
    """
    Outputs a tuple of the distribution for the longest node and a matrix
    """
    longest_key = sorted(keys, key=len)[-1]
    p_generated = {}
    terminal_nodes = [node for node in G.nodes() if G.out_degree(node) == 1]

    for node in terminal_nodes:  # one hot terminal nodes
        oh_node = create_one_hot(keys, node)
        p_generated[node] = oh_node

    out_edges = {
        source: [target for _, target in G.out_edges(source) if source != target]
        for source in G.nodes()
    }

    while len(p_generated.keys()) < len(keys):
        solvables = get_solvable(out_edges, p_generated)
        for solvable in solvables:
            p_generated[solvable] = np.zeros(len(keys))
            for source, target in G.out_edges(solvable):
                if source == target:
                    p_target = create_one_hot(keys, source)
                else:
                    p_target = p_generated[target]
                w_source_target = G[source][target]["weight"]
                p_generated[source] += w_source_target * p_target
    guess = {keys[i]: p_generated[longest_key][i] for i in range(len(keys))}
    return guess


def create_one_hot(keys, key):
    one_hot = np.zeros(len(keys))
    one_hot[keys.index(key)] = 1
    return one_hot


def get_solvable(out_edges, p_generated):
    solvable = []
    for source, targets in out_edges.items():
        if (
            set(targets).issubset(set((p_generated.keys())))
            and source not in p_generated.keys()
        ):
            solvable.append(source)
    return solvable


def update_weights(G, P, p_generated, lr):
    P = normalize_dict(P)  # observerade fördelningen
    p_generated = normalize_dict(p_generated)  # gissningen
    for key in P.keys():
        out_edges = G.out_edges(key, data=True)
        for _, target, data in out_edges:
            source_copy_number = P[key]
            target_copy_number = P[target]
            generated_target_copy_number = p_generated[target]
            diff = target_copy_number - generated_target_copy_number

            add_to_weight = diff * lr

            new_weight = data["weight"] + add_to_weight
            new_weight = max(0, new_weight)

            nx.set_edge_attributes(G, {(key, target): {"weight": new_weight}})

        # Normalize
        out_edges = G.out_edges(key, data=True)
        total_out = sum([data["weight"] for s, t, data in out_edges])

        if total_out != 0:
            for key, target, data in out_edges:
                nx.set_edge_attributes(
                    G, {(key, target): {"weight": data["weight"] / total_out}}
                )
        else:
            weights = np.random.uniform(0.01, 1, len(out_edges))
            total_out = sum(weights)
            for edge, weight in zip(out_edges, weights):
                key, target, _ = edge
                nx.set_edge_attributes(
                    G, {(key, target): {"weight": weight / total_out}}
                )

    return G




"""

Old way of generating guess (stochastic)


def generate_guess(G, keys, N_T):
    longest_key = sorted(keys, key=len)[-1]
    p_generated = {key: 0 for key in keys}
    for _ in range(N_T):
        sequence = longest_key
        next_edge = None
        while True:
            out_edges = G.out_edges(sequence, data=True)
            weights = np.array([weight["weight"] for _, _, weight in out_edges])
            edges_to = [edge_to for _, edge_to, _ in out_edges]
            next_edge = np.random.choice(edges_to, p=weights)
            if next_edge == sequence:
                break
            sequence = next_edge

        p_generated[sequence] += 1
    return p_generated
"""
