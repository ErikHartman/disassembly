"""
Functions for estimating weights in disassembly graph.

Main function is esimate_weights().

"""

import networkx as nx
import numpy as np
from disassembly.simulate_proteolysis import amino_acids


def normalize_dict(d):
    s = sum(d.values())
    return {k: v / s for k, v in d.items()}


def KL(a, b):
    a = np.asarray(list(a))
    b = np.asarray(list(b))
    a = 1e-8 + a / np.sum(a)
    b = 1e-8 + b / np.sum(b)

    return np.sum(np.where(a != 0, a * np.log(a / b), 0))


def estimate_weights(
    P: dict,
    meta_enzyme: dict = {
        amino_acid: 1 / len(amino_acids) for amino_acid in amino_acids.values()
    },
    exo_mult_factor: float = 10,
    lr: float = 0.1,
    n_iterations: int = 100,
    N_T: int = 1000,
):
    keys = list(P.keys())
    values = list(P.values())
    P = {k: N_T * v / sum(values) for k, v in zip(keys, values)}
    G = nx.DiGraph()
    G.add_nodes_from([(k, {"layer": len(k)}) for k in keys])
    for key1 in keys:
        for key2 in keys:
            if key1 == key2:
                G.add_edge(key2, key1, weight=1)
            elif key2.startswith(key1) or key2.endswith(
                key1
            ):  # key 1 = ABC, key 2 = ABCD
                G.add_edge(key2, key1, weight=0)

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
    lr_cooldown = 50
    for i in range(n_iterations):
        lr_cooldown -= 1
        p_generated = generate_guess(G, keys, N_T)
        generated[i] = p_generated
        kl = KL(P.values(), p_generated.values())
        trend = get_trend(kls[-50:])
        print(
            f"\r {i} / {n_iterations} | {kl:.2f}, mean: {np.mean(kls[-25:]):.2f} | {trend}",
            end="",
        )
        if lr_cooldown <= 0 and trend == "Increasing" and lr > 0.0001:
            lr_cooldown = 100
            lr = lr / 2
            print(f"\nLearning rate decreased to {lr}")
        G = update_weights(
            G, kl, P, p_generated, lr, meta_enzyme, exo_mult_factor
        )
        kls.append(kl)
        weights[:, i] = [data["weight"] for _, _, data in G.edges(data=True)]

    return G, kls, generated, weights


def generate_guess(G, keys, N_T):
    longest_key = sorted(keys, key=len)[-1]
    p_generated = {key: 0 for key in keys}
    """
    Starting from the longest key
    get random choice from out edges (weighted)
    if random choice == longest key
        p_generated[longest key] += 1
    else:
        get new random choice from out edges of random choice.
        if new random choice == random choice
            p_generated[random_choice] += 1
        repeat
    """
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


def update_weights(G, kl, P, p_generated, lr, meta_enzyme, exo_mult_factor):
    p_hat = normalize_dict(P)
    p_generated = normalize_dict(p_generated)
    for key in P.keys():
        out_edges = G.out_edges(key, data=True)
        for _, target, data in out_edges:
            #source_copy_number = p_hat[key]
            target_copy_number = p_hat[target]
            generated_target_copy_number = p_generated[target]
            diff = target_copy_number - generated_target_copy_number

            add_to_weight = diff * lr

            if len(key) - len(target) == 1:
                if diff > 0:
                    add_to_weight *= exo_mult_factor # idea to increase weight more if exoprotease
                mult_to_new_weight = 1

            elif key == target:
                mult_to_new_weight = 1

            elif key.startswith(target):
                p1 = key[len(target) - 1]
                mult_to_new_weight = 1 if meta_enzyme[p1] > 0 else 0

            else:
                p1 = key[-len(target) - 1]
                mult_to_new_weight = 1 if meta_enzyme[p1] > 0 else 0

            new_weight = max(0, data["weight"] + add_to_weight) * mult_to_new_weight

            nx.set_edge_attributes(G, {(key, target): {"weight": new_weight}})

        # Normalize
        out_edges = G.out_edges(key, data=True)
        total_out = sum([data["weight"] for s, t, data in out_edges])
        for key, target, data in out_edges:
            nx.set_edge_attributes(
                G, {(key, target): {"weight": data["weight"] / total_out}}
            )

    return G


def get_trend(data):
    diff = np.diff(data)
    mean_diff = np.mean(diff)
    std_diff = np.std(diff)
    trend_threshold = 1e-3
    stochastic_threshold = 1
    if mean_diff > trend_threshold:
        return "Increasing"
    elif mean_diff < -trend_threshold:
        return "Decreasing"
    elif std_diff > stochastic_threshold:
        return "Stochastic"
    else:
        return "Plateau"
