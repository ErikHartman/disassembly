import networkx as nx
import numpy as np

"""

TODO: implement initialization 

0 1 is bad

Maybe uniform? Or random uniform.

"""

def normalize_dict(d):
    s = sum(d.values())
    return {k: v / s for k, v in d.items()}

def KL(a, b):
    a = np.asarray(list(a))
    b = np.asarray(list(b))
    a = 1e-8 + a / np.sum(a)
    b = 1e-8 + b / np.sum(b)

    return np.sum(np.where(a != 0, a * np.log(a / b), 0))


def estimate_weights(P: dict, lr : float = 0.1, n_iterations: int = 100):
    keys = list(P.keys())
    values = list(P.values())
    N_T = int(sum(values))
    G = nx.DiGraph()
    G.add_nodes_from([(k, {"layer": len(k)}) for k in keys])
    for key1 in keys:
        for key2 in keys:
            if key2.startswith(key1) or key2.endswith(key1):  # key 1 = ABC, key 2 = ABCD
                if key1 == key2:
                    G.add_edge(key2, key1, weight=1)  # self prob = 1 first
                else:
                    G.add_edge(key2, key1, weight=0)
    generated = {}
    kls = []
    weights = np.zeros((len(G.edges()), n_iterations), dtype=float)
    for i in range(n_iterations):
       
        p_generated = generate_guess(G, keys, N_T)
        generated[i] = p_generated
        kl = KL(P.values(), p_generated.values())
        print(f"\r {i} / {n_iterations} | {kl:.2f}, mean: {np.mean(kls):.2f}", end="")
        G = update_weights(G, kl, P, p_generated, lr)
        kls.append(kl)  
        weights[:, i] = [data["weight"] for _,_, data in G.edges(data=True)]
    
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

def update_weights(G, kl, P, p_generated, lr):
    p_hat = normalize_dict(P)
    p_generated = normalize_dict(p_generated)
    for key in P.keys():
        out_edges = G.out_edges(key, data=True)
        for _, target, weight in out_edges:
            if key == target:
                continue
            source_copy_number = p_hat[key]
            target_copy_number = p_hat[target]
            generated_target_copy_number = p_generated[target]
            diff = (source_copy_number) * (
                (target_copy_number) / (generated_target_copy_number + 0.01)
            )

            add_to_weight = diff * lr * kl

            if np.abs(len(key) - len(target)) == 1:
                add_to_weight *= 2

            new_weight = max(
                0, weight["weight"] + add_to_weight
            )  # weight cannot be less than 0

            nx.set_edge_attributes(G, {(key, target): {"weight": new_weight}})
        total_out = sum([data["weight"] for s, t, data in out_edges])
        for key, target, data in out_edges:
            nx.set_edge_attributes(
                G, {(key, target): {"weight": data["weight"] / total_out}}
            )
    return G