import networkx as nx
import numpy as np
import pandas as pd
from disassembly.util import KL, normalize_dict, get_trend


def estimate_weights(
    true_dict: dict,
    lr: float = 0.1,
    n_iterations: int = 100,
):
    true_dict = normalize_dict(true_dict)
    true_dict_vals = list(true_dict.values())
    keys = list(true_dict.keys())

    graph = nx.DiGraph()
    graph.add_nodes_from([(k, {"layer": len(k)}) for k in keys])
    for key1 in keys:
        for key2 in keys:
            if (key1 in key2) and (key1 != key2):
                graph.add_edge(key2, key1, weight=np.random.uniform(0, 1))

    for node in graph.nodes():
        out_edges = graph.out_edges(node, data=True)
        total_out = sum(
            [data["weight"] for _, _, data in graph.out_edges(node, data=True)]
        )
        for key, target, data in out_edges:
            nx.set_edge_attributes(
                graph, {(key, target): {"weight": 0.5 * data["weight"] / total_out}}
            )

    generated = {}
    kls = []
    weights = np.zeros((len(graph.edges()), n_iterations), dtype=float)
    lr_cooldown = 100
    for i in range(n_iterations):
        lr_cooldown -= 1
        guess, df = generate_guess(graph, keys)
        generated[i] = guess
        kl = KL(true_dict.values(), guess.values())
        trend = get_trend(kls[-50:])
        print(
            f"\r {i} / {n_iterations} | {kl:.2f}, mean: {np.mean(kls[-25:]):.2f} | {trend} | nz: {len(np.nonzero(weights[:, i-1])[0])}",
            end="",
        )

        dp_dw = compute_dp_dw(graph, keys, df)
        dL_dp = compute_dL_dp(true_dict_vals, list(guess.values()))
        grad = compute_dL_dw(dL_dp, dp_dw)

        if (
            lr_cooldown <= 0
            and (trend == "Increasing" or trend == "Stochastic")
            and lr > 0.00001
        ):
            lr_cooldown = 50
            lr = lr / 2
            print(f"\nLearning rate decreased to {lr}")

        weights[:, i] = [data["weight"] for _, _, data in graph.edges(data=True)]
        graph = update_weights(graph, grad, lr)
        kls.append(kl)

        if np.mean(kls[-50:]) < 0.02:
            break

    return graph, kls, generated, weights


def generate_guess(graph: nx.DiGraph, keys):
    """
    Outputs a tuple of the distribution for the longest node and a matrix
    """
    longest_key = sorted(keys, key=len)[-1]
    p_generated = {}
    terminal_nodes = [node for node in graph.nodes() if graph.out_degree(node) == 0]

    for node in terminal_nodes:  # one hot terminal nodes
        oh_node = create_one_hot(keys, node)
        p_generated[node] = oh_node

    out_edges = {
        source: [target for _, target in graph.out_edges(source) if source != target]
        for source in graph.nodes()
    }

    while len(p_generated.keys()) < len(keys):
        solvables = get_solvable(out_edges, p_generated)
        for solvable in solvables:
            p_generated[solvable] = np.zeros(len(keys))

            for source, target in graph.out_edges(solvable):
                p_target = p_generated[target]
                w_source_target = graph[source][target]["weight"]
                p_generated[source] += w_source_target * p_target

            w_source_target = 1 - sum(
                [data["weight"] for _, _, data in graph.out_edges(source, data=True)]
            )
            p_target = create_one_hot(keys, source)
            p_generated[source] += w_source_target * p_target

    guess = {keys[i]: p_generated[longest_key][i] for i in range(len(keys))}
    return guess, pd.DataFrame(p_generated, index=keys)


def compute_dp_dw(graph, keys, df):
    """
    dP / dw
    Change of P based on w
    Sx1 vector
    """
    longest_key = sorted(keys, key=len)[-1]
    prob_traversed = {key: 0 for key in keys}
    prob_traversed[longest_key] = 1

    for sequence, n in prob_traversed.items():
        out_edges = [
            (source, target, data)
            for source, target, data in graph.out_edges(sequence, data=True)
        ]
        weights = np.array([weight["weight"] for _, _, weight in out_edges])
        edges_to = [edge_to for _, edge_to, _ in out_edges]
        for w, e in zip(weights, edges_to):
            prob_traversed[e] += w * n

    # TODO: spara matrisen
    dp_dw = {}

    for key in keys:
        out_edges = graph.out_edges(key)
        for (
            source,
            target,
        ) in out_edges:  # P(longest to source) * (P(target) - onehot(source))
            dp_dw[(source, target)] = prob_traversed[source] * (
                df[target].values - create_one_hot(keys, source)
            )

    return dp_dw


def compute_dL_dp(true, guess):
    """
    q = real
    p = guess
    """
    return -np.array(true) / (np.array(guess) + 1e-8)


def compute_dL_dw(dL_dp, dp_dw):
    """
    Gradient
    Sx1 * 1xS = 1x1
    """
    dL_dw = {}
    for edge, val in dp_dw.items():
        dL_dw[edge] = np.sum(val * dL_dp)
    return dL_dw


def get_l1(graph):
    return sum([abs(data["weight"]) for _, _, data in graph.edges(data=True)])


def get_l2(graph):
    return sum([data["weight"] ** 2 for _, _, data in graph.edges(data=True)])


def get_elastic_net(graph, lambda_1, lambda_2):
    return sum(
        [
            (lambda_1 * abs(data["weight"]))  # L1
            + (lambda_2 * data["weight"] ** 2)  # L2
            for _, _, data in graph.edges(data=True)
        ]
    )


def update_weights(graph, grad, lr):
    """
    Updates weights
    Makes sure that sum_new_weight < 1
    """
    for source in graph.nodes():
        sum_old_weight = sum(
            [data["weight"] for _, _, data in graph.out_edges(source, data=True)]
        )
        sum_diffs = 0
        diffs = {}
        for source, target in graph.out_edges(source):
            old_weight = graph[source][target]["weight"]
            new_weight = max(0, old_weight - lr * grad[(source, target)])
            diff = new_weight - old_weight
            sum_diffs += diff
            diffs[target] = diff
        k = 1
        while (sum_old_weight + k * sum_diffs) > 1:
            k = k / 2

        for source, target in graph.out_edges(source):
            nx.set_edge_attributes(
                graph,
                {
                    (source, target): {
                        "weight": max(
                            0, graph[source][target]["weight"] + diffs[target] * k
                        )
                    }
                },
            )

    return graph


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
