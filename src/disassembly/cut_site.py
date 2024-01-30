"""


Idea: 

Given a graph, get p1 specificity and endo/exo ratio

If exoprotease, don't add to specificity

if endo: add weight.


"""

import networkx as nx
import numpy as np
from disassembly.util import normalize_dict


def get_p1(graph: nx.DiGraph, N_particles: int = 10000):
    graph.remove_edges_from(nx.selfloop_edges(graph))

    p1_dict = {}  # endo p1
    p1_exo = {}
    n_endo = 0
    n_exo = 0

    for _ in range(N_particles):
        starting_sequence = [n for n in graph if graph.in_degree(n) == 0][0]
        sequence: str = starting_sequence
        while True:
            out_edges = graph.out_edges(sequence, data=True)
            if len(out_edges) == 0:
                break
            weights = np.array([weight["weight"] for _, _, weight in out_edges])
            sum_weights = sum(weights)
            weights = [w / sum_weights for w in weights]
            targets = [target for _, target, _ in out_edges]

            next_node = np.random.choice(targets, p=weights)

            if len(next_node) == len(sequence) - 1:
                n_exo += 1
                if sequence.startswith(next_node):
                    p1 = sequence[len(next_node) - 1]
                    if p1 in p1_exo.keys():
                        p1_exo[p1] = p1_exo[p1] + 1
                    else:
                        p1_exo[p1] = 1
                elif sequence.endswith(next_node):
                    p1 = sequence[-len(next_node) - 1]

                    if p1 in p1_exo.keys():
                        p1_exo[p1] = p1_exo[p1] + 1
                    else:
                        p1_exo[p1] = 1
            else:
                n_endo += 1
                if sequence.startswith(next_node):
                    p1 = sequence[len(next_node) - 1]
                    if p1 in p1_dict.keys():
                        p1_dict[p1] = p1_dict[p1] + 1
                    else:
                        p1_dict[p1] = 1

                elif sequence.endswith(next_node):
                    p1 = sequence[-len(next_node) - 1]

                    if p1 in p1_dict.keys():
                        p1_dict[p1] = p1_dict[p1] + 1
                    else:
                        p1_dict[p1] = 1

                else:
                    p1_left = sequence[sequence.find(next_node) - 1]
                    if p1_left in p1_dict.keys():
                        p1_dict[p1_left] = p1_dict[p1_left] + 1
                    else:
                        p1_dict[p1_left] = 1

                    p1_right = next_node[-1]
                    if p1_right in p1_dict.keys():
                        p1_dict[p1_right] = p1_dict[p1_right] + 1
                    else:
                        p1_dict[p1_right] = 1

            sequence = next_node
    p1_dict = normalize_dict(p1_dict)
    p1_exo = normalize_dict(p1_exo)
    return p1_dict, p1_exo, n_exo, n_endo


"""
TODO: Make function to plot logos!

Use LogoMaker

"""
