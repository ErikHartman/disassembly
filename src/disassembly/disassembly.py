import networkx as nx
import math

def get_disassembly(P : dict, disassembly_indexes : dict):
    """
    P is a dict of {object:copy_number}
    disassembly_indexes is a dict of {object:disassembly index}
    """
    disassembly = 0 
    n_t = sum(P.values()) # total number of copies in ensemble

    for sequence in P.keys():
        n_i_minus_one = P[sequence]-1
        if n_i_minus_one > 0:
            disassembly += (math.e**disassembly_indexes[sequence])*(n_i_minus_one/n_t)
    return disassembly


def get_disassembly_indexes(G : nx.DiGraph):
    disassembly_indexes = {sequence:0 for sequence in G.nodes()}
    for sequence, _ in disassembly_indexes.items():
        longest_object = sorted(list(disassembly_indexes.keys()), key=len)[-1]
        if sequence == longest_object:
            continue
        paths_from_sequence_to_longest_object = list(
            nx.all_simple_paths(G.reverse(), sequence, longest_object)
        )

        weighted_length = 0
        for path in paths_from_sequence_to_longest_object:
            prob = 1
            length = len(path) - 1
            for i in range(len(path) - 1):
                prob *= (G.get_edge_data(path[i+1], path[i])["weight"] / sum([data["weight"] for _,_, data in G.in_edges(path[i], data=True)]))

            weighted_length += length*prob
        disassembly_indexes[sequence] = weighted_length
        
    return disassembly_indexes

