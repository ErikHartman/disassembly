import networkx as nx
import math

def get_disassembly(objects : dict):
    """
    objects is a dict of {object:copy_number}
    """
    disassembly = 0 
    n_t = sum(objects.values()) # total number of copies in ensemble
    disassembly_indexes = get_disassembly_indexes(objects)

    for seq in objects.keys():
        n_i_minus_one = objects[seq]-1
        if n_i_minus_one > 0:
            disassembly += (math.e**disassembly_indexes[seq])*(n_i_minus_one/n_t)
    return disassembly


def get_disassembly_indexes(objects: dict):

    disassembly_indexes = {}
    for object in objects.keys():
        d_object = get_disassembly_index(objects, object)
        disassembly_indexes[object] = d_object
    return disassembly_indexes


def get_disassembly_index(objects, object):
    subset_G = subset_graph(objects, object)
    fundamental_object = get_largest_object(subset_G)  
    path_probs = get_path_prob(subset_G, fundamental_object, object)
    d = get_path_length(path_probs)
    return d

def get_largest_object(G):
    return sorted([x for x in  G.nodes()], key=len, reverse=True)[0]

def subset_graph(objects, object):
    G = nx.DiGraph()

    for u in objects.keys():
        for v in objects.keys():
            if v in u and u != v:
                G.add_edge(u, v, weight=objects[v])

    subsetted_nodes = [n for n in nx.bfs_successors(G.reverse(), object)][0][1] + [
        object
    ]
    new_dict = {}
    for seq, count in objects.items():
        if seq in subsetted_nodes:
            new_dict[seq] = count
    nodes = new_dict
    G = nx.DiGraph()
    for u in nodes.keys():
        G.add_node(u, count=nodes[u])

    for u in nodes.keys():
        for v in nodes.keys():
            if v in u and u != v:
                G.add_edge(u, v, weight=nodes[v])
    return G


def get_p(G, source, target):
    out_edges = sum([edge[2]["weight"] for edge in G.edges(source, data=True)])
    node_cnt = G.nodes(data=True)[target]["count"]
    p = node_cnt / out_edges  # nodes in current "layer" / all outgoing edges
    return p


def get_paths(G, source, target):
    return list(nx.all_simple_paths(G, source, target))


def get_path_prob(G, source, target):
    paths = get_paths(G, source, target)
    path_probs = []
    for path in paths:
        probs = []
        for i in range(len(path) - 1):
            probs.append(get_p(G, path[i], path[i + 1]))

        prob = 1
        for p in probs:
            prob *= p
        path_probs.append((path, prob))
    return path_probs


def get_path_length(path_probs):
    length = 0
    for path_prob in path_probs:
        l = len(path_prob[0]) - 1
        p = path_prob[1]
        length += l * p
    return length
