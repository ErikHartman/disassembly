import networkx as nx
import numpy as np

amino_acids = {
    "VAL": "V",
    "ILE": "I",
    "LEU": "L",
    "GLU": "E",
    "GLN": "Q",
    "ASP": "D",
    "ASN": "N",
    "HIS": "H",
    "TRP": "W",
    "PHE": "F",
    "TYR": "Y",
    "ARG": "R",
    "LYS": "K",
    "SER": "S",
    "THR": "T",
    "MET": "M",
    "ALA": "A",
    "GLY": "G",
    "PRO": "P",
    "CYS": "C",
}


class enzyme:
    """
    Enzyme with p1 specificity
    """

    def __init__(self, specificity: dict, name: str):
        self.specificity = specificity
        assert sum(specificity.values()) - 1 < 0.001
        self.aas = list(specificity.keys())
        self.name = name


class enzyme_set:
    def __init__(self, enzymes: list, activities: list):
        total_activity = sum(activities)
        activities = [a / total_activity for a in activities]
        assert len(enzymes) == len(activities)
        self.enzyme_dict = {e: a for e, a in zip(enzymes, activities)}

    def get_enzyme(self, name):
        for enzyme in self.enzyme_dict.keys():
            if enzyme.name == name:
                return enzyme
        print("enzyme not in dict")
        return None


def simulate_proteolysis(
    starting_sequence: str,
    enzymes: enzyme_set = enzyme_set(
        [
            enzyme({"K": 1}, "protease_iv"),
            enzyme({"K": 0.5, "R": 0.5}, "trypsin"),
            enzyme({"V": 0.5, "I": 0.25, "A": 0.15, "T": 0.1}, "elne"),
        ],
        [3, 2, 3],
    ),
    n_start: int = 100,
    n_iterations: int = 100,
) -> (dict, nx.DiGraph):
    sequence_dict = {starting_sequence: n_start}  # a dict with sequence:copy_number
    sequence_graph = (
        nx.DiGraph()
    )  # a weighted graph, telling us how degradation has occured
    sequence_graph.add_node(starting_sequence, len=len(starting_sequence))

    for _ in range(n_iterations):
        exo_or_endo = np.random.choice(["exo", "endo"], p=[0.5, 0.5])
        if exo_or_endo == "exo":
            sequences_longer_than_2 = [s for s in sequence_dict.keys() if len(s) > 2]

            sequence_to_chew = np.random.choice(
                sequences_longer_than_2,
                p=[
                    sequence_dict[v]
                    / sum([sequence_dict[s] for s in sequences_longer_than_2])
                    for v in sequences_longer_than_2
                ],
            )

            n_or_c_term = np.random.choice(["n", "c"], p=[0.5, 0.5])
            if n_or_c_term == "n":
                new_sequence = sequence_to_chew[1:]
            else:
                new_sequence = sequence_to_chew[:-1]
            sequence_dict = update_sequence_dict(
                sequence_dict, sequence_to_chew, new_sequence, endo_or_exo="exo"
            )
            sequence_graph = update_sequence_graph(
                sequence_graph, sequence_to_chew, new_sequence
            )

        elif exo_or_endo == "endo":
            cut_probabilities = get_cut_probability(sequence_dict, enzymes)
            enzyme_that_will_cut = np.random.choice(
                list(cut_probabilities.keys()), p=list(cut_probabilities.values())
            )
            e: enzyme = enzymes.get_enzyme(enzyme_that_will_cut)

            # Here I say that p(sequence) is prop. to n_cutsites * copy_number. This might not be true. Maybe only copy_number.
            sequence_frequencies = {}
            for sequence in sequence_dict.keys():
                n_cut_sites_in_sequence = 0
                for aminoacid in e.specificity.keys():
                    n_cut_sites_in_sequence += (
                        len(find_aminoacids_in_sequence(sequence, aminoacid))
                        * e.specificity[aminoacid]
                    )
                sequence_frequencies[sequence] = (
                    n_cut_sites_in_sequence * sequence_dict[sequence]
                )

            sequence_to_cut = np.random.choice(
                list(sequence_frequencies.keys()),
                p=[
                    p / sum(sequence_frequencies.values())
                    for p in sequence_frequencies.values()
                ],
            )

            index_to_cut = {}
            for aminoacid in e.specificity.keys():
                indices_for_aminoacid = find_aminoacids_in_sequence(
                    sequence_to_cut, aminoacid
                )
                for index in indices_for_aminoacid:
                    if index != len(sequence_to_cut):
                        index_to_cut[index] = e.specificity[aminoacid]

            cutting_index = np.random.choice(
                list(index_to_cut.keys()),
                p=[p / sum(index_to_cut.values()) for p in index_to_cut.values()],
            )
            left = sequence_to_cut[: cutting_index + 1]
            right = sequence_to_cut[cutting_index + 1 :]
            sequence_dict = update_sequence_dict(
                sequence_dict, sequence_to_cut, left, endo_or_exo="endo"
            )
            sequence_dict = update_sequence_dict(
                sequence_dict, sequence_to_cut, right, endo_or_exo="endo"
            )
            sequence_graph = update_sequence_graph(
                sequence_graph, sequence_to_cut, left
            )
            sequence_graph = update_sequence_graph(
                sequence_graph, sequence_to_cut, right
            )
    return sequence_dict, sequence_graph

def find_aminoacids_in_sequence(protein_sequence, target_aminoacid):
    return [i for i, aa in enumerate(protein_sequence) if aa == target_aminoacid]


def get_cut_probability(sequence_dict: dict, enzymes: enzyme_set) -> dict:
    """
    Returns probability for enzyme cut given enzymes and sequences
    """
    freqs = {}
    for enzyme, activity in enzymes.enzyme_dict.items():
        specificities = enzyme.specificity
        freq = 0
        for sequence, copy_number in sequence_dict.items():
            for aminoacid, specificity in specificities.items():
                n_cut_sites = len(find_aminoacids_in_sequence(sequence, aminoacid))
                freq += n_cut_sites * specificity * activity * copy_number
        freqs[enzyme.name] = freq
    probs = {k: v / sum(freqs.values()) for k, v in freqs.items()}
    return probs


def update_sequence_dict(
    sequence_dict: dict, source_sequence, target_sequence, endo_or_exo: str
):
    if endo_or_exo == "endo":
        mult_factor = 0.5  # This is here since we want to call this function twice for endoproteases.
    else:
        mult_factor = 1
    sequence_dict[source_sequence] -= 1 * mult_factor
    if sequence_dict[source_sequence] == 0:
        sequence_dict.pop(source_sequence)
    if target_sequence in sequence_dict.keys():
        sequence_dict[target_sequence] += 1
    else:
        sequence_dict[target_sequence] = 1
    return sequence_dict


def update_sequence_graph(sequence_graph: nx.DiGraph, source_sequence, target_sequence):
    if target_sequence not in sequence_graph.nodes():
        sequence_graph.add_node(target_sequence, len=len(target_sequence))

    if ~sequence_graph.has_edge(source_sequence, target_sequence):
        sequence_graph.add_edge(source_sequence, target_sequence, n=1)
    else:
        previous_n = sequence_graph.get_edge_data(source_sequence, target_sequence)["n"]
        new_n = previous_n + 1
        sequence_graph[source_sequence][target_sequence]["n"] = new_n
    return sequence_graph
