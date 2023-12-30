"""
Simulated proteolysis given a sequence and set of enzymes

An enzyme:
 - has a p1-specificity dict
 - has an activity
 - has an abundance

A protein is degraded by an endoprotease (enzyme) or exoprotease 
with probabilities endo_or_exo_probability (default [0.5, 0.5]).

If degraded by an exoprotease, it is either degraded 1 step at the 
c or n-terminal.

If degraded by an endoprotease, the probability that it is cut
by a given enzyme is computed by:

for each amino acid in specificity:
    n_available_cutsites * enzyme_abundance * specificity * substrate_copy_number

The output is both a sequence_dict, i.e., what we observe at the end
and a sequence_graph, which has weights for all degradation paths (ground truth).


"""

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
    def __init__(self, enzymes: list, activities: list, abundances: list):
        total_activity = sum(activities)
        activities = [a / total_activity for a in activities]
        assert len(enzymes) == len(activities)

        self.enzyme_dict = {
            e: (activity, abundance)
            for e, activity, abundance in zip(enzymes, activities, abundances)
        }

        meta_enzyme_dict = {aa: 0 for aa in amino_acids.values()}
        for enzyme, activity, abundance in zip(enzymes, activities, abundances):
            for amino_acid, value in enzyme.specificity.items():
                meta_enzyme_dict[amino_acid] += value * activity * abundance

        self.meta_enzyme = {
            aa: value / sum(meta_enzyme_dict.values())
            for aa, value in meta_enzyme_dict.items()
        }

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
        [3, 2, 3],  # activities
        [1, 1, 3],  # abundances
    ),
    n_start: int = 100,
    n_iterations: int = 100,
    endo_or_exo_probability: list = [0.5, 0.5],
) -> (dict, nx.DiGraph):
    sequence_dict = {starting_sequence: n_start}  # a dict with sequence:copy_number
    sequence_graph = (
        nx.DiGraph()
    )  # a weighted graph, telling us how degradation has occured
    sequence_graph.add_node(starting_sequence, len=len(starting_sequence))

    for i in range(n_iterations):
        print(f"\r {i} / {n_iterations}", end="")
        exo_or_endo = np.random.choice(["endo", "exo"], p=endo_or_exo_probability)
        if exo_or_endo == "exo":
            sequences_longer_than_4 = [s for s in sequence_dict.keys() if len(s) > 4]

            sequence_to_chew = np.random.choice(
                sequences_longer_than_4,
                p=[
                    sequence_dict[v]
                    / sum([sequence_dict[s] for s in sequences_longer_than_4])
                    for v in sequences_longer_than_4
                ],
            )

            n_or_c_term = np.random.choice(["weight", "c"], p=[0.5, 0.5])
            if n_or_c_term == "weight":
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
            sequences_longer_than_8 = {
                s: sequence_dict[s] for s in sequence_dict.keys() if len(s) > 8
            }

            sequence_frequencies = {}
            for sequence in sequences_longer_than_8.keys():
                n_cut_sites_in_sequence = 0
                for aminoacid in enzymes.meta_enzyme.keys():
                    n_cut_sites_in_sequence += (
                        len(find_aminoacids_in_sequence(sequence, aminoacid))
                        * enzymes.meta_enzyme[aminoacid]
                    )
                sequence_frequencies[sequence] = (
                    n_cut_sites_in_sequence * sequences_longer_than_8[sequence]
                )

            sequence_to_cut = np.random.choice(
                list(sequence_frequencies.keys()),
                p=[
                    p / sum(sequence_frequencies.values())
                    for p in sequence_frequencies.values()
                ],
            )

            index_to_cut = {}
            for aminoacid in enzymes.meta_enzyme.keys():
                indices_for_aminoacid = find_aminoacids_in_sequence(
                    sequence_to_cut, aminoacid
                )
                for index in indices_for_aminoacid:
                    if index != len(sequence_to_cut):
                        index_to_cut[index] = enzymes.meta_enzyme[aminoacid]

            cutting_index = np.random.choice(
                list(index_to_cut.keys()),
                p=[p / sum(index_to_cut.values()) for p in index_to_cut.values()],
            )
            left = sequence_to_cut[: cutting_index + 1]
            right = sequence_to_cut[cutting_index + 1 :]
            if len(left) > 3 and len(right) > 3:
                sequence_dict = update_sequence_dict(
                    sequence_dict, sequence_to_cut, left, endo_or_exo="endo"
                )
                sequence_graph = update_sequence_graph(
                    sequence_graph, sequence_to_cut, left
                )
                sequence_dict = update_sequence_dict(
                    sequence_dict, sequence_to_cut, right, endo_or_exo="endo"
                )
                sequence_graph = update_sequence_graph(
                    sequence_graph, sequence_to_cut, right
                )

    for node in sequence_graph.nodes():
        sequence_graph.add_edge(node, node, weight=sequence_dict[node])
    return sequence_dict, sequence_graph


def find_aminoacids_in_sequence(protein_sequence, target_aminoacid):
    indexes = [i for i, aa in enumerate(protein_sequence) if aa == target_aminoacid]
    if len(protein_sequence) - 1 in indexes:
        indexes.remove(len(protein_sequence) - 1)
    return indexes


def get_cut_probability(sequence_dict: dict, enzymes: enzyme_set) -> dict:
    """
    Returns probability for enzyme cut given enzymes and sequences
    """
    freqs = {}
    for enzyme, (activity, enzyme_abundance) in enzymes.enzyme_dict.items():
        specificities = enzyme.specificity
        freq = 0
        for sequence, copy_number in sequence_dict.items():
            for aminoacid, specificity in specificities.items():
                n_cut_sites = len(find_aminoacids_in_sequence(sequence, aminoacid))
                freq += (
                    n_cut_sites
                    * specificity
                    * activity
                    * copy_number
                    * enzyme_abundance
                )
        freqs[enzyme.name] = freq
    probs = {k: v / sum(freqs.values()) for k, v in freqs.items()}
    return probs


def update_sequence_dict(
    sequence_dict: dict, source_sequence, target_sequence, endo_or_exo: str
):
    # Removed this for assumption-reasons.
    # if endo_or_exo == "endo":
    #     sequence_dict[source_sequence] -= .5  # This is here since we want to call this function twice for endoproteases.
    # else:
    #     sequence_dict[source_sequence] -= 1

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
        sequence_graph.add_edge(source_sequence, target_sequence, weight=1)
    else:
        previous_n = sequence_graph.get_edge_data(source_sequence, target_sequence)[
            "weight"
        ]
        new_n = previous_n + 1
        sequence_graph[source_sequence][target_sequence]["weight"] = new_n
    return sequence_graph
