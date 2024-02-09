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
from disassembly.util import amino_acids
import random


from scipy.stats import gamma

a, shape, scale = 10, 1, 5
g = gamma(a=a, scale=scale, loc=shape)
x = np.linspace(0, 100)
s = g.pdf(x)
us = np.random.uniform(0, 1, size=10000)


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
    n_start: int = 10,
    n_generate: int = 100,
    endo_or_exo_probability: list = [0.5, 0.5],
    verbose: bool = True,
    graph: bool = True,
    accept_condition: bool = True,
) -> (dict, nx.DiGraph):
    sequence_dict = {starting_sequence: n_start}

    if graph:
        sequence_graph = (
            nx.DiGraph()
        )  # a weighted graph, telling us how degradation has occured
        for sequence in sequence_dict.keys():
            sequence_graph.add_node(sequence, len=len(sequence))  # Adding nodes

    total_iterations = 0
    n_generated_peptides = 0

    exo_or_endos = np.random.choice(
        ["endo", "exo"], size=n_generate * 100, p=endo_or_exo_probability
    )
    n_or_c_terms = np.random.choice(["n", "c"], size=n_generate * 100, p=[0.5, 0.5])

    while n_generated_peptides < n_generate:
        total_iterations += 1
        if verbose:
            print(
                f"\r {n_generated_peptides} / {n_generate} ({total_iterations})",
                end="",
                flush=True,
            )
        exo_or_endo = exo_or_endos[total_iterations]

        if exo_or_endo == "exo":
            seq_dict_keys = list(sequence_dict.keys())
            sequence_to_chew = random.choices(
                seq_dict_keys, weights=sequence_dict.values()
            )[0]

            accept = accept_addition(len(sequence_to_chew) - 1, n_generated_peptides)

            if accept:
                n_generated_peptides += 1
                n_or_c_term = n_or_c_terms[n_generate]
                if n_or_c_term == "n":
                    new_sequence = sequence_to_chew[1:]
                else:
                    new_sequence = sequence_to_chew[:-1]
                sequence_dict = update_sequence_dict(
                    sequence_dict, sequence_to_chew, new_sequence, endo_or_exo="exo"
                )
                if graph:
                    sequence_graph = update_sequence_graph(
                        sequence_graph, sequence_to_chew, new_sequence
                    )

        elif exo_or_endo == "endo":
            sequences_longer_than_12 = {
                s: sequence_dict[s] for s in sequence_dict.keys() if len(s) > 12
            }

            sequence_frequencies = {}
            for sequence in sequences_longer_than_12.keys():
                n_cut_sites_in_sequence = 0
                for aminoacid in enzymes.meta_enzyme.keys():
                    if enzymes.meta_enzyme[aminoacid] != 0:
                        n_cut_sites_in_sequence += (
                            len(find_aminoacids_in_sequence(sequence, aminoacid))
                            * enzymes.meta_enzyme[aminoacid]
                        )

                    sequence_frequencies[sequence] = (
                        n_cut_sites_in_sequence * sequences_longer_than_12[sequence]
                    )

            sequence_to_cut = random.choices(
                list(sequence_frequencies.keys()),
                weights=[
                    p / sum(sequence_frequencies.values())
                    for p in sequence_frequencies.values()
                ],
            )[0]

            index_to_cut = {}
            for aminoacid in enzymes.meta_enzyme.keys():
                indices_for_aminoacid = find_aminoacids_in_sequence(
                    sequence_to_cut, aminoacid
                )
                for index in indices_for_aminoacid:
                    if index != len(sequence_to_cut):
                        index_to_cut[index] = enzymes.meta_enzyme[aminoacid]

            # Perform two cuts
            # Perform two cuts
            cutting_index1 = int(
                random.choices(
                    list(index_to_cut.keys()),
                    weights=[
                        p / sum(index_to_cut.values()) for p in index_to_cut.values()
                    ],
                )[0]
            )

            a, shape, scale = 5.5, 8, 4.5
            g = gamma(a=a, scale=scale, loc=shape)
            for index in index_to_cut.keys():
                index_to_cut[index] = (
                    index_to_cut[index] * g.pdf(abs(index - cutting_index1)) + 1e-8
                )

            cutting_index2 = int(
                random.choices(
                    list(index_to_cut.keys()),
                    weights=[
                        p / sum(index_to_cut.values()) for p in index_to_cut.values()
                    ],
                )[0]
            )

            left = sequence_to_cut[: min(cutting_index1, cutting_index2) + 1]
            middle = sequence_to_cut[
                min(cutting_index1, cutting_index2) : max(
                    cutting_index1, cutting_index2
                )
            ]
            right = sequence_to_cut[max(cutting_index1, cutting_index2) + 1 :]

            # Accept middle
            n_generated_peptides += 1
            sequence_dict = update_sequence_dict(
                sequence_dict, sequence_to_cut, middle, endo_or_exo="endo"
            )
            if graph:
                sequence_graph = update_sequence_graph(
                    sequence_graph, sequence_to_cut, middle
                )
            # Check if accept others
            for sequence in [left, right]:

                accept = accept_addition(len(sequence) - 1, n_generated_peptides) and (
                    (not starting_sequence.startswith(sequence))
                    and (not starting_sequence.endswith(sequence))
                )

                if accept:
                    n_generated_peptides += 1
                    sequence_dict = update_sequence_dict(
                        sequence_dict, sequence_to_cut, sequence, endo_or_exo="endo"
                    )
                    if graph:
                        sequence_graph = update_sequence_graph(
                            sequence_graph, sequence_to_cut, sequence
                        )

    if graph:
        for node in sequence_graph.nodes():
            sequence_graph.add_edge(node, node, weight=sequence_dict[node])
    if verbose:
        print(
            f"\n{len(sequence_dict)} unique peptides. {sum(sequence_dict.values())} total"
        )
    if graph:
        return sequence_dict, sequence_graph
    else:
        return sequence_dict


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

    # if sequence_dict[source_sequence] == 0:
    #    sequence_dict.pop(source_sequence)

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


def accept_addition(length, iteration, min_length=4):
    if length < min_length:
        return False

    a = g.pdf(length) / max(s)

    if us[iteration % 10000] < a:
        return True
    return False
