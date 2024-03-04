import networkx as nx
import numpy as np
import random

import re
import numpy as np

from scipy.stats import gamma


class Enzyme:
    def __init__(self, name: str, cleavage_rules: dict):
        self.name = name
        self.cleavage_rules = cleavage_rules

    def cleave(self, protein_sequence: str):
        cleavage_probabilities = np.zeros(len(protein_sequence))
        protein_sequence = "X" + protein_sequence + "X"
        n_cuts = 0
        for pattern, amount in self.cleavage_rules.items():
            regex_pattern = re.compile(pattern)
            matches = regex_pattern.finditer(protein_sequence)
            cleavage_indices = [match.start() + 2 for match in matches]
            for cleavage_index in cleavage_indices:
                cleavage_probabilities[cleavage_index] += amount
                n_cuts += 1
        sum_cleavage_probs = sum(cleavage_probabilities)
        if sum_cleavage_probs == 0:
            cleavage_probabilities = np.ones(len(protein_sequence)-2)
            sum_cleavage_probs = len(protein_sequence)
        cleavage_probabilities = cleavage_probabilities / sum_cleavage_probs
        return cleavage_probabilities, n_cuts


class ProteolysisSimulator:
    def __init__(
        self,
        verbose: bool = True,
    ):
        self.verbose = verbose
        self.in_vitro_params = (1.438694550365602, 6.776824404926209, 4.564749111402035)
        self.in_vivo_params = (1.9155386766281643, 6.463266107715064, 4.557523852970956)

    def simulate_proteolysis(
        self,
        starting_sequence: str,
        enzyme: Enzyme = Enzyme("trypsin", {"(.)(.)([R|K])([^P])(.)(.)": 1}),
        n_start: int = 10,
        n_generate: int = 100,
        endo_or_exo_probability: list = [0.5, 0.5],
        graph: bool = True,
        length_params="vitro",
        noise: float = 1e-8,
    ):

        self.enzyme = enzyme
        if length_params == "vitro":
            a, shape, scale = self.in_vitro_params
        else:
            a, shape, scale = self.in_vivo_params

        self.g = gamma(a=a, scale=scale, loc=shape)
        self.max_gamma_pdf = max(self.g.pdf(np.linspace(0, 100)))
        self.gamma_us = np.random.uniform(0, 1, size=10000)

        self.sequence_dict = {starting_sequence: n_start}

        if graph:
            self.sequence_graph = (
                nx.DiGraph()
            )  # a weighted graph, telling us how degradation has occured
            for sequence in self.sequence_dict.keys():
                self.sequence_graph.add_node(
                    sequence, len=len(sequence)
                )  # Adding nodes

        total_iterations = 0
        self.n_generated_peptides = 0

        exo_or_endos = np.random.choice(
            ["endo", "exo"], size=n_generate * 100, p=endo_or_exo_probability
        )
        n_or_c_terms = np.random.choice(["n", "c"], size=n_generate * 100, p=[0.5, 0.5])

        while self.n_generated_peptides < n_generate:
            total_iterations += 1
            if self.verbose:
                print(
                    f"\r {self.n_generated_peptides} / {n_generate} ({total_iterations})",
                    end="",
                    flush=True,
                )
            exo_or_endo = exo_or_endos[total_iterations]

            if exo_or_endo == "exo":
                seq_dict_keys = list(self.sequence_dict.keys())
                sequence_to_chew = random.choices(
                    seq_dict_keys, weights=self.sequence_dict.values()
                )[0]

                accept = self.accept_addition(len(sequence_to_chew) - 1)

                if accept:
                    self.n_generated_peptides += 1
                    n_or_c_term = n_or_c_terms[n_generate]
                    if n_or_c_term == "n":
                        new_sequence = sequence_to_chew[1:]
                    else:
                        new_sequence = sequence_to_chew[:-1]
                    self.sequence_dict = self.update_sequence_dict(new_sequence)
                    if graph:
                        self.sequence_graph = self.update_sequence_graph(
                            sequence_to_chew, new_sequence
                        )

            elif exo_or_endo == "endo":
                sequences_longer_than_10 = {
                    s: self.sequence_dict[s]
                    for s in self.sequence_dict.keys()
                    if len(s) > 10
                }
                sequence_frequencies = {}
                sequence_probabilities = {}
                for sequence in sequences_longer_than_10.keys():
                    cut_probabilities, n_cut_sites_in_sequence = self.enzyme.cleave(
                        sequence
                    )
                    sequence_probabilities[sequence] = cut_probabilities
                    sequence_frequencies[sequence] = n_cut_sites_in_sequence

                sum_freq = sum(sequence_frequencies.values())
                if sum_freq == 0:
                    _add = 1
                    sum_freq = len(sequence_frequencies.keys())
                else:
                    _add = 0

                sequence_to_cut = random.choices(
                    list(sequence_frequencies.keys()),
                    weights=[
                        (p + _add) / sum_freq for p in sequence_frequencies.values()
                    ],
                )[0]

                # Perform two cuts
                try:
                    cutting_index1 = int(
                        random.choices(
                            list(range(len(sequence_to_cut))),
                            weights=sequence_probabilities[sequence_to_cut],
                        )[0]
                    )
                except:
                    print(sequence_probabilities[sequence_to_cut])
                    print(len(sequence_probabilities[sequence_to_cut]), len(sequence_to_cut))

                index_to_cut = sequence_probabilities[sequence_to_cut]

                if sum(index_to_cut) == 0:
                    for index in range(len(index_to_cut)):
                        index_to_cut[index] = self.g.pdf(abs(index - cutting_index1))
                else:
                    for index in range(len(index_to_cut)):
                        index_to_cut[index] = sequence_probabilities[sequence_to_cut][
                            index
                        ] * self.g.pdf(abs(index - cutting_index1))

                cutting_index2 = int(
                    random.choices(
                        list(range(len(sequence_to_cut))),
                        weights=index_to_cut + noise,
                    )[0]
                )

                left = sequence_to_cut[: min(cutting_index1, cutting_index2)]
                middle = sequence_to_cut[
                    min(cutting_index1, cutting_index2) : max(
                        cutting_index1, cutting_index2
                    )
                ]
                right = sequence_to_cut[max(cutting_index1, cutting_index2) :]

                # Accept middle
                if len(middle) > 5:
                    self.n_generated_peptides += 1
                    self.sequence_dict = self.update_sequence_dict(middle)
                    if graph:
                        self.sequence_graph = self.update_sequence_graph(
                            sequence_to_cut, middle
                        )
                # Check if accept others
                for sequence in [left, right]:

                    accept = self.accept_addition(len(sequence) - 1) and (
                        (not starting_sequence.startswith(sequence))
                        and (not starting_sequence.endswith(sequence))
                    )

                    if accept:
                        self.n_generated_peptides += 1
                        self.sequence_dict = self.update_sequence_dict(sequence)
                        if graph:
                            self.sequence_graph = self.update_sequence_graph(
                                sequence_to_cut, sequence
                            )

        if graph:
            for node in self.sequence_graph.nodes():
                self.sequence_graph.add_edge(
                    node, node, weight=self.sequence_dict[node]
                )
        if self.verbose:
            print(
                f"\n{len(self.sequence_dict)} unique peptides. {sum(self.sequence_dict.values())} total"
            )
        if graph:
            return self.sequence_dict, self.sequence_graph
        else:
            return self.sequence_dict

    def update_sequence_dict(self, target_sequence):
        if target_sequence in self.sequence_dict.keys():
            self.sequence_dict[target_sequence] += 1
        else:
            self.sequence_dict[target_sequence] = 1
        return self.sequence_dict

    def update_sequence_graph(self, source_sequence, target_sequence):
        if target_sequence not in self.sequence_graph.nodes():
            self.sequence_graph.add_node(target_sequence, len=len(target_sequence))

        if ~self.sequence_graph.has_edge(source_sequence, target_sequence):
            self.sequence_graph.add_edge(source_sequence, target_sequence, weight=1)
        else:
            previous_n = self.sequence_graph.get_edge_data(
                source_sequence, target_sequence
            )["weight"]
            new_n = previous_n + 1
            self.sequence_graph[source_sequence][target_sequence]["weight"] = new_n
        return self.sequence_graph

    def accept_addition(self, length, min_length=4):
        if length < min_length:
            return False

        a = self.g.pdf(length) / self.max_gamma_pdf

        if self.gamma_us[self.n_generated_peptides % 10000] < a:
            return True
        return False

    def format_graph(self):
        formatted_graph = nx.DiGraph()
        for source, target, data in self.sequence_graph.edges(data=True):
            sum_out_edges = sum(
                [
                    data["weight"]
                    for _, _, data in self.sequence_graph.out_edges(source, data=True)
                ]
            )
            if source != target:
                formatted_graph.add_edge(
                    source, target, weight=data["weight"] / (sum_out_edges)
                )
        return formatted_graph
