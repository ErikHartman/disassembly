"""
TODO: add code to simulate proteolysis given set of enzymes and substate

returns dict of seq:copy_number


Could also create a method that returns the weighted network depending on which paths were taken.
"""

import numpy as np
import networkx as nx


def find_aa(protein_sequence, target):
    indices = []
    for i, amino_acid in enumerate(protein_sequence):
        if amino_acid == target:
            indices.append(i)
    return indices

def get_probability_for_cut(sequences_dict, enzyme):
    # Probability should be n_cut_sites*copy_number
    total_prob = 0
    prob_sequences = []
    for seq, copy_nr  in sequences_dict.items():
        cut_sites = 0
        for aa in enzyme.aas:
            cut_sites += len(find_aa(seq, aa))
        prob_sequence = copy_nr*cut_sites
        prob_sequences.append(prob_sequence)
        total_prob += prob_sequence
    return [x/total_prob for x in prob_sequences]


def cleave(protein, enzyme=None, n_cleaves=10, n_starting_proteins=10):
    sequences_dict = {protein:n_starting_proteins}

    for _ in range(n_cleaves):
        sequence = np.random.choice(list(sequences_dict.keys()), p=get_probability_for_cut(sequences_dict, enzyme))
        sequences_dict[sequence] = sequences_dict[sequence] - 1

        indices = []
        probabilities = []
        if enzyme:
            for aa in enzyme.aas:
                new_indices= find_aa(sequence, aa)
                indices += new_indices
                probabilities += [enzyme.specificity[aa]]*len(new_indices)
        total_probabilities = sum(probabilities)
        probabilities = [p/total_probabilities for p in probabilities]

        # cleave at random aa in sequence
        index = np.random.choice(indices, p=probabilities)

        left = sequence[:index+1]
        if left in sequences_dict.keys():
            sequences_dict[left] = sequences_dict[left] + 1

        right = sequence[index+1:]

        if right in sequences_dict.keys():
            sequences_dict[right] = sequences_dict[right] + 1
            
        if len(right) > 6:
            if right not in sequences_dict.keys():
                sequences_dict[right] = 1
        if len(left) > 6:
            if left not in sequences_dict.keys():
                sequences_dict[left] = 1
        new_dict = {}
        for k, v in sequences_dict.items():
            if v > 0:
                new_dict[k] = v
        sequences_dict = new_dict

    return sequences_dict




class enzyme():
    """
    Enzyme with p1 specificity
    """
    def __init__(self, specificity:dict):
        self.specificity = specificity
        assert sum(specificity.values()) - 1 < 0.001
        self.aas = list(specificity.keys())






def cleaveAlternative(sequences, enzyme=None, n_cleaves=10):
    G = nx.Graph()
    G.add_node(sequences, copy_number=n_cleaves)
    if isinstance(sequences,str):
        sequences = [sequences]

    for _ in range(n_cleaves):
        sequence = np.random.choice(sequences, p=get_probability_for_cut(G, enzyme))
        G.nodes()[sequence]["copy_number"] = G.nodes()[sequence]["copy_number"] - 1

        indices = []
        probabilities = []
        if enzyme:
            for aa in enzyme.aas:
                new_indices= find_aa(sequence, aa)
                indices += new_indices
                probabilities += [enzyme.specificity[aa]]*len(new_indices)
        total_probabilities = sum(probabilities)
        probabilities = [p/total_probabilities for p in probabilities]

        # cleave at random aa in sequence
        index = np.random.choice(indices, p=probabilities)

        left = sequence[:index+1]
        if left in G.nodes():
            G.nodes()[left]["copy_number"] = G.nodes()[left]["copy_number"] + 1

        right = sequence[index+1:]

        if right in G.nodes():
            G.nodes()[right]["copy_number"] = G.nodes()[right]["copy_number"] + 1
            
        if len(right) > 6:
            if right not in sequences:
                sequences.append(right)
                G.add_node(right, copy_number=1)
                G.add_edge(sequence, right)
                
        if len(left) > 6:
            if left not in sequences:
                sequences.append(left)
                G.add_node(left, copy_number=1)
                G.add_edge(sequence, left)
    return G
