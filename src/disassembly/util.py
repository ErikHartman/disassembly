import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def plot_peptidome(protein : str, sequence_dict : dict, ax):
    sequence_dict = dict(
        sorted(sequence_dict.items(), key=lambda item: len(item[0]), reverse=True)
    )
    cmap = matplotlib.cm.coolwarm
    spaces = np.zeros(((len(sequence_dict.keys())), len(protein)))
    for sequence, copy_number in sequence_dict.items():
        start = protein.find(sequence)
        end = start + len(sequence)
        for height in range(spaces.shape[0]):
            position = spaces[height, start:end]
            if sum(position) == 0:
                spaces[height, start:end] = 1
                ax.plot(
                    [start + 1, end - 1],
                    [-height, -height],
                    linewidth=2,
                    color=cmap(copy_number),
        
                )
                break
