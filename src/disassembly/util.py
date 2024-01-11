import numpy as np
import matplotlib


def normalize_dict(d):
    s = sum(d.values())
    return {k: v / s for k, v in d.items()}


def KL(true, guess):
    true = np.asarray(list(true))
    guess = np.asarray(list(guess))
    true = 1e-8 + true / np.sum(true)
    guess = 1e-8 + guess / np.sum(guess)

    return np.sum(np.where(true != 0, true * np.log(true / guess), 0))


def get_nrmse(true, observed):
    true = np.array(true)
    observed = np.array(observed)
    rmse = np.sqrt(np.mean((true - observed) ** 2))
    nrmse = rmse / np.mean(true)
    return nrmse

def get_trend(data):
    diff = np.diff(data)
    mean_diff = np.mean(diff)
    std_diff = np.std(diff)
    trend_threshold = 1e-3
    stochastic_threshold = 1
    if mean_diff > trend_threshold:
        return "Increasing"
    elif mean_diff < -trend_threshold:
        return "Decreasing"
    elif std_diff > stochastic_threshold:
        return "Stochastic"
    else:
        return "Plateau"


def plot_peptidome(protein: str, sequence_dict: dict, ax):
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
