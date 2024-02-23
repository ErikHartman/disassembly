from disassembly.simulate_proteolysis import ProteolysisSimulator
from disassembly.estimate_weights_alg import estimate_weights
from disassembly.disassembly import get_disassembly, get_disassembly_indexes_mc
from disassembly.estimate_weights_gd import WeightEstimatorGD
from disassembly.estimate_parameters import ParameterEstimator
import matplotlib.pyplot as plt
from disassembly.util import get_nrmse
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
import math


class Benchmark:
    """
    Benchmarking disassembly methods


    """

    def __init__(self):
        self.results = {}
        self.simulated_peptidomes = {}
        self.simulated_graphs = {}
        self.generated_graphs = {}
        self.ps = ProteolysisSimulator()

    def simulate_degradation(
        self,
        protein,
        enzyme_sets: list,
        n_generate: int = 500,
        iterations: int = 5,
        n_start: int = 3,
        endo_or_exo_probability: list = [0.9, 0.1],
        enzyme_names: list = None,
        di_mc_n=10000,
    ):
        """
        Simulates proteolysis and saves to
        `simulated_graphs` and `simulated_peptidomes`
        """
        self.iterations = iterations
        self.enzyme_sets = enzyme_sets
        self.protein = protein
        self.n_generate = n_generate

        self.results["real"] = {}

        if not enzyme_names:
            self.enzyme_names = list(range(len(enzyme_sets)))
        else:
            self.enzyme_names = enzyme_names

        for enzyme_name, enzyme_set in zip(self.enzyme_names, self.enzyme_sets):
            self.results["real"][enzyme_name] = {}
            self.simulated_peptidomes[enzyme_name] = {}
            self.simulated_graphs[enzyme_name] = {}
            for iteration in range(iterations):
                print(f"Running {enzyme_name}, {iteration}")
                self.results["real"][enzyme_name][iteration] = {}
                simulated_peptidome, _ = self.ps.simulate_proteolysis(
                    protein,
                    enzyme_set,
                    n_start=n_start,
                    n_generate=self.n_generate,
                    endo_or_exo_probability=endo_or_exo_probability,
                )
                self.simulated_peptidomes[enzyme_name][iteration] = simulated_peptidome
                self.simulated_graphs[enzyme_name][iteration] = self.ps.format_graph()

                self.results["real"][enzyme_name][iteration]["di"] = (
                    get_disassembly_indexes_mc(
                        self.simulated_graphs[enzyme_name][iteration], di_mc_n
                    )
                )
                self.results["real"][enzyme_name][iteration]["d"] = get_disassembly(
                    self.simulated_peptidomes[enzyme_name][iteration],
                    self.results["real"][enzyme_name][iteration]["di"],
                )

    def estimate_weights(
        self,
        method: str = "gd",  # gd, alg
        lam1: float = 0,
        lam2: float = 0,
        n_iterations=100,
        lr=0.1,
        parameter_estimator=False,
        n_iterations_endo=1,
        n_iterations_exo=10,
        di_mc_n=10000,
        exo=0.2,
        method_name=None,
        lr_scheduler={},
    ):

        if not method_name:
            method_name = method

        self.results[method_name] = {}
        self.generated_graphs[method_name] = {}

        if parameter_estimator:
            self.results["param"]= {}
            self.generated_graphs["param"] = {}

        if method == "gd":
            wegd = WeightEstimatorGD(
                lr=lr,
                n_iterations=n_iterations,
                lam1=lam1,
                lam2=lam2,
                lr_scheduler=lr_scheduler,
            )
        for enzyme_name in self.enzyme_names:
            print(f"---{enzyme_name}---")
            self.results[method_name][enzyme_name] = {}
            self.generated_graphs[method_name][enzyme_name] = {}
            if parameter_estimator:
                self.results["param"][enzyme_name] = {}
                self.generated_graphs["param"][enzyme_name] = {}
            for iteration in range(self.iterations):
                self.results[method_name][enzyme_name][iteration] = {}
                if parameter_estimator:
                    self.results["param"][enzyme_name][iteration] = {}
                if method == "alg":
                    G, losses, _, _ = estimate_weights(
                        P=self.simulated_peptidomes[enzyme_name][iteration],
                        lr=lr,
                        n_iterations=n_iterations,
                    )
                elif method == "gd":
                    if parameter_estimator:
                        pe = ParameterEstimator(exo=exo)
                        pe.estimate(
                            self.protein,
                            self.simulated_peptidomes[enzyme_name][iteration],
                            n_iterations_endo=n_iterations_endo,
                            n_iterations_exo=n_iterations_exo,
                        )

                        parameters = pe.parameters
                        losses = pe.best_losses
                        # Save results before gd
                        wegd.parameters = parameters
                        G = wegd.create_graph_from_parameters(
                            self.simulated_peptidomes[enzyme_name][iteration]
                        )
                        self.generated_graphs["param"][enzyme_name][iteration] = G
                        self.results["param"][enzyme_name][iteration]["loss"] = losses

                        self.results["param"][enzyme_name][iteration]["di"] = (
                            get_disassembly_indexes_mc(G, di_mc_n)
                        )
                        self.results["param"][enzyme_name][iteration]["d"] = (
                            get_disassembly(
                                self.simulated_peptidomes[enzyme_name][iteration],
                                self.results["param"][enzyme_name][iteration]["di"],
                            )
                        )
                    else:
                        parameters = None

                    G = wegd.run(
                        self.simulated_peptidomes[enzyme_name][iteration],
                        verbose=True,
                        parameters=parameters,
                    )
                    losses = wegd.losses
                else:
                    raise ValueError("method must be either gd or alg")

                self.generated_graphs[method_name][enzyme_name][iteration] = G
                self.results[method_name][enzyme_name][iteration]["loss"] = losses

                self.results[method_name][enzyme_name][iteration]["di"] = (
                    get_disassembly_indexes_mc(G, di_mc_n)
                )
                self.results[method_name][enzyme_name][iteration]["d"] = (
                    get_disassembly(
                        self.simulated_peptidomes[enzyme_name][iteration],
                        self.results[method_name][enzyme_name][iteration]["di"],
                    )
                )

    def plot(
        self,
        ax,
        ptype="loss",
        method_name="gd",
        colors=["darkblue", "purple", "pink", "gray"],
    ):
        if len(colors) < len(self.enzyme_names):
            raise ValueError(
                "Length of colors must be equal to or longer than length of enzymes"
            )

        if ptype == "loss":
            for i, test_name in enumerate(self.enzyme_names):
                for iteration in range(self.iterations):
                    ax.plot(
                        self.results[method_name][test_name][iteration]["loss"],
                        label=test_name,
                        c=colors[i],
                        alpha=0.25,
                    )
        elif ptype == "corr_di":
            nrmse = {}
            for i, test_name in enumerate(self.enzyme_names):
                nrmse[test_name] = {}
                for iteration in range(self.iterations):
                    real_di = self.results["real"][test_name][iteration]["di"]
                    estimated_di = self.results[method_name][test_name][iteration]["di"]

                    r_di = []
                    e_di = []
                    for sequence in real_di.keys():
                        r_di.append(real_di[sequence])
                        e_di.append(estimated_di[sequence])

                    ax.scatter(r_di, e_di, c=colors[i], alpha=0.4, s=10)
                    nrmse[test_name][iteration] = get_nrmse(r_di, e_di)
        elif ptype == "d":
            for i, test_name in enumerate(self.enzyme_names):
                for iteration in range(self.iterations):
                    real_d = self.results["real"][test_name][iteration]["d"]
                    estimated_d = self.results[method_name][test_name][iteration]["d"]

                    ax.scatter(real_d, estimated_d, color=colors[i])

        patches = []
        for test_name, color in zip(self.enzyme_names, colors):
            patches.append(mpatches.Patch(color=color, label=test_name))

        plt.legend(
            handles=patches,
            bbox_to_anchor=(-0.1, 1.25),
            ncol=4,
            title="Enzyme complexity",
        )

    def plot_d_error(self):
        df = {"alg": [], "enzyme": [], "d": []}
        for alg in self.results.keys():
            if alg == "real":
                continue
            for enzyme_name in self.enzyme_names:
                for iteration in range(5):
                    df["alg"].append(alg)
                    df["enzyme"].append(enzyme_name)
                    df["d"].append(
                        math.abs(
                            self.results["real"][enzyme_name][iteration]["d"]
                            - self.results[alg][enzyme_name][iteration]["d"]
                        )
                    )
        df = pd.DataFrame(df)
        sns.boxplot(df, x="enzyme", y="d", hue="alg")

    def plot_weight_corr(self, alg_name="gd"):
        fig = plt.figure(
            layout="constrained",
            figsize=(len(self.enzyme_names) * 3, self.iterations * 3),
        )
        subfigs = fig.subfigures(
            self.iterations,
            len(self.enzyme_names),
        )
        for enzyme_name in self.enzyme_names:
            for iteration in range(3):
                in_both = 0
                in_estimate = 0
                in_real = 0
                real_g = self.simulated_graphs[enzyme_name][iteration]
                g = self.generated_graphs[alg_name][enzyme_name][iteration]
                real_vs_estimated_weights = []
                for node in real_g.nodes():
                    sum_out_edges = sum(
                        [
                            data["weight"]
                            for _, _, data in real_g.out_edges(node, data=True)
                        ]
                    )
                    # TODO: Should P[source] be added to sum_out_edges?
                    for source, target, _ in g.out_edges(node, data=True):
                        estimated_weight = g[source][target]["weight"]
                        if real_g.has_edge(source, target):
                            real_weight = (
                                real_g[source][target]["weight"] / sum_out_edges
                            )

                            real_vs_estimated_weights.append(
                                (real_weight, estimated_weight)
                            )
                            in_both += estimated_weight
                        else:
                            in_estimate += estimated_weight

                    for source, target, _ in real_g.out_edges(node, data=True):
                        if ~g.has_edge(source, target):
                            in_real += real_g[source][target]["weight"] / sum_out_edges

                subfigs[list(range(self.iterations)).index(iteration)][
                    self.enzyme_names.index(enzyme_name)
                ].suptitle(f"{enzyme_name} {iteration}")
                axs = subfigs[list(range(self.iterations)).index(iteration)][
                    self.enzyme_names.index(enzyme_name)
                ].subplots(1, 2, width_ratios=[1, 3])
                axs[0].bar(
                    x=["In both", "In est.", "In real"],
                    height=[in_both, in_estimate, in_real],
                    color=["blue", "darkorange", "red"],
                )
                axs[0].set_xticks(
                    ["In both", "In est.", "In real"],
                    labels=["In both", "In est.", "In real"],
                    rotation=90,
                )
                axs[0].set_ylabel("Sum of weights")
                axs[1].scatter(
                    x=[r for r, e in real_vs_estimated_weights],
                    y=[e for r, e in real_vs_estimated_weights],
                    color="black",
                    alpha=0.1,
                )
                axs[1].set_xlabel("real")
                axs[1].set_ylabel("est.")
