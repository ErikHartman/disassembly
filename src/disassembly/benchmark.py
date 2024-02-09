from disassembly.simulate_proteolysis import simulate_proteolysis
from disassembly.estimate_weights import estimate_weights
from disassembly.disassembly import get_disassembly, get_disassembly_indexes_mc
from disassembly.estimate_weights_gd import WeightEstimatorGD
from disassembly.estimate_parameters import ParameterEstimator
import matplotlib.pyplot as plt
from disassembly.util import get_nrmse
import matplotlib.patches as mpatches


class Benchmark:

    def __init__(self):
        self.results = {
            "real": {},
            "gd": {},
            "alg": {},
            "param_gd": {},
            "param_alg": {},
        }
        self.simulated_peptidomes = {}
        self.simulated_graphs = {}
        self.generated_graphs = {
            "real": {},
            "gd": {},
            "alg": {},
            "param_gd": {},
            "param_alg": {},
        }

    def simulate_degradation(
        self,
        protein,
        enzyme_sets: list,
        n_generate: int = 500,
        iterations: int = 5,
        n_start: int = 3,
        endo_or_exo_probability: list = [0.5, 0.5],
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
                simulated_peptidome, simulated_graph = simulate_proteolysis(
                    protein,
                    enzyme_set,
                    n_start=n_start,
                    n_generate=n_generate,
                    endo_or_exo_probability=endo_or_exo_probability,
                )
                self.simulated_peptidomes[enzyme_name][iteration] = simulated_peptidome
                self.simulated_graphs[enzyme_name][iteration] = simulated_graph

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
        method: str = "gd",
        lam: float = 0,
        n_iterations=100,
        lr=0.1,  # gd, alg,
        parameter_estimator=False,
        n_iterations_endo=1,
        n_iterations_exo=10,
        n_generate=500,
        di_mc_n=10000,
        exo=0.2,
    ):
        method_name = method
        if parameter_estimator:
            method_name += "_param"

        if method == "gd":
            wegd = WeightEstimatorGD(lr=lr, n_iterations=n_iterations, lam=lam)
        for enzyme_name in self.enzyme_names:
            self.results[method_name][enzyme_name] = {}
            self.generated_graphs[method_name][enzyme_name] = {}
            for iteration in range(self.iterations):
                self.results[method_name][enzyme_name][iteration] = {}
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
                            n_generate=n_generate,
                        )

                        parameters = pe.parameters
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
        if not len(colors) == len(self.enzyme_names):
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
