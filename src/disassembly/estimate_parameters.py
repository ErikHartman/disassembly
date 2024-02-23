from disassembly.util import normalize_dict, amino_acids, KL
from disassembly.simulate_proteolysis import ProteolysisSimulator, enzyme, enzyme_set
import random


class ParameterEstimator:

    def __init__(
        self,
        parameters: dict = None,
        exo: float = 0.25,
    ) -> None:
        self.ps = ProteolysisSimulator(verbose=False)
        if parameters:
            self.parameters = parameters
        else:
            self.parameters = {
                "endo": {
                    aa: 1 / len(amino_acids.values()) for aa in amino_acids.values()
                },
                "exo": exo,
            }

    def estimate(
        self,
        protein: str,
        true_dict: dict,
        n_iterations_endo=3,
        n_iterations_exo=20,
        lr_endo=0.25,
        lr_exo=0.05,
    ):
        """
        Main run-method to estimate parameters
        """
        self.n_generate = int(sum(true_dict.values()))
        print(self.n_generate)
        self.true_dict = true_dict
        self.protein = protein
        self.best_losses = []

        for i in range(n_iterations_endo):
            print(f"Endo iteration: {i}")
            guess = self.generate_guess()
            p, q = compare(self.true_dict, guess)
            self.loss_to_beat = KL(p, q) + KL(q, p)  # baseline
            for aa in self.parameters["endo"].keys():
                _, new_loss = self.update_parameter(aa, lr_endo, verbose=True)
                while new_loss < self.loss_to_beat:
                    diff = self.loss_to_beat - new_loss
                    print(f"{aa} better!")
                    self.loss_to_beat = new_loss
                    self.best_losses.append(new_loss)
                    self.parameters, new_loss = self.update_parameter(
                        aa, lr_endo * diff, verbose=True
                    )
                self.parameters["endo"][aa] -= lr_endo  # this resets the initial guess

        new_guess = self.generate_guess()
        p, q = compare(self.true_dict, new_guess)
        self.loss_to_beat = KL(p, q) + KL(q, p)  # baseline
        for i in range(n_iterations_exo):
            exo_diff = lr_exo * random.choice([-1,1])
            self.parameters["exo"] = self.parameters["exo"] + exo_diff # update parameter
            new_guess = self.generate_guess()
            p, q = compare(self.true_dict, new_guess)
            new_loss = KL(p, q) + KL(q, p)
            if new_loss > self.loss_to_beat:
                self.parameters["exo"] -= exo_diff #revert
            else:
                self.loss_to_beat = new_loss
                self.best_losses.append(new_loss)

            print(f" exo: {new_loss:.2f} | {self.loss_to_beat:.2f} ({self.parameters['exo']:.2f})")
        
        self.parameters["endo"] = normalize_dict(self.parameters["endo"])
        self.parameters["exo"] = max(0, self.parameters["exo"])
        return self.parameters

    def update_parameter(self, aa: str, e: float, verbose: bool = False):
        self.parameters["endo"][aa] += e
        new_guess = self.generate_guess()
        p, q = compare(self.true_dict, new_guess)
        new_loss = KL(p, q) + KL(q, p)
        if verbose:
            print(f"\t{aa}: {new_loss:.2f} | {self.loss_to_beat:.2f}")
        return self.parameters, new_loss

    def generate_guess(self) -> dict:
        """
        Generates a guess from parameters
        """
        norm_params = {"endo":{}, "exo":0}
        for aa in self.parameters["endo"]:
            norm_params["endo"][aa] = self.parameters["endo"][aa]
        norm_params["exo"] = max(0, self.parameters["exo"])
        norm_params["endo"] = normalize_dict(norm_params["endo"])
        parameter_enzyme = enzyme_set([enzyme(norm_params["endo"], "")], [1], [1])
        guess = self.ps.simulate_proteolysis(
            self.protein,
            parameter_enzyme,
            n_start=1,
            n_generate=self.n_generate,
            endo_or_exo_probability=[
                1 - norm_params["exo"],
                norm_params["exo"],
            ],
            graph=False,
        )
        return guess


def compare(P, generated):
    """
    Aligns two dicts and outputs the aligned value-vectors
    """
    P = normalize_dict(P)
    generated = normalize_dict(generated)
    P_vec = []
    generated_vec = []
    for key in P.keys():
        P_vec.append(P[key])
        if key in generated.keys():
            generated_vec.append(generated[key])
        else:
            generated_vec.append(0)
    for key in generated.keys():
        if key not in P.keys():
            P_vec.append(0)
            generated_vec.append(generated[key])
    return P_vec, generated_vec
