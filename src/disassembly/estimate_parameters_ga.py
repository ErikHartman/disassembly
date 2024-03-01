from disassembly.simulate_proteolysis_regex import Enzyme
from disassembly.util import amino_acids
import random
import re

amino_acids = list(amino_acids.values())


class Individual:
    """
    Class that represents an individual
    Has checks to assure that the regex doesn't break any rules
    """

    def __init__(self, rules: str = None) -> None:
        if rules:
            self.rules = rules

    def get_length_of_match(self):
        return len(self.split())

    def split(self):
        components = []
        for i, char in enumerate(self.rules):
            if char == "(":
                k = i
            if char == ")":
                components.append(self.rules[k + 1 : i])
        return components

    def compile(self, list_of_operators):
        self.rules = ""
        for operator in list_of_operators:
            self.rules += f"({operator})"
        try:
            re.compile(self.rules)
        except re.error:
            print(self.rules, " is not a valid regex pattern")

    def mutate(self):
        position = random.randint(0, 5)  # index to mutate
        add_or_delete = random.choice(["add", "delete"])
        type = random.choice(["inclusion", "exclusion"])
        list_of_operators = self.split()
        operator: str = list_of_operators[position]
        aa = random.choice(amino_acids)

        if add_or_delete == "add":

            if type == "inclusion":
                print(f"Add inclusion {aa} at {position}")
                if operator == ".":
                    new_operator = f"[{aa}]"
                else:
                    new_operator = (
                        "[" + operator.removeprefix("[").removesuffix("]") + f"|{aa}]"
                    )

            elif type == "exclusion":
                print(f"Add exclusion {aa} at {position}")
                if operator == ".":
                    new_operator = f"[?!{aa}]"
                else:
                    new_operator = (
                        "[" + operator.removeprefix("[").removesuffix("]") + f"|?!{aa}]"
                    )

        elif add_or_delete == "delete":
            if operator == ".":
                print("Do nothing")
                self.compile(list_of_operators)
                return
            print(f"Delete {position}")
            ops = operator.split("|")
            if len(ops) == 1:
                new_operator = "."
            else:
                n_ops = len(ops)
                i_del = random.randint(0, n_ops - 1)
                ops.remove(ops[i_del])
                new_operator = "|".join(ops)
                new_operator = (
                    "[" + new_operator.removeprefix("[").removesuffix("]") + "]"
                )

        print(new_operator)
        list_of_operators[position] = new_operator
        print(list_of_operators)
        self.compile(list_of_operators)


class ParameterEstimatorGA:
    def __init__(
        self, true_distribution, mutation_rate: float = 0.05, n_individuals: int = 10
    ) -> None:
        """
        Initiate a random population
        """
        self.true_distribution = true_distribution
        self.mutation_rate = mutation_rate
        self.n_individuals = n_individuals

        self.population = {}  # dict of Enzyme:fitness
        for _ in n_individuals:
            pass

    def run(self, n_generations: int = 10):
        """
        Main run-method
        """
        for generation in range(n_generations):
            pass

    def evaluate_fitness(self):
        """
        Evaluate fitness of individuals:

        - simulate proteolysis with given rules (individuals)
        - Compute similarity to true distribution
        """
        pass

    def kill(self):
        """
        Kill off a fraction of the population
        """
        pass

    def reproduce(self):
        """
        Creates new individuals from individuals in population
        """
        pass

    def mutate(self):
        """
        Mutates individuals randomly
        """
        pass
