from disassembly.simulate_proteolysis_regex import Enzyme, ProteolysisSimulator
from disassembly.util import amino_acids, normalize_dict, KL
import random
import math
import re

amino_acids = list(amino_acids.values())


class Individual:
    """
    Class that represents an individual
    Has checks to assure that the regex doesn't break any rules
    """

    def __init__(self, rules: str = None) -> None:
        if rules:
            if isinstance(rules, str):
                self.rules = rules
            elif isinstance(rules, list):
                self.compile(rules)

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
        if operator == ".":
            add_or_delete = "add"
        aa = random.choice(amino_acids)

        if add_or_delete == "add":

            if type == "inclusion":

                if operator == ".":
                    new_operator = f"[{aa}]"
                else:
                    new_operator = (
                        "[" + operator.removeprefix("[").removesuffix("]") + f"|{aa}]"
                    )

            elif type == "exclusion":

                if operator == ".":
                    new_operator = f"[^{aa}]"
                else:
                    new_operator = (
                        "[" + operator.removeprefix("[").removesuffix("]") + f"|^{aa}]"
                    )

        elif add_or_delete == "delete":
            if operator == ".":

                self.compile(list_of_operators)
                return

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

        list_of_operators[position] = new_operator
        self.compile(list_of_operators)

    def __repr__(self) -> str:
        return self.rules


class ParameterEstimatorGA:
    def __init__(
        self,
        true_distribution: dict,
        protein_sequence: str,
        mutation_rate: float = 0.1,
        n_individuals: int = 10,
        pattern_len: int = 6,
        kill_fraction: int = 0.5,
        n_generate : int = 200
    ) -> None:
        """
        Initiate a random population
        """
        self.n_generate = int(sum(true_distribution.values()))
        self.true_distribution = normalize_dict(true_distribution)
        self.mutation_rate = mutation_rate
        self.n_individuals = n_individuals
        self.protein_sequence = protein_sequence
        self.pattern_len = pattern_len
        self.n_kill = math.floor(kill_fraction * n_individuals)
        print("N Kill ", self.n_kill)
        self.n_generate=n_generate
        self.ps = ProteolysisSimulator(verbose=False)

        self.population = {}  # dict of Enzyme:fitness
        for _ in range(n_individuals):
            random_regex = self._generate_random_regex()
            individual = Individual(random_regex)
            self.population[individual] = 0  # fitness

        self.evaluate_fitness()

    def run(self, n_generations: int = 10):
        """
        Main run-method
        """
        self.fitness = {}
        self.fitness[-1] = list(self.population.values())
        for generation in range(n_generations):
            self.kill_reproduce()
            self.mutate()
            self.evaluate_fitness()
            print("Best ", list(self.population.keys())[0])
            self.fitness[generation] = list(self.population.values())

    def evaluate_fitness(self):
        """
        Evaluate fitness of individuals:

        - simulate proteolysis with given rules (individuals)
        - Compute similarity to true distribution
        """
        for individual in self.population.keys():
            e = Enzyme("_", cleavage_rules={individual.rules: 1})
            t_hat = self.ps.simulate_proteolysis(
                self.protein_sequence,
                enzyme=e,
                n_generate=self.n_generate,
                graph=False,
                length_params="vitro",
            )
            t_hat = normalize_dict(t_hat)
            t1, t2 = self._make_comparable(t_hat, self.true_distribution)
            loss = KL(t1, t2) + KL(t2, t1)
            self.population[individual] = loss

        self.population = {
            k: v for k, v in sorted(self.population.items(), key=lambda item: item[1])
        }
        print("Sorted ", self.population)

    def kill_reproduce(self):
        """
        Combine kill and reroduce into one method
        """
        keys = list(self.population.keys())
        not_kill = keys[: self.n_individuals - self.n_kill]
        kill = keys[self.n_individuals - self.n_kill:]

        for key in kill:
            self.population.pop(key)  # kill individual
            parents = random.choices(not_kill, k=2)
            offspring = self.mate(parents[0], parents[1])  # reproduce
            self.population[offspring] = 0

    def mate(self, ind1: Individual, ind2: Individual):
        """
        Pair 2 individuals to create a new
        """
        position = random.randint(0, 5)
        ind1_operators = ind1.split()
        ind2_operators = ind2.split()

        part_from1 = ind1_operators[:position]
        part_from2 = ind2_operators[position:]
        new_operators = part_from1 + part_from2
        offspring = Individual(new_operators)
        return offspring

    def mutate(self):
        """
        Mutates individuals randomly
        """
        keys = list(self.population.keys())
        for individual in keys[1:]:  # elitism
            u = random.random()
            if u <= self.mutation_rate:
                self.population.pop(individual)
                individual.mutate()
                self.population[individual] = 0

    def _generate_random_regex(self):
        """
        Generate a random regex
        """
        random_regex = ""
        for _ in range(self.pattern_len):
            u = random.randint(0, 1)
            if u == 0:
                random_regex += "(.)"
            else:
                u = random.randint(0, 1)
                aa = random.choice(amino_acids)
                if u == 0:  # exclusion
                    random_regex += f"([^{aa}])"
                else:
                    random_regex += f"([{aa}])"
        return random_regex

    def _make_comparable(self, t1, t2):
        """
        Aligns two dicts and outputs the aligned value-vectors
        """
        t1_vec = []
        t2_vec = []
        for key in t1.keys():
            t1_vec.append(t1[key])
            if key in t2.keys():
                t2_vec.append(t2[key])
            else:
                t2_vec.append(0)
        for key in t2.keys():
            if key not in t1.keys():
                t1_vec.append(0)
                t2_vec.append(t2[key])
        return t1_vec, t2_vec