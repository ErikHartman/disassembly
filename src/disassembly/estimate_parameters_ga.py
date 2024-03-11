from disassembly.simulate_proteolysis_regex import Enzyme, ProteolysisSimulator
from disassembly.util import amino_acids, normalize_dict, KL
import random
import math
import re
import numpy as np
import pandas as pd
import logomaker

amino_acids = list(amino_acids.values())

class Individual:
    """
    Class that represents an individual
    Has checks to assure that the regex doesn't break any rules
    """
    id = 0

    def __init__(self, rules: list = None) -> None:
        self.rules = rules
        Individual.id += 1
        self.id = Individual.id

    def get_length_of_match(self):
        return len(self.split())

    def split(self):
        list_of_components = []
        for rule, amount in self.rules:
            components = []
            for i, char in enumerate(rule):
                if char == "(":
                    k = i
                if char == ")":
                    components.append(rule[k + 1 : i])
            list_of_components.append((components, amount))
        return list_of_components

    def compile(self, list_of_operators):
        rule = ""
        for operator in list_of_operators:
            rule += f"({operator})"
        try:
            re.compile(rule)
        except re.error:
            print(rule, " is not a valid regex pattern")
        return rule

    def mutate_rule(self):
        rule_index = random.randint(0, len(self.rules) - 1)
        position = random.randint(0, 5)  # index to mutate
        add_or_delete = random.choice(["add", "delete"])
        type = random.choice(["inclusion", "exclusion"])
        list_of_operators, amount = self.split()[rule_index]
        operator: str = list_of_operators[position]

        if operator == ".":
            add_or_delete = "add"
        aa = random.choice(amino_acids)

        if aa in operator:
            add_or_delete = "delete"

        if add_or_delete == "add":

            if type == "inclusion":

                if operator == ".":
                    new_operator = f"[{aa}]"
                else:
                    if "^" in operator:
                        new_operator = (
                            "["
                            + operator.removeprefix("[^").removesuffix("]")
                            + f"|{aa}]"
                        )
                    else:
                        new_operator = (
                            "["
                            + operator.removeprefix("[").removesuffix("]")
                            + f"|{aa}]"
                        )

            elif type == "exclusion":

                if operator == ".":
                    new_operator = f"[^{aa}]"
                else:
                    if "^" in operator:
                        new_operator = (
                            "["
                            + operator.removeprefix("[").removesuffix("]")
                            + f"|{aa}]"
                        )
                    else:
                        new_operator = (
                            "[^"
                            + operator.removeprefix("[").removesuffix("]")
                            + f"|{aa}]"
                        )

        elif add_or_delete == "delete":
            if operator == ".":

                new_rule = self.compile(list_of_operators)
                self.rules[rule_index] = new_rule
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
        new_rule = self.compile(list_of_operators)
        self.rules[rule_index] = (new_rule, amount)

    def mutate_amount(self, rule):
        regex, amount = rule
        self.rules.remove(rule)
        add_or_dec = random.choice(["add", "dec"])
        if add_or_dec == "add":
            amount += 1
        else:
            amount -= 1
        if amount > 0:
            self.rules.append((regex, amount))
        elif len(self.rules) == 0:
            self.rules.append(("(.)(.)(.)(.)(.)(.)", 1))

    def __repr__(self) -> str:
        s = str(self.id) + ":"
        for rule, amount in self.rules:
            s +=  rule + " " + str(amount) + " "
        return s


class ParameterEstimatorGA:
    def __init__(
        self,
        true_distributions: list,
        protein_sequences: list,
        mutation_rate: float = 0.1,
        n_individuals: int = 10,
        pattern_len: int = 6,
        kill_fraction: int = 0.5,
        n_generate: int = 200,
        length_penalty: float = 0.1,
    ) -> None:
        """
        Initiate a random population
        """
        self.true_distributions = [
            normalize_dict(true_distribution)
            for true_distribution in true_distributions
        ]
        self.mutation_rate = mutation_rate
        self.n_individuals = n_individuals
        self.protein_sequences = protein_sequences
        self.pattern_len = pattern_len
        self.n_generate = n_generate
        self.n_kill = math.floor(kill_fraction * n_individuals)
        print("N Kill ", self.n_kill)

        self.length_penalty = length_penalty
        self.ps = ProteolysisSimulator(verbose=False)

        self.population = {}  # dict of Enzyme:fitness
        for n_ind in range(n_individuals):
            ind_regex = self._initialize_aminoacids(
                n_ind
            )  # self._generate_random_regex()
            individual = Individual([(ind_regex, 1)])
            self.population[individual] = 0  # fitness

        self.evaluate_fitness(0)

    def run(self, n_generations: int = 10):
        """
        Main run-method
        """
        self.all_results = {}
        self.best_ever = None
        self.best_ever_fitness = 10e3
        self.fitness = {}
        self.fitness[-1] = list(self.population.values())
        
        for generation in range(n_generations):
            self.kill_reproduce()
            self.mutate()
            self.evaluate_fitness(generation)

            for individual, fitness in self.population.items():
                self.all_results[individual] = fitness

            best_individual = list(self.population.keys())[0]
            best_fitness = self.population[best_individual]

            if best_fitness < self.best_ever_fitness:
                self.best_ever = best_individual
                self.best_ever_fitness = best_fitness

            print(best_individual, best_fitness)

            self.fitness[generation] = list(self.population.values())

    def evaluate_fitness(self, generation):
        """
        Evaluate fitness of individuals:

        - simulate proteolysis with given rules (individuals)
        - Compute similarity to true distribution
        """
        n_protein_sequences = len(self.protein_sequences)
        protein_sequence = self.protein_sequences[generation % n_protein_sequences]
        true_distribution = self.true_distributions[generation % n_protein_sequences]
        for individual in self.population.keys():
            e = Enzyme(name="_", cleavage_rules=individual.rules)
            t_hat = self.ps.simulate_proteolysis(
                protein_sequence,
                enzyme=e,
                n_generate=self.n_generate,
                graph=False,
                length_params="vitro",
            )
            penalty = self.length_penalty * sum(
                [len(re.findall("[A-Z]", rule)) for rule, _ in individual.rules]
            )
            t_hat = normalize_dict(t_hat)
            t1, t2 = self._make_comparable(t_hat, true_distribution)
            loss = (
                (KL(t1, t2) + KL(t2, t1)) / KL(np.ones(len(t1)) / len(t1), t1)
            ) + penalty
            self.population[individual] = loss

        self.population = {
            k: v for k, v in sorted(self.population.items(), key=lambda item: item[1])
        }

    def kill_reproduce(self):
        """
        Combine kill and reroduce into one method
        """
        keys = list(self.population.keys())
        not_kill = keys[: self.n_individuals - self.n_kill]
        kill = keys[self.n_individuals - self.n_kill :]

        for key in kill:
            self.population.pop(key)  # kill individual
            parents = random.choices(not_kill, k=2)
            offspring = self.mate(parents[0], parents[1])  # reproduce
            self.population[offspring] = 0

    def mate(self, ind1: Individual, ind2: Individual):
        """
        Pair 2 individuals to create a new
        """
        offspring_rules = ind1.rules + ind2.rules
        offspring = Individual(offspring_rules)
        return offspring

    def mutate(self):
        """
        Mutates individuals randomly
        """
        keys = list(self.population.keys())
        for individual in keys[3:]:  # elitism
            individual: Individual
            u = random.random()
            if u <= self.mutation_rate:
                self.population.pop(individual)
                individual.mutate_rule()
                self.population[individual] = 0

        for individual in keys[3:]:
            for rule in individual.rules:
                u = random.random()
                if u <= self.mutation_rate:
                    self.population.pop(individual)
                    individual.mutate_amount(rule)
                    self.population[individual] = 0

    def _initialize_aminoacids(self, i):
        n_amino_acids = len(amino_acids)
        i = i % n_amino_acids
        aa = amino_acids[i]
        return f"(.)(.)([{aa}])(.)(.)(.)"

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

    def plot(self, topn : int = 10):
        d = {}
        top_population = {
            k: v for k, v in sorted(self.all_results.items(), key=lambda item: item[1])[0:topn]
        }
        for pos in range(self.pattern_len):
            d[pos] = {}
            for aa in amino_acids:
                d[pos][aa] = 0

        max_fitness = max(top_population.values())
        for seqs, fitness in top_population.items():
            seqs = seqs.split()
            
            for seq, amount in seqs:
                for pos in range(len(seq)):
                    if seq[pos].startswith("[^"):
                        s = seq[pos].removeprefix("[^").removesuffix("]")
                        s = s.split("|")
                        for aa in s:
                            d[pos][aa] -= amount * abs(fitness - max_fitness)
                    elif seq[pos].startswith("["):
                        s = seq[pos].removeprefix("[").removesuffix("]")
                        s = s.split("|")
                        for aa in s:
                            d[pos][aa] += amount * abs(fitness - max_fitness)
        df = pd.DataFrame(d).T
        df = df.apply(abs)
        row_sums = df.sum(axis=1) + 1e-5
        df = df.div(row_sums, axis=0)
        ginis = []
        for _, row in df.iterrows():
            gini = np.sum(row**2)
            ginis.append(gini)
        df = pd.DataFrame(d).T.mul(ginis, axis=0)
        logomaker.Logo(df, fade_below=0.75)
