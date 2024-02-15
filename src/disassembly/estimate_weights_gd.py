import networkx as nx
import numpy as np
import pandas as pd
from disassembly.util import KL


class WeightEstimatorGD:
    """
    Class to estimate weights using gradient descent.

    ```
    wegd = WeightEstimatorGD(lr, n_iterations, lam)
    generated_graph = wegd(true_dict, verbose=True)
    ```
    """

    def __init__(
        self,
        lr: float,
        n_iterations: int,
        lam1: float = 0,
        lam2: float = 0,
    ) -> None:
        self.lr = lr
        self.n_iterations = n_iterations
        self.lam1 = lam1
        self.lam2 = lam2

    def run(self, true_dict: dict, verbose: bool, parameters: dict = None):
        self.true_dict = true_dict
        self.keys = list(true_dict.keys())
        self.true_dict_vals = list(true_dict.values())
        if parameters:
            self.parameters = parameters
            self.graph = (
                self.create_graph_from_parameters()
            )  # creates the graph from params
        else:
            self.graph = self.create_graph()

        self.generated = {}
        self.losses = []
        self.weights = {}

        for iteration in range(self.n_iterations):
            guess, guess_df = self.generate_output(self.graph)
            self.generated[iteration] = guess
            self.weights[iteration] = np.array(
                [data["weight"] for _, _, data in self.graph.edges(data=True)]
            )
            # Compute loss
            kl = KL(self.true_dict_vals, guess.values()) + KL(guess.values(), self.true_dict_vals) 
            reg = get_l1(self.graph) * self.lam1
            reg += get_l2(self.graph) * self.lam2
            loss = kl + reg
            self.losses.append(loss)

            if verbose:
                print(
                    f"\r {iteration} / {self.n_iterations} | {loss:.2f}, kl: {kl:.2f}, reg: {reg:.2f}  | nz: { np.sum( self.weights[iteration] > 0.0 )} | ",
                    end="",
                    flush=True,
                )

            # Compute gradient
            dp_dw = self.compute_dp_dw(guess_df)
            dL_dp = self.compute_dL_dp(self.true_dict_vals, list(guess.values()))
            gradient = self.compute_dL_dw(dL_dp, dp_dw)
            grad_reg = self.get_grad_reg(self.graph)

            # Update graph
            self.graph = self.update_weights(gradient, grad_reg)

            if loss < 0.01:
                break

        return self.graph

    def generate_output(self, graph: nx.DiGraph):
        """
        Generates an output dict from a graph
        """
        longest_key = sorted(self.keys, key=len)[-1]
        p_generated = {}
        terminal_nodes = [node for node in graph.nodes() if graph.out_degree(node) == 0]

        for node in terminal_nodes:  # one hot terminal nodes
            oh_node = create_one_hot(self.keys, node)
            p_generated[node] = oh_node

        out_edges = {
            source: [
                target for _, target in graph.out_edges(source) if source != target
            ]
            for source in graph.nodes()
        }

        while len(p_generated.keys()) < len(self.keys):
            solvables = get_solvable(out_edges, p_generated)
            for solvable in solvables:
                p_generated[solvable] = np.zeros(len(self.keys))

                for source, target in graph.out_edges(solvable):
                    p_target = p_generated[target]
                    w_source_target = graph[source][target]["weight"]
                    p_generated[source] += w_source_target * p_target

                w_source_target = 1 - sum(
                    [
                        data["weight"]
                        for _, _, data in graph.out_edges(source, data=True)
                    ]
                )
                p_target = create_one_hot(self.keys, source)
                p_generated[source] += w_source_target * p_target

        guess = {
            self.keys[i]: p_generated[longest_key][i] for i in range(len(self.keys))
        }
        return guess, pd.DataFrame(p_generated, index=self.keys)

    def create_graph(self):
        graph = nx.DiGraph()
        graph.add_nodes_from([(k, {"layer": len(k)}) for k in self.keys])
        for key1 in self.keys:
            for key2 in self.keys:
                if (key1 in key2) and (key1 != key2):
                    graph.add_edge(key2, key1, weight=np.random.uniform(0, 1))
        # normalize
        for node in graph.nodes():
            out_edges = graph.out_edges(node, data=True)
            total_out = sum(
                [data["weight"] for _, _, data in graph.out_edges(node, data=True)]
            )
            for key, target, data in out_edges:
                nx.set_edge_attributes(
                    graph,
                    {(key, target): {"weight": 0.75 * data["weight"] / total_out}},
                )
        return graph

    def create_graph_from_parameters(self):
        graph = nx.DiGraph()
        graph.add_nodes_from([(k, {"layer": len(k)}) for k in self.keys])
        for key1 in self.keys:
            for key2 in self.keys:
                if (key1 in key2) and (key1 != key2):
                    if len(key1) == len(key2) - 1:
                        p1_left = None
                        p1_right = None
                        w = self.parameters["exo"]
                    elif key2.startswith(key1):
                        p1_left = None
                        p1_right = key1[-1]
                        w = self.parameters["endo"][p1_right]
                    elif key2.endswith(key1):
                        p1_left = key2[-len(key1) - 1]
                        p1_right = None
                        w = self.parameters["endo"][p1_left]
                    else:  # middle
                        p1_left = key2[key2.find(key1) - 1]
                        p1_right = key1[-1]
                        w = (
                            self.parameters["endo"][p1_left]
                            * self.parameters["endo"][p1_right]
                        ) ** 0.5
                  #  if w > 0.01:
                    graph.add_edge(
                        key2,
                        key1,
                        weight=w,
                    )
                        
        # normalize
        for node in graph.nodes():
            out_edges = graph.out_edges(node, data=True)
            total_out = sum(
                [data["weight"] for _, _, data in graph.out_edges(node, data=True)]
            )
            for key, target, data in out_edges:
                if total_out == 0:
                    w = 0
                else:
                    w = 0.75 * data["weight"] / total_out
                nx.set_edge_attributes(
                    graph,
                    {(key, target): {"weight": w}},
                )
        return graph

    def compute_dp_dw(self, guess_df: pd.DataFrame) -> dict:
        """
        dP / dw
        Change of P based on w
        Sx1 vector
        """
        longest_key = sorted(self.keys, key=len)[-1]
        prob_traversed = {key: 0 for key in self.keys}
        prob_traversed[longest_key] = 1

        for sequence, n in prob_traversed.items():
            out_edges = [
                (source, target, data)
                for source, target, data in self.graph.out_edges(sequence, data=True)
            ]
            weights = np.array([weight["weight"] for _, _, weight in out_edges])
            edges_to = [edge_to for _, edge_to, _ in out_edges]
            for w, e in zip(weights, edges_to):
                prob_traversed[e] += w * n

        dp_dw = {}

        for key in self.keys:
            out_edges = self.graph.out_edges(key)
            for (
                source,
                target,
            ) in out_edges:  # P(longest to source) * (P(target) - onehot(source))
                dp_dw[(source, target)] = prob_traversed[source] * (
                    guess_df[target].values - create_one_hot(self.keys, source)
                )

        return dp_dw

    def update_weights(self, grad, grad_reg=None) -> nx.DiGraph:
        diffs = {}
        k = 1

        old_graph = nx.DiGraph()
        for source, target, data in self.graph.edges(data=True):
            old_graph.add_edge(source, target, weight=data["weight"])

        old_loss = self.losses[-1]

        for source in self.graph.nodes():
            sum_old_weight = sum(
                [
                    data["weight"]
                    for _, _, data in self.graph.out_edges(source, data=True)
                ]
            )
            sum_diffs = 0

            for source, target in self.graph.out_edges(source):
                old_weight = self.graph[source][target]["weight"]
                grad_weight = grad[(source, target)]

                if grad_reg:  # if we regularize
                    grad_weight += grad_reg[(source, target)]

                new_weight = max(0, old_weight - self.lr * grad_weight)
                diff = new_weight - old_weight  # diff is -lr*grad
                sum_diffs += diff
                diffs[(source, target)] = diff

            while (sum_old_weight + k * sum_diffs) >= 1:
                k = k / 2

        new_graph = self.graph.copy()

        while True:
            # Update graph

            for source, target in new_graph.edges():

                
                nx.set_edge_attributes(
                    new_graph,
                    {
                        (source, target): {
                            "weight": max(
                                0,
                                old_graph[source][target]["weight"]
                                + diffs[(source, target)] * k,
                            )
                        }
                    },
                )



            # Get new KL
            new_guess, _ = self.generate_output(new_graph)

            new_loss = KL(self.true_dict_vals, list(new_guess.values())) + KL(list(new_guess.values()), self.true_dict_vals) + (
                get_l2(new_graph) * self.lam2 + get_l1(new_graph) * self.lam1
            )

            if new_loss <= old_loss:
                return new_graph

            if k < 1e-8:
                return old_graph

            k = k / 2
            new_graph = self.graph  # resets the new_graph to graph

    def compute_dL_dp(self, true, guess):
        return -np.array(true) / (np.array(guess) + 1e-8)

    def compute_dL_dw(self, dL_dp, dp_dw):
        """
        Gradient
        """
        dL_dw = {}
        for edge, val in dp_dw.items():
            dL_dw[edge] = np.sum(val * dL_dp)
        return dL_dw

    def get_grad_reg(self, graph):
        grad_reg = {}
        for source in graph.nodes():
            for _, target, data in graph.out_edges(source, data=True):
                grad_reg[(source, target)] = 2 * data["weight"] * self.lam2 + self.lam1
        return grad_reg

    def drop_weights(self, threshold: float = 0.01):
        """
        Idea, drop edges that are very small
        """
        new_graph = nx.DiGraph()
        for source, target, data in self.graph.edges(data=True):
            if data["weight"] > threshold:
                new_graph.add_edge(source, target, weight=data["weight"])
        return new_graph


def create_one_hot(keys, key):
    one_hot = np.zeros(len(keys))
    one_hot[keys.index(key)] = 1
    return one_hot


def get_solvable(out_edges, p_generated):
    solvable = []
    for source, targets in out_edges.items():
        if (
            set(targets).issubset(set((p_generated.keys())))
            and source not in p_generated.keys()
        ):
            solvable.append(source)
    return solvable


#


def get_l1(graph):
    return sum([abs(data["weight"]) for _, _, data in graph.edges(data=True)])


def get_l2(graph):
    return sum([data["weight"] ** 2 for _, _, data in graph.edges(data=True)])


def get_elastic_net(graph, lambda_1, lambda_2):
    return sum(
        [
            (lambda_1 * abs(data["weight"]))  # L1
            + (lambda_2 * data["weight"] ** 2)  # L2
            for _, _, data in graph.edges(data=True)
        ]
    )
