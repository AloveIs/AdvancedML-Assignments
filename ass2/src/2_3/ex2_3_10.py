import numpy as np
import pickle
from ex_2_3_tree_helper import Tree, Node


def compute_probability_observation(tree):
    root = tree.root
    root.marginal = root.cat[0]
    leaves_results = []

    # compute marginal distribution across the tree
    recursive_marginal(root, leaves_results)

    # multiply the found leaves probability

    probability = np.exp(np.sum(np.log(leaves_results))) # np.prod(leaves_results)

    return probability


def recursive_marginal(node, leaves_results):
    # if the node is a leaf
    if len(node.descendants) == 0:
        # add the probability into the list
        leaves_results.append(node.marginal[node.sample])
        return

    # recurse on children
    for child in node.descendants:
        compute_marginal(child)
        recursive_marginal(child, leaves_results)


def compute_marginal(node):
    K = len(node.cat)
    parent = node.ancestor
    node.marginal = []
    # values of the current node
    for j in range(K):
        marginal_p_j = []
        # values of the parent node
        for i in range(K):
            marginal_p_j.append(node.cat[i][j] * parent.marginal[i])
        node.marginal.append(np.sum(marginal_p_j))

def custom_tree():
    tree = Tree()
    tree.root = Node("1", [[0.9, 0.1]])
    child1 = Node("2", [[0.15, 0.85], [0.1, 0.9]])
    child2 = Node("3", [[0.5, 0.5], [0.1, 0.9]])
    child3 = Node("3", [[0.3, 0.7], [0.4, 0.6]])
    child1.sample = 0
    child2.sample = 0
    child3.sample = 0
    child1.ancestor = tree.root
    child2.ancestor = tree.root
    child3.ancestor = tree.root
    tree.root.descendants.append(child1)
    tree.root.descendants.append(child2)
    tree.root.descendants.append(child3)
    tree.root.sample = 1
    return tree


if __name__ == "__main__":
    t = Tree()

    my_data_path = ''

    # get data to load into the tree
    params = np.load(my_data_path + 'tree_params.npy').tolist()
    samples = np.load(my_data_path + 'tree_samples.npy').tolist()

    # choose a set of param
    params_name = params.keys()[0]
    parameters = params[params_name]

    # choose samples
    samples_name = params_name + '_sample_2'
    sample = samples[samples_name]

    # Load params into tree and samples
    t.load_params(parameters)
    t.load_sample(sample)
    print("### printing loaded tree")
    t.print_tree(True, True)

    t = custom_tree()
    t.print_tree(True, True)
    print(compute_probability_observation(t))
