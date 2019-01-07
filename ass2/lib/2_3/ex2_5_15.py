import numpy as np
import pickle
from ex_2_3_tree_helper import Tree, Node
from ex_2_3 import load_sample, load_params

def compute_posterior(tree):
    root = tree.root
    compute_s(root)
    compute_t(root)


def compute_factor_s(child, father_value):
    K = len(child.cat[0])
    results = [0] * K

    for i in range(K):
        results[i] = child.cat[father_value][i] * child.s_i[i]

    return np.sum(results)


def compute_s(node):
    # if this is a leaf
    K = len(node.cat[0])
    s = [0] * K

    if len(node.descendants) == 0:
        s[node.sample] = 1
        # stop recursion
        node.s_i = s
        return

    # recurse to compute s in subtrees
    for child in node.descendants:
        compute_s(child)

    # once s is computed then compute it
    # for the current node.

    for i in range(K):
        children_sub = []
        for child in node.descendants:
            s_uv_i = compute_factor_s(child, i)
            children_sub.append(s_uv_i)
        # update s distribution

        s[i] = np.exp(np.sum(np.log(children_sub)))  # np.prod(children_sub)

    # save in the tree structure
    node.s_i = s
    return


def compute_siblings_s(node, father, j):
    result = 1.0
    for c in range(len(father.descendants)):
        child = father.descendants[c]
        if child == node:
            continue
        child_res = 0.0
        for k in range(len(child.cat[0])):
            child_res = child_res + child.cat[j][k] * child.s_i[k]
        result = result * child_res
    return result


def compute_factor_t(node):
    father = node.ancestor
    K = len(node.cat[0])
    F = len(father.cat[0])
    results = [0] * K

    for i in range(K):
        for j in range(F):
            s_siblings = compute_siblings_s(node, father, j)
            results[i] = results[i] + node.cat[j][i] * s_siblings * father.t[j]
    return results


def compute_t(node):
    # if the node is the root
    if node.ancestor is None:
        node.t = list(node.cat[0])
    else:
        # if it's a non root node
        node.t = compute_factor_t(node)

    node.posterior = np.multiply(node.s_i, node.t)
    # normalize the distribution
    node.posterior = node.posterior / np.sum(node.posterior)
    for child in node.descendants:
        compute_t(child)


def compute_node_posterior(node):
    # if node is the root
    if node.ancestor is None:
        node.posterior = np.multiply(node.s_i, node.cat[0])
        node.posterior = node.posterior / np.sum(node.posterior)
    else:
        # if it is a inner node
        node.posterior = np.multiply(node.cat[node.ancestor.sample], node.s_i)
        node.posterior = node.posterior / np.sum(node.posterior)



def sample_node(node):
    # leaves
    if len(node.descendants) == 0:
        return

    compute_node_posterior(node)
    node.sample = np.random.choice(len(node.posterior), p=node.posterior)

    for child in node.descendants:
        sample_node(child)


def sample_tree(tree):
    compute_s(tree.root)
    sample_node(tree.root)
    print("8======D sampled the tree")


def custom_tree():
    tree = Tree()
    tree.root = Node("1", [[0.5, 0.5]])
    child1 = Node("2", [[0.9, 0.1], [0.1, 0.9]])  # [[0.15, 0.85], [0.1, 0.9]])
    child2 = Node("3", [[0.9, 0.1], [0.1, 0.9]])  # [[0.5, 0.5], [0.1, 0.9]])
    child3 = Node("4", [[0.9, 0.1], [0.1, 0.9]])  # [[0.3, 0.7], [0.4, 0.6]])
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


def compute_leaves_probability(tree):
    root = tree.root
    total_sum = []
    for i in range(len(root.cat[0])):
        total_sum.append(root.s_i[i] * root.cat[0][i])

    return np.sum(total_sum)


def compute_observation_prob(params, samples):
    # choose a set of param
    for key in params.keys():
        print("#### TREE : " + key)
        parameters = params[key]
        for j in [1, 2, 3, 4, 5]:
            t = Tree()
            t.load_params(parameters)
            samples_name = key + '_sample_' + str(j)
            try:
                sample = samples[samples_name]
            except:
                break
            t.load_sample(sample)
            compute_posterior(t)
            print("\tsample_" + str(j) + " : " + str(compute_leaves_probability(t)))


def sample_given_tree(root):
    # choose a set of param
    tree = Tree()
    tree.root = root
    sample_tree(tree)
    tree.print_tree(True, True, True)


if __name__ == "__main__":
    t = Tree()

    my_data_path = ''

    # with open(my_data_path + 'tree_params.pickle', 'rb') as handle:
    #     params = pickle.load(handle)
    #
    # with open(my_data_path + 'tree_samples.pickle', 'rb') as handle:
    #     samples = pickle.load(handle)

    # get data to load into the tree
    params = np.load(my_data_path + 'tree_params.npy').tolist()
    samples = np.load(my_data_path + 'tree_samples.npy').tolist()

    if False:
        compute_observation_prob(params, samples)
    elif True:
        with open('../2_5/tree_proto2.pkl', 'rb') as handle2:
            samples2 = pickle.load(handle2)
        with open('../2_5/tree_proto2_cpd.pkl', 'rb') as handle2:
            params2 = pickle.load(handle2)

        root2 = load_params(params2)
        load_sample(root2, samples2)
        sample_given_tree(root2)
        exit(0)
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

    print ("### custom tree")
    #t = custom_tree()
    t.print_tree(True, True)

    # compute_posterior(t)
    sample_tree(t)
    print("root : " + str(t.root.posterior) + "\t|\tsample : " + str(t.root.sample))
    for d in t.root.descendants:
        print(d.name + " : " + str(d.posterior))

    print(compute_leaves_probability(t))
