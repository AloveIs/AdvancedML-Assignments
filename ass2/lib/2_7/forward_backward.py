import numpy as np
from generator import define_HMMs, generate_data, nr_vehicles, nr_classes, nr_rows, nr_columns, class_prob, \
    start_prob, transition_prob, emission_prob, targets, data


def get_transition(r1, m1, r2, param):
    # calculate A((r1, m1), (r2, m1+1)) (for test purpose we set below)
    return param[r1][r2]


def get_emission(r, o, param):
    # calculate O(m, o) (for test purpose we set below)
    return param[r][o]


def get_init(param):
    # provide an array containing the initial state probability having size R (for test purpose we set below)
    return np.array(param), nr_rows


def forward(get_init, get_transition, get_emission, observations, trans, emiss, start_prob):
    pi, R = get_init(start_prob)
    M = len(observations)
    alpha = np.zeros((M, R))

    # base case
    O = []
    for r in range(R):
        O.append(get_emission(r, observations[0], emiss))
    alpha[0, :] = pi * O[:]

    # recursive case
    for m in range(1, M):
        for r2 in range(R):
            for r1 in range(R):
                transition = get_transition(r1, m, r2, trans)
                emission = get_emission(r2, observations[m], emiss)
                alpha[m, r2] += alpha[m - 1, r1] * transition * emission

    return (alpha, np.sum(alpha[M - 1, :]))


def backward(get_init, get_transition, get_emission, observations, trans, emiss, start_prob):
    pi, R = get_init(start_prob)
    M = len(observations)
    beta = np.zeros((M, R))

    # base case
    beta[M - 1, :] = 1

    # recursive case
    for m in range(M - 2, -1, -1):
        for r1 in range(R):
            for r2 in range(R):
                transition = get_transition(r1, m, r2, trans)
                emission = get_emission(r2, observations[m + 1], emiss)
                beta[m, r1] += beta[m + 1, r2] * transition * emission

    O = []
    for r in range(R):
        O.append(get_emission(r, observations[0], emiss))

    return beta, np.sum(pi * O[:] * beta[0, :])


def generate_parameters(observations):
    rand_class_prob = np.random.dirichlet([1] * nr_classes)

    rand_start_prob = np.matrix([[1.0 / float(nr_rows)] * nr_rows]*nr_classes)

    rand_transition_prob = np.zeros((nr_classes, nr_rows, nr_rows))
    rand_emission_prob = np.zeros((nr_classes, nr_rows, nr_rows))

    for c in range(nr_classes):
        for r in range(nr_rows):
            hyperpar = [1] * nr_rows
            rand_transition_prob[c, r, :] = np.random.dirichlet(hyperpar)
            hyperpar[r] = nr_rows
            rand_emission_prob[c, r, :] = np.random.dirichlet(hyperpar)

    return rand_class_prob, rand_start_prob, rand_transition_prob, rand_emission_prob


def compute_class_responsibility(params, obs):
    _class_prob, _start_prob, _transition_prob, _emission_prob = params
    # matrix C (classes) x N (samples)
    class_resp = np.zeros([nr_classes, nr_vehicles])

    for c in range(nr_classes):
        for n in range(nr_vehicles):
            o = obs[n, :]
            b, L = forward(get_init, get_transition, get_emission, o,
                           _transition_prob[c, :, :], _emission_prob[c, :, :], _start_prob[c,:])
            class_resp[c, n] = L * _class_prob[c]

    class_resp = class_resp / np.sum(class_resp, axis=0)

    return class_resp


def compute_gamma_resp(params, obs):
    _class_prob, _start_prob, _transition_prob, _emission_prob = params

    gamma = dict()

    for n in range(nr_vehicles):
        gamma[n] = np.zeros([nr_columns, nr_rows, nr_classes])

        for c in range(nr_classes):
            # compute responsibility for column
            # m of class c of the observation of
            # car n
            o = obs[n, :]
            alpha, A = forward(get_init, get_transition, get_emission, o,
                               _transition_prob[c, :, :], _emission_prob[c, :, :],_start_prob[c,:])
            beta, B = backward(get_init, get_transition, get_emission, o,
                               _transition_prob[c, :, :], _emission_prob[c, :, :],_start_prob[c,:])

            for m in range(nr_columns):
                gamma[n][m, :, c] = np.multiply(alpha[m, :], beta[m, :])
                # normalize
                gamma[n][m, :, c] = gamma[n][m, :, c] / np.sum(gamma[n][m, :, c])

    return gamma


def compute_xi_resp(params, obs):
    # size : dict(car){dic(class){[row x row x columns]}}
    _class_prob, _start_prob, _transition_prob, _emission_prob = params
    xi = dict()
    for n in range(nr_vehicles):
        xi[n] = dict()
        o = obs[n, :]
        for c in range(nr_classes):
            xi[n][c] = np.zeros([nr_rows, nr_rows, nr_columns])

            alpha, A = forward(get_init, get_transition, get_emission, o,
                               _transition_prob[c, :, :], _emission_prob[c, :, :],_start_prob[c,:])
            beta, B = backward(get_init, get_transition, get_emission, o,
                               _transition_prob[c, :, :], _emission_prob[c, :, :],_start_prob[c,:])

            transition_v = _transition_prob[c, :, :]
            # now compute through time
            for m in range(1, nr_columns):
                emission_v = _emission_prob[c, :, o[m]]
                temp = np.multiply(emission_v, np.transpose(beta[m, :]))
                xi[n][c][:, :, m] = np.multiply(transition_v, np.outer(alpha[m - 1, :], temp))
    return xi


def E_step(old_parameters, obs):
    # compute class resp for HMM (nu in assignment)
    # size : classes x cars
    class_resp = compute_class_responsibility(old_parameters, obs)

    # compute gamma values
    # size : dict(car) = column x row_value x classes
    gamma_resp = compute_gamma_resp(old_parameters, obs)

    # compute xi
    # size : dict(car){dic(class){[row x row x columns]}}
    xi_resp = compute_xi_resp(old_parameters, obs)

    return class_resp, gamma_resp, xi_resp


def M_step(responsibilities, obs):
    class_resp, gamma_resp, xi_resp = responsibilities

    # new class prob
    new_class_prob = np.sum(class_resp, axis=1) / np.sum(class_resp)

    # transition prob
    new_transition_prob = np.zeros([nr_classes, nr_rows, nr_rows])
    for c in range(nr_classes):
        for n in range(nr_vehicles):
            new_transition_prob[c, :, :] += class_resp[c, n] * np.sum(xi_resp[n][c], axis=2)
        # normalization
        for r in range(nr_rows):
            new_transition_prob[c, r, :] = new_transition_prob[c, r, :] / (np.sum(new_transition_prob[c, r, :]))

    # emission prob
    new_emission_prob = np.zeros([nr_classes, nr_rows, nr_rows])

    for c in range(nr_classes):
        for n in range(nr_vehicles):

            resp = gamma_resp[n][:, :, c]
            for m in range(nr_columns):
                o = obs[n, m]
                for r in range(nr_rows):
                    new_emission_prob[c, :, r] += class_resp[c, n] * resp[m, :] * (o == r)
        # normalization
        for r in range(nr_rows):
            new_emission_prob[c, r, :] = new_emission_prob[c, r, :] / (np.sum(new_emission_prob[c, r, :]))

    new_start_prob = np.zeros([nr_classes,nr_rows])

    for n in range(nr_vehicles):
        for c in range(nr_classes):
            new_start_prob[c, :] += class_resp[c, n] * gamma_resp[n][1, :, c]

    # normalizing
    for r in range(nr_classes):
        new_start_prob[r,:] = new_start_prob[r,:] / np.sum(new_start_prob[r,:])

    return new_class_prob, new_start_prob, new_transition_prob, new_emission_prob


def EM(obs, threshold=0.1):
    # (means, precision)
    old_parameters = generate_parameters(obs)
    print("Computed parameters:")
    for e in old_parameters:
        print(e)
    converged = False
    rounds = 0

    while not converged:
        responsibilities = E_step(old_parameters, obs)
        parameters = M_step(responsibilities, obs)
        rounds += 1

        distance = 0
        for i in range(len(parameters)):
            distance += np.linalg.norm(parameters[i] - old_parameters[i])

        if np.isnan(distance):
            raise Exception("division by zero somewhere")

        print(str(rounds) + ") distance : " + str(distance))
        if distance < threshold or rounds > 15:
            converged = True
        else:
            old_parameters = parameters
    return parameters


def ML_assignment(params, obs):
    _class_prob, _start_prob, _transition_prob, _emission_prob = params
    M = np.zeros([nr_vehicles, nr_classes])

    for n in range(nr_vehicles):
        for c in range(nr_classes):
            a, likelihood = forward(get_init, get_transition, get_emission, obs[n, :],
                                    _transition_prob[c, :, :], _emission_prob[c, :, :],_start_prob[c,:])
            M[n, c] = _class_prob[c] * likelihood

    return np.argmax(M,axis=1)


if __name__ == "__main__":
    # test examples
    # alpha, A = forward(get_init, get_transition, get_emission, [0, 0, 1, 1, 1, 1])
    # print(backward(get_init, get_transition, get_emission, [0, 0, 1, 1, 1, 1]))
    flag = True
    while flag:
        try:
            print("#\n"*20)
            param = EM(data)
            flag = False
        except Exception:
            continue

    print("finished EM")
    e_class_prob, e_start_prob, e_transition_prob, e_emission_prob = param
    cluster = ML_assignment(param, data)
    print(cluster)
    print(targets)

    # print("Class probabilities\n" + str(class_prob))
    # print("Estimated Class probabilities\n" + str(e_class_prob))
    # print("\nStart probabilities\n" + str(start_prob))
    # print("\nEstimated Start probabilities\n" + str(e_start_prob))
    # print("\nTransition probabilities\n" + str(transition_prob))
    # print("\nEstimated Transition probabilities\n" + str(e_transition_prob))
    # print("\nEmission probabilities\n" + str(emission_prob))
    # print("\nEstimated Emission probabilities\n" + str(e_emission_prob))
