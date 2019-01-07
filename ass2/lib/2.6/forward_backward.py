import numpy as np
import scipy
from scipy.stats import norm

from generator import generate_data, M, N, mu as original_means, beta as original_beta


def get_transition(r1, m1, r2):
    # calculate A((r1, m1), (r2, m1+1)) (for test purpose we set below)
    if r1 == r2:
        return 0.25
    else:
        return 0.75


def get_emission(r, o, m, means, std_dev):
    # calculate O((m,r), o) (for test purpose we set below)
    mean = means[r][m]

    result = norm.pdf(o, loc=mean, scale=np.sqrt(2) * std_dev)

    return result


def get_init():
    # provide an array containing the initial state probability having size R (for test purpose we set below)
    pi = np.array([0.5, 0.5])
    # number of rows
    R = pi.shape[0]
    return pi, R


def forward(get_init, get_transition, get_emission, observations, means, std_dev):
    pi, R = get_init()
    M = len(observations)
    alpha = np.zeros((M, R))

    # base case
    O = []
    for r in range(R):
        O.append(get_emission(r, observations[0], 0, means, std_dev))
    alpha[0, :] = pi * O[:]

    # recursive case
    for m in range(1, M):
        for r2 in range(R):
            for r1 in range(R):
                transition = get_transition(r1, m, r2)
                emission = get_emission(r2, observations[m], m, means, std_dev)
                alpha[m, r2] += alpha[m - 1, r1] * transition * emission

    return (alpha, np.sum(alpha[M - 1, :]))


def backward(get_init, get_transition, get_emission, observations, means, std_dev):
    pi, R = get_init()
    M = len(observations)
    beta = np.zeros((M, R))

    # base case
    beta[M - 1, :] = 1

    # recursive case
    for m in range(M - 2, -1, -1):
        for r1 in range(R):
            for r2 in range(R):
                transition = get_transition(r1, m, r2)
                emission = get_emission(r2, observations[m + 1], m + 1, means, std_dev)
                beta[m, r1] += beta[m + 1, r2] * transition * emission

    O = []
    for r in range(R):
        O.append(get_emission(r, observations[0], 0, means, std_dev))

    return beta, np.sum(pi * O[:] * beta[0, :])


def compute_obs_responsibilities(observation, mean, std_dev):
    alpha, A = forward(get_init, get_transition, get_emission, observation, mean, std_dev)
    beta, B = backward(get_init, get_transition, get_emission, observation, mean, std_dev)

    R = np.multiply(alpha, beta)
    # normalize
    R = R / np.sum(R, axis=1)[0]

    return R


def E_step(parameters, observations):
    means = parameters[0]
    std_dev = parameters[1]
    responsibilities = dict()

    for couple_p in observations.keys():
        responsibilities[couple_p] = []

        # summing up the means of the two players
        couple_means = [means[:, couple_p[0] - 1, 0] + means[:, couple_p[1] - 1, 0],
                        means[:, couple_p[0] - 1, 1] + means[:, couple_p[1] - 1, 1]]

        for instance in observations[couple_p].values():
            instance_resp = compute_obs_responsibilities(instance, couple_means, std_dev)
            responsibilities[couple_p].append(instance_resp)
    return responsibilities


def M_step(responsibility, observ):
    new_means = np.zeros((M, N, 2))
    new_variance = 0
    for m in range(M):
        for k in range(2):

            # compute system of equation
            linear_system = []
            linear_known = []
            for p in range(1, N + 1):
                coefficients = [0] * N
                known_term = 0

                # compute known term
                for pair in observ.keys():
                    if p not in pair:
                        continue
                    for r in observ[pair].keys():
                        known_term = known_term + observ[pair][r][m] * responsibility[pair][r - 1][m, k]

                    for r in observ[pair].keys():
                        p2 = pair[0] if p == pair[1] else pair[1]
                        coefficients[p - 1] += responsibility[pair][r - 1][m, k]
                        coefficients[p2 - 1] += responsibility[pair][r - 1][m, k]

                linear_system.append(coefficients)
                linear_known.append(known_term)

            param_m_k = np.linalg.solve(linear_system, linear_known)
            new_means[m, :, k] = param_m_k

    # compute new variance
    denominator = 0
    numerator = 0
    for pair in observation.keys():
        for m in range(M):
            for k in range(2):
                pair_mean = new_means[m, pair[0] - 1, k] + new_means[m, pair[1] - 1, k]
                for r in observation[pair].keys():
                    numerator += responsibility[pair][r - 1][m, k] * np.power(observation[pair][r][m] - pair_mean, 2)
                    denominator += responsibility[pair][r - 1][m, k]

    new_variance = numerator / (2.0 * denominator)

    return new_means, np.sqrt(new_variance)


def decompose_means(M, N, means):
    combinations = int(scipy.misc.comb(N, 2))
    contributions = np.zeros([combinations, N])

    for idx, keys in enumerate(means.keys(), 0):
        contributions[idx, keys[0] - 1] = 1
        contributions[idx, keys[1] - 1] = 1

    values = means.values()
    players_mean_s = []

    for i in range(M):
        results = []

        for j in range(combinations):
            results.append(values[j][i])

        players_mean_s.append(np.dot(np.linalg.pinv(contributions), results))

    return players_mean_s


def generate_parameters(M, N, obs):
    means = dict()
    variances = []
    for key, seq in obs.items():
        means[key] = np.mean(seq.values(), axis=0)
        variances.append(np.mean(np.var(seq.values(), axis=0)) / 2)

    starting_variance = np.mean(variances)
    starting_std_dev = np.sqrt(starting_variance)

    decomposed_m = np.matrix(decompose_means(M, N, means))
    starting_means = np.zeros([M, N, 2])
    starting_means[:, :, 0] = decomposed_m + np.random.normal(0, starting_std_dev)
    starting_means[:, :, 1] = decomposed_m + np.random.normal(0, starting_std_dev)

    return starting_means, starting_std_dev


def EM(obs, M, N, threshold=1):
    # (means, precision)
    old_parameters = generate_parameters(M, N, obs)
    converged = False
    rounds = 0

    while not converged:
        responsibilities = E_step(old_parameters, obs)
        parameters = M_step(responsibilities, obs)
        rounds += 1

        distance = np.linalg.norm(old_parameters[0][:, :, 0] - parameters[0][:, :, 0]) \
                   + np.linalg.norm(old_parameters[0][:, :, 1] - parameters[0][:, :, 1]) \
                   + np.linalg.norm(old_parameters[1] - parameters[0])
        print(str(rounds) + ") distance : " + str(distance))
        if distance < 1e-6 or rounds > 200:
            converged = True
        else:
            old_parameters = parameters
    return parameters


if __name__ == "__main__":
    # test examples
    # print(forward(get_init, get_transition, get_emission, [0, 0, 1, 1, 1, 1]))
    # print(backward(get_init, get_transition, get_emission, [0, 0, 1, 1, 1, 1]))

    observation = generate_data()
    estimated_parameters = EM(observation, M, N)
    print("Estimated parameters")
    print(estimated_parameters[0][:, :, 0])
    print(estimated_parameters[0][:, :, 1])
    print("Original parameters")
    print(original_means)

    print("##### variance")

    print(estimated_parameters[1])
    print(original_beta)
