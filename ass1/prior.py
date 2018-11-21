import itertools as it
from math import exp , sqrt , pi
from scipy.stats import logistic
import numpy as np
import matplotlib.pyplot as plt




def create_dataset():
    d = np.matrix(list(it.product([-1 , 1], repeat=9)))
    x = np.matrix(list(it.product([-1 ,0, 1], repeat=2)))
    return (x,d)

def view_dataset(x,y):
    fig = plt.figure()

    marker_mapping = {-1 : 's',1: '^'}
    color_mapping = {-1 : 'b', 1: 'r'}

    for i in range(9):

        plt.scatter(x[i,0],x[i,1], marker=marker_mapping[y[0,i]], color=color_mapping[y[0,i]], s=60**2)


    plt.show()


def likelyhood_M0():
    return (1./2.)**9

def likelyhood_M1(x,y,w):
    result = 1.0
    for i in range(9):
        result = result * logistic(y[0,i]*(w[0]*x[i,0]))

    return result

def likelyhood_M2(x,y,w):
    result = 1.0
    for i in range(9):
        result = result * logistic(y[0,i]*(w[0]*x[i,0]+w[1]*x[i,1]))

    return result

def likelyhood_M3(x,y,w):
    result = 1.0
    for i in range(9):
        result = result * logistic(y[0,i]*(w[0]*x[i,0]+w[1]*x[i,1]+w[2]))

    return result

def sample_prior(size):
    if size == 0:
        return 1
    return np.random.normal(0,10,size)

def montecarlo_integral(likelihood, weight_lenght, iteration=1e8):
    result = 0

    for i in range(1e8):
        w = sample_prior(weight_lenght)
        result = result + likelihood(w) * # here compute prior value 



def compute_evidence(x,y):

    models = [likelyhood_M0,likelyhood_M1,likelyhood_M2,likelyhood_M3]

    datasets_n = y.shape[0]
    models_n = len(models)

    result = np.zeroes(datasets_n,models_n)

    for i,j in enumerate(models):
        for d in range(datasets_n):
            result[d, i] = montecarlo_integral(i,j)


if __name__ == '__main__':
    x,y = create_dataset()
    print(x)
    print(y.shape)
    view_dataset(x, y[57,:])
    print(sample_prior(1))
    print(np.var(sample_prior(10000)))
