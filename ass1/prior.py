import itertools as it
from math import exp , sqrt , pi
from scipy.stats import logistic
import numpy as np
import matplotlib.pyplot as plt
from threading import Thread

class Evidence(Thread):
    def __init__(self, likelihood,model_n,x,y, result):
        Thread.__init__(self)
        self.x = x
        self.y = y
        self.likelihood = likelihood
        self.model_n = model_n
        self.result = result

    def run(self):
        d = 0
        datasets_n = self.y.shape[0]
        while d < datasets_n:
            if d % 10 == 0 :
                print("\t" + str(self.model_n) + ":" + str(d))

            self.result[d, self.model_n] = montecarlo_integral(self.x,self.y[d,:],self.likelihood,self.model_n,1000)
            d = d + 1




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


def likelyhood_M0(x,y,w):
    return (1./2.)**9

def likelyhood_M1(x,y,w):
    result = 1.0

    for i in range(9):
        result = result * logistic.cdf(y[0,i]*(w[0]*x[i,0]))

    return result

def likelyhood_M2(x,y,w):
    result = 1.0
    for i in range(9):
        result = result * logistic.cdf(y[0,i]*(w[0]*x[i,0]+w[1]*x[i,1]))

    return result

def likelyhood_M3(x,y,w):
    result = 1.0
    for i in range(9):
        result = result * logistic.cdf(y[0,i]*(w[0]*x[i,0]+w[1]*x[i,1]+w[2]))

    return result

def sample_prior(size):
    if size == 0:
        return 1
    return np.random.normal(0,10,size)

def montecarlo_integral(x,y,likelihood, weight_lenght, iteration=1e8):
    result = 0

    i = 0
    while i < iteration:
        w = sample_prior(weight_lenght)
        result = result + likelihood(x,y,w)
        i = i + 1
    return result


def compute_evidence(x,y):

    models = [likelyhood_M0,likelyhood_M1,likelyhood_M2,likelyhood_M3]

    datasets_n = y.shape[0]
    print("dataset of size " + str(datasets_n))
    models_n = len(models)
    result_size = (datasets_n,models_n)
    result = np.zeros(result_size)

    for i,j in enumerate(models):
        print("Calculating for model " + str(i))
        for d in range(datasets_n):
            print("\t" + str(i) + ":" + str(d))
            result[d, i] = montecarlo_integral(x,y[d,:],j,i,10)
    return result

def compute_evidence_multithread(x,y):

    models_n = 4
    datasets_n = y.shape[0]
    result_size = (datasets_n,models_n)
    result = np.zeros(result_size)

    threads = [Evidence(likelyhood_M0,0,x.copy(),y.copy(),result),Evidence(likelyhood_M1,1,x.copy(),y.copy(),result),Evidence(likelyhood_M2,2,x.copy(),y.copy(),result),Evidence(likelyhood_M3,3,x.copy(),y.copy(),result)]

    for t in threads:
        t.start()
        print("Started")

    for t in threads:
        t.join()
        print("Joined")

    return result



if __name__ == '__main__':
    x,y = create_dataset()
    print(x)
    print(y.shape)
    view_dataset(x, y[57,:])
    print(sample_prior(1))
    print(np.var(sample_prior(10000)))
    if True:
        evidence = compute_evidence_multithread(x,y)
    else:
        evidence = compute_evidence(x,y)
    np.save("evidence.mat", evidence)
