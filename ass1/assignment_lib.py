import numpy as np
from math import sqrt
from scipy.spatial.distance import cdist, pdist

class GP_prior:

    def __init__(self, sigma=1, lenghtscale=1):

        self.sigma = sigma
        self.lenghtscale = lenghtscale
        self.covariance = None
        self.mean = None

    def exp_kernel(self, X, Y):
        X = X[:,None]
        Y = Y[:,None]
        return (self.sigma**2) * np.exp(-cdist(X, Y, 'sqeuclidean')/(self.lenghtscale*self.lenghtscale))
        #x = np.matrix(X)
        #x = np.tile(x.transpose(),(1,Y.size))
        #m = np.subtract(x,Y)
        #M = np.multiply(m,m)
        #e = 1.0/(self.lenghtscale**2)
        #return (self.sigma**2) * np.exp(-e * M)


    def sample_function(self,inputs):
        # compute covariance
        self.covariance = self.exp_kernel(inputs, inputs)
        self.mean = np.zeros(inputs.shape)
        y = self.sample()
        return y


    def sample(self):
        if self.covariance is None or self.mean is None:
            return None
        return np.random.multivariate_normal(self.mean, self.covariance)


class GP_posterior:

    def __init__(self, var = 1,sigma=1, lenghtscale=1,diag_covar=0):
        self.err_dev = sqrt(var)
        self.err_var = var
        self.sigma = sigma
        self.lenghtscale = lenghtscale
        self.diag_covariance = diag_covar
        self.data_x = None
        self.data_t = None


    def compute_prior(self):
        self.prior = GP_prior(self.sigma,self.lenghtscale)
        #self.prior.sample_function(self.data_x)
        self.prior.covariance = self.compute_kernel(self.data_x,self.data_x) + self.err_var * np.eye(self.data_x.size)
        self.prior.mean = np.zeros(self.data_x.size)
        self.prior_covariance_inv = np.linalg.inv(self.prior.covariance)


    def compute_kernel(self,X,Y=None):
        if Y is None:
            Y = self.data_x

        X = X[:,None]
        Y = Y[:,None]
        return (self.sigma**2) * np.exp(-cdist(X, Y, 'sqeuclidean')/(self.lenghtscale*self.lenghtscale)) + self.diag_covariance * np.eye(X.size, Y.size)

        #x = np.matrix(X)
        #x = np.tile(x.transpose(),(1,Y.size))
        #m = np.subtract(x,Y)

        #M = np.multiply(m,m)
        #e = 1.0/(self.lenghtscale**2)
        #return (self.sigma**2) * np.exp(-e * M) + self.diag_covariance * np.eye(X.size, Y.size)



    def sample_predictive_mean(self,x):
        K = self.compute_kernel(x,self.data_x)
        T = np.matrix(self.data_t)
        R = np.dot(K,np.dot(self.prior_covariance_inv,T.T))
        return R

    def sample_predictive_variance_errorfree(self,x):
        K = self.compute_kernel(x,self.data_x)
        C = self.compute_kernel(x,x) #+ self.err_var * np.eye(x.size)
        return C - np.dot(K,np.dot(self.prior_covariance_inv,K.T))


    def sample_predictive_variance(self,x):
        K = self.compute_kernel(x,self.data_x)
        C = self.compute_kernel(x,x) + self.err_var * np.eye(x.size)
        return C - np.dot(K,np.dot(self.prior_covariance_inv,K.T))

    def generate_data_f(self,det=True):
        if det:
            self.data_x = np.linspace(-4,6,9)
        else:
            self.data_x = np.random.uniform(-4,6,9)

        self.data_t = (2+np.power(0.5*self.data_x-1,2)) * np.sin(3*self.data_x) + np.random.normal(0, self.err_dev, len(self.data_x))
        self.compute_prior()
        return (self.data_x, self.data_t)
