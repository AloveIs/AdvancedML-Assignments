import numpy as np

class BayesLinearRegression:
	def __init__(self, prior_mean,prior_cov, precision):
		self.prior_mean = prior_mean
		self.prior_cov = prior_cov

		self.posterior_mean = prior_mean
		self.posterior_cov = prior_cov

		self.precision = precision

	def sample_posterior():
		 return []

	def put_data(self,x,y):
		self.posterior_cov = np.linalg.inv(np.linalg.inv(self.prior_cov) + self.precision* np.outer(x,x))
		self.posterior_mean = np.dot(self.posterior_cov, ( np.dot(np.linalg.inv(self.prior_cov),self.prior_mean)) + (self.precision*x*y) )
		
		self.prior_cov = self.posterior_cov
		self.prior_mean = self.posterior_mean
