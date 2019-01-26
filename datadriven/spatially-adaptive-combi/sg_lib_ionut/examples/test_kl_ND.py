import numpy as np
import chaospy as cp
from scipy.optimize import minimize
from itertools import product
from matplotlib.pyplot import *

def model(a, b, x):
	return np.exp(-b*x)*np.sin(a**2*b) + a*x + b 

def compute_evidence(func, quad_deg):
	nodes, weights 	= cp.generate_quadrature(quad_deg, mvar_prior, rule='G')
	res 			= np.sum([func(n[0], n[1])*w for n, w in zip(nodes.T, weights)])

	return res

def compute_mean(func, quad_deg):
	nodes, weights 	= cp.generate_quadrature(quad_deg, mvar_prior, rule='G')
	res 			= np.array([np.sum([n[0]*func(n[0], n[1])*w for n, w in zip(nodes.T, weights)]), \
						np.sum([n[1]*func(n[0], n[1])*w for n, w in zip(nodes.T, weights)])])

	return res

def compute_var(func, mean, quad_deg):
	nodes, weights 	= cp.generate_quadrature(quad_deg, mvar_prior, rule='G')
	second_moment 	= np.array([np.sum([n[0]**2*func(n[0], n[1])*w for n, w in zip(nodes.T, weights)]), \
						np.sum([n[1]**2*func(n[0], n[1])*w for n, w in zip(nodes.T, weights)])])

	res = second_moment - mean**2

	return res

def compute_cov(func, mean, quad_deg):
	nodes, weights 	= cp.generate_quadrature(quad_deg, mvar_prior, rule='G')
	res 			= np.sum([n[0]*n[1]*func(n[0], n[1])*w for n, w in zip(nodes.T, weights)]) - mean[0]*mean[1]

	return res

def multivar_dens(m1, m2, v1, cov, v2, x1, x2):
	multivar_mean 		= np.array([m1, m2])
	multivar_cov 		= np.array([[v1, cov], [cov, v2]])

	ll 				= 0.5*np.dot(np.dot((np.array([x1, x2]) - multivar_mean).T, np.linalg.inv(multivar_cov)),  \
														(np.array([x1, x2]) - multivar_mean))
	norm_constant 	= np.sqrt((2.0*np.pi)**2*np.linalg.det(multivar_cov))

	return np.exp(-ll)/norm_constant

def log_dens(m1, m2, v1, cov, v2, x1, x2):
	multivar_mean 		= np.array([m1, m2])
	multivar_cov 		= np.array([[v1, cov], [cov, v2]])

	ll 				= 0.5*np.dot(np.dot((np.array([x1, x2]) - multivar_mean).T, np.linalg.inv(multivar_cov)),  \
														(np.array([x1, x2]) - multivar_mean))
	norm_constant 	= np.sqrt((2.0*np.pi)**2*np.linalg.det(multivar_cov))

	return  ll + np.log(norm_constant)

def KL_div(x):
	m1_  	= x[0]
	m2_  	= x[1] 
	v1_ 	= x[2] 
	cov_ 	= x[3] 
	v2_ 	= x[4]

	res = np.sum(np.array([likelihood(n[0], n[1])*log_dens(m1_, m2_, v1_, cov_, v2_, n[0], n[1])*w \
									for n, w in zip(nodes.T, weights)]) )

	return res

def get_std_points(deg_1D=5):
	nodes, weights 	= np.polynomial.legendre.leggauss(deg_1D)
	std_nodes 		= np.array(list(product(nodes, repeat=2)))

	return std_nodes

def map_nodes(std_nodes, mean, E):
	new_nodes = np.array([mean + np.dot(E, node) for node in std_nodes])

	return new_nodes

if __name__ == '__main__':
	quad_deg 		= 30
	quad_deg_nodes 	= 5

	true_param_a 	= 0.4
	true_param_b 	= 0.6
	sigma_noise 	= 0.05
	loc_interest 	= np.arange(0.1, 1.0, 0.1)

	mvar_prior = cp.J(cp.Uniform(), cp.Uniform())
	
	prior 	= lambda a, b: 1.0
	data 	= np.array([model(true_param_a, true_param_b, loc) + np.random.normal(loc=0., scale=sigma_noise) for loc in loc_interest])

	ll 			= lambda a, b: np.sum([(model(a, b, loc) - d)**2 for loc, d in zip(loc_interest, data)])
	likelihood 	= lambda a, b: np.exp(-ll(a, b)/(2.0*sigma_noise**2))
	evidence 	= compute_evidence(likelihood, quad_deg)
	posterior 	= lambda a, b: prior(a, b)*likelihood(a, b)/evidence

	mean 	= compute_mean(posterior, quad_deg)
	var 	= compute_var(posterior, mean, quad_deg)
	cov 	= compute_cov(posterior, mean, quad_deg)

	nodes, weights 		= cp.generate_quadrature(quad_deg, mvar_prior, rule='G')
	x0 					= np.array([mean[0], mean[1], var[0], cov, var[1]])
	var_Bayes_params 	= minimize(KL_div, x0, method='nelder-mead').x

	param_cov_matrix 	= np.array([[var_Bayes_params[2], var_Bayes_params[3]], [var_Bayes_params[3], var_Bayes_params[4]]])
	param_var_matrix 	= np.array([[var_Bayes_params[2], 0.], [0., var_Bayes_params[4]]])
	param_mean 			= np.array([var_Bayes_params[0], var_Bayes_params[1]])
	L 					= np.linalg.cholesky(param_cov_matrix)
	V 					= np.linalg.cholesky(param_var_matrix)

	std_nodes 				= get_std_points(quad_deg_nodes)
	mapped_nodes_correct 	= map_nodes(std_nodes, param_mean, L)
	mapped_nodes_uncorr 	= map_nodes(std_nodes, param_mean, V)

	print mean[0], mean[1], var[0], cov, var[1]
	print var_Bayes_params[0], var_Bayes_params[1], var_Bayes_params[2], var_Bayes_params[3], var_Bayes_params[4]

	delta 	= 0.01
	xx 		= np.arange(0.0, 1.0 + delta, delta)
	yy 		= np.arange(0.0, 1.0 + delta, delta)

	X, Y 				= np.meshgrid(xx, yy)
	num_rows, num_cols 	= X.shape

	true_posterior_eval 	= np.zeros((num_rows, num_cols))
	approx_posterior_eval 	= np.zeros((num_rows, num_cols))

	for i in xrange(num_rows):
		for j in xrange(num_cols):
			true_posterior_eval[i, j] 	= posterior(X[i][j], Y[i][j])
			approx_posterior_eval[i, j]	= multivar_dens(var_Bayes_params[0], var_Bayes_params[1], var_Bayes_params[2], \
											var_Bayes_params[3], var_Bayes_params[4], X[i][j], Y[i][j])

	fig = figure()

	c1 = contour(X, Y, true_posterior_eval, cmap=cm.coolwarm, linewidths=2.0)
	clabel(c1, inline=0, fontsize=0)
	c1.collections[0].set_label('true posterior')

	c2 = contour(X, Y, approx_posterior_eval, cmap=cm.Spectral, linewidths=2.0)	
	clabel(c2, inline=0, fontsize=0)
	c2.collections[0].set_label('approx posterior variational Bayes')

	legend(loc='best', fontsize=20)

	
	fig = figure()

	subplot(1, 2, 1)

	c = contour(X, Y, approx_posterior_eval, cmap=cm.Spectral, linewidths=2.0)
	clabel(c, inline=0, fontsize=0)
	c.collections[0].set_label('approx posterior variational Bayes')
	plot(mapped_nodes_correct.T[0], mapped_nodes_correct.T[1], 'ok', markersize=10, label='"correlated" quad nodes')
	legend(loc='best', fontsize=20)
	
	
	subplot(1, 2, 2)

	c = contour(X, Y, approx_posterior_eval, cmap=cm.Spectral, linewidths=2.0)
	clabel(c, inline=0, fontsize=0)
	c.collections[0].set_label('approx posterior variational Bayes')
	plot(mapped_nodes_uncorr.T[0], mapped_nodes_uncorr.T[1], 'sg', markersize=10, label='"uncorrelated" quad nodes')
	legend(loc='best', fontsize=20)

	show()