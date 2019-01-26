import numpy as np
import chaospy as cp
from scipy.optimize import minimize
from scipy.stats import norm
from matplotlib.pyplot import *

def model(a, x):
	return a*x + 1.0 + np.sin(a**2)

def compute_evidence(func, deg=40):
	nodes, weights 	= np.polynomial.legendre.leggauss(deg)
	res 			= np.sum([func(n)*w for n, w in zip(nodes, weights)])

	return res

def compute_mean(func, deg=40):
	nodes, weights 	= np.polynomial.legendre.leggauss(deg)
	res 			= np.sum([n*func(n)*w for n, w in zip(nodes, weights)])

	return res

def compute_std_dev(func, mean, deg=40):
	nodes, weights 	= np.polynomial.legendre.leggauss(deg)
	second_mom 		= np.sum([n**2*func(n)*w for n, w in zip(nodes, weights)])

	res = np.sqrt(second_mom - mean**2)

	return res

def log_parametric_dens(m, sigma, x):
	 ll 			= (x - m)**2/(2.0*sigma**2)
	 norm_constant 	= np.sqrt(2.0*np.pi)*sigma 

	 res = -ll - np.log(norm_constant)
	 
	 return res

def KL_div(x):
	m_ = x[0]
	s_ = x[1]

	res = np.sum(np.array([-likelihood(n)*log_parametric_dens(m_, s_, n)*w for n, w in zip(nodes, weights)]))

	return res

if __name__ == '__main__':
	true_param 		= 0.53
	sigma_noise 	= 0.1
	loc_interest 	= np.arange(0.1, 1.0, 0.1)
	
	prior 	= lambda a: 1.0
	data 	= np.array([model(true_param, loc) + np.random.normal(loc=0., scale=sigma_noise) for loc in loc_interest])

	ll 			= lambda a: np.sum([(model(a, loc) - d)**2 for loc, d in zip(loc_interest, data)])
	likelihood 	= lambda a: np.exp(-ll(a)/(2.0*sigma_noise**2))
	evidence 	= compute_evidence(likelihood)
	posterior 	= lambda a: prior(a)*likelihood(a)/evidence

	mean 	= compute_mean(posterior)
	std_dev = compute_std_dev(posterior, mean) 

	nodes, weights 		= np.polynomial.legendre.leggauss(40)
	x0 					= np.array([mean, std_dev])
	var_Bayes_params 	= minimize(KL_div, x0, method='nelder-mead', tol=1e-14).x
	parametric_dens 	= lambda x: norm.pdf(x, loc=var_Bayes_params[0], scale=var_Bayes_params[1])

	print var_Bayes_params[0], var_Bayes_params[1] 
	print mean, std_dev

	x_axis = np.linspace(0, 1, 100)

	prior_eval 		= [prior(x_) for x_ in x_axis]
	post_eval 		= [posterior(x_) for x_ in x_axis]
	post_var_Bayes 	= [parametric_dens(x_) for x_ in x_axis] 

	plot(x_axis, prior_eval, 'r-', linewidth=3.0, label='prior')
	plot(x_axis, post_eval, 'g-', linewidth=3.0, label='true posterior')
	plot(x_axis, post_var_Bayes, 'b--', linewidth=3.0, label='variatial Bayes posterior')

	legend()
	show()