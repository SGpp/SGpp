import numpy as np
import chaospy as cp

def test_function(x):
	a = 7
	b = 0.05

	x1 = x[0]
	x2 = x[1]
	x3 = x[2]

	test = np.sin(x1) + a*np.sin(x2)**2 + b*x3**2*np.sin(x1)

	return test

if __name__ == '__main__':
	poly_deg = 6
	quad_deg = 4
	
	distr 				= cp.J(cp.Uniform(), cp.Uniform(), cp.Uniform())
	nodes, weights 		= cp.generate_quadrature(quad_deg, distr, rule='G')
	poly 				= cp.orth_ttr(poly_deg, distr, normed=True)
	approx 				= np.array([test_function(node) for node in nodes.T])
	gpc_approx, coeff 	= cp.fit_quadrature(poly, nodes, weights, approx, retall=True)

	sobol_indices1 = cp.Sens_m(gpc_approx, distr)
	sobol_indices2 = cp.Sens_t(gpc_approx, distr)

	print len(nodes.T)

	print coeff