import numpy as np
from scipy.special import eval_hermitenorm, eval_sh_legendre
from math import factorial

norm_hermite = {0: 1.0, 1: 1.0, 2: 1.4142135623730951, 3: 2.4494897427831783, 4: 4.898979485566357, 5: 10.954451150103324, 6: 26.832815729997478, 7: 70.9929573971954, 8: 200.79840636817815, 9: 602.3952191045345, 10: 1904.9409439665055, 11: 6317.974358922329, 12: 21886.105181141756, 13: 78911.47445080469, 14: 295259.7012800765, 15: 1143535.9058639132, 16: 4574143.623455653, 17: 18859677.306253154, 18: 80014834.28544986, 19: 348776576.6344295, 20: 1559776268.6284986, 21: 7147792818.185868, 22: 33526120082.371723, 23: 160785623545.4059, 24: 787685471322.9384, 25: 3938427356614.692, 26: 20082117944245.96, 27: 104349745809073.98, 28: 552166953567228.56, 29: 2973510046012911.0, 30: 1.6286585271694958e+16, 31: 9.067986906793549e+16, 32: 5.129628026803635e+17, 33: 2.9467469553410734e+18, 34: 1.7182339742875652e+19, 35: 1.016520927791757e+20, 36: 6.099125566750542e+20, 37: 3.709953246501409e+21, 38: 2.28696877430935e+22, 39: 1.4282115417961528e+23, 40: 9.032802905233224e+23}


def get_1D_orth_poly(deg, left, right, x_1D):

	poly_eval = 0.

	if left == 0. and right == 1.:
		poly_eval = eval_sh_legendre(deg, x_1D)*np.sqrt(2*deg + 1)
	elif left == -5. and right == 5.:
		poly_eval = eval_hermitenorm(deg, x_1D)/norm_hermite[deg]
	else:
		print 'not implemented!'
		exit(0)

	return poly_eval

def get_1D_barycentric_weights(grid_1D):
		
	size    = len(grid_1D)
	w       = np.zeros(size)

	w[0] = 1.0
	for j in xrange(1, size):
		for k in xrange(j):
			w[k] *= (grid_1D[k] - grid_1D[j])

		w[j] = np.prod([(grid_1D[j] - grid_1D[k]) for k in xrange(j)])			

	for j in xrange(size):
	    w[j] = 1./w[j]

	return w

def eval_1D_barycentric_interpolant(grid_1D, barycentric_weights, j, x_1D):
	
	l = np.prod([(x_1D - x_i) for x_i in grid_1D])

	barycentric_interpolant = 1.0

	if np.abs(x_1D - grid_1D[j]) > 1e-14:
		barycentric_interpolant = l * barycentric_weights[j]/(x_1D - grid_1D[j])
		
	return barycentric_interpolant

def eval_1D_orth_poly(deg, left, right, weight, x_1D):
	
	poly_eval= get_1D_orth_poly(deg, left, right, x_1D)

	return poly_eval

def compute_1D_quad_weights(grid_1D, left, right, weight):
		
	N = len(grid_1D)
	V = np.zeros((N, N))

	for i in xrange(N):
		for j in xrange(N):
			V[i, j] = get_1D_orth_poly(j, left, right, grid_1D[i])

	weights = np.linalg.inv(V)[0, :]

	return weights
