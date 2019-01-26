from sg_lib.grid.grid import *
from sg_lib.algebraic.multiindex import *
from sg_lib.operation.spectral_projection import *

def test_function(x):
	x1 = x[0]
	x2 = x[1]

	test = x1*np.sin(x2) + x2*np.cos(x1)

	return test
    
if __name__ == '__main__':
	dim 			= 5
	left_bounds 	= np.zeros(dim)
	right_bounds 	= np.ones(dim)
	grid_level 		= 10
	level_to_nodes 	= 'symmetric'


	weights = [lambda x: 1.0 for i in xrange(dim)]
	
	Grid_obj 				= Grid(dim, grid_level, level_to_nodes, left_bounds, right_bounds, weights)	
	Multiindex_obj 			= Multiindex(dim)
	SpectralProjection_obj 	= SpectralProjection(dim, left_bounds, right_bounds, weights)

	multiindex_set 	= Multiindex_obj.get_std_total_degree_mindex(grid_level)
	sg_points 		= Grid_obj.get_std_sg_surplus_points(multiindex_set) 

	print len(sg_points)

	exit(0)

	for sg_point in sg_points:
		sg_val = test_function(sg_point)
		SpectralProjection_obj.update_sg_evals_lut(sg_val)

	SpectralProjection_obj.update_func_evals(Grid_obj, multiindex_set)

	print SpectralProjection_obj.get_spectral_coeff_sg(multiindex_set)

	exit(0)

	x = [0.32, 0.42]

	approx = lambda x_: SpectralProjection_obj.eval_operation_sg(multiindex_set, x_)

	kk = approx(x)

	print kk	
	print test_function(x)