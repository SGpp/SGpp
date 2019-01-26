from sg_lib.grid.grid import *
from sg_lib.algebraic.multiindex import *
from sg_lib.operation.interpolation_to_spectral import *

def test_function(x):
	x1 = x[0]
	x2 = x[1]

	test = x1*np.sin(x2) + x2*np.cos(x1)

	return test
    
if __name__ == '__main__':
	dim 			= 2
	left_bounds 	= np.zeros(dim)
	right_bounds 	= np.ones(dim)
	grid_level 		= 3
	level_to_nodes 	= 'symmetric'

	weight_dim1 = lambda x: 1.0
	weight_dim2 = lambda x: 1.0

	weights = []
	weights.append(weight_dim1)
	weights.append(weight_dim2)

	Grid_obj 					= Grid(dim, grid_level, level_to_nodes, left_bounds, right_bounds, weights)	
	Multiindex_obj 				= Multiindex(dim)
	InterpolationToSpectral_obj = InterpolationToSpectral(dim, left_bounds, right_bounds, weights)

	#multiindex_set 	= Multiindex_obj.get_std_total_degree_mindex(grid_level)

	multiindex_set = np.array([[1, 1], [1, 2], [2, 1], [1, 3], [2, 2]])

	sg_points 		= Grid_obj.get_std_sg_surplus_points(multiindex_set) 

	for sg_point in sg_points:
		sg_val = test_function(sg_point)
		InterpolationToSpectral_obj.update_sg_evals_lut(sg_val)

	InterpolationToSpectral_obj.update_func_evals(Grid_obj, multiindex_set)

	coeff = InterpolationToSpectral_obj.get_spectral_coeff_sg(multiindex_set)

	#print InterpolationToSpectral_obj.global_basis(multiindex_set)

	print coeff

	x = [0.32, 0.42]

	approx = lambda x_: InterpolationToSpectral_obj.eval_operation_sg(multiindex_set, x_)

	print approx(x)
	print test_function(x)