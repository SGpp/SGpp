from sg_op.grid.grid import *
from sg_op.algebraic.multiindex import *
from sg_op.operation.quadrature import *

def test_function(x):
	x1 = x[0]
	x2 = x[1]

	test = x1*np.sin(x2) + x2*np.cos(x1)

	return test
    
if __name__ == '__main__':
	dim 			= 2
	left_bounds 	= np.zeros(dim)
	right_bounds 	= np.ones(dim)
	grid_level 		= 5
	growth_factor 	= 2

	weight_dim1 = lambda x: 1.0
	weight_dim2 = lambda x: 1.0

	weights 	= []
	weights.append(weight_dim1)
	weights.append(weight_dim2)

	a = [0.2, 0.4]
	b = [1.1, 2.2]

	Grid_obj 		= Grid(dim, grid_level, growth_factor, left_bounds, right_bounds, weights)	
	Multiindex_obj 	= Multiindex(dim)
	Quadrature_obj 	= Quadrature(dim, left_bounds, right_bounds, weights)

	multiindex_set 	= Multiindex_obj.get_std_total_degree_mindex(grid_level)
	sg_points 		= Grid_obj.get_std_sg_surplus_points(multiindex_set)

	for sg_point in sg_points:

		eval_point = np.zeros(dim)
		for d in xrange(dim):
			eval_point[d] = a[d] + (b[d] - a[d])*sg_point[d]

		sg_val = test_function(eval_point)
		Quadrature_obj.update_sg_evals_lut(sg_val)

	Quadrature_obj.update_func_evals(Grid_obj, multiindex_set)
	res = Quadrature_obj.eval_operation_sg(multiindex_set)

	print res