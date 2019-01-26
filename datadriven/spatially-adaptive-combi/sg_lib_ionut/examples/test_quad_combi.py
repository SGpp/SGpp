from sg_op.grid.grid_combi import *
from sg_op.algebraic.multiindex_combi import *
from sg_op.operation.quadrature_combi import *
from matplotlib.pyplot import *

def test_function(x):
	x1 = x[0]
	x2 = x[1]

	test = np.sin(x1)*x2 + np.cos(x2)*x1

	return test

if __name__ == '__main__':
	dim 			= 2
	left_bounds 	= np.zeros(dim)
	right_bounds 	= np.ones(dim)
	grid_level 		= 4
	growth_factor 	= 2

	weight_dim1 = lambda x: 1.0
	weight_dim2 = lambda x: 1.0

	weights 	= []
	weights.append(weight_dim1)
	weights.append(weight_dim2)

	Grid_obj 		= Grid(dim, grid_level, growth_factor, left_bounds, right_bounds, weights)	
	Multiindex_obj 	= Multiindex(dim)
	Quadrature_obj 	= Quadrature(dim, left_bounds, right_bounds, weights)

	multiindex_set 	= Multiindex_obj.get_std_mindex(grid_level)
	combi_coeffs 	= Multiindex_obj.get_std_combi_coeff(grid_level, multiindex_set)

	fg_res = np.zeros(len(multiindex_set))
	for index, multiindex in enumerate(multiindex_set):
		grid_points_fg, grid_points_1D = Grid_obj.get_multiindex_grid_points(multiindex)

		func_eval = np.zeros(len(grid_points_fg))
		for i, grid_point in enumerate(grid_points_fg):
			func_eval[i] = test_function(grid_point)

		fg_res[index] = Quadrature_obj.eval_full_grid_quadrature(func_eval, grid_points_1D)

	sg_res = Quadrature_obj.eval_combi_quadrature(fg_res, combi_coeffs)

	print sg_res