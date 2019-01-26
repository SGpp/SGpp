from sg_op_combi.grid.grid_combi import *
from sg_op_combi.algebraic.multiindex_combi import *
from sg_op_combi.operation.spectral_projection_combi import *

def test_function(x):
	x1 = x[0]
	x2 = x[1]

	test = x1*np.sin(x2) + x2*np.cos(x1)

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

	Grid_obj 				= Grid(dim, grid_level, growth_factor, left_bounds, right_bounds, weights)	
	Multiindex_obj 			= Multiindex(dim)
	SpectralProjection_obj 	= SpectralProjection(dim, left_bounds, right_bounds, weights)

	multiindex_set 	= Multiindex_obj.get_std_mindex(grid_level)
	combi_coeffs 	= Multiindex_obj.get_std_combi_coeff(grid_level, multiindex_set)

	fg_res = []
	for index, multiindex in enumerate(multiindex_set):
		grid_points_fg, grid_points_1D = Grid_obj.get_multiindex_grid_points(multiindex)

		func_evals = np.zeros(len(grid_points_fg))
		for i, grid_point in enumerate(grid_points_fg):
			func_evals[i] = test_function(grid_point)

		spectral_coeff = SpectralProjection_obj.get_spectral_coeff_fg_full(func_evals, grid_points_1D, multiindex_set, growth_factor)
		fg_res.append(spectral_coeff)

	spectral_coeff = SpectralProjection_obj.get_combi_spectral_coeff(fg_res, combi_coeffs)

	my_E 			= SpectralProjection_obj.get_mean(spectral_coeff)
	my_Var 			= SpectralProjection_obj.get_variance(spectral_coeff)
	my_Sobol_first 	= SpectralProjection_obj.get_first_order_sobol_indices(spectral_coeff, growth_factor, multiindex_set)
	my_Sobol_total 	= SpectralProjection_obj.get_total_sobol_indices(spectral_coeff, growth_factor, multiindex_set)

	print SpectralProjection_obj.get_spectral_approx_combi(growth_factor, multiindex_set, spectral_coeff, [0.32, 0.42])
	print test_function([0.32, 0.42])