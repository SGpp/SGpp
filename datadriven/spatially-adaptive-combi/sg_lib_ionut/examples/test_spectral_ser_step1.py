from sg_op.grid.grid import *
from sg_op.algebraic.multiindex import *
from sg_op.operation.spectral_projection import *
# from sg_op.adaptivity.spectral_error_work import *
from sg_op.adaptivity.spectral_error_work_dir_var import *
from sg_op.util.serialize import *
from setup import *

def test_function(x):
	x1 = x[0]
	x2 = x[1]
	x3 = x[2]

	test = (1.0 + 0.33/(2*x1 + x2 + 3.5))*np.exp(-(0.5*(x2 - 0.2)*(x1 + 1))**2) + x1*x2*x3

	return test
    
if __name__ == '__main__':

	Grid_obj 				= Grid(dim, grid_level, growth_factor, left_bounds, right_bounds, weights)	
	Multiindex_obj 			= Multiindex(dim)
	SpectralProjection_obj 	= SpectralProjection(dim, growth_factor, left_bounds, right_bounds, weights)
	Adaptivity_obj 			= SpectralErrorWorkDirVar(dim, tol, tol_dims, init_multiindex, max_level)
	# Adaptivity_obj 			= SpectralErrorWork(dim, tol, init_multiindex, max_level)
	Serialize_obj 			= SerializeRefInfo()

	init_multiindex_set = Multiindex_obj.get_std_total_degree_mindex(grid_level)
	init_grid_points 	= Grid_obj.get_std_sg_surplus_points(init_multiindex_set)
	init_no_points 		= Grid_obj.get_no_fg_grid_points(init_multiindex_set)

	for sg_point in init_grid_points:
		sg_val = test_function(sg_point)
		SpectralProjection_obj.update_sg_evals_lut(sg_val)

	SpectralProjection_obj.update_func_evals(Grid_obj, init_multiindex_set)
	
	init_delta = SpectralProjection_obj.eval_operation_delta(init_multiindex, init_multiindex_set, eval_point)
	init_coeff = SpectralProjection_obj.get_spectral_coeff_delta(init_multiindex, init_multiindex_set)

	Adaptivity_obj.init_adaption(init_delta, init_coeff, init_no_points)
	
	total_no_gp += init_no_points
	SpectralProjection_obj.serialize_data(data_file)
	Adaptivity_obj.serialize_data(ref_file)
	Serialize_obj.serialize_ref_info(ref_step, init_no_points, total_no_gp, info_file)