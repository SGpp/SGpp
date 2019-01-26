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

	ref_step, prev_len, total_no_gp = Serialize_obj.unserialize_ref_info(info_file)
	SpectralProjection_obj.unserialize_data(data_file)
	Adaptivity_obj.unserialize_data(ref_file)

	new_multiindices = Adaptivity_obj.do_one_adaption_step_preproc()

	curr_multiindex_set = Adaptivity_obj.multiindex_set
	curr_grid_points 	= Grid_obj.get_std_sg_surplus_points(curr_multiindex_set)

	curr_len = len(curr_grid_points)
	for i in xrange(prev_len, curr_len):
		sg_val = test_function(curr_grid_points[i])
		SpectralProjection_obj.update_sg_evals_lut(sg_val)

	SpectralProjection_obj.update_func_evals(Grid_obj, curr_multiindex_set)

	SpectralProjection_obj.serialize_data(data_file)
	Adaptivity_obj.serialize_multiindices_adapt(new_multiindices, mindex_file)
	Adaptivity_obj.serialize_data(ref_file)

	ref_step 	+= 1
	total_no_gp += curr_len

	Serialize_obj.serialize_ref_info(ref_step, curr_len, total_no_gp, info_file)