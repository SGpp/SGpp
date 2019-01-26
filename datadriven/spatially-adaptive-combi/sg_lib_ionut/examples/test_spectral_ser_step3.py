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
	new_multiindices = Adaptivity_obj.unserialize_multiindices_adapt(mindex_file)

	curr_multiindex_set = Adaptivity_obj.multiindex_set

	SpectralProjection_obj.update_func_evals(Grid_obj, curr_multiindex_set)

	curr_no_points 	= Grid_obj.get_no_surplus_grid_points(new_multiindices)
	curr_delta 		= np.zeros(len(new_multiindices))
	curr_coeffs 	= [] 

	for i, multiindex in enumerate(new_multiindices):
		delta 			= SpectralProjection_obj.eval_operation_delta(multiindex, curr_multiindex_set, eval_point)
		curr_delta[i] 	= delta

		delta_ceoff = SpectralProjection_obj.get_spectral_coeff_delta(multiindex, curr_multiindex_set)
		curr_coeffs.append(delta_ceoff)

	curr_coeffs = np.array(curr_coeffs)
	
	# Adaptivity_obj.do_one_adaption_step_postproc(curr_delta, curr_coeffs, curr_no_points)
	Adaptivity_obj.do_one_adaption_step_postproc(curr_delta, curr_coeffs, curr_no_points, SpectralProjection_obj)
	Adaptivity_obj.check_termination_criterion()

	Adaptivity_obj.serialize_data(ref_file)

	if Adaptivity_obj.stop_adaption:
		print 'REFINEMENT FINISHED'
		print 'FINAL RESULTS'
		coeff = SpectralProjection_obj.get_spectral_coeff_sg(Adaptivity_obj.multiindex_set)

		print "Sobol' indices", SpectralProjection_obj.get_total_sobol_indices(coeff, Adaptivity_obj.multiindex_set)
		print 'number of refinement steps', ref_step
		print 'total number of grid points', total_no_gp