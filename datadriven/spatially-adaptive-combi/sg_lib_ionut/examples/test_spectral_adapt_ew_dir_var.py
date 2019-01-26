from sg_lib.grid.grid import *
from sg_lib.algebraic.multiindex import *
from sg_lib.operation.interpolation_to_spectral import *
from sg_lib.adaptivity.spectral_error_work_dir_var import *
from matplotlib.pyplot import *
from matplotlib.patches import *

def test_function(x):

	test = 1.0
	for i in xrange(dim):
		test *= np.exp(-(i + 1)*x[i])*np.sin(x[i])

	return test
    
if __name__ == '__main__':
	dim 			= 3
	left_bounds 	= np.zeros(dim)
	right_bounds 	= np.ones(dim)
	grid_level 		= 1
	level_to_nodes 	= 2

	tol 		= 1e-4
	tol_dims 	= 1e-20*np.ones(dim)
	max_level 	= 20

	weights = [lambda x: 1. for i in xrange(dim)]

	init_multiindex = np.ones(dim, dtype=int)

	eval_point = 0.33*np.ones(dim)

	Grid_obj 				= Grid(dim, grid_level, level_to_nodes, left_bounds, right_bounds, weights)	
	Multiindex_obj 			= Multiindex(dim)
	InterpToSpectral_obj 	= InterpolationToSpectral(dim, level_to_nodes, left_bounds, right_bounds, weights)
	Adaptivity_obj 			= SpectralErrorWorkDirVar(dim, tol, tol_dims, init_multiindex, max_level)

	init_multiindex_set = Multiindex_obj.get_std_total_degree_mindex(grid_level)
	init_grid_points 	= Grid_obj.get_std_sg_surplus_points(init_multiindex_set)
	init_no_points 		= Grid_obj.get_no_fg_grid_points(init_multiindex_set)

	for sg_point in init_grid_points:
		sg_val = test_function(sg_point)
		InterpToSpectral_obj.update_sg_evals_lut(sg_val)

	InterpToSpectral_obj.update_func_evals(Grid_obj, init_multiindex_set)
	
	init_coeff = InterpToSpectral_obj.get_spectral_coeff_delta(init_multiindex, init_multiindex_set)
	Adaptivity_obj.init_adaption(init_coeff, init_no_points)

	prev_len 		= len(init_no_points)
	no_adapt_steps 	= 0
	total_len 		= init_no_points
	while not Adaptivity_obj.stop_adaption:
		no_adapt_steps += 1

		new_multiindices = Adaptivity_obj.do_one_adaption_step_preproc()
	
		curr_multiindex_set = Adaptivity_obj.multiindex_set
		curr_grid_points 	= Grid_obj.get_std_sg_surplus_points(curr_multiindex_set)

		for multiindex in new_multiindices:
			new_grid_points = Grid_obj.get_sg_surplus_points_multiindex(multiindex)
			total_len 		+= len(new_grid_points)

			for sg_point in new_grid_points:
				sg_val = test_function(sg_point)
				InterpToSpectral_obj.update_sg_evals_lut(sg_val)

		InterpToSpectral_obj.update_func_evals(Grid_obj, curr_multiindex_set)

		print 'multiindices active set'
		for multiindex in new_multiindices:
			print InterpToSpectral_obj.get_var_level_active_set(curr_multiindex_set, multiindex)

		curr_no_points 	= Grid_obj.get_no_surplus_grid_points(new_multiindices)
		curr_coeffs 	= [] 

		for multiindex in new_multiindices:
			delta_coeff = InterpToSpectral_obj.get_spectral_coeff_delta(multiindex, curr_multiindex_set)
			curr_coeffs.append(delta_coeff)

		curr_coeffs = np.array(curr_coeffs)
		
		Adaptivity_obj.do_one_adaption_step_postproc(curr_coeffs, curr_no_points, InterpToSpectral_obj)
		Adaptivity_obj.check_termination_criterion()

		coeff = InterpToSpectral_obj.get_spectral_coeff_sg(Adaptivity_obj.multiindex_set)

		print 'active set dir var'
		print InterpToSpectral_obj.get_directional_var_active_set(Adaptivity_obj.multiindex_set, Adaptivity_obj.A.values())

		print 'total Sobol indices'
		print InterpToSpectral_obj.get_total_sobol_indices(coeff, Adaptivity_obj.multiindex_set)

	print 'adaptivity done'

	# print 'adaptivity terminated after', no_adapt_steps, 'steps'
	# print Adaptivity_obj.eta
	# print Adaptivity_obj.approx

	print InterpToSpectral_obj.eval_operation_sg(Adaptivity_obj.multiindex_set, eval_point)
	print test_function(eval_point)

	print Adaptivity_obj.multiindex_set

	coeff = InterpToSpectral_obj.get_spectral_coeff_sg(Adaptivity_obj.multiindex_set)

	print InterpToSpectral_obj.get_total_sobol_indices(coeff, Adaptivity_obj.multiindex_set)

	print total_len[0]