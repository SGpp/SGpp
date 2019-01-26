from sg_lib.grid.grid import *
from sg_lib.algebraic.multiindex import *
from sg_lib.operation.quadrature import *
from sg_lib.adaptivity.error_work_ratio import *

def test_function(x):
	test = 1.

	for i in xrange(len(x)):
		test *= np.sin(x[i] * np.cos(x[i]))

	return test
    
if __name__ == '__main__':
	dim 			= 2
	left_bounds 	= np.zeros(dim)
	right_bounds 	= np.ones(dim)
	grid_level 		= 1
	growth_factor 	= 2

	tol 		= 1e-5
	max_level 	= 10

	weight 	= lambda x: 1.0
	weights = []

	for d in xrange(dim):
		weights.append(weight)

	init_multiindex = np.ones(dim, dtype=int)

	Grid_obj 		= Grid(dim, grid_level, growth_factor, left_bounds, right_bounds, weights)	
	Multiindex_obj 	= Multiindex(dim)
	Quadrature_obj 	= Quadrature(dim, left_bounds, right_bounds, weights)
	Adaptivity_obj 	= ErrorWorkRatio(dim, tol, init_multiindex, max_level)

	init_multiindex_set = Multiindex_obj.get_std_total_degree_mindex(grid_level)
	init_grid_points 	= Grid_obj.get_std_sg_surplus_points(init_multiindex_set)
	init_no_points 		= Grid_obj.get_no_fg_grid_points(init_multiindex_set)

	for sg_point in init_grid_points:
		sg_val = test_function(sg_point)
		Quadrature_obj.update_sg_evals_lut(sg_val)

	Quadrature_obj.update_func_evals(Grid_obj, init_multiindex_set)
	
	init_delta = Quadrature_obj.eval_operation_delta(init_multiindex, init_multiindex_set)
	Adaptivity_obj.init_adaption(init_delta, init_no_points)

	prev_len 		= len(init_no_points)
	no_adapt_steps 	= 0
	while not Adaptivity_obj.stop_adaption:

		new_multiindices = Adaptivity_obj.do_one_adaption_step_preproc()
	
		curr_multiindex_set = Adaptivity_obj.multiindex_set
		curr_grid_points 	= Grid_obj.get_std_sg_surplus_points(curr_multiindex_set)

		curr_len = len(curr_grid_points)
		for i in xrange(prev_len, curr_len):
			sg_val = test_function(curr_grid_points[i])
			Quadrature_obj.update_sg_evals_lut(sg_val)

		Quadrature_obj.update_func_evals(Grid_obj, curr_multiindex_set)

		curr_no_points 	= Grid_obj.get_no_fg_grid_points(new_multiindices)
		curr_delta 		= np.zeros(len(new_multiindices))

		for i, multiindex in enumerate(new_multiindices):
			delta 			= Quadrature_obj.eval_operation_delta(multiindex, curr_multiindex_set)
			curr_delta[i] 	= delta

		Adaptivity_obj.do_one_adaption_step_postproc(curr_delta, curr_no_points)
		Adaptivity_obj.check_termination_criterion()

		prev_len = curr_len

		no_adapt_steps += 1

		print Adaptivity_obj.eta

	print 'adaptivity terminated after', no_adapt_steps, 'steps'
	print Adaptivity_obj.eta
	print Adaptivity_obj.approx

	print Adaptivity_obj.multiindex_set
