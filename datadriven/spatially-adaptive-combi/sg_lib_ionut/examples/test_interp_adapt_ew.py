from sg_op.grid.grid import *
from sg_op.algebraic.multiindex import *
from sg_op.operation.interpolation import *
from sg_op.adaptivity.error_work_ratio import *

def test_function(x):
	x1 = x[0]
	x2 = x[1]

	test = x1*np.sin(x2) + x2*np.cos(x1)

	test = 0.75 * np.exp(-(9*x1 - 2)**2/4 - (9*x2 - 2)**2/4) + \
	      0.75 * np.exp(-(9*x1 + 1)**2/49 - (9*x2 + 1)/10) + \
	      0.5 * np.exp(-(9*x1 - 7)**2/4 - (9*x2 - 3)**2/4) - \
	      -0.2 * np.exp(-(9*x1 - 4)**2 - (9*x2 - 7)**2)

	return test
    
if __name__ == '__main__':
	dim 			= 2
	left_bounds 	= np.zeros(dim)
	right_bounds 	= np.ones(dim)
	grid_level 		= 1
	growth_factor 	= 2

	tol 		= 1e-5
	max_level 	= 10

	weight_dim1 = lambda x: 1.0
	weight_dim2 = lambda x: 1.0

	weights 	= []
	weights.append(weight_dim1)
	weights.append(weight_dim2)
	weights.append(weight_dim2)

	init_multiindex = np.ones(dim, dtype=int)

	eval_point = [0.32, 0.42]

	Grid_obj 			= Grid(dim, grid_level, growth_factor, left_bounds, right_bounds, weights)	
	Multiindex_obj 		= Multiindex(dim)
	Interpolation_obj 	= Interpolation(dim)
	Adaptivity_obj 		= ErrorWorkRatio(dim, tol, init_multiindex, max_level)

	init_multiindex_set = Multiindex_obj.get_std_total_degree_mindex(grid_level)
	init_grid_points 	= Grid_obj.get_std_sg_surplus_points(init_multiindex_set)
	init_no_points 		= Grid_obj.get_no_fg_grid_points(init_multiindex_set)

	for sg_point in init_grid_points:
		sg_val = test_function(sg_point)
		Interpolation_obj.update_sg_evals_lut(sg_val)

	Interpolation_obj.update_func_evals(Grid_obj, init_multiindex_set)
	
	init_delta = Interpolation_obj.eval_operation_delta(init_multiindex, init_multiindex_set, eval_point)
	Adaptivity_obj.init_adaption(init_delta, init_no_points)

	prev_len 		= len(init_no_points)
	no_adapt_steps 	= 0
	while not Adaptivity_obj.stop_adaption:
		no_adapt_steps += 1

		new_multiindices = Adaptivity_obj.do_one_adaption_step_preproc()
	
		curr_multiindex_set = Adaptivity_obj.multiindex_set
		curr_grid_points 	= Grid_obj.get_std_sg_surplus_points(curr_multiindex_set)

		curr_len = len(curr_grid_points)
		for i in xrange(prev_len, curr_len):
			sg_val = test_function(curr_grid_points[i])
			Interpolation_obj.update_sg_evals_lut(sg_val)

		Interpolation_obj.update_func_evals(Grid_obj, curr_multiindex_set)

		curr_no_points 	= Grid_obj.get_no_fg_grid_points(new_multiindices)
		curr_delta 		= np.zeros(len(new_multiindices))

		for i, multiindex in enumerate(new_multiindices):
			delta 			= Interpolation_obj.eval_operation_delta(multiindex, curr_multiindex_set, eval_point)
			curr_delta[i] 	= delta

		Adaptivity_obj.do_one_adaption_step_postproc(curr_delta, curr_no_points)
		Adaptivity_obj.check_termination_criterion()

		prev_len = curr_len

		print Adaptivity_obj.eta

	print 'adaptivity terminated after', no_adapt_steps, 'steps'
	print Adaptivity_obj.eta
	print Adaptivity_obj.approx
	print test_function(eval_point)

	print Adaptivity_obj.multiindex_set