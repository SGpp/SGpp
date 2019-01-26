from sg_op.grid.grid import *
from sg_op.algebraic.multiindex import *
from sg_op.operation.spectral_projection import *
from sg_op.adaptivity.spectral_error_work import *
from matplotlib.pyplot import *
from matplotlib.patches import *

def test_function(x):
	x1 = x[0]
	x2 = x[1]
	x3 = x[2]

	test = (1.0 + 0.33/(2*x1 + x2 + 3.5))*np.exp(-(0.5*(x2 - 0.2)*(x1 + 1))**2) + x1*x2*x3

	return test
    
if __name__ == '__main__':
	dim 			= 3
	left_bounds 	= np.zeros(dim)
	right_bounds 	= np.ones(dim)
	grid_level 		= 1
	growth_factor 	= 2

	tol 		= 1e-6
	max_level 	= 10

	weight = lambda x: 1.0

	weights 	= []
	for d in xrange(dim):
		weights.append(weight)

	init_multiindex = np.ones(dim, dtype=int)

	eval_point = 0.33*np.ones(dim)

	Grid_obj 				= Grid(dim, grid_level, growth_factor, left_bounds, right_bounds, weights)	
	Multiindex_obj 			= Multiindex(dim)
	SpectralProjection_obj 	= SpectralProjection(dim, growth_factor, left_bounds, right_bounds, weights)
	Adaptivity_obj 			= SpectralErrorWork(dim, tol, init_multiindex, max_level)

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

	prev_len 		= len(init_no_points)
	no_adapt_steps 	= 0
	total_len 		= prev_len 
	while not Adaptivity_obj.stop_adaption:
		no_adapt_steps += 1

		new_multiindices = Adaptivity_obj.do_one_adaption_step_preproc()
	
		curr_multiindex_set = Adaptivity_obj.multiindex_set
		curr_grid_points 	= Grid_obj.get_std_sg_surplus_points(curr_multiindex_set)

		curr_len = len(curr_grid_points)
		for i in xrange(prev_len, curr_len):
			sg_val = test_function(curr_grid_points[i])
			SpectralProjection_obj.update_sg_evals_lut(sg_val)

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
		
		Adaptivity_obj.do_one_adaption_step_postproc(curr_delta, curr_coeffs, curr_no_points)
		Adaptivity_obj.check_termination_criterion()

		prev_len 	= curr_len
		total_len 	+= curr_len

	# print 'adaptivity terminated after', no_adapt_steps, 'steps'
	# print Adaptivity_obj.eta
	# print Adaptivity_obj.approx
	# print test_function(eval_point)

	# print Adaptivity_obj.multiindex_set

	coeff = SpectralProjection_obj.get_spectral_coeff_sg(Adaptivity_obj.multiindex_set)

	print SpectralProjection_obj.get_total_sobol_indices(coeff, Adaptivity_obj.multiindex_set)

	print total_len

	# plotting
	fig = figure()
	ax 	= fig.add_subplot(111, aspect='equal')

	old_index_set 	= Adaptivity_obj.O
	active_set 		= Adaptivity_obj.A

	temp 				= np.ones(dim, dtype=int)
	max_reached_level 	= np.amax(Adaptivity_obj.multiindex_set)

	for p in [
		Rectangle(
	    tuple(multiindex - temp),
	    1.0,
	    1.0,
	    hatch='\\',
	    facecolor='lightblue',
	    edgecolor='black', 
	    linewidth=4.0,
	    label='old index set' ) for multiindex in old_index_set.values()]:
		ax.add_patch(p)
		old_index_set_handle = p

	for multiindex in old_index_set.values():
		x_coord = (2*tuple(multiindex - temp)[0] + 0.5)/2.
		y_coord = (2*tuple(multiindex - temp)[1] + 0.9)/2.

		ax.text(x_coord, y_coord, str(tuple(multiindex)), fontsize=20)

	if len(active_set) >= 1:
		for p in [
			Rectangle(
		    tuple(multiindex - temp),
		    1.0,
		    1.0,
	        facecolor='orangered',
		    edgecolor='black', 
		    linewidth=4.0,
		    label='active set' ) for multiindex in active_set.values()]:
			ax.add_patch(p)
			active_set_handle = p

		for multiindex in active_set.values():
			x_coord = (2*tuple(multiindex - temp)[0] + 0.5)/2.
			y_coord = (2*tuple(multiindex - temp)[1] + 0.9)/2.

			ax.text(x_coord, y_coord, str(tuple(multiindex)), fontsize=20)

		legend(handles=[old_index_set_handle, active_set_handle], loc='best', fontsize=25)
	else:
		legend(handles=[old_index_set_handle], loc='best', fontsize=25)

	x_coord = (2*max_reached_level - 3)/2.
	y_coord = (2*max_reached_level - 2)/2.	
	fig.suptitle(str(total_len) + ' grid points', fontsize=20)

	ax.set_xticks([])
	ax.set_yticks([])
	ax.set_xlim([0, max_reached_level])
	ax.set_ylim([0, max_reached_level])

	show()