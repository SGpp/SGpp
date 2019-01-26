from sg_op.grid.grid import *
from sg_op.algebraic.multiindex import *
import chaospy as cp

if __name__ == '__main__':
	dim 			= 3
	left_bounds 	= np.zeros(dim)
	right_bounds 	= np.ones(dim)
	grid_level 		= 6
	growth_factor 	= 2

	weight 	= lambda x: 1.0
	weights = []

	for d in xrange(dim):
		weights.append(weight)

	Grid_obj 		= Grid(dim, grid_level, growth_factor, left_bounds, right_bounds, weights)	
	Multiindex_obj 	= Multiindex(dim)

	multiindex_set 	= Multiindex_obj.get_std_total_degree_mindex(grid_level)
	sg_points 		= Grid_obj.get_std_sg_surplus_points(multiindex_set)

	print len(sg_points)


	distr = cp.J(cp.Uniform(1.22, 2.43), cp.Uniform(1.22, 2.43), cp.Uniform(0.78, 1.25))
	nodes, weights = cp.generate_quadrature(7, distr, rule='G')

	print len(nodes.T)
	for node in nodes.T:
		print node