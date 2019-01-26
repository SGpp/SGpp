import numpy as np
from scipy.optimize import fmin
from itertools import product
from collections import OrderedDict

class Grid(object):
	def __init__(self, dim, grid_level, linear_growth_factor, left_bounds, right_bounds, weights):
		self.__dim 					= dim
		self.__grid_level 			= grid_level
		self.__linear_growth_factor = linear_growth_factor
		self.__left_bounds 			= left_bounds
		self.__right_bounds 		= right_bounds
		self.__weights 				= weights

		self.__init_level 	= 1
		self.__machine_eps 	= 7./3 - 4./3 - 1.
		self.__test_eval_x 	= (left_bounds + right_bounds)/2.0

	def __get_no_1D_grid_points(self, level):

		no_points = 0

		if self.__linear_growth_factor == 'sym':
			no_points = 2*level - 1
		else:
			if self.__linear_growth_factor == 1:
				no_points = level
			elif self.__linear_growth_factor >= 2:
				if level == 1:
					no_points = 1
				else:
					no_points = self.__linear_growth_factor*level - 1
			else:
				raise NotImplementedError

		return no_points

	def __minimize_function(self, func, a, b):

		guess = a + (b - a)/2.0
		f_min = fmin(func, guess, xtol=self.__machine_eps, maxiter=10000, disp=False)[-1]

		return f_min 

	def __get_lleja_poly(self, x, sorted_points, a, b, weight = lambda x: 1.0):

		if (x < a or x > b):
		    return -1

		poly = 1.0
		for i in xrange(len(sorted_points)):
		    poly *= np.abs(x - sorted_points[i])

		poly *= weight(x)

		return poly

	def __get_neg_lleja_poly(self, x, sorted_points, a, b, weight = lambda x: 1.0):

		return (-1.0) * self.__get_lleja_poly(x, sorted_points, a, b, weight)

	def __get_starting_point(self, a, b, weight = lambda x: 1.0):

		neg_weight 		= lambda x: -weight(x)
		starting_point 	= self.__minimize_function(neg_weight, a, b)

		return starting_point

	def __get_1D_level_points(self, curr_level, left_bound, right_bound, weight = lambda x: 1.0, eps=1e-14):

		sorted_points 	= []
		unsorted_points = []
		no_points 		= self.__get_no_1D_grid_points(curr_level)

		starting_point = self.__get_starting_point(left_bound, right_bound, weight)
		sorted_points.append(starting_point)
		unsorted_points.append(starting_point)

		for point in xrange(1, no_points):
		    x_val = []
		    y_val = []

		    a = 0.
		    b = left_bound
		    for i in xrange(len(sorted_points) + 1):
		        sorted_points = sorted(sorted_points)

		        a = b
		        if i < len(sorted_points):
		        	b = sorted_points[i]
		        else:
		        	b = right_bound

		        x_min = (a + b)/2.0
		        y_min = 0.0

		        if np.abs(b - a) > eps:
			        wlleja_func = lambda x : self.__get_neg_lleja_poly(x, sorted_points, a, b, weight)
			        x_min 		= self.__minimize_function(wlleja_func, a, b)

		        	x_val.append(x_min)
		        	y_min = wlleja_func(x_val[-1])
		        	y_val.append(y_min)
		        else:
		        	x_val.append(x_min)
		        	y_val.append(y_min)

			wlleja_point = x_val[y_val.index(np.min(y_val))]
			sorted_points.append(wlleja_point)
		    unsorted_points.append(wlleja_point)

		unsorted_points = np.array(unsorted_points, dtype=np.float64)

		return unsorted_points

	def __get_1D_surplus_and_level_points(self, unsorted_points_prev_level, next_level, left_bound, right_bound, \
													weight = lambda x: 1.0, eps=1e-14):

		sorted_points 	= unsorted_points_prev_level.tolist()[:]
		unsorted_points = unsorted_points_prev_level.tolist()[:]
		surplus_points 	= []
		no_points 		= self.__get_no_1D_grid_points(next_level)

		x0 				= unsorted_points[0]
		starting_point 	= self.__get_starting_point(left_bound, right_bound, weight)

		for point in xrange(len(sorted_points), no_points):

			if self.__linear_growth_factor == 'sym' and point >= 3 and np.mod(point, 2) == 1:
				x_val = []
				y_val = []

				a = 0.
				b = left_bound
				for i in xrange(len(sorted_points) + 1):
				    sorted_points = sorted(sorted_points)

				    a = b
				    if i < len(sorted_points):
				    	b = sorted_points[i]
				    else:
				    	b = right_bound

				    x_min = (a + b)/2.0
				    y_min = eps

				    if np.abs(b - a) > eps:
				        wlleja_func = lambda x : self.__get_neg_lleja_poly(x, sorted_points, a, b, weight)
				        x_min 		= self.__minimize_function(wlleja_func, a, b)

				    	x_val.append(x_min)
				    	y_min = wlleja_func(x_val[-1])
				    	y_val.append(y_min)
				    else:
				    	x_val.append(x_min)
				    	y_val.append(y_min)

				wlleja_point = x_val[y_val.index(np.min(y_val))]
				sorted_points.append(wlleja_point)
				unsorted_points.append(wlleja_point)
				surplus_points.append(wlleja_point)

			elif self.__linear_growth_factor == 'sym' and point >= 3 and np.mod(point, 2) == 0:
				wlleja_point = 2*starting_point - unsorted_points[point - 1]
				sorted_points.append(wlleja_point)
				unsorted_points.append(wlleja_point)
				surplus_points.append(wlleja_point)
			else:
				x_val = []
				y_val = []

				a = 0.
				b = left_bound
				for i in xrange(len(sorted_points) + 1):
				    sorted_points = sorted(sorted_points)

				    a = b
				    if i < len(sorted_points):
				    	b = sorted_points[i]
				    else:
				    	b = right_bound

				    x_min = (a + b)/2.0
				    y_min = eps

				    if np.abs(b - a) > eps:
				        wlleja_func = lambda x : self.__get_neg_lleja_poly(x, sorted_points, a, b, weight)
				        x_min 		= self.__minimize_function(wlleja_func, a, b)

				    	x_val.append(x_min)
				    	y_min = wlleja_func(x_val[-1])
				    	y_val.append(y_min)
				    else:
				    	x_val.append(x_min)
				    	y_val.append(y_min)

				wlleja_point = x_val[y_val.index(np.min(y_val))]
				sorted_points.append(wlleja_point)
				unsorted_points.append(wlleja_point)
				surplus_points.append(wlleja_point)

		unsorted_points = np.array(unsorted_points, dtype=np.float64)
		surplus_points 	= np.array(surplus_points, dtype=np.float64)

		return unsorted_points, surplus_points

	def __get_indices_all(self, level):

		no_points 	= self.__get_no_1D_grid_points(level)
		indices_all = np.array([i for i in xrange(no_points)])

		return indices_all

	def __get_indices_surpluses_all(self, no_surplus_points, level):

		no_points 			= self.__get_no_1D_grid_points(level)
		indices_surpluses 	= np.array([i for i in xrange(no_points - no_surplus_points, no_points)])

		return indices_surpluses

	def __tensorize(self, univariate_list):
		
		tensorization = np.array(list(product(*univariate_list)))

		return tensorization

	def __get_1D_full_points_per_dim(self, max_level, d, eps=1e-14):
		
		all_1D_grid_points_dim_d 	= []
		all_1D_indices_dim_d 		= []
		unique_levels 				= np.array(range(self.__init_level, max_level + 1), dtype=int)

		unsorted_prev_points = self.__get_1D_level_points(unique_levels[0], self.__left_bounds[d], \
															self.__right_bounds[d], self.__weights[d], eps)

		all_1D_grid_points_dim_d.append(unsorted_prev_points)
		for i in xrange(1, len(unique_levels)):
			unsorted_prev_points, surplus_level_points = self.__get_1D_surplus_and_level_points(unsorted_prev_points, \
								unique_levels[i], self.__left_bounds[d], self.__right_bounds[d], self.__weights[d], eps)

			all_1D_grid_points_dim_d.append(unsorted_prev_points)

		return all_1D_grid_points_dim_d

	def __get_1D_full_indices_per_dim(self, max_level, d):

		all_1D_indices_dim_d 		= []
		unique_levels 				= np.array(range(self.__init_level, max_level + 1), dtype=int)

		indices_all = self.__get_indices_all(unique_levels[0])

		all_1D_indices_dim_d.append(indices_all)
		for i in xrange(1, len(unique_levels)):
			indices_all = self.__get_indices_all(unique_levels[i])
			
			all_1D_indices_dim_d.append(indices_all)

		return all_1D_indices_dim_d

	def __get_1D_surplus_points_per_dim(self, max_level, d, eps=1e-14):
		
		all_1D_grid_points_dim_d 	= []
		unique_levels 				= np.array(range(self.__init_level, max_level + 1), dtype=int)

		unsorted_prev_points = self.__get_1D_level_points(unique_levels[0], self.__left_bounds[d], \
															self.__right_bounds[d], self.__weights[d], eps)

		all_1D_grid_points_dim_d.append(unsorted_prev_points)
		for i in xrange(1, len(unique_levels)):
			unsorted_prev_points, surplus_level_points = self.__get_1D_surplus_and_level_points(unsorted_prev_points, \
								unique_levels[i], self.__left_bounds[d], self.__right_bounds[d], self.__weights[d], eps)

			all_1D_grid_points_dim_d.append(surplus_level_points)

		return all_1D_grid_points_dim_d

	def __get_1D_surplus_indices_per_dim(self, max_level, d):
		
		all_1D_indices_dim_d 		= []
		unique_levels 				= np.array(range(self.__init_level, max_level + 1), dtype=int)

		prev_indices_all = self.__get_indices_all(unique_levels[0])
		all_1D_indices_dim_d.append(prev_indices_all)
		for i in xrange(1, len(unique_levels)):
			curr_indices_all 	= self.__get_indices_all(unique_levels[i])
			no_surplus_points 	= len(curr_indices_all) - len(prev_indices_all)
			indices_surplus 	= self.__get_indices_surpluses_all(no_surplus_points, unique_levels[i])

			all_1D_indices_dim_d.append(indices_surplus)
			prev_indices_all = curr_indices_all[:]

		return all_1D_indices_dim_d

	def __get_all_1D_indices(self, max_level):
		
		all_1D_indices 		= []

		all_1D_indices_dim = self.__get_1D_full_indices_per_dim(max_level, 0)
		all_1D_indices.append(all_1D_indices_dim)

		for d1 in xrange(1, self.__dim):
			for d2 in xrange(self.__dim):
				if d1 != d2:
					if self.__weights[d1](self.__test_eval_x[d1]) - self.__weights[d2](self.__test_eval_x[d2]) == 0.:
						if d1 < d2:
							all_1D_indices_dim = all_1D_indices[d1]
						else:
							all_1D_indices_dim = all_1D_indices[d2]

						all_1D_indices.append(all_1D_indices_dim)
						break
					else:
						all_1D_indices_dim = self.__get_1D_full_indices_per_dim(max_level, d1)
						all_1D_indices.append(all_1D_indices_dim)

		return all_1D_indices

	def __get_all_1D_surplus_indices(self, max_level):
		
		all_1D_indices 		= []

		all_1D_indices_dim = self.__get_1D_surplus_indices_per_dim(max_level, 0)
		all_1D_indices.append(all_1D_indices_dim)

		for d1 in xrange(1, self.__dim):
			for d2 in xrange(self.__dim):
				if d1 != d2:
					if self.__weights[d1](self.__test_eval_x[d1]) - self.__weights[d2](self.__test_eval_x[d2]) == 0.:
						if d1 < d2:
							all_1D_indices_dim = all_1D_indices[d1]
						else:
							all_1D_indices_dim = all_1D_indices[d2]

						all_1D_indices.append(all_1D_indices_dim)
						break
					else:
						all_1D_indices_dim = self.__get_1D_surplus_indices_per_dim(max_level, d1)
						all_1D_indices.append(all_1D_indices_dim)

		return all_1D_indices

	def __get_all_1D_surplus_points(self, max_level, eps=1e-14):
		
		all_1D_grid_points 	= []

		all_1D_grid_points_dim = self.__get_1D_surplus_points_per_dim(max_level, 0, eps)
		all_1D_grid_points.append(all_1D_grid_points_dim)

		for d1 in xrange(1, self.__dim):
			for d2 in xrange(self.__dim):
				if d1 != d2:
					if self.__weights[d1](self.__test_eval_x[d1]) - self.__weights[d2](self.__test_eval_x[d2]) == 0.:
						if d1 < d2:
							all_1D_grid_points_dim 	= all_1D_grid_points[d1]
						else:
							all_1D_grid_points_dim 	= all_1D_grid_points[d2]

						all_1D_grid_points.append(all_1D_grid_points_dim)
						break
					else:
						all_1D_grid_points_dim = self.__get_1D_surplus_points_per_dim(max_level, d1, eps)
						all_1D_grid_points.append(all_1D_grid_points_dim)

		return all_1D_grid_points

	def __get_all_indices(self, multiindex_set):

		indices = []

		max_level 		= np.max(multiindex_set)
		all_indices_1D 	= self.__get_all_1D_indices(max_level)

		for mindex in multiindex_set:

			level_indices = []
			for d in xrange(self.__dim):
				indices_1D 	= all_indices_1D[d][mindex[d] - 1]
				level_indices.append(indices_1D)

			fg_level_indices = self.__tensorize(level_indices)
			indices.extend(fg_level_indices)

		indices = np.array(indices, dtype=int)

		return indices

	def __get_surplus_indices(self, multiindex_set):

		indices = []

		max_level 		= np.max(multiindex_set)
		all_indices_1D 	= self.__get_all_1D_surplus_indices(max_level)
		for mindex in multiindex_set:

			level_indices = []
			for d in xrange(self.__dim):
				indices_1D 	= all_indices_1D[d][mindex[d] - 1]
				level_indices.append(indices_1D)

			fg_level_indices = self.__tensorize(level_indices)
			indices.extend(fg_level_indices)

		indices = np.array(indices, dtype=int)

		return indices

	def __get_no_fg_grid_points_mindex(self, multiindex):

		no_grid_points = 1
		for d in xrange(self.__dim):
			no_grid_points *= self.__get_no_1D_grid_points(multiindex[d])

		return no_grid_points

	def __get_no_surplus_grid_points_mindex(self, multiindex):

		no_grid_points = 1
		for d in xrange(self.__dim):
			if multiindex[d] != 1:
				no_surplus_points = self.__get_no_1D_grid_points(multiindex[d]) - \
										self.__get_no_1D_grid_points(multiindex[d] - 1) 
				no_grid_points *= no_surplus_points

		return no_grid_points

	def get_all_1D_points(self, multiindex_set, eps=1e-14):
		
		all_1D_grid_points 	= []

		max_level 				= np.max(multiindex_set)
		all_1D_grid_points_dim 	= self.__get_1D_full_points_per_dim(max_level, 0, eps)
		all_1D_grid_points.append(all_1D_grid_points_dim)

		for d1 in xrange(1, self.__dim):
			for d2 in xrange(self.__dim):
				if d1 != d2:
					if self.__weights[d1](self.__test_eval_x[d1]) - self.__weights[d2](self.__test_eval_x[d2]) == 0.:
						if d1 < d2:
							all_1D_grid_points_dim 	= all_1D_grid_points[d1]
						else:
							all_1D_grid_points_dim 	= all_1D_grid_points[d2]
							
						all_1D_grid_points.append(all_1D_grid_points_dim)
						break
					else:
						all_1D_grid_points_dim = self.__get_1D_full_points_per_dim(max_level, d1, eps)
						all_1D_grid_points.append(all_1D_grid_points_dim)

		return all_1D_grid_points

	def get_sg_surplus_points_multiindex(self, multiindex, eps=1e-14):

		max_level 			= np.max(multiindex)
		all_1D_grid_points 	= self.__get_all_1D_surplus_points(max_level, eps)	

		level_points 	= []
		for d in xrange(self.__dim):
			points_1D 	= all_1D_grid_points[d][multiindex[d] - 1]
			level_points.append(points_1D)

		sg_points = self.__tensorize(level_points)

		return sg_points

	def get_std_sg_surplus_points(self, multiindex_set, eps=1e-14):

		sg_points 	= []

		max_level 			= np.max(multiindex_set)
		all_1D_grid_points 	= self.__get_all_1D_surplus_points(max_level, eps)	
		for i, mindex in enumerate(multiindex_set):

			level_points 	= []
			for d in xrange(self.__dim):
				points_1D 	= all_1D_grid_points[d][mindex[d] - 1]
				level_points.append(points_1D)

			fg_level_points = self.__tensorize(level_points)
			sg_points.extend(fg_level_points)

		sg_points = np.array(sg_points, dtype=np.float64)

		return sg_points

	def map_std_sg_surplus_points(self, std_sg_points, left_stoch_boundary, right_stoch_boundary):

		mapped_sg_points 	= np.zeros(std_sg_points.shape)
		map_0_1_a_b 		= lambda x, d: left_stoch_boundary[d] + (right_stoch_boundary[d] - left_stoch_boundary[d])*x

		for i, sg_point in enumerate(std_sg_points):
			for j in xrange(self.__dim):
				mapped_sg_points[i, j] = map_0_1_a_b(sg_point[j], j)

		return mapped_sg_points

	def get_no_fg_grid_points(self, multiindex_set):

	    no_grid_points = np.zeros(len(multiindex_set), dtype=int)
	    
	    for i, index in enumerate(multiindex_set):
	        no_grid_points[i] = self.__get_no_fg_grid_points_mindex(index)

	    return no_grid_points

	def get_no_surplus_grid_points(self, multiindex_set):
		
		no_grid_points = np.zeros(len(multiindex_set), dtype=int)

		for i, index in enumerate(multiindex_set):
		    no_grid_points[i] = self.__get_no_surplus_grid_points_mindex(index)

		return no_grid_points

	def get_local_global_indices(self, multiindex_set):

		global_indices_dict = OrderedDict()

		indices_all 	= self.__get_all_indices(multiindex_set)
		indices_surplus = self.__get_surplus_indices(multiindex_set)

		for i, index_all in enumerate(indices_all):
			for j, index_surplus in enumerate(indices_surplus):

				if index_surplus.tolist() == index_all.tolist():
					global_indices_dict[i] = j
					break
				else:
					global_indices_dict[i] = i

		return global_indices_dict
