import numpy as np
from collections import OrderedDict

class Multiindex(object):
	
	def __init__(self, dim):
		self._dim = dim

		self.__poly_mindex_dict 	= OrderedDict()
		self.__poly_mindex_dict_key = 0
		self.__poly_prev_pos 		= 0

		self.__poly_mindex_dict[self.__poly_mindex_dict_key] = [0 for d in xrange(dim)]

	def __get_l1_norm(self, vec):
		
		return np.sum(vec)

	def get_successors(self, multiindex):
		
		successors = np.zeros((self._dim, self._dim), dtype=int)

		for i in xrange(self._dim):
			temp 				= np.zeros(self._dim, dtype=int)
			temp[i] 			= 1
			successors[i, :]	= multiindex + temp

		return successors

	def get_successors_right_left(self, multiindex):
		
		successors = np.zeros((self._dim, self._dim), dtype=int)

		for i in xrange(self._dim):
			temp 					= np.zeros(self._dim, dtype=int)
			temp[self._dim - i - 1] = 1
			successors[i, :]		= multiindex + temp

		return successors

	def get_std_total_degree_mindex(self, level):
		
		multiindex_set = []

		init_multiindex = np.ones(self._dim).tolist()
		multiindex_set.append(init_multiindex)

		prev_pos = 0
		curr_pos = len(multiindex_set)

		not_finished = True

		if level == 1:
			not_finished = False

		while not_finished:
			for i in xrange(prev_pos, curr_pos):

				successors = self.get_successors(multiindex_set[i]).tolist()
				for successor in successors:
					if np.sum(successor) <= self._dim + level - 1 and successor not in multiindex_set:
						multiindex_set.append(successor)

					if np.sum(successor) > self._dim + level - 1:
						not_finished = False

			prev_pos = curr_pos
			curr_pos = len(multiindex_set)

		multiindex_set = np.array(multiindex_set, dtype=int)

		return multiindex_set

	def get_poly_mindex(self, level):

		multiindex_set = []

		init_multiindex = np.zeros(self._dim).tolist()
		multiindex_set.append(init_multiindex)

		prev_pos = 0
		curr_pos = len(multiindex_set)

		not_finished = True
		
		while not_finished:
			for i in xrange(prev_pos, curr_pos):

				successors = self.get_successors(multiindex_set[i]).tolist()
				for successor in successors:
					if np.sum(successor) <= level and successor not in multiindex_set:
						multiindex_set.append(successor)

					if np.sum(successor) > level:
						not_finished = False

			prev_pos = curr_pos
			curr_pos = len(multiindex_set)

		multiindex_set = np.array(multiindex_set, dtype=int)

		return multiindex_set

	def get_poly_degs(self, max_degs):

		multiindex_set = []

		init_multiindex = np.zeros(self._dim).tolist()
		multiindex_set.append(init_multiindex)

		prev_pos = 0
		curr_pos = len(multiindex_set)

		not_finished = True
		
		while not_finished:
			for i in xrange(prev_pos, curr_pos):

				successors = self.get_successors(multiindex_set[i]).tolist()
				for successor in successors:
					if successor not in multiindex_set:
						truth = 1

						for d in xrange(self._dim):
							if successor[d] > max_degs[d]:
								truth = 0
								not_finished = False
						if truth:
							multiindex_set.append(successor)						

			prev_pos = curr_pos
			curr_pos = len(multiindex_set)

		multiindex_set = np.array(multiindex_set, dtype=int)

		return multiindex_set

	def get_poly_mindex_experiment(self, level):

		not_finished = True

		curr_pos = len(self.__poly_mindex_dict.values())

		print 'curr values'
		print self.__poly_mindex_dict.values()
		
		while not_finished:

			for i in xrange(self.__poly_prev_pos, curr_pos):

				curr_mindex = self.__poly_mindex_dict[i]
				successors 	= self.get_successors(curr_mindex).tolist()

				for successor in successors:
					if np.sum(successor) <= level and successor not in self.__poly_mindex_dict.values():

						self.__poly_mindex_dict_key += 1
						self.__poly_mindex_dict[self.__poly_mindex_dict_key] = successor

					if np.sum(successor) > level:
						not_finished = False

			self.__poly_prev_pos 	= curr_pos
			curr_pos 				= len(self.__poly_mindex_dict.values())

		multiindex_set = np.array(self.__poly_mindex_dict.values(), dtype=int)

		return multiindex_set