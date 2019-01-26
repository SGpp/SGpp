import numpy as np
from abc import ABCMeta, abstractmethod
from collections import OrderedDict
from pickle import dump, load
from sg_lib.algebraic.multiindex import *

class DimensionAdaptivity(object):

	__metaclass__ = ABCMeta

	_dim 				= 0
	_init_multiindex 	= None

	@property
	def dim(self):

	    return self._dim

	@property
	def approx(self):

	    return self._approx

	@property
	def eta(self):

	    return np.asscalar(self._eta)

	@property
	def A(self):
		
		return self._A

	@property
	def O(self):
		
		return self._O

	@property
	def local_error(self):
		
		return self._local_error

	@property
	def multiindex_set(self):

		return self._multiindex_set

	@property
	def stop_adaption(self):

		return self._stop_adaption

	@property
	def local_basis_local(self):
		
		return self._local_basis_local

	@property
	def local_basis_global(self):
		
		return self._local_basis_global

	def _is_O_admissible(self, multiindex):

		admisability 	= 1
		old_index_set 	= [list(item) for item in self._O.values()]

		for q in xrange(self._dim):
			temp  	= np.zeros(self._dim, dtype=int)
			temp[q] = 1

			if multiindex[q] >= 2:
				diff = (multiindex - temp).tolist()
				
				if diff in old_index_set:
					admisability = 1
				else:
					admisability = 0
					break

		return admisability

	def _select_highest_priority_index(self):
		
		max_index = max(self._local_error, key=self._local_error.get)

		return max_index

	def _get_no_1D_grid_points(self, level):

		no_points = 0

		if self._level_to_nodes == 'sym':
			no_points = 2*level - 1
		else:
			if self._level_to_nodes == 1:
				no_points = level
			elif self._level_to_nodes >= 2:
				if level == 1:
					no_points = 1
				else:
					no_points = self._level_to_nodes*level - 1
			else:
				raise NotImplementedError

		return no_points

	def _get_local_hierarchical_basis(self, largest_basis):

		temp 			= np.eye(self._dim, dtype=int)
		first_basis 	= [0 for i in xrange(self._dim)]

		local_basis = [largest_basis.tolist()]
	 	 		
		for d in xrange(self._dim):
			curr_multiindex = local_basis[-1]

			temp_multiindex = np.array(curr_multiindex, dtype=int) - temp[d, :]
			temp_multiindex = temp_multiindex.tolist()

			if temp_multiindex[d] >= 0 and temp_multiindex not in local_basis:
				
				local_basis.append(temp_multiindex)

			curr_multiindex = local_basis[-1]

			if sum(curr_multiindex) == 0:
					break

		local_basis = local_basis[::-1]

		return local_basis

	def _update_local_basis(self, neighbor, local_basis_neighbor):

		neigh_basis 			= self._get_local_hierarchical_basis(local_basis_neighbor)
		local_basis_prev_step 	= self._local_basis_global[:]

		local_basis_curr_step = []
		for basis in local_basis_prev_step:
			
			basis_temp 	= np.array(basis, dtype=int)
			diff 		= local_basis_neighbor - basis_temp
			
			if ((diff >= 0).sum() == diff.size).astype(int):
				local_basis_curr_step.append(basis)

		for basis in neigh_basis:
			if basis not in self._local_basis_global:
				local_basis_curr_step.append(basis)
				self._local_basis_global.append(basis)

		self._local_basis_local[repr(neighbor)] = local_basis_curr_step

	def serialize_multiindices_adapt(self, new_multiindices, serialization_file):

		with open(serialization_file, "wb") as output_file:
			dump(new_multiindices, output_file)

		output_file.close()

	def unserialize_multiindices_adapt(self, serialization_file):

		with open(serialization_file, "rb") as input_file:
			new_multiindices = load(input_file)

		input_file.close()

		return new_multiindices

	@abstractmethod
	def init_adaption(self):

		return 

	@abstractmethod
	def do_one_adaption_step_preproc(self):

		return
		
	@abstractmethod
	def do_one_adaption_step_postproc(self):

		return

	@abstractmethod
	def check_termination_criterion(self):

		return

	@abstractmethod
	def serialize_data(self):
		
		return

	@abstractmethod
	def unserialize_data(self):

		return