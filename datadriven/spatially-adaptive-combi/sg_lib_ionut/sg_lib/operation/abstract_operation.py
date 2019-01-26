import numpy as np
from abc import ABCMeta, abstractmethod
from itertools import product
from collections import OrderedDict
from pickle import dump, load
from sg_lib.algebraic.multiindex import *

class AbstractOperation(object):
    __metaclass__ = ABCMeta

    _dim = 0

    @property
    def dim(self):

        return self._dim

    def _get_differences_sign(self, multiindex):

        differences_indices = OrderedDict()
        differences_signs   = OrderedDict()

        possible_diff_indices = np.zeros((self._dim, 2))

        for d in xrange(self._dim):
            possible_diff_indices[d][0] = multiindex[d]
            possible_diff_indices[d][1] = multiindex[d] - 1

        differences_indices_temp = np.array(list(product(*possible_diff_indices)), dtype=int)

        key_indices = 0
        for d in differences_indices_temp:

            if all(i >= 1 for i in d):
                differences_indices[key_indices] = d
                key_indices                     += 1                
                
        possible_differences_signs = np.zeros((self._dim, 2))

        for d in xrange(self._dim):

            if multiindex[d] == 1:
                possible_differences_signs[d][0] = 1
                possible_differences_signs[d][1] = 0
            else:
                possible_differences_signs[d][0] = 1
                possible_differences_signs[d][1] = -1

        differences_signs_temp = np.array(list(product(*possible_differences_signs)), dtype=int)

        key_signs = 0
        for element in differences_signs_temp:

            if np.prod(element) != 0:   
                differences_signs[key_signs] = np.prod(element)
                key_signs                   += 1

        return differences_indices, differences_signs

    def _get_multiindex_dict(self, multiindex_set):

        multiindex_dict = OrderedDict()

        for index, multiindex in enumerate(multiindex_set):
            multiindex_dict[repr(multiindex.tolist())] = index

        return multiindex_dict

    def update_sg_evals_lut(self, func_eval):

        self._sg_func_evals_lut[self._sg_func_evals_key] = func_eval
        self._sg_func_evals_key += 1

    def update_func_evals(self, grid_object, multiindex_set):

        self._fg_func_evals = []

        self._sg_func_evals         = self._sg_func_evals_lut.values()
        self._grid_points_1D        = grid_object.get_all_1D_points(multiindex_set)
        self._global_indices_dict   = grid_object.get_local_global_indices(multiindex_set)
        self._no_fg_grid_points     = grid_object.get_no_fg_grid_points(multiindex_set)

        global_index = 0
        for index, multiindex in enumerate(multiindex_set):

            no_points           = self._no_fg_grid_points[index]
            level_func_evals    = np.zeros(no_points)

            for i in xrange(no_points): 
                level_func_evals[i] = self._sg_func_evals[self._global_indices_dict[global_index]]
                global_index        += 1
            
            self._fg_func_evals.append(level_func_evals)

    def reset_datastructures(self):

        self._fg_func_evals     = []
        self._sg_func_evals_key = 0
        self._sg_func_evals_lut = OrderedDict()  

    def serialize_data(self, serialization_file):
        
        with open(serialization_file, "wb") as output_file:
            dump(self._sg_func_evals_lut, output_file)

        output_file.close()

    def unserialize_data(self, serialization_file):

        data = []
        with open(serialization_file, "rb") as input_file:
            while True:
                try: 
                    data.append(load(input_file))
                except EOFError: 
                    break
        
        input_file.close() 

        self._sg_func_evals_lut = data[-1]
        self._sg_func_evals_key = len(self._sg_func_evals_lut)
        
    @abstractmethod
    def _eval_operation_fg(self):
        
        return

    @abstractmethod
    def eval_operation_delta(self):
        
        return 
 
    @abstractmethod
    def eval_operation_sg(self):
        
        return  