
import math
import numpy as np
import ErrorCalculator
import abc,logging
from combiScheme import CombiScheme


# This class implements a general container that can be filled with refinementObjects (typically specified by the refinement strategy)
# In addition it stores accumulated values over all refinementObjects (like integral, numberOfEvaluations)
class RefinementContainer(object):
    def __init__(self, initial_objects, dim, error_estimator):
        self.refinementObjects = initial_objects
        self.dim = dim
        self.evaluationstotal = 0
        self.integral = 0
        self.popArray = []
        self.startNewObjects = 0
        self.errorEstimator = error_estimator
        self.searchPosition = 0

    # returns the error that is associated with the specified refinementObject
    def get_error(self, object_id):
        return self.refinementObjects[object_id].error

    # refines the specified refinementObject
    def refine(self, object_id):
        # at first refinement in current refinement round we have
        # to save where the new RefinementObjects start
        if self.startNewObjects == 0:
            self.startNewObjects = len(self.refinementObjects)
        # refine RefinementObject;
        # returns the new RefinementObjects a possible update to lmax and update information for other RefinementObjects
        new_objects, lmax_update, update_information = self.refinementObjects[object_id].refine()
        # update other RefinementObjects if necessary
        if update_information is not None:
            for r in self.refinementObjects:
                r.update(update_information)
        # remove refined (and now outdated) RefinementObject
        self.prepare_remove(object_id)
        # add new RefinementObjects
        self.add(new_objects)
        return lmax_update

    def refinement_postprocessing(self):
        self.searchPosition = 0

    # if strategy decides from outside to update elements this function can be used
    def update_objects(self, update_info):
        for r in self.refinementObjects:
            r.update(update_info)

    # reset everything so that all RefinementObjects will be iterated
    def reinit_new_objects(self):
        self.startNewObjects = 0
        self.integral = 0
        self.evaluationstotal = 0

    # return the maximal error among all RefinementObjects
    def get_max_error(self):
        max_error = 0
        for i in self.refinementObjects:
            if i.error > max_error:
                max_error = i.error
        return max_error

    # return the maximal error among all RefinementObjects
    def get_total_error(self):
        total_error = 0
        for i in self.refinementObjects:
            total_error += i.error
        return total_error

    # indicate that all objects have been processed and new RefinementObjects will be added at the end
    def clear_new_objects(self):
        self.startNewObjects = len(self.refinementObjects)

        # returns only newly added RefinementObjects

    def get_new_objects(self):
        return self.refinementObjects[self.startNewObjects:]

    # returns amount of newly added RefinementObjects
    def new_objects_size(self):
        return len(self.refinementObjects) - self.startNewObjects

    # prepares removing a RefinementObject (will be removed after refinement round)
    def prepare_remove(self, objectID):
        self.popArray.append(objectID)

    # remove all RefinementObjects that are outdated from container
    def apply_remove(self):
        for position in reversed(sorted(self.popArray)):
            self.integral -= self.refinementObjects[position].integral
            self.evaluationstotal -= self.refinementObjects[position].evaluations
            self.refinementObjects.pop(position)
            if self.startNewObjects != 0:
                self.startNewObjects -= 1
        self.popArray = []

    # add new RefinementObjects to the container
    def add(self, new_refinement_objects):
        self.refinementObjects.extend(new_refinement_objects)

    # calculate the error according to the error estimator for specified RefinementObjects
    def calc_error(self, object_id, f):
        refine_object = self.refinementObjects[object_id]
        refine_object.set_error(self.errorEstimator.calc_error(f, refine_object))

    # returns all RefinementObjects in the container
    def get_objects(self):
        return self.refinementObjects

    def get_next_object_for_refinement(self, tolerance):
        if self.startNewObjects == 0:
            end = self.size()
        else:
            end = self.startNewObjects
        for i in range(self.searchPosition, end):
            if self.refinementObjects[i].error >= tolerance:
                self.searchPosition = i + 1
                return True, i, self.refinementObjects[i]
        return False, None, None

    # returns the specified RefinementObject from container
    def getObject(self, object_id):
        return self.refinementObjects[object_id]

    # returns amount of RefinementObjects in container
    def size(self):
        return len(self.refinementObjects)

    # sets the number of evaluations associated with specified RefinementObject
    def set_evaluations(self, object_id, evaluations):
        # add evaluations also to total number of evaluations
        self.evaluationstotal += evaluations
        self.refinementObjects[object_id].evaluations = evaluations

    # sets the integral for area associated with specified RefinementObject
    def set_integral(self, object_id, integral):
        # also add integral to global integral value
        self.integral += integral
        self.refinementObjects[object_id].set_integral(integral)


# this class defines a container of refinement containers for each dimension in the single dimension test case
# it delegates methods to subcontainers and coordinates everything
class MetaRefinementContainer(object):
    def __init__(self, refinement_containers):
        self.refinementContainers = refinement_containers
        self.evaluations = 0
        self.integral = None

    # return the maximal error among all RefinementContainers
    def get_max_error(self):
        max_error = 0
        for c in self.refinementContainers:
            error = c.get_max_error()
            if max_error < error:
                max_error = error
        return max_error

    # sets the integral for area associated with whole meta container
    def set_integral(self, objectID, integral):
        self.integral = integral

    # sets the number of evaluations associated with whole meta container
    def set_evaluations(self, objectID, evaluations):
        self.evaluations = evaluations

    # delegate to containers
    def reinit_new_objects(self):
        for c in self.refinementContainers:
            c.reinit_new_objects()

    def size(self):
        return 1

    def new_objects_size(self):
        return 1

    def clear_new_objects(self):
        pass

    def get_next_object_for_refinement(self, tolerance):
        pass
        # toDo

    # delegate to containers
    def apply_remove(self):
        for c in self.refinementContainers:
            c.apply_remove()

    def get_refinement_container_for_dim(self, d):
        return self.refinementContainers[d]

    # apply refinement
    def refine(self, position):
        lmax_change = self.refinementContainers[position[0]].refine(position[1])
        if lmax_change is not None:
            for d, c in enumerate(self.refinementContainers):
                if d != position[0]:
                    c.update(1)
        return lmax_change

'''
class RefinementContainerCell(RefinementContainer):
    def __init__(self, initial_objects, dim, error_estimator, a, b):
        super().__init__(self, initial_objects, dim, error_estimator)
        self.grid_dictionary = {}
        index_set = CombiScheme.get_index_set()
        for object in initial_objects:
            if tuple(object.levelvec) in self.grid_dictionary:
                self.grid_dictionary[tuple(object.levelvec)].append((object.start,object.end))
            else:
                self.grid_dictionary[tuple(object.levelvec)] = [(object.start,object.end)]

    def refinement_postprocessing(self, new_objects, object_id):
        for object in new_objects:
            previous_levelvec = [max(CombiScheme.lmin,object.levelvec - 1) for d in range(self.dim)]

            self.grid_dictionary[tuple(object.levelvec)] = self.grid_dictionary[tuple(previous_levelvec)]

            area_info_old = (self.grid_dictionary[tuple(self.getObject(object_id).levelvec)].start, self.grid_dictionary[tuple(self.getObject(object_id).levelvec)].end)
            self.grid_dictionary[tuple(object.levelvec)].remove(area_info_old)

        for object in new_objects:
            self.grid_dictionary[tuple(object.levelvec)].append((object.start, object.end))

    def set_integral(self, object_id, integral):
        self.integral += integral

    def set_evaluations(self, object_id, evaluations):
        self.evaluationstotal += evaluations

    def calc_error(self, object_id, f):
        pass
'''


