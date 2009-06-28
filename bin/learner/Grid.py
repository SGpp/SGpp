#!/usr/bin/env python
# -*- coding: utf-8 -*-


class Grid(object):
    """ generated source for Grid

    """
    level = Integer()
    size = Integer()
    polynom = Integer()
    border = BorderTypes()
    points = DataVector()
    values = DataVector()
    myCheckpointController = CheckpointController()
    grid = Learner()
    myCheckpointController = CheckpointController()
    adapter = GridAdapter()

    def save(self, iteration):
        return

    @classmethod
    def refine(cls, numOfPoints):
        return

    def eval(self, point):
        return

    def load(self, source):
        return

    def addPoint(self, point):
        return

    def deletePoint(self, point):
        return

    def Grid(self):
        return

    def setLevel(self, level):
        return

    def setPolynom(self, deg):
        return

    def setBorder(self, type):
        return

    def setAdapter(self, adapter):
        return

    def getNumberOfRefinablePoints(self):
        return

    def createCOperator(self):
        return

    def createBOperator(self):
        return


