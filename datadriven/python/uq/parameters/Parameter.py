from pysgpp.extensions.datadriven.uq.dists import Dist

import pysgpp.extensions.datadriven.uq.jsonLib as ju
from pysgpp.extensions.datadriven.uq.transformation.Transformation import Transformation


class Parameter(object):
    """
    Super-class for parameters
    """

    def __init__(self, name):
        self._name = name
        self._value = None
        self._trans = None

    def getName(self):
        return self._name

    def setName(self, name):
        self._name = name

    def getCount(self):
        raise NotImplementedError()

    def hasValue(self, value):
        self._value = value

    def getProbabilisticValue(self):
        return self._value

    def getUnitValue(self):
        return self._trans.probabilisticToUnit(self._value)

    def getTransformation(self):
        return self._trans

    def isUncertain(self):
        raise NotImplementedError()

    def isActive(self):
        return self.isUncertain()

    def toJson(self):
        """
        Returns a string that represents the object
        """
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        # serialize
        for attrName in ('_name', '_value', '_trans'):
            attrValue = self.__getattribute__(attrName)
            serializationString += ju.parseAttribute(attrValue, attrName)

        # serialize dist
        if self.isUncertain():
            attrName = "_dist"
            attrValue = self.__getattribute__(attrName)
            serializationString += ju.parseAttribute(attrValue, attrName)

        return "{" + serializationString.rstrip(",\n") + "}"

    @classmethod
    def fromJson(cls, jsonObject):
        """
        Restores the Parameter object from the json object
        with its attributes.
        @param jsonObject: json object
        @return: the restored Parameter object
        """
        key = '_name'
        if key in jsonObject:
            name = jsonObject[key]

        key = '_value'
        if key in jsonObject and jsonObject[key] != "null":
            value = jsonObject[key]
        else:
            value = None

        key = '_trans'
        if key in jsonObject:
            trans = Transformation.fromJson(jsonObject[key])

        key = '_dist'
        if key in jsonObject:
            dist = Dist.fromJson(jsonObject[key])

        if 'parameters.UncertainParameter' in jsonObject['module']:
            from pysgpp.extensions.datadriven.uq.parameters.UncertainParameter import UncertainParameter
            return UncertainParameter(name, dist, trans, value)
        elif 'parameters.DeterministicParameter' in jsonObject['module']:
            from pysgpp.extensions.datadriven.uq.parameters.DeterministicParameter import \
                DeterministicParameter
            return DeterministicParameter(name, value)
        else:
            raise TypeError("Unknown parameter '%s' => Please register it in fromJson function in Parameters.py" % jsonObject["module"])
