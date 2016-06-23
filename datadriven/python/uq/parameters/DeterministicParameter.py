from Parameter import Parameter


class DeterministicParameter(Parameter):
    """
    Deterministic parameter
    """

    def __init__(self, name, value, trans):
        super(DeterministicParameter, self).__init__(name)

        self._value = value
        self._trans = trans

    def getCount(self):
        return 1

    def isUncertain(self):
        return False

    def getSample(self):
        return [self._value]

    def __str__(self):
        return "%s = %s" % (self._name, self._value)

    # @classmethod
    # def fromJson(cls, jsonObject):
    #     """
    #     Restores the DeterministicParameter object from the json object
    #     with its attributes.

    #     Arguments:
    #     jsonObject -- json object

    #     Return the restored DeterministicParameter object
    #     """
    #     # restore surplusses
    #     key = '_name'
    #     if jsonObject.isContaining(key):
    #         name = jsonObject[key]

    #     key = '_value'
    #     if jsonObject.isContaining(key):
    #         value = jsonObject[key]

    #     return DeterministicParameter(name, value)


    # def toJson(self):
    #     """
    #     Returns a string that represents the object
    #     """
    #     serializationString = '"module" : "' + \
    #                           self.__module__ + '",\n'

    #     # serialize
    #     for attrName in ('_name', '_value'):
    #         attrValue = self.__getattribute__(attrName)
    #         serializationString += ju.parseAttribute(attrValue, \
    #                                                  attrName)


    #     s = serializationString.rstrip(",\n").replace("'",'"')

    #     return "{" + s + "}"
