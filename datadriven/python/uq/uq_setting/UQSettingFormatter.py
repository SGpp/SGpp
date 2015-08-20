#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    UQSettingFormatter.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Fri Jul 19 13:08:22 2013

@brief   Provides functionality for the runtime serialization of the
         UQSetting subclasses

         This design intends to separate the binary object
         representation and its business logic from the text
         representation that can be saved into file.  The class is a
         part of <a
         href='http://en.wikipedia.org/wiki/Memento_pattern'
         target='new'>Memento design pattern</a> described in details
         in @link
         bin.controller.CheckpointController.CheckpointController
         CheckpointController @endlink.

         Currently, the UQSetting memento object is a dictionary with
         attributes describing the UQSetting object.

@version  0.1

"""

from pysgpp.extensions.datadriven.utils.GzipSerializer import GzipSerializer
import json
import types
import pysgpp.extensions.datadriven.utils.json as json_manual


class UQSettingFormatter(GzipSerializer):
    """
    UQSettingFormatter class.
    """

    def deserialize(self, serializationStream):
        """
        Deserializes the UQSetting memento object from the stream.

        Arguments:
        serializationStream -- the stream to deserialize

        Return The UQSetting object
        """
        # text = GzipSerializer.deserialize(self, serializationStream)
        # jsonObject = json.JsonReader().read(text)
        jsonObject = json.load(serializationStream, encoding='utf8')
        return jsonObject

    def deserializeFromFile(self, filename):
        """
        Deserializes the UQSettting memento object from the file.
        The file may or may be not gzip compressed.

        Arguments:
        filename -- The name of the file with serialized UQSetting

        Return The UQSetting memento object
        """
        if not isinstance(filename, types.StringType):
            raise AttributeError("Filename as destination expected")
        serializationStream = self.gzOpen(filename, "r")
        try:
            memento = self.deserialize(serializationStream)
        finally:
            serializationStream.close()

        return memento

    def serialize(self, memento, serializationStream):
        """
        Serializes the UQSetting memento to the stream

        Arguments:
        memento -- The UQSetting memento structure
        serializationStream -- output stream where memento should
                               be serialized to

        Return None
        """
        # text = json_manual.JsonWriter().write(memento)
        text = json.dumps(memento, indent=None)
        serializationStream.write(text)

    def serializeToFile(self, memento, filename):
        """
        Serializes the UQSetting memento to file

        Arguments:
        memento -- The UQSetting memento structure
        filename -- The name of the file where the memento should be
                    serialized to

        Return None
        """
        fstream = self.gzOpen(filename, "w")
        try:
            self.serialize(memento, fstream)
        finally:
            fstream.close()

    def toString(self, memento):
        """
        Returns a string that represents the UQSetting memento object.

        Arguments:
        memento -- The UQSetting memento object

        Return String representing the UQSetting memento
        """
        return json.dumps(memento, indent=4)


