#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    ASGCSamplerFormatter.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Fri Jul 19 13:08:22 2013

@brief   Provides functionality for the runtime serialization of the
         ASGC subclasses

         This design intends to separate the binary object
         representation and its business logic from the text
         representation that can be saved into file.  The class is a
         part of <a
         href='http://en.wikipedia.org/wiki/Memento_pattern'
         target='new'>Memento design pattern</a> described in details
         in @link
         bin.controller.CheckpointController.CheckpointController
         CheckpointController @endlink.

         Currently, the ASGC memento object is a dictionary with
         attributes describing the ASGC object.

@version  0.1

"""
from pysgpp.extensions.datadriven.utils.GzipSerializer import GzipSerializer
import types

import pysgpp.extensions.datadriven.utils.json as json


class ASGCSamplerFormatter(GzipSerializer):
    """
    ASGCSamplerFormatter class.
    """

    def deserialize(self, serializationStream):
        """
        Deserializes the ASGC memento object from the stream.

        Arguments:
        serializationStream -- the stream to deserialize

        Return The ASGC object
        """
        text = GzipSerializer.deserialize(self, serializationStream)
        jsonObject = json.JsonReader().read(text)
        return jsonObject

    def deserializeFromFile(self, filename):
        """
        Deserializes the UQSettting memento object from the file.
        The file may or may be not gzip compressed.
        @param filename: the name of the file with serialized ASGC
        @return: the ASGC memento object
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
        Serializes the ASGC memento to the stream
        @param memento: the ASGC memento structure
        @param serializationStream: output stream where memento should
                                     be serialized to
        """
        text = json.JsonWriter().write(memento)
        serializationStream.write(text)

    def serializeToFile(self, memento, filename):
        """
        Serializes the ASGC memento to file
        @param memento: The ASGC memento structure
        @param filename: The name of the file where the memento should be
                          serialized to
        """
        fstream = self.gzOpen(filename, "w")
        try:
            self.serialize(memento, fstream)
        finally:
            fstream.close()

    def toString(self, memento):
        """
        Returns a string that represents the ASGC memento object.
        @param memento: The ASGC memento object
        @return: string representing the ASGC memento
        """
        return json.JsonWriter().write(memento)
