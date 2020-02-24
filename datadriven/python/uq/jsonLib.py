# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp import DataVector
import types

import numpy as np


def reprVal(attrValue):
    """
    Computes a string representation of attrValue that is
    valid json. Supported are:
    1) boolean
    2) integer
    3) dictionary
    4) float
    5) list, np.ndarray
    6) tuple
    7) string
    8) DataVector
    9) objects with ".toJson" method
    10) objects with ".serialize" method
    @param attrValue: arbitrary object
    """
    if isinstance(attrValue, bool):
        return reprBool(attrValue)
    elif isinstance(attrValue, int):
        return reprInt(attrValue)
    elif isinstance(attrValue, dict):
        return reprDict(attrValue)
    elif type(attrValue) in (float, np.float64, np.float32, np.float):
        return reprFloat(attrValue)
    elif isinstance(attrValue, list) or \
            isinstance(attrValue, np.ndarray):
        return reprList(attrValue)
    elif isinstance(attrValue, tuple):
        return reprTuple(attrValue)
    elif isinstance(attrValue, bytes) or \
            isinstance(attrValue, str):
        return reprString(attrValue)
    elif isinstance(attrValue, DataVector):
        return reprList(attrValue.array())
    elif 'toJson' in dir(attrValue):
        return attrValue.toJson()
    elif 'serialize' in dir(attrValue):
        text = attrValue.serialize()
        # hack to make the text json compatible
        return reprString(text.replace('\n', '__'))
    elif attrValue is None:
        return "null"
    else:
        raise AttributeError('jsonLib: reprVal - Unknown type "%s" of "%s"' %
                             (type(attrValue), attrValue))


def reprString(attrValue):
    return '"' + attrValue + '"'


def reprList(attrValue):
    s = [None] * len(attrValue)
    for i, val in enumerate(attrValue):
        s[i] = reprVal(val)
    return '[%s]' % ','.join(s)


def reprTuple(attrValue):
    s = [None] * len(attrValue)
    for i, val in enumerate(attrValue):
        s[i] = reprVal(val)
    return '"(%s)"' % ','.join(s)


def reprFloat(attrValue):
    attrValue = float(attrValue)
    if np.isnan(attrValue) or np.isinf(attrValue):
        return '"%g"' % attrValue
    else:
        return repr(attrValue)


def reprBool(attrValue):
    return '"%s"' % str(attrValue).lower()


def reprInt(attrValue):
    return str(attrValue)


def reprKey(attrValue):
    if attrValue.find('"') > -1:
        return attrValue
    else:
        return '"%s"' % attrValue


def reprDict(attrValue):
    s = [None] * len(attrValue)
    for i, (key, value) in enumerate(attrValue.items()):
        r = reprVal(key)
        s[i] = "%s: %s" % (reprKey(r), reprVal(value))

    return '{%s}' % ", ".join(s)


def parseAttribute(attrValue, attrName):
    """
    Parses a single class attribute to the correspondent string
    representation in json.
    @param attrValue: value of the attribute
    @param attrName: name of the attribute
    @return: json string equivalence
    """
    if attrName.find('__') != 0:
        try:
            s = reprVal(attrValue)
            return '%s: %s,\n' % (reprString(attrName), s)
        except AttributeError as e:
            print( e )
            return ''
    else:
        return ''


def parseKeyAsTuple(d):
    """
    Converts the string representation of a dict to an actual dict
    @param d: dictionary
    """
    ans = {}
    for key, value in list(d.items()):
        if (isinstance(key, bytes) or isinstance(key, str)) and \
                key[0] == "(" and key[-1] == ")":
            kt = stringToTupleOfFloats(key)

        if (isinstance(value, bytes) or isinstance(value, str)) and \
                value[0] == "(" and value[-1] == ")":
            vt = stringToTupleOfFloats(value)
        else:
            vt = value

        ans[kt] = vt

    return ans


def stringToTupleOfFloats(s):
    """
    Converts s to a tuple
    @param s: string
    @return: tuple represented by s
    """
    ans = []
    for i in s.strip("()").split(","):
        if i.strip() != "":
            if i == "null":
                ans.append(None)
            else:
                ans.append(float(i))
    return tuple(ans)


def stringToList(s, f=float):
    return list(map(f, s[1:-1].split(',')))


def stringToListOfLists(s, f=float):
    strs = s.replace('[', '').split('],')
    lists = [list(map(f, st.replace(']', '').split(','))) for st in strs]
    return lists
