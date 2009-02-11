# This file is part of sg++, a program package making use of spatially adaptive sparse grids to solve numerical problems.
#
# Copyright (C) 2007  Joerg Blank (blankj@in.tum.de)
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with pyclass. If not, see <http://www.gnu.org/licenses/>.
#

import unittest

import test_GridIndex
import test_GridStorage
import test_algorithms
import test_laplace
import test_hierarchisation

import test_GridFactory
import test_DataVector

if __name__ == '__main__':
    alltests = unittest.TestSuite([
         unittest.defaultTestLoader.loadTestsFromModule(test_GridIndex),
         unittest.defaultTestLoader.loadTestsFromModule(test_GridStorage),
         unittest.defaultTestLoader.loadTestsFromModule(test_algorithms),
         unittest.defaultTestLoader.loadTestsFromModule(test_laplace),
         unittest.defaultTestLoader.loadTestsFromModule(test_GridFactory),
         unittest.defaultTestLoader.loadTestsFromModule(test_DataVector),
         unittest.defaultTestLoader.loadTestsFromModule(test_hierarchisation),
        ])

    unittest.TextTestRunner().run(alltests)


    
