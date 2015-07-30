# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+'/src/python')

from src.python import learner as learner
sys.modules['datadriven.learner'] = learner

from src.python import controller as controller
sys.modules['datadriven.controller'] = controller

from src.python import data as data
sys.modules['datadriven.data'] = data

from src.python import uq as uq
sys.modules['datadriven.uq'] = uq

from src.python import tools as tools
sys.modules['datadriven.tools'] = tools
