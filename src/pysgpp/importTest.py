#!/usr/bin/python

print "executing pysgpp import test (pysgpp/importTest.py)"
print "trying to import python lib (pysgpp) ...",

try:
    import pysgpp
except Exception as e:
    print "failed, error:"
    print str(e)
print "success"

# import os
# print os.environ['LD_LIBRARY_PATH']
# print os.environ['PYTHONPATH']

#content = dict([(name, cls) for name, cls in pysgpp.__dict__.items() if isinstance(cls, type)])
