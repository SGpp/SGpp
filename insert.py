#!/usr/bin/python

import re
import os

print "hello :*)"

files = os.listdir(".")

files = [i for i in files if i.endswith("_tuned.cfg")]

# files = ["testing_tuned.cfg"]

for fileName in files:
    f = open(fileName, "r+")
    fixedString = ""
    for line in f:
        # print line,
        m = re.search(r"\"StreamingOCLMultiPlatform\": {", line)
        if m:
            print "has match"
            line = re.sub(r"\"StreamingOCLMultiPlatform\": {",
                              "\"StreamingOCLMultiPlatform\": {\n" +
                              "                     " +
                              "\"KERNEL_SCHEDULE_SIZE\": 1024000,", line)
            
        m = re.search(r"\"StreamingModOCLMaskMultiPlatform\": {", line)
        if m:
            print "has match (mod)"
            line = re.sub(r"\"StreamingModOCLMaskMultiPlatform\": {",
                              "\"StreamingModOCLMaskMultiPlatform\": {\n" +
                              "                     " +
                              "\"KERNEL_SCHEDULE_SIZE\": 1024000,", line)
        m = re.search(r"\"OPTIMIZATION_FLAGS\": -cl-strict-aliasing -cl-fast-relaxed-math",
                      line)
        if m:
            print "fixing optimization flags"
            line = re.sub(r"\"OPTIMIZATION_FLAGS\": -cl-strict-aliasing -cl-fast-relaxed-math",
                          "\"OPTIMIZATION_FLAGS\": \"-cl-strict-aliasing -cl-fast-relaxed-math\"",
                          line)
        fixedString += line                
    f.close()

    f = open(fileName, "w")
    f.write(fixedString)
    f.close()
