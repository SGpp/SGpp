#!/bin/bash
# parameters are: nodeName floatOnly

echo "executing:"
echo "./datadriven/performanceTests/test_datadriven_boost --run_test=AutoTuningPaper --log_level=all > $1.log 2>&1"
echo "./datadriven/examplesOCL/evaluateConsistency $1 $2 > evaluateConsistency_$1.log 2>&1"

./datadriven/performanceTests/test_datadriven_boost --run_test=AutoTuningPaper --log_level=all > $1.log 2>&1
./datadriven/examplesOCL/evaluateConsistency $1 $2 > evaluateConsistency_$1.log 2>&1