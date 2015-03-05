#!/bin/sh

for i in `find ./../src/ -name "*.?pp"`; do astyle --options=astylerc $i ; done
