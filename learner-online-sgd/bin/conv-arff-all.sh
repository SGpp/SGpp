#!/bin/bash

for f in ../data/*.arff
do
  echo "$f -> ${f}.mod"
  ./conv-arff-file.sed "$f" > "${f}.mod"
done
