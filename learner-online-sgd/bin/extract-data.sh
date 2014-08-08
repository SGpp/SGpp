#!/bin/bash

for f in ../data/*.gz
do
  echo $f
  gunzip "$f"
done

rm ../data/*.gz
