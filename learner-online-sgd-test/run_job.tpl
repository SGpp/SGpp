#!/bin/bash

DIR=../results/cond$cond-num$num-ref_type$ref_type

if [ -d "$DIR" ]; then 
  echo "Directory $DIR already exists"
else 
  mkdir -p "$DIR"
fi

cd $DIR

MAIN=../../main

$MAIN $cond $num $ref_type | tee log 2>&1 
