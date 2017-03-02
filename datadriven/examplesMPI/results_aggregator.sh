#!/bin/bash
echo "dataset: strong_c10_dim2.arff"
./get-ergs.sh strong-scaling-run-strong_c10_dim2
echo "dataset: strong_c10_dim5.arff"
./get-ergs.sh strong-scaling-run-strong_c10_dim5
echo "dataset: strong_c10_dim8.arff"
./get-ergs.sh strong-scaling-run-strong_c10_dim8
echo "dataset: strong_c10_dim10.arff"
./get-ergs.sh strong-scaling-run-strong_c10_dim10
