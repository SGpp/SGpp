#!/bin/bash
set -e
source jobconfig.sh

datasets=(strong_c10_dim2.arff strong_c10_dim5.arff strong_c10_dim8.arff strong_c10_dim10.arff)
walltimes_single_node=(1000 2500 7000 37000)
packetsize=10000
max_nodes=1024

# overwrite possibly existing older result aggregator
echo "#!/bin/bash" > results_aggregator.sh

# overwrite possibly existing older trigger all script
echo "#!/bin/bash" > run-all-jobs.sh

for i in `seq 0 $((${#walltimes_single_node[@]}-1))`; do
    echo "echo \"dataset: ${datasets[$i]}\"" >> results_aggregator.sh
    ./create-strong-scaling-test.sh $packetsize $max_nodes strong-scaling ${datasetFolder} ${datasets[$i]} ${walltimes_single_node[$i]}
done

walltimes_weak_scaling_single_node=800
dataset_size_single_node=20000
echo "echo \"weak scaling:\"" >> results_aggregator.sh
./create-weak-scaling-test.sh $packetsize $max_nodes weak-scaling ${datasetFolder} "weak_N%i_c10_size%i_dim10.arff" $walltimes_weak_scaling_single_node $dataset_size_single_node
