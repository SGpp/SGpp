#!/bin/bash
set -e
source jobconfig.sh

# hazelhen or pizdaint
machine=$1

if [[ ! ($machine == "pizdaint" || $machine == "hazelhen") ]]; then
    echo "error: specified machine \"$machine\" invalid";
    exit 1
fi;


echo "Create tests for $machine"
read -p "Continue? (y/n)" -n 1 -r
echo # newline
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    exit 1
fi

datasets=(strong_c10_dim2.arff strong_c10_dim5.arff strong_c10_dim8.arff strong_c10_dim10.arff)
walltimes_single_node=(240 480 3840 15360)
max_nodes=(4 8 64 128)
packetsize=50000


# overwrite possibly existing older result aggregator
echo "#!/bin/bash" > results_aggregator.sh

# overwrite possibly existing older trigger all script
echo "#!/bin/bash" > run-all-jobs.sh

for i in `seq 0 $((${#walltimes_single_node[@]}-1))`; do
    echo "echo \"dataset: ${datasets[$i]}\"" >> results_aggregator.sh
    ./create-strong-scaling-test.sh $machine $packetsize ${max_nodes[$i]} strong-scaling ${datasetFolder} ${datasets[$i]} ${walltimes_single_node[$i]}
done

# walltimes_weak_scaling_single_node=800
# dataset_size_single_node=20000
# echo "echo \"weak scaling:\"" >> results_aggregator.sh
# ./create-weak-scaling-test.sh $machine $packetsize ${max_nodes[$i]} weak-scaling ${datasetFolder} "weak_N%i_c10_size%i_dim10.arff" $walltimes_weak_scaling_single_node $dataset_size_single_node
