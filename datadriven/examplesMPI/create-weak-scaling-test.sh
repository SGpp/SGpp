#!/bin/bash
set -e
# Create folder for the testrun
# today=`date +%Y-%m-%d_%H:%M`
machine=$1
packagesize=$2
max_nodes=$3
name=$4
datasetFolder=$5
datasetPattern=$6
walltime_single_node=$7
dataset_size_single_node=$8
test_folder="weak-scaling"
echo "Creating testrun for log_2(${max_nodes}) cases into folder ${test_folder}..."
read -p "Continue? (y/n)" -n 1 -r
echo # newline
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    exit 1
fi

mkdir ${test_folder}
nodes=1
dataset_size=$dataset_size_single_node
walltime=${walltime_single_node}
while [ ${nodes} -le ${max_nodes} ]; do
# for i in `seq 1 $2`; do
    echo Creating config for ${nodes} nodes
		dataset=`printf "$datasetPattern" $nodes $dataset_size`
		echo Dataset: $dataset
    # Create folder
    mkdir ${test_folder}/Nodes-${nodes}
    # Create MPI config file
    ./generateConfigFiles $packagesize ${nodes} ${test_folder}/Nodes-${nodes}/MPIConf.cfg
    # Copy OCL config file
    cp MyOCLConf.cfg ${test_folder}/Nodes-${nodes}/MyOCLConf.cfg
    # Copy executable
    cp ./mpi_examples ${test_folder}/Nodes-${nodes}/mpi_examples
    # Load runscript template (double qoutes preserver newline characters)
    if [[ $machine == "pizdaint" ]]; then
        templatefile="$(cat template-pizdaint)"
    elif [[ $machine == "hazelhen" ]]; then
        templatefile="$(cat template-hazelhen)"
    fi
    # Replace placeholders
    procs=$((${nodes}+1))

    sec=$((walltime%60))
    walltime_in_min=${walltime}
    walltime_in_min=$((walltime_in_min/60))
    min=$((walltime_in_min%60))
    hrs=$((walltime_in_min/60))
    
    formatted_walltime=$(printf "%02d:%02d:%02d" $hrs $min $sec)
    printf -v templatefile "$templatefile" "Clustering-Weak" "$procs" "$formatted_walltime" "${datasetFolder}/${dataset}"
    # Print and save to file
    echo "$templatefile"
    if [[ $machine == "pizdaint" ]]; then
        echo "$templatefile" > ${test_folder}/Nodes-${nodes}/runscript.sh
    elif [[ $machine == "hazelhen" ]]; then
        echo "$templatefile" > ${test_folder}/Nodes-${nodes}/runscript.pbs
    fi
    let nodes=nodes*2
		let dataset_size=dataset_size*2
    walltime=`echo "s=$walltime; s/1.8" | bc`
    # let walltime=walltime/2
done

if [[ $machine == "pizdaint" ]]; then
    runtemplatefile="$(cat run-template-pizdaint)"
elif [[ $machine == "hazelhen" ]]; then
    runtemplatefile="$(cat run-template-hazelhen)"
fi

# Replace placeholders
printf -v runtemplatefile "$runtemplatefile" "$2"
# Print and save to file
echo "$runtemplatefile"
echo "$runtemplatefile" > ${test_folder}/run.sh

# add entry in result aggregator
echo "./get-ergs.sh ${test_folder}" >> results_aggregator.sh

# add entry in run-all script
echo "cd ${test_folder}" >> run-all-jobs.sh
echo "bash run.sh" >> run-all-jobs.sh
echo "cd .." >> run-all-jobs.sh
