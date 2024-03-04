#!/bin/bash
# Modify according to your data - array job is defined on line 60. 

# Set the number of CPU cores and memory (modify as needed)
ncpu=$(nproc)  # Get the number of available CPU cores
mem="10g"     # Set the memory limit

# Function to process a specific id
process_id() {
    local id=$1

    # Get list of species (modify the directory path as needed)
    indirs=($(ls /home/MitoFinder/2_cleanedreads_PMC))

    # Get the specific sample
    sample="${indirs[$id]}"
    # Load data
    r1=$(echo "/home/MitoFinder/2_cleanedreads_PMC/$sample/split-adapter-quality-trimmed/$sample-READ1.fastq.gz")
    r2=$(echo "/home/MitoFinder/2_cleanedreads_PMC/$sample/split-adapter-quality-trimmed/$sample-READ2.fastq.gz")


    # echo $(date)
    # echo "------------------------------------------------------"
    # echo "This job is allocated $ncpu CPU core(s) and $mem of memory"
    # echo "------------------------------------------------------"
    # echo "Running on host $(hostname)"
    # echo "Working directory is $(pwd)"
    # echo "Sample index is $id"
    # echo "Sample name is $sample"
    # echo "------------------------------------------------------"

    # Load necessary modules (modify as needed)
    module load java/openjdk8
    module load blast/2.12.0
    module load spades/3.15.5
    #source activate megahit
    #source activate spades-3.15.2

    # Get into the directory where you want to store your results
    cd /scratch/MitoFinder/MF_RUN1

    # Run mitofinder
    mitofinder --metaspades \
    -j $sample \
    -1 $r1 \
    -2 $r2 \
    -r PMC_mtRefSeq.gb \
    -o 2 \
    -p 20 \
    -t mitfi \
    -e 0.000001 \
    --min-contig-size 500 \
    --max-memory 100 \
    --new-genes
}

export -f process_id

# Loop through the range of values from 0 to 67 and run jobs in parallel
for id in {0..67}; do
    process_id "$id" &
done

# Wait for all background processes to finish
wait

mkdir logs
mv *.log logs
wait

echo "All samples have been processed!"
