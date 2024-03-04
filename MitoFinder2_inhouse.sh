#!/bin/bash
# Modify according to your data - array job is defined on line 55

# Set the number of CPU cores and memory (modify as needed)
ncpu=$(nproc)  # Get the number of available CPU cores
mem="10g"     # Set the memory limit

# Function to process a specific id
process_id() {
    local id=$1

    # get list of contigs from first mitofinder run
    indirs=($(ls /scratch/MitoFinder/PMC_contigs_MFRUN1))
    # get infile
    sample="${indirs[$id]}"
    # MFRUN1 contigs
    contigs="/scratch/MitoFinder/PMC_contigs_MFRUN1/$sample/${sample}_MFRUN1_mtDNA_contigs.fasta"
    
    # echo $(date)
    # echo "------------------------------------------------------"
    # echo "This job is allocated $ncpu CPU core(s) and $mem of memory"
    # echo "------------------------------------------------------"
    # echo "Running on host $(hostname)"
    # echo "Working directory is $(pwd)"
    # echo "Sample index is $id"
    # echo "Sample name is $sample"
    # echo "------------------------------------------------------"

    # Get into the directory where you want to store your results
    cd /scratch/MitoFinder/MF_RUN2
    # Load necessary modules (modify as needed)
    module load java/openjdk8
    module load blast/2.12.0
    module load spades/3.15.5

    # run mitofinder to assemble with metaspades
    mitofinder -j $sample \
        -a $contigs \
        -r PMC_mtRefSeq.gb \
        -o 2 \
        -p 20 \
        -m 40 \
        -t mitfi \
        -e 0.000001 \
        --new-genes \
        --allow-intron \
        --intron-size 12000 \
        --cds-merge \
        --adjust-direction
}

export -f process_id

# Loop through the range of values run jobs in parallel
for id in {0..58}; do
    process_id "$id" &
done

# Wait for all background processes to finish
wait

mkdir logs
mv *.log logs
wait

echo "All samples have been processed!"
