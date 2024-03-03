#!/bin/bash
# This script was written to conduct two runs of MitoFinder to extract mitogenomes from UCE captured data. Mitofinder is a pipeline to assemble and annotate mitochondrial DNA from trimmed reads. Please check instruction on MF's GitHub to make sure that your options are adequate for your data and study organism - Instructions on https://github.com/RemiAllio/MitoFinder
# Written by lauriane.baraf@my.jcu.edu.au on 12/06/2023 - last modified 30/10/2023
# Associated GitHub repo includes - 1 annotated bash script (MitoFinder.sh) + 1 PBS script (run_MitoFinder.pbs) + 2 bash scripts to run in parallel on your local  machine (mitofinder*_inhouse.sh) + 1 mt genes text file (genes_list.txt)


--------------------------------------------------- REQUIREMENTS & STEPS ---------------------------------------------------
#REQUIREMENTS - varies depending on the method, here we used pair-end reads
# Reference_file.gb containing at least one mitochondrial genome of reference - available on NCBI
# left_reads.fastq.gz containing the left reads of paired-end sequencing
# right_reads.fastq.gz containing the right reads of paired-end sequencing
# SE_reads.fastq.gz containing the reads of single-end sequencing
# assembly.fasta containing the assembly on which MitoFinder have to find and annotate mitochondrial contig.s - can be obtained after the first run of MitoFinder

# STEPS
#1 Download your reference genomes in GenBank format
#2 If multiple ref genomes - make sure all annotations are consistent across references - e.g. D-loop/control region, COX1/COI
#3 First MitoFinder run - assemble and/or annotate mitochondrial-like genes
#4 Prepare input data for second run
#5 Second MitoFinder run - annotate mt genes containing intron(s)
#6 Store results
#7 Plot results
-----------------------------------------------------------------------------------------------------------------------------

################################################################## 2 Consistent Annotations #################################################################
# It's easier to visualise differences in annotations across sequences in genome visualisation softwares like Geneious, if you have a license you can change the annotations in there directly. Or you can do it using bash commands

sed -i 's/control\ region\ D-loop/D-loop/g' /home/MitoFinder/RefSeq_mtDNA_Pomacanthidae.gb
sed -i 's/control\ region/D-loop/g' /home/MitoFinder/RefSeq_mtDNA_Pomacanthidae.gb


############################################# 3 First MitoFinder run - MitoFinder.pbs or mitofinder_inhouse.sh #############################################
# Note: If you are using captured data (e.g. UCE), you'll need to specify --min-contig-size option based on your contig coverage. The default value is 1,000 bp which might be too high for mitochondrial contigs assembled from off-target reads. The same applies for --blast-size so I also reduced the e-value to 1e-06 (default 0.00001, 30%) following Allio et al (2019). This will allow you to get more hits and detect more genes. You can run multiple test runs with different values for these options to check how your output assemblies look like.
# There are two scripts available, one to submit as a pbs job (.pbs) but if you don't have access to an HPC or if the queue in there is longer than for a newly opened Disney ride, there are two other scripts (.sh) for run1 and run2 that will run your samples in parallel on your local machine.

# Get array index
id=${PBS_ARRAY_INDEX}
# Get list of species
cleaned_reads=($(ls /home/MitoFinder/2_cleanedreads_PMC))
# Get individual samples
sample="${cleaned_reads[$id]}"
r1=$(echo "/home/MitoFinder/2_cleanedreads_PMC/$sample/split-adapter-quality-trimmed/$sample-READ1.fastq.gz")
r2=$(echo "/home/2/MitoFinder/2_cleanedreads_PMC/$sample/split-adapter-quality-trimmed/$sample-READ2.fastq.gz")

mitofinder --metaspades \
    -j $sample \
    -1 $r1 \
    -2 $r2 \
    -r RefSeq_mtDNA_Pomacanthidae.gb \
    -o 2 \
    -p 20 \
    -t mitfi \
    -e 0.000001 \
    --min-contig-size 500 \
    --max-memory 100 \
    --new-genes


############################################################# 4 Prepare input data for second run #############################################################
# Make a new directory to store mt contigs from MF run 1
mkdir -p /scratch/MitoFinder/PMC_contigs_MFRUN1
# Concatenate all the mt contigs from the first MitoFinder run for each species - you'll have to tweak these loops a little depending on how your directories are organised
cd /scratch/MitoFinder/MF_RUN1
for sp in "${cleaned_reads[@]}"; do
    mkdir -p ../PMC_contigs_MFRUN1/$sp;
    outdir="../PMC_contigs_MFRUN1/$sp"
    if [ -f "$sp/${sp}_MitoFinder_metaspades_mitfi_Final_Results/${sp}_mtDNA_contig.fasta" ]; then
        cat $sp/${sp}_MitoFinder_metaspades_mitfi_Final_Results/${sp}_mtDNA_contig.fasta >> $outdir/${sp}_MFRUN1_mtDNA_contigs.fasta
    else
        # this bit does the same but for multiple contigs
        for i in {1..9}; do
            cat $sp/${sp}_MitoFinder_metaspades_mitfi_Final_Results/${sp}_mtDNA_contig_${i}.fasta >> $outdir/${sp}_MFRUN1_mtDNA_contigs.fasta
        done 2>/dev/null # this prevents your terminal from printing bunch of "file not found" when there's no contig for i=8 for example
    fi
done

# You can clean it up by only keeping directories with mt contigs and see for how many species MitoFinder has found mt contigs for
cd /scratch/MitoFinder/PMC_contigs_MFRUN1
count=0
for file in */*_MFRUN1_mtDNA_contigs.fasta; do
    sp=($(echo $file | cut -d/ -f1));
    # Check which sp MF didn't find contigs for and delete directories.
    if [ $(wc -l < "$file") -eq 0 ]; then
            rm -r "$sp" # Delete the file - comment out if you want to keep it and just know how many sp have data
            echo "Mitofinder didn't find contigs for $sp. Directory was removed."
    else
        ((count++))
    fi
done
echo "Number of species with mtDNA contigs: $count"


########################################## 5 Second MitoFinder run - SCRIPT MitoFinder.pbs or mitofinder2_inhouse.sh ##########################################
# Only run this one if you have a close reference available
# Note: modify intron size based on your study species - here it's set for fishes. Have a look at the different options using -h
/home/MitoFinder/MitoFinderPackage/mitofinder -j $sample \
        -a $contigs \
        -r PMC_mtRefSeq.gb \
        -o 2 \
        -p 20 \
        -m 40 \
        -t mitfi \
        -e 0.000001 \
        --min-contig-size 500 \
        --new-genes \
        --allow-intron \
        --intron-size 12000 \
        --cds-merge \
        --adjust-direction


#Note: sometimes genes can be found more than once in different contigs - MitoFinder will keep the longest sequence - this suggests either fragmentation, NUMT annotations, or potential contamination of the sequencing data. Different contigs may be part of different organisms thus "sample_final_genes_NT.fasta" and "sample_final_genes_AA.fasta" could be erroneous - It's important to check contigs and associated genes separately
# Because MitoFinder generates a file [Seq_ID]_link_[assembler].scafSeq which contains assembled UCE contigs, they can be annotated and fed directly to Phyluce from the phyluce_assembly_match_contigs_to_probes step (Finding UCE)


###################################################################### 6 Store & Counts ######################################################################
# Make a new directory to store mt contigs
mkdir -p /scratch/MitoFinder/PMC_contigs_MFRUN2
# Concatenate all the mt contigs from the second MitoFinder run for each species
cd /scratch/MitoFinder/MF2_out
for dir in */*_MitoFinder_mitfi_Final_Results; do
    sp=($(echo $dir | cut -d/ -f1))
    if [ -z "$(ls -A "$dir")" ]; then
        echo "MitoFinder didn't find mt contigs for $sp."
    else
        mkdir -p /scratch/MitoFinder/PMC_contigs_MFRUN2/$sp
        outdir="/scratch/MitoFinder/PMC_contigs_MFRUN2/$sp"
        if [ -f $dir/${sp}_mtDNA_contig.fasta ]; then
                cat $sp/${sp}_MitoFinder_mitfi_Final_Results/${sp}_mtDNA_contig.fasta >> $outdir/${sp}_MFRUN2_mtDNA_contigs.fasta
            else
                echo "Fragmented or partial mt genome for $sp. Will still appear in the out directory, comment off if not desired"
                cat $sp/${sp}_MitoFinder_mitfi_Final_Results/${sp}_mtDNA_contig_[1-9].fasta >> $outdir/${sp}_MFRUN2_mtDNA_contigs.fasta
        fi
    fi
done

# Maybe you want to check how long the contigs are in total to see how many potential whole genomes you got
for contigs in PMC_contigs_MFRUN2/*/*.fasta; do
    sp=$(echo $contigs | cut -d/ -f3 | cut -d_ -f1,2)
    length=$(awk '/^>/ {next} {letters += length} END {print letters}' $contigs)
    echo "$sp: $length"
done | sort -t' ' -k2 -rn

# Or get assembly coverage - you'll need that info if you decide to submit your sequences on GenBank for example. Here we assembled the reads with MetaSpades which has the k-mer coverage for the largest k value used for each contig in the metaspades/contigs.fasta output file - note that k-mer coverage is always lower than the read per-base coverage. You can still extract and average those values to get estimated k-mer coverage of your assembly. To get per-base coverage, map your assembly back to your trimmed raw reads using bwa, bbmap or minimap for example 

# K-mer coverage
awk -F'_' '/^>/{sum+=$6; count++} END{if(count>0) print "Sum:", sum, "\nCount:", count, "\nAverage:", sum/count; else print "No values found."}' /scratch/MitoFinder/MF_RUN1/sp1/sp1_metaspades/contigs.fasta

# Per-base coverage
# map assemblies back to trimmed reads (the same that you used to run MitoFinder)
bbmap.sh in=./sp1-READ1.fastq in2=./sp1-READ2.fastq ref=./sp1_mtDNA_contig.fasta covstats=covstats.txt 1>screen.txt 2>&1


######################################################################## 7 Plot away ########################################################################
# Create text file with genes you want to plot, here it's named genes_list.txt (creative I know) - make sure that annotations match those of MitoFinder
# The following bits of code is going to count the number of genes found by Mitofinder after round1 and round2. Run it twice on the output folders by replacing the $file and $csv_file - don't forget to reset count!

## First initialize an associative array to store counts for each string (which are the genes in the gene list)
declare -A gene_counts
# Loop through the lines in genes_list.txt
while read -r line; do
    # Use grep to search for the genes in the original output files from MitoFinder (MF_RUN*) - not the contigs ones (name_contigs_MFRUN*)
    file="/scratch/MitoFinder/test/*/*_MitoFinder_metaspades_mitfi_Final_Results/*_final_genes_NT.fasta"
    count=$(grep -l "$line" $file | wc -l)
    # Store the count in the associative array
    ((gene_counts["$line"] += count))
done < "genes_list.txt"

# Create the CSV file
csv_file="MF_gene_counts_R1.csv"
echo "Genes,Count" > $csv_file 
# Display the counts for each string
for gene in $(<genes_list.txt); do
    echo "$gene, ${gene_counts[$gene]}" >> $csv_file
done

# /!\ reset count in between each iteration /!\
for gene in "${!gene_counts[@]}"; do
    gene_counts["$gene"]=0
done

## If you want gene count across both runs with species names
output_file="total_MF_genes.csv"
# Header for the output file
echo -e "Species\tGenes\tCount" > "$output_file"

# Make a function to process run files
process_directory() {
    local directory="$1"

    for fasta in "$directory"/*/*_MitoFinder*_Final_Results/*_final_genes_NT.fasta; do
        sp=$(basename "${fasta%_final_genes_NT.fasta}")
        while read -r line; do
            present=$(grep -qwF "$line" "$fasta" && echo "1" || echo "0")
            echo -e "$sp\t$line\t$present" >> "$output_file"
        done < "genes_list.txt"
    done
}

# Process files in MF_RUN1 directory
process_directory "/scratch/MitoFinder/MF_RUN1"
# Process files in MF_RUN2 directory
process_directory "/scratch/MitoFinder/MF_RUN2"

# Sort the output file so you keep only unique occurrence of genes across runs per species
sort -u -k1,2 "$output_file" -o "$output_file"
# See how many genes in total
awk '{sum += $3} END {print sum}' "$output_file"

## You can then use the csv files to plot results in R
