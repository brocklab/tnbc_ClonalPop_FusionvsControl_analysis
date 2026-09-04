#!/bin/bash

#Use environment_CNV.yml

# to make file executable type this in command line: chmod +x 2_WGSfilepipeline_Trim.sh
# then type ./2_WGSfilepipeline_Trim.sh to acually run the script

# ------------------------------
# CONFIGURATION
# ------------------------------
ADAPTER_P7="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
THREADS=4  # Adjust based on your system

# ------------------------------
# FUNCTION TO PROCESS A PAIR OF FILES
# ------------------------------
trim_pair() {
    R1="$1"
    R2="$2"
    out_dir=$(dirname "$R1")
    base_R1=$(basename "$R1")
    base_R2=$(basename "$R2")

    echo "Trimming: $base_R1 and $base_R2"

    cutadapt \
        -a "$ADAPTER_P7" -A "$ADAPTER_P7" \
        -a "G{20}" -A "G{20}" \
        -q 20 \
        --minimum-length 30 \
        -j "$THREADS" \
        -o "${out_dir}/trimmed_${base_R1}" \
        -p "${out_dir}/trimmed_${base_R2}" \
        "$R1" "$R2"
}

# ------------------------------
# MAIN LOOP: Find and Trim Paired Reads
# ------------------------------
echo "Searching for paired-end FASTQ files..."

find . -type f \( -name "*_1.f*q.gz" -o -name "*R1*.f*q.gz" \) | while read -r R1; do
    # Replace only the trailing _1.fq.gz or _R1.fq.gz to get the R2 file
    R2=$(echo "$R1" \
        | sed -E 's/(_L[0-9]+)_1\.fq\.gz$/\1_2.fq.gz/' \
        | sed -E 's/(_L[0-9]+)_R1\.fq\.gz$/\1_R2.fq.gz/')

    if [[ -f "$R2" ]]; then
        trim_pair "$R1" "$R2"
    else
        echo "No matching R2 file for: $R1"
    fi
done

echo "Trimming complete."