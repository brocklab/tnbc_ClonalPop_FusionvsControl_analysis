#!/bin/bash

# Use environment_CNV.yml

# Base directories to search for samples
BASE_DIRS=(
    "/stor/work/Brock/kennedy/CNV_analysis/CMDuoPopCNV/01.RawData"
)

# Target output base directory
FREEC_BASE="/stor/work/Brock/kennedy/CNV_analysis/CMDuoPopCNV/new_freec_outputs"

# Static paths in config
CHR_LEN="/stor/work/Brock/kennedy/CNV_analysis/Pipeline_creation/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
CHR_FILES="/stor/work/Brock/kennedy/CNV_analysis/Pipeline_creation/ref/chr_files"
FASTA_FILE="/stor/work/Brock/kennedy/CNV_analysis/Pipeline_creation/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Search each base directory
for BASE in "${BASE_DIRS[@]}"; do
    find "$BASE" -type d -name "aligned_reads" | while read aligned_dir; do
        sample_dir=$(dirname "$aligned_dir")
        sample_name=$(basename "$sample_dir")
        bam_file="$aligned_dir/${sample_name}_paired.dedup.bam"

        if [[ ! -f "$bam_file" ]]; then
            echo "BAM file not found for $sample_name. Skipping."
            continue
        fi

        output_dir="$FREEC_BASE/$sample_name"
        mkdir -p "$output_dir"

        config_path="$output_dir/freec_config_${sample_name}"

        echo "Generating config for $sample_name"

        cat > "$config_path" <<EOF
[general]
bedtools=bedtools
samtools=samtools
ploidy=3
chrLenFile=$CHR_LEN
chrFiles = $CHR_FILES
window=500000
outputDir=$output_dir
maxThreads=16
breakPointType=4
minimalSubclonePresence = 0.3
noisyData=FALSE
BedGraphOutput=TRUE

[sample]
mateFile=$bam_file
inputFormat=BAM
mateOrientation=FR

[GCcontent]
fastaFile=$FASTA_FILE
EOF

    done
done

echo "Config generation complete."