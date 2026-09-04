#!/bin/bash

# Use environment_CNV.yml

# Reference genome path
REF="/stor/work/Brock/kennedy/CNV_analysis/Pipeline_creation/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
THREADS=8
ROOT_DIR="/stor/work/Brock/kennedy/CNV_analysis/CMDuoPopCNV/01.RawData"

for sample_dir in "$ROOT_DIR"/*; do
    [ -d "$sample_dir" ] || continue
    echo "Processing sample: $(basename "$sample_dir")"

    R1_LIST=($(find "$sample_dir" -maxdepth 1 -name "trimmed*_1.fq.gz" | sort))
    R2_LIST=($(find "$sample_dir" -maxdepth 1 -name "trimmed*_2.fq.gz" | sort))

    if [[ ${#R1_LIST[@]} -eq 0 || ${#R2_LIST[@]} -eq 0 || ${#R1_LIST[@]} -ne ${#R2_LIST[@]} ]]; then
        echo "Missing or mismatched read pairs in $sample_dir"
        continue
    fi

    echo "Found ${#R1_LIST[@]} R1/R2 files."

    # Decide whether to concatenate or not
    if [[ ${#R1_LIST[@]} -gt 1 ]]; then
        echo "Concatenating reads from multiple lanes..."
        MERGED_R1="$sample_dir/merged_R1.fq.gz"
        MERGED_R2="$sample_dir/merged_R2.fq.gz"
        cat "${R1_LIST[@]}" > "$MERGED_R1"
        cat "${R2_LIST[@]}" > "$MERGED_R2"
    else
        MERGED_R1="${R1_LIST[0]}"
        MERGED_R2="${R2_LIST[0]}"
    fi

    mkdir -p "$sample_dir/aligned_reads"

    BAM="$sample_dir/aligned_reads/$(basename "$sample_dir")_paired.bam"
    SORTED_BAM="$sample_dir/aligned_reads/$(basename "$sample_dir")_paired.sorted.bam"
    RG_BAM="$sample_dir/aligned_reads/$(basename "$sample_dir")_paired.sorted.rg.bam"
    DEDUP_BAM="$sample_dir/aligned_reads/$(basename "$sample_dir")_paired.dedup.bam"
    METRICS="$sample_dir/aligned_reads/metrics.txt"

    echo "Running BWA MEM..."
    bwa mem -t $THREADS "$REF" "$MERGED_R1" "$MERGED_R2" | samtools view -Sb - > "$BAM"

    echo "Sorting BAM..."
    samtools sort -@ $THREADS -o "$SORTED_BAM" "$BAM"

    echo "Attempting MarkDuplicates without read groups..."
    if ! picard MarkDuplicates I="$SORTED_BAM" O="$DEDUP_BAM" M="$METRICS"; then
        echo "MarkDuplicates failed — adding read groups and retrying..."
        picard AddOrReplaceReadGroups \
            I="$SORTED_BAM" \
            O="$RG_BAM" \
            RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$(basename "$sample_dir")

        if picard MarkDuplicates I="$RG_BAM" O="$DEDUP_BAM" M="$METRICS"; then
            echo "MarkDuplicates succeeded after adding read groups."
        else
            echo "MarkDuplicates failed again. Skipping sample: $(basename "$sample_dir")"
            continue
        fi
    else
        echo "MarkDuplicates succeeded without read groups."
    fi

    echo "Indexing final BAM..."
    samtools index "$DEDUP_BAM"

    echo "Finished processing: $(basename "$sample_dir")"
    echo ""
done

echo "All samples processed!"