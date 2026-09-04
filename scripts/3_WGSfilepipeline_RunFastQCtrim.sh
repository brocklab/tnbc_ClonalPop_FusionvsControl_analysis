# Use environment_CNV.yml

# Run FastQC only on files that start with 'trimmed_' and end with '.fq.gz' in all subdirectories
find . -type f -name "trimmed_*.fq.gz" | parallel -j 4 "fastqc {} -o trim_qc_reports"

# Generate a combined QC report from those results
multiqc trim_qc_reports -o trim_qc_reports