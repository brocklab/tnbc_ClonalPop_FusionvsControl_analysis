# Use environment_CNV.yml

# Run FastQC on all fastq.gz files in subfolders using 4 parallel jobs
find . -name "*.fq.gz" | parallel -j 5 "fastqc {} -o qc_reports"

# Generate a combined QC report
multiqc qc_reports -o qc_reports