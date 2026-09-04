# Processing and analyzing TagSeq and low-pass WGS data from clonally-derived homotypic cell-cell fusion triple-negative breast cancer cell populations and corresponding controls

## Introduction

### Summary
To investigate the genomic and transcriptomic impact of homotypic cell-cell fusion, we compared copy number variation (CNV) of eight single-cell isolated clonal cell populations (controls) to eight fusion-derived clonal populations (fusions) via low-pass whole-genome sequencing (WGS) in two triple-negative breast cancer (TNBC) cell lines (HCC1806 and MDA-MB-231). Additionally, we co-cultured pairs of control clonal populations and single-cell isolated spontaneously arising fusion events to generate clonal fusion populations with matched parental controls, referred to as matched fusion-parent trios. All fusion clonal populations from matched trios were also profiled using WGS to analyze CNV versus matched parents, and all surviving fusion and control clonal populations from the matched trios were profiled in biological triplicate by TagSeq (3' tag-based RNA-seq) for transcriptomic analysis.

### Data Availability
The raw data analyzed in this study is available via SRA at [PRJNA1522362](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1522362) and the processed TagSeq data are available via GEO at [GSE346002](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE346002).

## Environments
All environments were handled with conda and are saved as yml files in misc/. Each respective notebook or script used throughout the analysis describes which environment to use in a line near the top.

## Workflow
Generally, analysis was performed using jupyter notebooks available in `notebooks/` and named intuitively. Some notebooks have cells near the top which allow you to run the data for one cell line or the other. Alternatively, these notebooks can be duplicated and one can be used for each cell line. Data necessary to replicate results is available in the `data/` folder under intuitive descriptions. Alternatively, analysis can be replicated from raw data available on SRA and GEO using scripts located in `scripts/` which will be described more thoroughly below.

### Processing raw data files
- To process low-pass WGS data starting with the raw files, first download the fq.gz files for all populations and replicates from SRA. Next, utilize the numbered WGS scripts in `scripts/` one by one, sequentially, to process the raw files. Scripts 1, 3, and 6 contain command line prompts which can simply be pasted into the terminal after setting up `misc/environment_CNV.yml` which all scripts utilize. Scripts 2, 4, and 5 can be executed as shell scripts via the terminal by first making the file executable "chmod u+x filename.sh" then running with "./filename.sh". Additionally, to run the final script `6_WGSfilepipeline_Run_freec_configs.sh` Control-FREEC v11.6b was downloaded from the [official Control-FREEC GitHub releases page](https://github.com/BoevaLab/FREEC) and installed as a standalone Linux command-line tool. After unpacking/compiling the source distribution, the compiled executable located in FREEC-11.6b/src/ was added to $PATH after which the command lines in script 6 can be executed. Note, scripts should be located in the same directory as files during execution.
- To process TagSeq files starting with the fq.gz files available via SRA, deduplicated and trimmed FASTQ files were aligned to the human reference genome (GRCh38) and quantified using the nf-core/rnaseq Nextflow pipeline (v3.26.0) with STAR for alignment and Salmon for transcript-level quantification. This process will produce the `salmon.merged.gene_counts.tsv` (among others) which should be used for transcriptomic analysis and can be downloaded directly from GEO or found in `data/TranscriptomicData`.

### FACS fusion rate analysis
This analysis can be replicated using `notebooks/FACSFusionRateAnalysis.ipynb` and data located in `data/FACSFusionRateAnalysis`.

### ImageBasedFusionRateAnalysis
This analysis first involved running images for each sample from each day of the investigation (representative images found here `data/ImageBasedFusionRateAnalysis/representative_HCC1806_images` or here `data/ImageBasedFusionRateAnalysis/representative_MDA-MB-231_images`) through custom segmentation code (found here `scripts/segment_cells_for_ImageBasedFusionRate.py`) to better visualize fusion cells within each well. Then manual counting was conducted of fusion events over the course of the experimental time frame. Next, analysis was conducted using `notebooks/ImageBasedFusionRateAnalysis.ipynb` and count data located in `data/ImageBasedFusionRateAnalysis`.

### Growth rate analysis
This analysis can be replicated using `notebooks/GrowthRate_Analysis_HCC1806.ipynb` or `notebooks/GrowthRate_Analysis_MDA-MB-231.ipynb` and confluence vs. time data located in `data/GrowthRate`.

### Morphology analysis
Morphology measurements of MDA-MB-231 samples utilized automatic segmentation. This process can be replicated starting with the image segmentation notebook `notebooks/SegmentationForMorphology_MDA-MB-231.ipynb` which should be located in the same folder as `scripts/incucyte_pipeline.py` when running the ipynb. All parameters in the config cell can be left alone except for the `use_selection_manifest_path` parameter which should be set to the respective population's manifest such as `data/Morphology/MDA-MB-231_ImageSegmentationForMorph/Outputs/C1_231/qc/C1_selection_manifest.json`. The control C1 population's results can be entirely replicated from its manifest and the image files located at `data/Morphology/MDA-MB-231_ImageSegmentationForMorph/Input_Representative_ControlSample`. Images for all other samples are available upon request, but are omitted here due to limited data storage. However, all outputs of the analysis for each respective sample group are available in their respective output folders such as `data/Morphology/MDA-MB-231_ImageSegmentationForMorph/Outputs/C1_231`. Following the segmentation for 231s, morphology analysis can be replicated using `notebooks/Morphology_MDA-MB-231.ipynb` and each population's respective morphology outputs produced in the image segmentation analysis above, for example: `data/Morphology/MDA-MB-231_ImageSegmentationForMorph/Outputs/C1_231/cell_morphology_C1_20260408_132607.csv` and `data/Morphology/MDA-MB-231_ImageSegmentationForMorph/Outputs/C1_231/nuclear_morphology_C1_20260408_132607.csv`. HCC1806 cell area measurements were conducted manually in ImageJ and nuclear segmentation data was exported from the automatic IncuCyte segmentation. Analysis can be replicated using `notebooks/Morphology_HCC1806.ipynb` and data located here `data/Morphology/HCC1806_Morphology`.

### Ploidy analysis
This analysis can be replicated using each cell line's respective notebook (`notebooks/Ploidy_HCC1806.ipynb` or `notebooks/Ploidy_MDA-MB-231.ipynb`) and the data located here `data/Ploidy`.

### CNV analysis
This analysis can be replicated using `notebooks/CNV_Analysis_AllPopulations.ipynb`. The re-created Control-FREEC `_ratio.txt` files or the provided ones located at `data/WGS_CNV_Analysis` can be used to run the notebook.

### Matched trio additive CNV analysis
This analysis can be replicated using `notebooks/MatchedTrioAdditiveCNVAnalysis.ipynb`. The re-created Control-FREEC `_ratio.txt` files or the provided ones located at `data/MatchedTrioAdditiveCNVAnalysis` can be used to run the notebook.

### Gene dosage effect analysis
This analysis can be replicated using the `notebooks/GeneDosageAnalysis.ipynb` notebook, the DEseq data for each sample located at `data/GeneDosageAnalysis/deseq_data` (or recreated), and the CNV data located at `data/GeneDosageAnalysis/CNV_data`.

### Transcriptomic analysis
This analysis can be replicated using the `notebooks/TranscriptomicAnalysis.ipynb` notebook and the re-created `salmon.merged.gene_counts.tsv` file, or the downloaded version from GEO or provided here `data/TranscriptomicData`.

### Transcriptomic IPA analysis
This analysis can be replicated using the DEseq .tsv files re-created for each population during transcriptomic analysis above, or from the DEseq data provided at `data/GeneDosageAnalysis/deseq_data`. These files were fed into IPA software which produced an `_summary.txt` file for each sample group which are located at `data/IPA_Analysis`. All of these files are necessary to run the `notebooks/IPA_Analysis.ipynb`.

