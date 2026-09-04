# Use environment_CNV.yml

#single file
freec -conf /stor/work/Brock/kennedy/CNV_analysis/CMDuoPopCNV/new_freec_outputs/C6C7_1806/freec_config_C6C7_1806


#all config files
find /stor/work/Brock/kennedy/CNV_analysis/CMDuoPopCNV/testing_FREEC_ploidyassignment/FREEC_configs -type f -name "freec_config_*" -exec freec -conf {} \;

