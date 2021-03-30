## Run all scripts

start <- Sys.time()
message(" \n Running primary analyses... \n ")

setwd(input.yaml$prime_dir)

scripts_path = input.yaml$prime_dir

## SUMMARIZE DATASET

# Create base and prop files
source(paste0(scripts_path, "/1_stx_base_prop_ic.R"))
rm(list=setdiff(ls(), c("vs", "scripts_path")))

# Summary output
source(paste0(scripts_path, "2_stx_summary.R"))
rm(list=setdiff(ls(), c("vs", "scripts_path")))

# Prop change in IC frequency
# source(paste0(scripts_path, "3_stx_hpo_prop_ic_change.R"))
# rm(list=setdiff(ls(), c("vs", "scripts_path")))

## MAIN ANALYSES

setwd(paste0("/Volumes/helbig_lab/projects/STXBP1/v", vs))

# Syndromes freq and HPO associations
source(paste0(scripts_path, "4_stx_combined_syndromes_hpo_assoc.R"))
rm(list=setdiff(ls(), c("vs", "scripts_path")))

# PTV vs missense HPO associations
source(paste0(scripts_path, "5_stx_ptv_missense_hpo_assoc.R"))
rm(list=setdiff(ls(), c("vs", "scripts_path")))

# PTV vs missense (common recurrent variants filtered out) HPO associations 
source(paste0(scripts_path, "6_stx_ptv_missense_filtered_hpo_assoc.R"))
rm(list=setdiff(ls(), c("vs", "scripts_path")))

# Recurrent variants HPO associations
source(paste0(scripts_path, "7_stx_recurrent_var_hpo_assoc.R"))
rm(list=setdiff(ls(), c("vs", "scripts_path")))

# Recurrent variants 406, 292, and 551 "individually collapsed" HPO associations
source(paste0(scripts_path, "8_stx_recurrent_combined_var_hpo_assoc.R"))
rm(list=setdiff(ls(), c("vs", "scripts_path")))

# Variants grouped into 406 recurrent variants, PTV, and other missense HPO associations
# source(paste0(scripts_path, "9_stx_variant_broad_hpo_assoc.R"))
# rm(list=setdiff(ls(), c("vs", "scripts_path")))



## PHENOGRAMS
source(paste0(scripts_path, "10_stx_phenograms.R"))
rm(list=setdiff(ls(), c("vs", "scripts_path")))

## Seizure onset cumulative
source(paste0(scripts_path, "11_stx_seizure_onset.R"))
rm(list=setdiff(ls(), c("vs", "scripts_path")))



## AED comparative effectiveness analysis
source("/Volumes/helbig_lab/projects/STXBP1/v_scripts/STXBP1_seizure_freq_aed_binning.R")
source("/Volumes/helbig_lab/projects/STXBP1/v_scripts/STXBP1_seizure_freq_aed_binning_seiz_free.R")

source("/Volumes/helbig_lab/projects/STXBP1/v_scripts/STXBP1_seizure_freq_aed_comp_effect.R")


## FIGURES
source("/Volumes/helbig_lab/projects/STXBP1/v_scripts/STXBP1_FIGURES_RUN.R")
rm(list=setdiff(ls(), c("vs")))

## Seizure severity natural histories- runs in figures scripts




