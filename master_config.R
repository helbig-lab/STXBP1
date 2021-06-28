library(tidyverse)
library(memoise)

start <- Sys.time()
message(" \n Starting config file... \n ")


message("\n  ...Checking for analyses to run... \n ")

#General yaml files needed for most analyses
capture <- commandArgs(trailingOnly = TRUE)

opt1 = list(make_option(c("--input"), type = "character", default = "input.yml", dest = "input"))


user_input <- function(name, argv) {
  return(name %in% names(argv))
}

argv <- parse_args(OptionParser(option_list = opt1))

if (user_input("input", argv)) {
  input = argv$input 
    if(file.exists(input) == T){
      input.yaml <- yaml::read_yaml(input)
    }else{
      message('\n Input YAML not found \n')
      break;
    }
  } else {
  message('Cannot proceed without input yaml file. Please use "--input" flag .\n')
}

#Main analyses directories for input/output
if(is.null(input.yaml$output_dir) == T){
  message('\n Please mention the Field output_dir in input config file - Cant Proceed without that \n')
  break;
}

if(is.null(input.yaml$file_path) == T){
  message('\n Please mention the Field file_path in input config file - Cant Proceed without that \n')
  break;
}

if(is.null(input.yaml$prime_dir) == T){
  message('\n Please mention the Field prime_dir in input config file - Cant Proceed without that \n')
  break;
}

if(is.null(input.yaml$prime_dir) == T){
  message('\n Please mention the Field prime_dir in input config file - Cant Proceed without that \n')
  break;
}

#Term propagation
if(is.null(input.yaml$pos_ic) == F ){
  pos_ic <- read_csv(paste0(input.yaml$file_path, "pos_IC.csv"))
}
else{
    message("\n  Manually propagating information content... \n ")
    source(paste0(input.yaml$secondary_dir,"compose_base_prop_ic.R"))
    input.yaml$pos_ic <- read_csv(paste0(input.yaml$file_path,"pos_IC_manual.csv"))
  if(is.null(input.yaml$pos_ic) == T){
    message("\n  Must use default propagation file or run compose_base_prop_ic.R to continue analyses... \n ")
    break;
  }
}


#Similarity analyses
if(is.null(input.yaml$sim_dir) == F ){
  message("\n  Running similarity analysis... \n ")
  source(sim_config.R)
}
else{
    message("\n  Checking for similarity analysis directory source... \n ")
    source(sim_config.R)
    if(is.null(input.yaml$sim_dir) == F){
      message("\n  Running similarity analyses... \n ")
      source(sim_config.R)
    }
    else {
      next;
    }
}

#Comparative effectiveness analyses
if(is.null(input.yaml$comp_dir) == F ){
  message("\n  Running comparative effectiveness analysis... \n ")
  source(comp_config.R)
}
else{
    message("\n  Checking for comparative effectiveness analysis directory source... \n ")
    source(comp_config.R)
    if(is.null(input.yaml$comp_dir) == F){
      message("\n  Running comparative effectiveness analysis... \n ")
      source(comp_config.R)
    }
    else {
      next;
    }
}


message("\n  ...all selected analyses complete \n ")
stop = Sys.time()
stop - start
