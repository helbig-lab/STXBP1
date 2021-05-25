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

#Manual term propagation and other analyses
if(is.null(input.yaml$secondary_dir) == T){
  message('\n Please mention the Field secondary_dir in input config file to complete manual term propagation...skipping this analysis \n')
} 
else {
      next;
    }


message("\n  ...all selected analyses complete \n ")
stop = Sys.time()
stop - start
