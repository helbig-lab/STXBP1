library(tidyverse)
library(memoise)

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

if(is.null(input.yaml$cluster) == T){
  message('\n Please determine processing method in input config file - Cant Proceed without that \n')
  break;
}if(input.yaml$cluster) == 'c'){
  cluster_loc <- readline(prompt="Enter cluster location: ")
  if(cluster_loc == "" | cluster_loc == " " )
    { 
    message('\n  Cluster name invalid, please try a different location \n')
    break;
    }
  else{
    input.yaml$sim_cluster <- as.character(cluster_loc)
    message('\n  Cluster specified \n')
  } 
}if(input.yaml$cluster) == 'm'){
  input.yaml$sim_cluster <- as.character(paste0(file_path, "/sim_analyses/"))
  message('\n  Manual processing specified - Note that this may take longer \n')
}else{
input.yaml$sim_cluster <- as.character(paste0(file_path, "/sim_analyses/"))
message('\n  Default to manual processing - Note that this may take longer \n')
}

if(is.null(input.yaml$sim_algorithm) == T){
  message('\n Please determine algorithm in input config file - Cant Proceed without that \n')
  break;
}if(input.yaml$sim_algorithm) == 1){
  input.yaml$sim_dir <- as.character("/sim_analyses/resnik/")
  message('\n  Similarity analysis algorithm selected - Resnik \n')
}if(input.yaml$sim_algorithm) == 2){
 input.yaml$sim_dir <- as.character("/sim_analyses/cube/")
 message('\n  Similarity analysis algorithm selected - Cube \n')
}else{
input.yaml$sim_dir <- as.character("/sim_analyses/resnik/")
message('\n  Default to Resnik Similarity analysis algorithm \n')
}

#Run sim_analysis files
if(input.yaml$sim_dir == as.character("/sim_analyses/cube/")){
  source(paste0(input.yaml$sim_dir, cube_sim_auto_chunks.R))
  source(paste0(input.yaml$sim_dir, gene_count_cube_auto.R))
}
if(input.yaml$sim_dir == as.character("/sim_analyses/resnik/")){
  source(paste0(input.yaml$sim_dir, res_mod_auto_chunks.R))
  source(paste0(input.yaml$sim_dir, gene_count_resnik_mod_auto.R))
}
else{
  message('\n Please determine directory in input config file - Cant Proceed without that \n')
  break;
}
