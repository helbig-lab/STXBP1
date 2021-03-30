# Compose Base and Prop

start <- Sys.time()
message(" \n Running base propogation... \n ")

library(readr)
# library(readxl)
library(tidyverse)
library(Hmisc)

hpo_ancs <- read_csv(paste0(input.yaml$file_path,"/HPO_ancestors_v0.1.2_dl-2019-02-05_rl_2018-12-21.csv")) %>% 
  rename(HPO = term)

# hpo_child <- read_csv(paste0(input.yaml$file_path,"HPO_child_paths_v0.1.2_dl-2019-02-05_rl_2018-12-21.csv")) %>% 
#   rename(HPO = term) %>% 
#   mutate(children = paste0(HPO,";", children))

# hpo_neg_child <- hpo_child %>% 
#   mutate(HPO = gsub("H","N",HPO)) %>% 
#   mutate(children = gsub("H","N", children))
  
tot_base <- read_csv(paste0(input.yaml$file_path,"STXBP1_full_base_v.csv")) %>% 
  # filter(grepl("HP:[[:digit:]]", HPO) | grepl("NP:[[:digit:]]", HPO) ) %>% 
  filter(grepl("HP:[[:digit:]]", HPO)) %>%
  select(famID, HPO) %>% 
  unique()

write_csv(tot_base, paste0(input.yaml$output_dir,"full_base_v.csv"))

########
# Base Files
########

# neg_base <- tot_base %>% 
#   filter(grepl("NP:[[:digit:]]", HPO)) %>% 
#   select(famID, HPO) %>% 
#   unique

pos_base <- tot_base %>%
  filter(grepl("HP:[[:digit:]]", HPO)) %>%
  unique

# write_csv(neg_base, paste0(input.yaml$output_dir,"neg_base_v.csv"))
# write_csv(pos_base, paste0(input.yaml$output_dir,"pos_base_v.csv"))

########
# Prop Files
########

# Remove modifier terms
# "HP:0012823" ('clinical modifier'); "HP:0003674" ('Onset'); "HP:0031797" ('Clinical course')
redundant_terms <- c("HP:0003674","HP:0031797", "HP:0012823")

pos_prop <- pos_base %>% 
  left_join(hpo_ancs) %>% 
  select(famID, ancs) %>% 
  separate_rows(ancs, sep=";") %>% 
  rename(HPO = ancs) %>% 
  filter(!is.na(HPO)) %>% 
  filter(HPO %nin% redundant_terms) %>%
  unique

# neg_prop <- neg_base %>% 
#   left_join(hpo_neg_child) %>% 
#   select(famID, children) %>% 
#   separate_rows(children, sep=";") %>% 
#   rename(HPO = children) %>% 
#   filter(!is.na(HPO), HPO !="NA") %>% 
#   unique

# write_csv(neg_prop, paste0(input.yaml$output_dir,"neg_prop_v.csv"))
# write_csv(pos_prop, paste0(input.yaml$output_dir,"pos_prop_v.csv"))


# full_prop <- neg_prop %>% rbind(pos_prop)
# write_csv(full_prop, paste0(input.yaml$output_dir,"full_prop_v.csv"))

write_csv(pos_prop, paste0(input.yaml$output_dir,"full_prop_v.csv"))

########################
# Create IC
########################

########################################################################
########################################################################
########################################################################
########################################################################


# neg_base_ct <- neg_base %>% 
#   count(HPO) %>% 
#   mutate(freq = n/length(unique(tot_base$famID)))

pos_base_ct <- pos_base %>% 
  count(HPO) %>% 
  mutate(freq_base = n/length(unique(tot_base$famID))) %>% 
  mutate(base.IC = -log2(freq_base)) %>% 
  select(-n)

pos_prop_ct <- pos_prop %>% 
  count(HPO) %>% 
  mutate(freq_prop = n/length(unique(tot_base$famID))) %>% 
  mutate(prop.IC = -log2(freq_prop)) %>% 
  select(-n)

pos_IC <- pos_base_ct %>% 
  full_join(pos_prop_ct) %>%
  right_join(hpo_ancs %>% select(HPO, def) %>% unique)

pos_IC[is.na(pos_IC)] <- 0

# tot_prop = pos_prop %>% rbind(neg_prop) %>% unique

length(unique(pos_prop$HPO))

write_csv(pos_IC, paste0(input.yaml$output_dir,"pos_IC_v.csv"))

message("\n  ...term propogation complete \n ")
stop = Sys.time()
stop - start
