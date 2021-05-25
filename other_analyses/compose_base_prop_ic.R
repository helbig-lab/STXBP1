#Compose Base and Prop

library(readr)
library(tidyverse)
library(Hmisc)

start <- Sys.time()
message(" \n Begin propagation... \n ")

hpo_ancs <- read_csv(input.yaml$hpo_ancestor) %>% 
  rename(HPO = term)

hpo_child <- read_csv(input.yaml$hpo_path) %>% 
  rename(HPO = term) %>% 
  mutate(children = paste0(HPO,";", children)) %>% 
  select(-X1)

hpo_neg_child <- hpo_child %>% 
  mutate(HPO = gsub("H","N",HPO)) %>% 
  mutate(children = gsub("H","N", children))
  

tot_base <- read_csv(input.yaml$variant_file) %>% 
  # rename(famID = ID, HPO = term_id) %>% 
  # mutate(famID = paste0("fam",famID)) %>% 
  filter(var_type_2 != "exclude") %>% 
  filter(grepl("HP:[[:digit:]]", HPO) | grepl("NP:[[:digit:]]", HPO) ) %>% 
  select(famID, HPO) %>% 
  unique()



write_csv(tot_base, paste0(input.yaml$file_path,"full_base.csv"))

########
# Base Files
########

neg_base <- tot_base %>% 
  filter(grepl("NP:[[:digit:]]", HPO)) %>% 
  select(famID, HPO) %>% 
  unique

pos_base <- tot_base %>% 
  filter(grepl("HP:[[:digit:]]", HPO)) %>% 
  unique

write_csv(neg_base, paste0(input.yaml$file_path,"neg_base.csv"))
write_csv(pos_base, paste0(input.yaml$file_path,"pos_base.csv"))

########
# Prop Files
########

# Remove modifier terms
redundant_terms <- c("HP:0003674","HP:0031797", "HP:0012823")

# "HP:0012823" ('clinical modifier'); "HP:0003674" ('Onset'); "HP:0031797" ('Clinical course')

pos_prop <- pos_base %>% 
  left_join(hpo_ancs) %>% 
  select(famID, ancs) %>% 
  separate_rows(ancs, sep=";") %>% 
  rename(HPO = ancs) %>% 
  filter(!is.na(HPO)) %>% 
  filter(HPO %nin% redundant_terms) %>%
  unique

neg_prop <- neg_base %>% 
  left_join(hpo_neg_child) %>% 
  select(famID, children) %>% 
  separate_rows(children, sep=";") %>% 
  rename(HPO = children) %>% 
  filter(!is.na(HPO), HPO !="NA") %>% 
  unique

write_csv(neg_prop, paste0(input.yaml$file_path,"neg_prop.csv"))
write_csv(pos_prop, paste0(input.yaml$file_path,"pos_prop.csv"))


full_prop <- neg_prop %>% rbind(pos_prop)

write_csv(full_prop, paste0(input.yaml$file_path,"full_prop.csv"))



#Negative Prop Pruned

neg_filter <- read_csv(paste0(input.yaml$file_path,"pruned_neg_filter_terms.csv"))

neg_prop_filt <- neg_prop %>% 
  filter(HPO %in% neg_filter$term)

write_csv(neg_prop_filt, paste0(input.yaml$file_path,"neg_prop_pruned.csv"))

########################
# Create IC
########################

########################################################################
########################################################################
########################################################################
########################################################################


neg_base_ct <- neg_base %>% 
  count(HPO) %>% 
  mutate(freq = n/length(unique(tot_base$famID)))

# pos_base = read_csv(paste0(input.yaml$file_path,"stxbp1_full_base_v.csv")) %>% 
#   rename(famID = ID, HPO = term_id) %>% 
#   # filter(!grepl("NP", HPO), grepl("HP:[[:digit:]]", HPO)) %>% 
#   filter(grepl("HP:[[:digit:]]", HPO)) %>% 
#   mutate(famID = paste0("fam_",famID)) %>% 
#   select(famID, HPO) %>% 
#   unique

pos_base_ct <- pos_base %>% 
  count(HPO) %>% 
  mutate(freq_base = n/length(unique(tot_base$famID))) %>% 
  mutate(base.IC = -log2(freq_base)) %>% 
  select(-n)

# hpo_ancs <- read_csv(input.yaml$hpo_ancestor) %>% 
#   rename(HPO = term) 
# 
# hpo_child <- read_csv("input.yaml$hpo_path) %>% 
#   rename(HPO = term) %>% 
#   mutate(children = paste0(HPO,";",children)) %>% 
#   mutate(HPO = gsub("H","N", HPO)) %>% 
#   mutate(children = gsub("H","N", children)) %>% 
#   select(-X1)

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

write_csv(pos_IC,paste0(input.yaml$file_path,"pos_IC.csv"))

message("\n  Propagation complete \n ")
stop = Sys.time()
message("Time for propagation:",stop - start)
