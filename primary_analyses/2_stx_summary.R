library(tidyverse)
library(Hmisc)

stx_fams <- read_csv(paste0("STXBP1_full_base_v", vs, ".csv"))

stx_base <- read_csv(paste0("full_base_v", vs, ".csv"))
stx_prop <- read_csv(paste0("full_prop_v", vs, ".csv"))
ic <- read_csv(paste0("pos_IC_v", vs, ".csv"))

# Create output
fileConn <- file(paste0("output/STXBP1_v", vs, "_output.txt"), open = "wt")
writeLines(c(paste0("version = ", vs)), fileConn)

# Individuals
writeLines(c(paste0("\nIndividuals = ", length(unique(stx_fams$famID)))), fileConn)

writeLines(c(paste0("\nVariants = ", length(unique(stx_fams$variant)))), fileConn)


# EGRP individuals
writeLines(c(paste0("EGRP = ", length(unique(stx_fams %>% filter(PMID == "EGRP") %>% pull(famID))))), fileConn)


# Collaborators and literature
collab <- stx_fams %>% 
  select(famID, Year, Journal, PMID, Notes) %>% unique() %>% 
  filter(Journal %in% c("GB_pats", "Referring physician: Emma Palmer", "Referring physician: Lewis-Smith, Rhys Thomas",
                        "Referring physician: Karl Martin Klein, Felix Rosenow, Philip Reif (from Epi25)",
                        "Referring physician: Sarah/Hannah", "Referring physician: C. Schropp", "Colorado_patients",
                        "Rikke_Elena", "Boston", "EGRP", "Seoul_pat", "Referring physician: Christan Hensbach, Yvonne Weber",
                        "Referring physician: Mochel Fanny, Camille Giron", "Referring physician: Hana Krijtova - froms NLES study H Mefford)",
                        "Referring physician I. Scheffer, A. Jane", "Amsterdam", "Cambridge_pats", "Spain_Luiza", "Zhang_pats",
                        "G_pats", "Cambridge_pats and OBrien paper"))

lit <- stx_fams %>% 
  select(famID, Year, Journal, PMID, Notes) %>% unique() %>% 
  filter(famID %nin% collab$famID)

writeLines(c(paste0("Previously unpublished in literature (collab only) = ", collab %>% filter(!grepl("published", Notes)) %>% 
         select(-Notes) %>% unique() %>% nrow())), fileConn)

writeLines(c(paste0("Patients from literature only = ", nrow(lit %>% filter(!grepl("Also GB's patients", Notes) & 
                                                                 !grepl("Obrien_KB", Notes) & 
                                                                 !grepl("Boston", Notes)) %>% 
                                                  select(-Notes) %>% unique()))), fileConn)

writeLines(c(paste0("Patients seen by collaborators (inc. literature) = ", nrow(collab %>% select(-Notes) %>% unique()) + 
         nrow(lit %>% filter(grepl("Also GB's patients", Notes) | grepl("Obrien_KB", Notes) | grepl("Boston", Notes))%>% select(-Notes) %>% unique()))), fileConn)

writeLines(c(paste0("Patients from literature reported more than once = ", nrow(lit %>% select(-Notes) %>% unique()) +
         nrow(collab %>% filter(grepl("published", Notes)) %>% select(-Notes) %>% unique()))), fileConn)


# Sex (M or F)
writeLines(c("\nSex:"), fileConn)
x <- stx_fams %>% 
  select(famID, Sex) %>% 
  unique() %>% 
  dplyr::group_by(Sex) %>% 
  dplyr::summarize(n=n())

writeLines(as.character(x), fileConn)

# Combined syndromes
writeLines(c("\nSyndromes:"), fileConn)
x <- stx_fams %>% 
  select(famID, combined_syndromes) %>% 
  unique() %>% 
  dplyr::group_by(combined_syndromes) %>% 
  dplyr::summarize(n=n()) %>% 
  arrange(desc(n))

# writeLines(as.character(x), fileConn)
for (i in x$combined_syndromes) {
  writeLines(paste0(i, " (n=", x[x$combined_syndromes == i,]$n, ")"), fileConn)
}

# Variant type
writeLines(c("\nVariant Type:"), fileConn)
x <- stx_fams %>% 
  select(famID, var_type) %>% 
  unique() %>% 
  dplyr::group_by(var_type) %>% 
  dplyr::summarize(n=n()) %>% 
  arrange(desc(n)) %>% 
  filter(!is.na(var_type))

# writeLines(as.character(x), fileConn)
for (i in x$var_type) {
  writeLines(paste0(i, " (n=", x[x$var_type == i,]$n, ")"), fileConn)
}


# Denovo missense variants
stx_fams %>% filter(var_type_2 == "missense") %>% count(Inheritance)



writeLines(c("\nPTV vs Missense:"), fileConn)
x <- stx_fams %>% 
  select(famID, var_type_2) %>% 
  unique() %>% 
  dplyr::group_by(var_type_2) %>% 
  dplyr::summarize(n=n()) %>% 
  arrange(desc(n)) %>% 
  filter(!is.na(var_type_2))

writeLines(as.character(x), fileConn)

# Recurrent hotspots
writeLines(c("\nRecurrent Hotspots:"), fileConn)
x <- stx_fams %>% 
  select(famID, variant_combined) %>% 
  unique() %>% 
  dplyr::group_by(variant_combined) %>% 
  dplyr::summarize(n=n()) %>% 
  arrange(desc(n)) %>% 
  head(10) %>% 
  filter(!is.na(variant_combined))

# writeLines(as.character(x), fileConn)
for (i in x$variant_combined) {
  writeLines(paste0(i, "\t(n=", x[x$variant_combined == i,]$n, ")"), fileConn)
}

# Recurrent variants
writeLines(c("\nRecurrent Variants:"), fileConn)
x <- stx_fams %>% 
  select(famID, variant) %>% 
  unique() %>% 
  dplyr::group_by(variant) %>% 
  dplyr::summarize(n=n()) %>% 
  arrange(desc(n)) %>% 
  filter(n>1) %>% 
  filter(!is.na(variant))

writeLines(paste0(NROW(x), " total recurrent variants (n>1); ", 
                  NROW(x %>% filter(n>=3)), " recurrent variants (n>=3); ",
                  NROW(x %>% filter(n>=5)), " recurrent variants (n>=5)"), 
           fileConn)

# writeLines(as.character(x), fileConn)
for (i in x$variant) {
  writeLines(paste0(i, "\t(n=", x[x$variant == i,]$n, ")"), fileConn)
}

# In total, recurrent variants identified in five or more individuals accounted 
# for 35% of all individuals with STXBP1-related disorders and the p.Arg406Cys/His 
# alone accounted for 8.0% of all individuals with STXBP1-related disorders. 

n_5 <- length(unique(stx_fams %>% filter(variant %in% x[x$n >=5,]$variant) %>% pull(famID)))
writeLines(paste0("Individuals with recurrent variants in five or more individuals = ", n_5, " (", round(n_5/length(unique(stx_fams$famID)), 2), " %)"), fileConn)

n_406 <- length(unique(stx_fams %>% filter(grepl("406", variant)) %>% pull(famID)))
writeLines(paste0("Individuals with recurrent variants in five or more individuals = ", n_406, " (", round(n_406/length(unique(stx_fams$famID)), 2), " %)"), fileConn)



# HPO terms
writeLines(c("\nHPO Base:"), fileConn)

writeLines(paste0("Terms = ", length(stx_base$HPO), " (", 
                  length(unique(stx_base$HPO)), " unique terms)"), fileConn)

stx_base %>% 
  dplyr::group_by(famID) %>% 
  dplyr::summarize(n=n()) %>% 
  summary(n) -> x
writeLines(as.character(x), fileConn)

term_dist <- stx_base %>% dplyr::group_by(famID) %>% dplyr::summarize(n=n()) %>% arrange(desc(n))
nrow(term_dist %>% filter(n>=20))
nrow(term_dist %>% filter(n<5))

writeLines(c("\nHPO Prop:"), fileConn)

writeLines(paste0("Terms = ", length(stx_prop$HPO), " (", 
                  length(unique(stx_prop$HPO)), " unique terms)"), fileConn)

stx_prop %>% 
  dplyr::group_by(famID) %>% 
  dplyr::summarize(n=n()) %>% 
  summary(n) -> x
writeLines(as.character(x), fileConn)



# Seizure types
(seizures <- ic %>% filter(grepl("seizure", def) | grepl("spasm", def)) %>% filter(freq_prop > 0) %>% arrange(desc(freq_prop)))
write_csv(seizures, paste0("output/STXBP1_Seizures_freq_v", vs, ".csv"))




# Close output
close(fileConn)

