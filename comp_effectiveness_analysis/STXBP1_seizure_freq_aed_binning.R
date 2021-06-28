library(tidyverse)
library(Hmisc)

stx_aed_raw <- read_csv(paste0(input.yaml$file_path,"/STXBP1_AED_clean_FINAL_v8.csv")) %>% 
  filter(!(RECORD_ID == "EG0086" & MED_FINAL %in% c("PREDNISOLONE", "METHYLPREDNISOLONE"))) # %>% 
  # mutate(RECORD_ID = paste0(RECORD_ID, "P"))

stx_aed_keto_ukiss <- read_csv(paste0(input.yaml$file_path,"/STXBP1_UKISS_KETO.csv")) %>% 
  mutate(AGE_ENCOUNTER = Age_m/12) %>% 
  rename(MED_FINAL = AED) %>% 
  mutate(AED = "AED") %>% 
  mutate(TITLE = MED_FINAL) %>% 
  select(colnames(stx_aed_raw))

stx_aed_raw <- stx_aed_raw %>% rbind(stx_aed_keto_ukiss)

egrp_seizures_raw <- read_csv(paste0(input.yaml$file_path,"/STXBP1_EGRP_seizure_nat_hist_base_v8.csv")) %>%
  filter(Age_y == 'freq') %>%
  filter(HPO_term == "HP:0001250")

# Infantile spasms
# egrp_seizures_raw <- read_csv(paste0(input.yaml$file_path,"/STXBP1_EGRP_seizure_nat_hist_base_v8.csv")) %>%
#   filter(Age_y == 'freq') %>%
#   filter(HPO_term == "HP:0012469")


# # Focal impaired awareness
# egrp_seizures_raw <- read_csv(paste0(input.yaml$file_path,"/STXBP1_EGRP_seizure_nat_hist_base_v8.csv")) %>%
#   filter(Age_y == 'freq') %>%
#   filter(HPO_term == "HP:0002384")


# # Neonatal seizures
# egrp_seizures_raw <- read_csv(paste0(input.yaml$file_path,"/STXBP1_EGRP_seizure_nat_hist_base_v8.csv")) %>%
#   filter(Age_y == 'freq') %>%
#   filter(HPO_term == "HP:0032807")



#####
#####

# Seizures with unknown frequency?
# for (i in 1:length(unique(egrp_seizures_raw$RECORD_ID))) {
#   
#   p_tab <- egrp_seizures_raw %>% filter(RECORD_ID == unique(egrp_seizures_raw$RECORD_ID)[i]) 
#   t_na <- which(p_tab[p_tab$Age_y == "hpo",] == "yes" & p_tab[p_tab$Age_y == "freq",] == "NAX")
#   egrp_seizures_raw[egrp_seizures_raw$RECORD_ID == unique(egrp_seizures_raw$RECORD_ID)[i] & p_tab$Age_y == "freq", t_na] <- "1"
#   
# }

#####
#####

seizure_freq <- egrp_seizures_raw %>% select(c(EGRP_ID, which(grepl("t_", colnames(egrp_seizures_raw)))))

pats <- seizure_freq$EGRP_ID %>% unique()
aeds <- stx_aed_raw$MED_FINAL %>% unique()

#####
#####

# Binning - monthly

interval_x = round(1/12, 3) # in years
max_age = 20

timex = seq(0, max_age, by = interval_x)

# Create AR matrix for seizure frequencies
# ar_freq <- array(dim = c(length(pats), length(timex)))

# Transform into selected bin length 
ar_freq_bin <- data.frame(array(dim = c(length(pats), length(timex))))
rownames(ar_freq_bin) = pats
colnames(ar_freq_bin) = paste0("t_", timex)

for (i in 1:length(timex)) {
  t = timex[i]
  f = which(as.numeric(gsub("t_", "", colnames(seizure_freq))) < round((t+interval_x), 3) & as.numeric(gsub("t_", "", colnames(seizure_freq))) >= t)
  foo <- seizure_freq[f]
  
  for (p in 1:nrow(foo)) {
    ar_freq_bin[p,i] = median(as.numeric(foo[p,]), na.rm = T)
  }
}

  

# Change in seizure severity between intervals - includes individuals with no seizures
ar_freq_change <- data.frame(array(dim = c(length(pats), length(timex))))
rownames(ar_freq_change) = pats
colnames(ar_freq_change) = paste0("t_", timex)

for (i in 1:(NCOL(ar_freq_bin)-1)) {
  
  ar_freq_change[,i+1] <- -ar_freq_bin[,i+1] + ar_freq_bin[,i]
  
}


# Create AR matrix for aeds
ar_aed <- array(dim = c(length(pats), 
                        length(timex),
                        length(aeds)))

dimnames(ar_aed)[[1]] = pats
dimnames(ar_aed)[[2]] = paste0("t_", timex)
dimnames(ar_aed)[[3]] = aeds

for (a in 1:length(aeds)) {
  
  aed_x <- aeds[a]
  
  pat_timex <- matrix(ncol=length(timex),nrow=length(pats)) %>% as.data.frame()
  rownames(pat_timex) <- pats
  names(pat_timex) <- paste0("t_", timex)
  

  for(p in 1:NROW(pats)) {
    
    patID <- rownames(pat_timex)[p]
    
    asm <- stx_aed_raw %>% filter(RECORD_ID == patID, MED_FINAL == aed_x, AGE_ENCOUNTER <= max_age)
    asm <- asm %>% 
      dplyr::group_by(RECORD_ID, MED_FINAL) %>% 
      dplyr::summarize(min = min(AGE_ENCOUNTER), max = max(AGE_ENCOUNTER))
    
    # Medication usage, assign 1
    if (NROW(asm) > 0) {
      aed_lower = which(timex < asm$min)
      if (length(aed_lower) == 0) {
        aed_lower = 0
      }
      
      aed_upper = which(timex > asm$max)
      if (length(aed_upper) == 0) {
        aed_upper = length(timex)
      }
      aed_ix <- c(max(aed_lower):(min(aed_upper)-1))
      pat_timex[patID,aed_ix] <- 1
    }
    
  }
  
  ar_aed[,,a] <- pat_timex %>% as.matrix()

}


# Seizure frequency changes across AEDs and time intervals

ar_merged <- data.frame(patID = NA, AED = NA, time_int = NA, freq_change = NA, seiz_freq = NA)

for (a in 1:length(aeds)) {
  
  ar_sub <- ar_aed[,,a]
  
  for (p in 1:length(pats)) {

    # Change in seizure frequency
    freq_yes <- ar_freq_change[pats[p], c(which(ar_sub[pats[p],] == 1))] %>% t() %>% 
      as.data.frame()
    
    # Seizure frequency
    freq_yes_1 <- ar_freq_bin[pats[p], c(which(ar_sub[pats[p],] == 1))] %>% t() %>% 
      as.data.frame()
    
    if (NROW(freq_yes) > 0) {
      
      rownames(freq_yes) <- paste0("t_", timex[c(which(ar_sub[pats[p],] == 1))])
      freq_yes$time_int <- rownames(freq_yes)
      
      rownames(freq_yes_1) <- paste0("t_", timex[c(which(ar_sub[pats[p],] == 1))])
      freq_yes_1$time_int <- rownames(freq_yes_1)
      
      for (i in 1:NROW(freq_yes)) {
        foo <- c(pats[p], aeds[a], freq_yes[i,2], freq_yes[i,1], freq_yes_1[i,1])
        ar_merged <- ar_merged %>% rbind(foo)
      }
    }
  }
  
}

ar_merged <- ar_merged %>% 
  filter(!is.na(patID)) %>% 
  unique() %>% 
  mutate(freq_change = as.numeric(freq_change),
         seiz_freq = as.numeric(seiz_freq),
         time_int = as.numeric(gsub("t_", "", time_int)))

ar_merged$patID %>% unique() %>% length()


# FILTER OUT FOR SELECT MEDS
# vigabatrin, topiramate, phenobarbital, levetiracetam, clobazam, UKISS/prednisolone/methylprednisolone, 
# cannabidiol, adrenocorticotropic hormone (ACTH), and oxcarbazepine 

ar_merged_clean <- ar_merged %>% filter(AED %in% c("VIGABATRIN", "TOPIRAMATE", "PHENOBARBITAL", "LEVETIRACETAM",
                                                   "CLOBAZAM", "PREDNISOLONE", "UKISS", "METHYLPREDNISOLONE",
                                                   "CANNABIDIOL", "ACTH", "OXCARBAZEPINE", "KETO", "UKISS PRED")) %>% 
  mutate(AED = case_when(AED %in% c("PREDNISOLONE", "UKISS", "METHYLPREDNISOLONE") ~ "UKISS PRED", TRUE ~ AED)) %>% 
  filter(time_int <= 20) %>% 
  unique()


# write_csv(ar_merged_clean, paste0(input.yaml$file_path,"/STX_EGRP_seizure_aed_1_month_ar_merged_FINAL_v8.csv"))

# write_csv(ar_merged_clean, paste0(input.yaml$file_path,"/STX_EGRP_seizure_aed_1_month_ar_merged_Infantile_Spasms_FINAL_v8.csv"))
# write_csv(ar_merged_clean, paste0(input.yaml$file_path,"/STX_EGRP_seizure_aed_1_month_ar_merged_Focal_Impaired_Awareness_FINAL_v8.csv"))
# write_csv(ar_merged_clean, paste0(input.yaml$file_path,"/STX_EGRP_seizure_aed_1_month_ar_merged_Neonatal_Seizures_FINAL_v8.csv"))


fileConn <- file(paste0(input.yaml$output_dir, "AED_summary.txt"), open = "wt")

# ar_merged_clean <- ar_merged %>% filter(AED %in% c("LEVETIRACETAM", "PHENOBARBITAL", "ACTH", "TOPIRAMATE", "VIGABATRIN", "OXCARBAZEPINE",
#                                                    "CLOBAZAM", "PREDNISOLONE", "CANNABIDIOL", "KETO", "UKISS")) %>% 
#   mutate(AED = case_when(AED %in% c("PREDNISOLONE", "UKISS") ~ "UKISS PRED", TRUE ~ AED)) %>% 
#   filter(time_int <= 20) %>% 
#   unique()

ar_merged_clean$patID %>% unique()

ar_merged_clean %>% unique %>% filter(AED != "KETO") %>% nrow()

ar_merged_clean %>% select(patID, AED) %>% unique() %>% dplyr::group_by(AED) %>% dplyr::summarize(n=n()) %>% arrange(desc(n))

writeLines(paste0("\nTotal month intervals = ", nrow(ar_merged_clean)), fileConn)

writeLines(paste0("\nTotal individuals = ", length(unique(ar_merged_clean$patID))), fileConn)

writeLines(paste0("\nClobazam (months)= ", nrow(ar_merged_clean[ar_merged_clean$AED == "CLOBAZAM",])), fileConn)
writeLines(paste0("\nClobazam (individuals)= ", length(unique(ar_merged_clean[ar_merged_clean$AED == "CLOBAZAM",]$patID))), fileConn)

writeLines(paste0("\nMidazolam (months)= ", nrow(ar_merged_clean[ar_merged_clean$AED == "MIDAZOLAM",])), fileConn)
writeLines(paste0("\nMidazolam (individuals)= ", length(unique(ar_merged_clean[ar_merged_clean$AED == "MIDAZOLAM",]$patID))), fileConn)

writeLines(paste0("\nTopiramate (months)= ", nrow(ar_merged_clean[ar_merged_clean$AED == "TOPIRAMATE",])), fileConn)
writeLines(paste0("\nTopiramate (individuals)= ", length(unique(ar_merged_clean[ar_merged_clean$AED == "TOPIRAMATE",]$patID))), fileConn)

writeLines(paste0("\nKETO (months)= ", nrow(ar_merged_clean[ar_merged_clean$AED == "KETO",])), fileConn)
writeLines(paste0("\nKETO (individuals)= ", length(unique(ar_merged_clean[ar_merged_clean$AED == "KETO",]$patID))), fileConn)

close(fileConn)
