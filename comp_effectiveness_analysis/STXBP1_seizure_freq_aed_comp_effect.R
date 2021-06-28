start <- Sys.time()
message(" \n Starting comp effectiveness analysis script... \n ")

library(tidyverse)

#####
# Comparative function
#####


# In years
comp_effect <- function(aedx, start_int, stop_int) {
  
  x <- ar_merged %>% filter(time_int < stop_int & time_int >= start_int)
  
  x_1 <- x %>% filter(AED == aedx)
  
  aed_pos_improv <- x_1 %>% filter(!is.na(freq_change)) %>% filter(freq_change > 0) %>% nrow()
  aed_pos_nimprov <- x_1 %>% filter(!is.na(freq_change)) %>% filter(freq_change <= 0) %>% nrow()
  
  aed_neg_improv <- x %>% filter(!is.na(freq_change)) %>% filter(freq_change > 0) %>% nrow() - aed_pos_improv
  aed_neg_nimprov <- x %>% filter(!is.na(freq_change)) %>% filter(freq_change <= 0) %>% nrow() - aed_pos_nimprov
  
  fish <- matrix(ncol=2,nrow=2) %>% as.data.frame()
  
  row.names(fish) <- c('aed_pos','aed_neg')
  names(fish) <- c('improvement','not_improvement')
  
  fish['aed_pos','improvement'] <- aed_pos_improv
  fish['aed_pos','not_improvement'] <- aed_pos_nimprov
  fish['aed_neg','improvement'] <- aed_neg_improv
  fish['aed_neg','not_improvement'] <- aed_neg_nimprov
  
  test.res <- fisher.test(fish)
  pval <- test.res$p.value
  OR <- test.res$estimate
  OR_lower <- test.res$conf.int[1]
  OR_upper <- test.res$conf.int[2]
  
  return(c(aedx, paste0("t_", round(start_int, 2), "_", round(stop_int, 2)), pval, OR, OR_lower, OR_upper, length(unique(x_1$patID)),
           aed_pos_improv, aed_pos_nimprov, aed_neg_improv, aed_neg_nimprov))
  
}

comp_effect_maintain <- function(aedx, start_int, stop_int) {
  
  x <- ar_merged %>% filter(time_int < stop_int & time_int >= start_int)
  
  x_1 <- x %>% filter(AED == aedx)
  
  aed_pos_improv <- x_1 %>% filter(!is.na(freq_change)) %>% filter(freq_change == 99) %>% nrow()
  aed_pos_nimprov <- x_1 %>% filter(!is.na(freq_change)) %>% filter(freq_change <= 0) %>% nrow()
  
  aed_neg_improv <- x %>% filter(!is.na(freq_change)) %>% filter(freq_change == 99) %>% nrow() - aed_pos_improv
  aed_neg_nimprov <- x %>% filter(!is.na(freq_change)) %>% filter(freq_change <= 0) %>% nrow() - aed_pos_nimprov
  
  fish <- matrix(ncol=2,nrow=2) %>% as.data.frame()
  
  row.names(fish) <- c('aed_pos','aed_neg')
  names(fish) <- c('improvement','not_improvement')
  
  fish['aed_pos','improvement'] <- aed_pos_improv
  fish['aed_pos','not_improvement'] <- aed_pos_nimprov
  fish['aed_neg','improvement'] <- aed_neg_improv
  fish['aed_neg','not_improvement'] <- aed_neg_nimprov
  
  test.res <- fisher.test(fish)
  pval <- test.res$p.value
  OR <- test.res$estimate
  OR_lower <- test.res$conf.int[1]
  OR_upper <- test.res$conf.int[2]
  
  return(c(aedx, paste0("t_", round(start_int, 2), "_", round(stop_int, 2)), pval, OR, OR_lower, OR_upper, length(unique(x_1$patID)),
           aed_pos_improv, aed_pos_nimprov, aed_neg_improv, aed_neg_nimprov))
  
}

#####
#####

for (s in c("seizure_improve", "seizure_improve_freedom", "seizure_freedom")) {
  
  if (s == "seizure_improve") {
    
    ar_merged <- read_csv(paste0(input.yaml$file_path,"/STX_EGRP_seizure_aed_1_month_ar_merged_FINAL_v8.csv"))
    
  } else if (s == "seizure_improve_freedom") {
    
    ar_merged <- read_csv(paste0(input.yaml$file_path,"/STX_EGRP_seizure_aed_1_month_ar_merged_seizure_free_FINAL_v8.csv"))
    
  } else if (s == "seizure_freedom") {
    
    ar_merged <- read_csv(paste0(input.yaml$file_path,"/STX_EGRP_seizure_aed_1_month_ar_merged_seizure_free_FINAL_v8.csv"))
    
  }
  
  ar_merged %>% select(patID, AED, time_int) %>% unique() %>% dplyr::group_by(AED) %>% dplyr::summarize(n=n()) %>% arrange(desc(n))
  ar_merged %>% select(patID, AED) %>% unique() %>% dplyr::group_by(AED) %>% dplyr::summarize(n=n()) %>% arrange(desc(n)) -> x
  
  
  # Empty data frame
  aed_effect <- data.frame(AED = NA, t = NA, pval = NA, OR = NA, OR_lower = NA, OR_upper = NA, 
                           n_pats = NA,
                           aed_pos_improve = NA, aed_pos_nimprove = NA, 
                           aed_neg_improve = NA, aed_neg_nimprove = NA)
  
  for (i in 1:length(unique(ar_merged$AED))) {
  
    if (s %in% c("seizure_improve", "seizure_improve_freedom")) {
      aed_effect <- aed_effect %>% rbind(comp_effect(unique(ar_merged$AED)[i], 0, 20))
    } else if (s == "seizure_freedom") {
      aed_effect <- aed_effect %>% rbind(comp_effect_maintain(unique(ar_merged$AED)[i], 0, 20))
    }
    
  }
  
  aed_effect <- aed_effect %>% 
    arrange(desc(OR)) %>% 
    filter(!is.na(AED)) %>%
    mutate(pval = as.numeric(pval),
           OR = as.numeric(OR),
           OR_lower = as.numeric(OR_lower),
           OR_upper = as.numeric(OR_upper),
           aed_pos_nimprove = as.numeric(aed_pos_nimprove),
           aed_neg_improve = as.numeric(aed_neg_improve),
           aed_neg_nimprove = as.numeric(aed_neg_nimprove)) %>% 
    mutate(pos_OR = case_when(OR > 1 ~ "yes", TRUE ~ "no")) %>% 
    filter(AED %in% c("LEVETIRACETAM", "PHENOBARBITAL", "ACTH", "TOPIRAMATE", "VIGABATRIN", "OXCARBAZEPINE",
                      "CLOBAZAM", "CANNABIDIOL", "KETO", "UKISS PRED"))
  
  write_csv(aed_effect, paste0(input.yaml$file_path,"/STXBP1_AED_", s, "_FINAL_v8.csv"))
  
}


## Supplementary (Ages 0-2 and 2-20)
for (s in c("seizure_improve", "seizure_improve_freedom", "seizure_freedom")) {
  
  if (s == "seizure_improve") {
    
    ar_merged <- read_csv(paste0(input.yaml$file_path,"/STX_EGRP_seizure_aed_1_month_ar_merged_FINAL_v8.csv")
    
  } else if (s == "seizure_improve_freedom") {
    
    ar_merged <- read_csv(paste0(input.yaml$file_path,"/STX_EGRP_seizure_aed_1_month_ar_merged_seizure_free_FINAL_v8.csv")
    
  } else if (s == "seizure_freedom") {
    
    ar_merged <- read_csv(paste0(input.yaml$file_path,"/STX_EGRP_seizure_aed_1_month_ar_merged_seizure_free_FINAL_v8.csv")
    
  }

  ar_merged %>% select(patID, AED, time_int) %>% unique() %>% dplyr::group_by(AED) %>% dplyr::summarize(n=n()) %>% arrange(desc(n))
  ar_merged %>% select(patID, AED) %>% unique() %>% dplyr::group_by(AED) %>% dplyr::summarize(n=n()) %>% arrange(desc(n)) -> x
  
  for (tf in c("early", "later")) {
    
    if (tf == "early") {
      si = 0; ei = 2
    } else if (tf == "later") {
      si = 2; ei = 20
    }
    
    # Empty data frame
    aed_effect <- data.frame(AED = NA, t = NA, pval = NA, OR = NA, OR_lower = NA, OR_upper = NA, 
                             n_pats = NA,
                             aed_pos_improve = NA, aed_pos_nimprove = NA, 
                             aed_neg_improve = NA, aed_neg_nimprove = NA)
    
    for (i in 1:length(unique(ar_merged$AED))) {
      
      if (s %in% c("seizure_improve", "seizure_improve_freedom")) {
        aed_effect <- aed_effect %>% rbind(comp_effect(unique(ar_merged$AED)[i], si, ei))
      } else if (s == "seizure_freedom") {
        aed_effect <- aed_effect %>% rbind(comp_effect_maintain(unique(ar_merged$AED)[i], si, ei))
      }
      
    }
    
    aed_effect <- aed_effect %>% 
      arrange(desc(OR)) %>% 
      filter(!is.na(AED)) %>%
      mutate(pval = as.numeric(pval),
             OR = as.numeric(OR),
             OR_lower = as.numeric(OR_lower),
             OR_upper = as.numeric(OR_upper),
             aed_pos_nimprove = as.numeric(aed_pos_nimprove),
             aed_neg_improve = as.numeric(aed_neg_improve),
             aed_neg_nimprove = as.numeric(aed_neg_nimprove)) %>% 
      mutate(pos_OR = case_when(OR > 1 ~ "yes", TRUE ~ "no")) %>% 
      filter(AED %in% c("LEVETIRACETAM", "PHENOBARBITAL", "ACTH", "TOPIRAMATE", "VIGABATRIN", "OXCARBAZEPINE",
                        "CLOBAZAM", "CANNABIDIOL", "KETO", "UKISS PRED"))
    
    write_csv(aed_effect, paste0(input.yaml$file_path,"/STXBP1_AED_", s, "_FINAL_v8_SUPP", si, "_", ei, ".csv"))
    
  }
  
}

## Supplementary (seizure types)
for (s in c("seizure_improve")) {
  
  for (tf in c("neonatal", "infantile", "focal")) {
    
    if (tf == "neonatal") {
      ar_merged <- read_csv(paste0(input.yaml$file_path,"/STX_EGRP_seizure_aed_1_month_ar_merged_Neonatal_Seizures_FINAL_v8.csv")
    } else if (tf == "infantile") {
      ar_merged <- read_csv(paste0(input.yaml$file_path,,"/STX_EGRP_seizure_aed_1_month_ar_merged_Infantile_Spasms_FINAL_v8.csv")
    } else if (tf == "focal") {
      ar_merged <- read_csv(paste0(input.yaml$file_path,"/STX_EGRP_seizure_aed_1_month_ar_merged_Focal_Impaired_Awareness_FINAL_v8.csv")
    }
  
    
    ar_merged %>% select(patID, AED, time_int) %>% unique() %>% dplyr::group_by(AED) %>% dplyr::summarize(n=n()) %>% arrange(desc(n))
    ar_merged %>% select(patID, AED) %>% unique() %>% dplyr::group_by(AED) %>% dplyr::summarize(n=n()) %>% arrange(desc(n)) -> x
    
    
    # Empty data frame
    aed_effect <- data.frame(AED = NA, t = NA, pval = NA, OR = NA, OR_lower = NA, OR_upper = NA, 
                             n_pats = NA,
                             aed_pos_improve = NA, aed_pos_nimprove = NA, 
                             aed_neg_improve = NA, aed_neg_nimprove = NA)
    
    for (i in 1:length(unique(ar_merged$AED))) {
      
      if (s %in% c("seizure_improve", "seizure_improve_freedom")) {
        aed_effect <- aed_effect %>% rbind(comp_effect(unique(ar_merged$AED)[i], 0, 20))
      } else if (s == "seizure_freedom") {
        aed_effect <- aed_effect %>% rbind(comp_effect_maintain(unique(ar_merged$AED)[i], 0, 20))
      }
      
    }
    
    aed_effect <- aed_effect %>% 
      arrange(desc(OR)) %>% 
      filter(!is.na(AED)) %>%
      mutate(pval = as.numeric(pval),
             OR = as.numeric(OR),
             OR_lower = as.numeric(OR_lower),
             OR_upper = as.numeric(OR_upper),
             aed_pos_nimprove = as.numeric(aed_pos_nimprove),
             aed_neg_improve = as.numeric(aed_neg_improve),
             aed_neg_nimprove = as.numeric(aed_neg_nimprove)) %>% 
      mutate(pos_OR = case_when(OR > 1 ~ "yes", TRUE ~ "no")) %>% 
      filter(AED %in% c("LEVETIRACETAM", "PHENOBARBITAL", "ACTH", "TOPIRAMATE", "VIGABATRIN", "OXCARBAZEPINE",
                        "CLOBAZAM", "CANNABIDIOL", "KETO", "UKISS PRED"))
    
    write_csv(aed_effect, paste0(input.yaml$file_path,"/STXBP1_AED_", s, "_FINAL_v8_SUPP", tf, ".csv"))
    
  }
  
}
                            
message("\n  ...comparative effectiveness analyses complete \n ")
stop = Sys.time()
stop - start                            


