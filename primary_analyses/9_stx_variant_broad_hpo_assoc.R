start <- Sys.time()
message(" \n Running variant broad association analysis... \n ")

library(tidyverse)
library(Hmisc)

stx_fams <- read_csv(paste0(input.yaml$file_path,"STXBP1_full_base_v.csv")) %>% 
  select(famID, variant_broad) %>% 
  unique()

stx_prop <- read_csv(paste0(input.yaml$file_path,"full_prop_v.csv"))

hpo_def <- read_csv(paste0(input.yaml$file_path, "pos_IC_v.csv")) %>% 
  select(HPO, def, freq_prop)

stx_merged <- stx_prop %>% left_join(stx_fams)

## Analysis

pheno_groups <- stx_merged$variant_broad %>% unique()

pheno_groups <- pheno_groups %>% na.exclude()

freq_all <- hpo_def %>% 
  filter(HPO %in% stx_merged$HPO) %>% 
  rename(ALL = freq_prop)

for (i in 1:length(pheno_groups)) {
  
  pg = pheno_groups[i]
  
  pg_terms <- stx_merged %>% filter(variant_broad == pg)
  pats = pg_terms$famID %>% unique()
  
  term_freq <- pg_terms %>% 
    dplyr::group_by(HPO) %>% 
    dplyr::summarize(freq=n()/length(pats))
  
  freq_all <- freq_all %>% left_join(term_freq)
  freq_all[is.na(freq_all)] = 0
  
  names(freq_all)[names(freq_all) == 'freq'] <- pg
  
}

write.csv(freq_all, paste0(input.yaml$output_dir,"stx_recurrent_variants_broad_freq_v.csv"), row.names = F)

## HPO associations

sig_test <- matrix(nrow=0,ncol=13) %>% as.data.frame()
names(sig_test) <- c('HPO','VAR_GROUP', 'pval', "OR", "OR_adjusted", "OR.lower", "OR.upper",
                     'present_group', 'absent_group', 'present_other', 'absent_other', "group_freq", "other_freq")

for (i in 1:length(pheno_groups)) {
  
  pg = pheno_groups[i]
  pg_terms <- stx_merged %>% filter(variant_broad == pg)
  pats <- pg_terms$famID %>% unique()
  yes_hpo <- stx_merged %>% filter(famID %in% pats)
  
  no_hpo <- stx_merged %>% filter(famID %nin% pats)
  no_pats <- no_hpo$famID %>% unique()
  
  y_hpo_count <- yes_hpo %>% count(HPO)
  n_hpo_count <- no_hpo %>% count(HPO)
  
  hpo_1 <- matrix(nrow=nrow(y_hpo_count), ncol=13) %>% as.data.frame()
  names(hpo_1) <- c('HPO','VAR_GROUP', 'pval', "OR", "OR_adjusted", "OR.lower", "OR.upper",
                    'present_group', 'absent_group', 'present_other', 'absent_other', "group_freq", "other_freq")
  hpo_1$HPO <- y_hpo_count$HPO
  hpo_1$VAR_GROUP <- pg
  
  for(h in 1:nrow(hpo_1)){
    
    hp <- hpo_1$HPO[h]
    
    fish <- matrix(ncol=2,nrow=2) %>% as.data.frame()
    names(fish) <- c('hpo_present','hpo_absent')
    row.names(fish) <- c('group_present','group_absent')
    
    fish['group_present','hpo_present'] <- y_hpo_count %>% filter(HPO == hp) %>% pull(n)
    
    fish['group_absent','hpo_present'] <- 0
    if (hp %in% n_hpo_count$HPO) {
      fish['group_absent','hpo_present'] <- n_hpo_count %>% filter(HPO == hp) %>% pull(n)
    }
    
    fish['group_present','hpo_absent'] <- length(pats) - fish['group_present','hpo_present']
    
    fish['group_absent','hpo_absent'] <- length(no_pats) - fish['group_absent','hpo_present']
    
    # Total counts
    hpo_1$present_group[h] = fish['group_present','hpo_present']
    hpo_1$absent_group[h] = fish['group_present','hpo_absent']
    hpo_1$present_other[h] = fish['group_absent','hpo_present']
    hpo_1$absent_other[h] = fish['group_absent','hpo_absent']
    hpo_1$group_freq[h] = fish['group_present','hpo_present']/length(pats)
    hpo_1$other_freq[h] = fish['group_absent','hpo_present']/length(no_pats)
    
    # Fisher's test
    test.res <- fisher.test(fish)
    hpo_1$pval[h] <- test.res$p.value
    hpo_1$OR[h] <- test.res$estimate
    
    # Odds ratio adjusted - get estimate for Inf odds ratio
    if (hpo_1$OR[h] == Inf) {
      fish.adj = fish
      fish.adj['group_absent','hpo_present'] <- 1
      fish.adj['group_absent','hpo_absent'] <- fish['group_absent','hpo_absent'] - 1
      hpo_1$OR_adjusted[h] <- fisher.test(fish.adj)$estimate
    } else {
      hpo_1$OR_adjusted[h] = hpo_1$OR[h]
    }
    
    # 95% CI for OR
    hpo_1$OR.lower[h] <- test.res$conf.int[1]
    hpo_1$OR.upper[h] <- test.res$conf.int[2]
    
  }
  
  sig_test <- sig_test %>% rbind(hpo_1)
  
}

hpo_sig <- sig_test %>% 
  left_join(hpo_def %>% select(HPO, def)) %>% 
  select(HPO, def, VAR_GROUP, pval, OR, OR_adjusted, OR.lower, OR.upper, group_freq, 
         other_freq, present_group, absent_group, present_other, absent_other) 

hpo_sig <- hpo_sig[order(sig_test$pval),]


## FDR

fdr_adjust_or = FALSE
source(paste0(input.yaml$prime_dir,"FDR.R"))

write.csv(fdr_res, paste0(input.yaml$output_dir,"stx_recurrent_variants_broad_hpo_assoc_v.csv"), 
          row.names = FALSE)
          
message("\n ...variant broad association analysis complete \n ")
stop = Sys.time()
stop - start

