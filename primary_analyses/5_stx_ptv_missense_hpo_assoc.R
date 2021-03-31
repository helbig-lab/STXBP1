library(tidyverse)
library(Hmisc)

# Required files

stx_fams <- read_csv(paste0("STXBP1_full_base_v", vs, ".csv"))

# stx_base <- read_csv(paste0("full_base_v", vs, ".csv"))
stx_prop <- read_csv(paste0("full_prop_v", vs, ".csv"))
ic <- read_csv(paste0("pos_IC_v", vs, ".csv"))


variants <- stx_prop %>%
  left_join(stx_fams %>% select(famID, var_type_2)) %>%
  filter(!is.na(var_type_2)) %>%
  unique()
  

## Analysis
hpo_sig <- matrix(nrow=0,ncol=13) %>% as.data.frame()
names(hpo_sig) <- c('HPO','var', 'yes_var', 'no_var', 'yes_cohort', 'no_cohort', 
                    'pval',"var_freq", "cohort_freq", "OR", "OR_adjusted", 
                    "ORlower", "OR.upper")

for(g in 1:length(unique(variants$var_type_2))){
  
  var_x = unique(variants$var_type_2)[g]
  
  pats_var <- variants %>% filter(var_type_2 == var_x)
  pats <- pats_var$famID %>% unique()
  yesg_hpo <- variants %>% filter(famID %in% pats)
  
  nog_hpo <- variants %>% filter(famID %nin% pats)
  nog_pats <- nog_hpo$famID %>% unique()
  
  # Hpo count
  y_hpo_count <- yesg_hpo %>% count(HPO)
  n_hpo_count <- nog_hpo %>% count(HPO)
  
  hpo_1 <- matrix(nrow=nrow(y_hpo_count),ncol=13) %>% as.data.frame()
  names(hpo_1) <- c('HPO','var', 'yes_var', 'no_var', 'yes_cohort', 'no_cohort', 
                    'pval',"var_freq", "cohort_freq", "OR", "OR_adjusted", 
                    "OR.lower", "OR.upper")
  hpo_1$HPO <- y_hpo_count$HPO
  hpo_1$var <- var_x
  
  for (h in 1:nrow(hpo_1)) {
    
    # Hpo term
    hp <- hpo_1$HPO[h]
    
    fish <- matrix(ncol=2,nrow=2) %>% as.data.frame()
    names(fish) <- c('hpo_present','hpo_absent')
    row.names(fish) <- c('var_present','var_absent')
    
    fish['var_present','hpo_present'] <- y_hpo_count %>% 
      filter(HPO == hp) %>% pull(n)
    
    fish['var_absent','hpo_present'] <- 0
    if (hp %in% n_hpo_count$HPO) {
      fish['var_absent','hpo_present'] <- n_hpo_count %>% filter(HPO == hp) %>% pull(n)
    }
    
    fish['var_present','hpo_absent'] <- length(pats) - fish['var_present','hpo_present']
    
    fish['var_absent','hpo_absent'] <- length(nog_pats) - fish['var_absent','hpo_present']
    
    # Counts
    hpo_1$yes_var[h] = fish['var_present','hpo_present']
    hpo_1$no_var[h] = fish['var_present','hpo_absent']
    hpo_1$yes_cohort[h] = fish['var_absent','hpo_present']
    hpo_1$no_cohort[h] = fish['var_absent','hpo_absent']
    hpo_1$var_freq[h] = fish['var_present','hpo_present']/length(pats)
    hpo_1$cohort_freq[h] = fish['var_absent','hpo_present']/length(nog_pats)
    
    # Fisher's test
    test.res <- fisher.test(fish)
    hpo_1$pval[h] <- test.res$p.value
    hpo_1$OR[h] <- test.res$estimate
    
    # Odds ratio adjusted - get estimate for Inf odds ratio
    if (hpo_1$OR[h] == Inf) {
      fish.adj = fish
      fish.adj['var_absent','hpo_present'] <- 1
      fish.adj['var_absent','hpo_absent'] <- fish['var_absent','hpo_absent'] - 1
      hpo_1$OR_adjusted[h] <- fisher.test(fish.adj)$estimate
    } else {
      hpo_1$OR_adjusted[h] = hpo_1$OR[h]
    }
    
    # 95% CI for OR
    hpo_1$OR.lower[h] <- test.res$conf.int[1]
    hpo_1$OR.upper[h] <- test.res$conf.int[2]
    
  }
  
  hpo_sig <- hpo_sig %>% rbind(hpo_1)
}

hpo_sig <- hpo_sig %>% left_join(ic %>% select(HPO, def))

hpo_sig <- hpo_sig %>% select(HPO, def, var, pval, OR, OR_adjusted, OR.lower, OR.upper,
                              var_freq, cohort_freq, yes_var, no_var, yes_cohort, no_cohort)


hpo_sig <- hpo_sig[order(hpo_sig$pval),]

# hpo_sig <- hpo_sig %>% filter(OR > 1)
# write.csv(hpo_sig, "ptv_missense_sig_assoc_base.csv", row.names = F)
# write.csv(hpo_sig, "ptv_missense_sig_assoc_prop.csv", row.names = F)

## FDR
fdr_adjust_or = FALSE
source(paste0("/Volumes/helbig_lab/projects/STXBP1/v_scripts/primary/FDR.R"))
write.csv(fdr_res, paste0("output/stx_ptv_missense_hpo_assoc_prop_v", 
                          vs, "_full.csv"), row.names = FALSE)
# fdr_adjust_or = TRUE
# source(paste0("/Volumes/helbig_lab/projects/STXBP1/v_scripts/primary/FDR.R"))
# write.csv(fdr_res, paste0("output/stx_ptv_missense_hpo_assoc_prop_v", 
#                           vs, ".csv"), row.names = FALSE)

