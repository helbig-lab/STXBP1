# FDR analysis

# Input file
infile <- hpo_sig

if (fdr_adjust_or) {
  infile <- infile %>% filter(OR_adjusted > 1)
}

# Rank pval ascending
infile <- arrange(infile, pval)
infile$n <- 1:NROW(infile)

## Parameters
infile$h = NROW(infile)

# Alpha values
infile$p_1 = 0.05
infile$p_2 = 0.10
infile$p_3 = 0.20

# Corrected pval

infile$pval_fdr_5 <- (infile$n/infile$h)*infile$p_1
infile <- infile %>% mutate(sig_fdr_5 = case_when(pval < pval_fdr_5 ~ "SIG", TRUE ~ "NOT"))

infile$pval_fdr_10 <- (infile$n/infile$h)*infile$p_2
infile <- infile %>% mutate(sig_fdr_10 = case_when(pval < pval_fdr_10 ~ "SIG", TRUE ~ "NOT"))

infile$pval_fdr_20 <- (infile$n/infile$h)*infile$p_3
infile <- infile %>% mutate(sig_fdr_20 = case_when(pval < pval_fdr_20 ~ "SIG", TRUE ~ "NOT"))

infile <- infile %>% select(-c(n, h, p_1, p_2, p_3, pval_fdr_5, pval_fdr_10, pval_fdr_20))

# Save file "fdr_res
fdr_res <- infile

