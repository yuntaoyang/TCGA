#---- Load libraries------------------------------------------------------------
# install.packages("survival")
# install.packages("survminer")
# install.packages("maxstat")
# BiocManager::install("edgeR")
library(dplyr)
library(survival)
library(survminer)
library(maxstat)
library(edgeR)
#---- Set up parameters --------------------------------------------------------
project <- 'TCGA-BLCA'
target <- 'SMARCB1'
dir_data <- './data'
dir_out <- './out'
if (!dir.exists(dir_out)) {
  dir.create(dir_out)
}
#---- Load TCGA data -----------------------------------------------------------
df_counts <- read.csv(file.path(dir_data, paste0(project, '_counts.csv')), 
                      check.names = FALSE)
df_clinical <- read.csv(file.path(dir_data, paste0(project,'_clinical.csv')))
df_biospecimen <- read.csv(file.path(dir_data, paste0(project,'_biospecimen.csv')))
#---- Data pre-processing ------------------------------------------------------
# Survival data
survival <- df_clinical %>%
  select(submitter_id, days_to_last_follow_up, days_to_death, vital_status) %>%
  filter(!is.na(days_to_last_follow_up) | !is.na(days_to_death), vital_status != "Not Reported") %>% # remove samples with missing values
  mutate(survival_days = pmax(days_to_last_follow_up, days_to_death, na.rm = TRUE), # use the larger days in days_to_last_follow_up and days_to_death
         survival_days = ifelse(survival_days < 0, 0, survival_days)) # days < 0 = 0
# Biospecimen data:  primary cancer samples
primary_cancer <- df_biospecimen %>%
  select(c(submitter_id ,sample_type)) %>%
  filter(sample_type == 'Primary Tumor')
# Gene expression data
expdat <- df_counts %>%
  filter(!is.na(gene_name))
expdat_for_logcpm <- expdat %>%
  select(-c(gene_id, gene_name))
# Log2 Count Per Million
log2cpm <- cpm(expdat_for_logcpm, log=TRUE)
rownames(log2cpm) <- expdat$gene_name
log2cpm_target <- log2cpm[rownames(log2cpm) == target, , drop = FALSE]
#---- Determine cutoff ---------------------------------------------------------
group <- data.frame(submitter_id = substring(colnames(log2cpm_target), 1, 16), 
                    log2cpm = as.numeric(log2cpm_target))
# Merge with primary cancer samples
group <- merge(group, primary_cancer)
group$submitter_id <- substring(group$submitter_id, 1, 12)
# Merge with survival data
group <- merge(survival, group)
# Get an average for patients with multiple primary cancer samples
group_average <- aggregate(log2cpm ~ submitter_id, group, mean)
group <- select(group, -c(log2cpm)) %>%
  distinct() %>%
  merge(group_average) %>%
  mutate(vital_status_code = if_else(vital_status == 'Alive', 1, 2))
# Convert days to years
group$survival_years <- group$survival_days/365
# maxstat estimate cutpoint
test <- maxstat.test(Surv(survival_years, vital_status_code) ~ log2cpm,
                     data=group, smethod="LogRank")
cutoff <- as.numeric(test$estimate)
#---- survival analysis --------------------------------------------------------
group <- group %>%
  mutate(sample_group = if_else(log2cpm > cutoff, paste0(target, '_High'), paste0(target, '_Low'))) %>%
  mutate(sample_group_code = if_else(sample_group == paste0(target,'_Low'), 1, 2))
# Save the survival data
write.csv(group, file.path(dir_out, paste0(project, '_survival.csv')), row.names = FALSE)
fit <- survfit(Surv(survival_years, vital_status_code) ~ sample_group_code, data = group)
print(fit)
survival_curve <- ggsurvplot(fit, pval = TRUE, 
                             conf.int = FALSE, 
                             conf.int.style = "step",
                             break.time.by = 2, 
                             risk.table = "absolute", 
                             risk.table.y.text.col = TRUE, 
                             risk.table.y.text = FALSE,
                             ncensor.plot = FALSE, 
                             surv.median.line = "none", 
                             fontsize = 2.5,
                             xlab = "Time in years", 
                             xlim = c(0, 10), 
                             legend.labs = c(paste0(target,'_Low'), paste0(target,'_High')), 
                             palette = c("#E7B800", "#2E9FDF"), 
                             ggtheme = theme_classic())
pdf(file = file.path(dir_out, paste0(project, '_survival.pdf')), onefile = FALSE,
    width = 4, height = 4)
print(survival_curve)
dev.off()
