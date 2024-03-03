listOfPackages <- c("data.table", "tidyverse", "ivpack")
for (i in listOfPackages) {
  if (!i %in% installed.packages()) {
    install.packages(i, dependencies = TRUE)
  }
}

library("data.table")
library("tidyverse")
library("ivpack")

setwd("/path/to/your/working/directory/")

source("parameters.R") # containing the definition of paths used in the scripts below
source("functions.R") # containing functions to perform IV analysis using individual level data

# Load PRS data
bmi10_prs_female <- fread(paste0(data_dir, "UKBB/prs/childhood_bmi/childhood_bmi_female_jackknifing_prs.profile"))
bmi10_prs_male <- fread(paste0(data_dir, "UKBB/prs/childhood_bmi/childhood_bmi_male_jackknifing_prs.profile"))
bmi10_prs_all <- rbind(bmi10_prs_female, bmi10_prs_male)

adult_bmi_prs_females <- fread(paste0(data_dir, "UKBB/prs/adult_bmi_categorical/adult_bmi_categorical_females_jackknifing_prs.profile"))
adult_bmi_prs_males <- fread(paste0(data_dir, "UKBB/prs/adult_bmi_categorical/adult_bmi_categorical_males_jackknifing_prs.profile"))
adult_bmi_prs_all <- rbind(adult_bmi_prs_females, adult_bmi_prs_males)

# Load phenotype data
adult_bmi <- fread(paste0(data_dir, "UKBB/phenotype/bmi_categorical_matching_bmi10_grp/ukbb_adult_bmi_categorical.txt"))
bmi10 <- fread(paste0(data_dir, "UKBB/phenotype/bmi10/ukbb_bmi10.txt"))
biomarkers <- fread(paste0(data_dir, "UKBB/blood_biochem_data_clean.csv"))

# Load covariates data
covar <- fread(paste0(data_dir, "UKBB/plink_covariates.txt")) # including age, sex, genotyping chip, genetic PCs

# Filter biomarkers based on sample size
to_remove <- names(which(colMeans(is.na(biomarkers)) > 0.5))
biomarkers <- biomarkers[, (to_remove):=NULL]

# Merge data
adult_bmi_data <- merge(adult_bmi_prs_all[,c("FID","IID","SCORE")], adult_bmi[,c("FID","IID","BMI_adult")], by = c("FID","IID"), all.x = T)
names(adult_bmi_data)[3] <- c("adult_bmi_prs")

bmi10_data <- merge(bmi10_prs_all[,c("FID","IID","SCORE")], bmi10[,c("FID","IID","bmi10")], by = c("FID","IID"), all.x = T)
names(bmi10_data)[3] <- c("bmi10_prs")

merge1 <- merge(adult_bmi_data, bmi10_data, by = c("FID","IID"), all = T)
merge2 <- merge(merge1, covar, by = c("FID","IID"), all.x = T)
merge_all <- merge(merge2, biomarkers, by = c("FID","IID"), all.x = T)

biomarker_ids <- names(biomarkers[,c(-1, -2)])
covar_names <- c(names(covar)[c(-1,-2,-5)])

data_m <- merge_all[sex == 1]
data_f <- merge_all[sex == 2]

# Univariable and multivariable MR
res_all <- mvmr_using_prs(exposure_list = c("bmi10", "BMI_adult"),
                          exposure_prs_list = c("bmi10_prs", "adult_bmi_prs"),
                          outcome_list = biomarker_ids, covar_list = covar_names,
                          data = merge_all)
res_all$sex <- "all"

res_m <- mvmr_using_prs(exposure_list = c("bmi10", "BMI_adult"),
                        exposure_prs_list = c("bmi10_prs", "adult_bmi_prs"),
                        outcome_list = biomarker_ids, covar_list = covar_names,
                        data = data_m)
res_m$sex <- "males"

res_f <- mvmr_using_prs(exposure_list = c("bmi10", "BMI_adult"),
                        exposure_prs_list = c("bmi10_prs", "adult_bmi_prs"),
                        outcome_list = biomarker_ids, covar_list = covar_names,
                        data = data_f)
res_f$sex <- "females"

res <- rbindlist(list(res_m, res_f, res_all))
fwrite(res, paste0(output_dir, "bmi10_biomarker_mvmr_results_", Sys.Date(), ".csv"))

