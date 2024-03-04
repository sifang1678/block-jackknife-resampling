# Sample scripts for conducting block jackknife resampling Mendelian randomization using individual-level data

This repo contains some scripts used in our paper: Fang, S., Hemani, G., Richardson, T. G., Gaunt, T. R., & Davey Smith, G. (2023) Evaluating and implementing block jackknife resampling Mendelian randomization to mitigate bias induced by overlapping samples. Human molecular genetics, 32(2), 192–203. https://doi.org/10.1093/hmg/ddac186

## File descriptions

- **script_for_JF_PGS.sh** contains scripts for performing clumping and PGS calculation based on block jackknifing GWAS (N<sub>block</sub>=10). JF GWASs based on every N-1 blocks of individuals need to be performed prior to using this script.
- **script_for_JF_MR.R** contains scripts for performing multivariable MR investigating the effects of childhood and adult body sizes on circulating biomarkers (Biochemistry) in the UK Biobank.
- **functions.R** contains some in-house functions for performing individual-level MR and multivariable MR used in our paper.
