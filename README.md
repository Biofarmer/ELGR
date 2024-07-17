1. The files in data fold are input for the script in script folder.
2. The description and the output for each script.
   
   1_preterm_ARGs_profile_host_eggnog.R
   This script is used to analyze the composition, host, and functions of ARGs in preterm infants, including Figure 1b, Figure 4d, Supplementary Figure 8, Figure 1c, Supplementary Figure 1d, Supplementary Figure 1f, Supplementary Figure 9a, Supplementary Figure 1e, Supplementary Figure 1g, Supplementary Figure 1h, Supplementary Figure 1b, Supplementary Figure 2b, Supplementary Figure 2a, Figure 1d, Supplementary Figure 3a, Supplementary Figure 3b, Supplementary Figure 3c, Supplementary Figure 3d, Supplementary Figure 4.

   2_term_ARGs_profile_host_comparisonWithpreterm.R
   This script is used to analyze the composition and host of ARGs in full-term infants, and also the comparisons between preterm and full-term infants, including Supplementary Figure 5d, Figure 2f, Supplementary Figure 6f, Supplementary Figure 5e, Figure 2g, Supplementary Figure 5a, Figure 2a, Figure 2c, Supplementary Figure 5b, Figure 2b, Figure 2e, Figure 2d, Supplementary Figure 6a, Supplementary Figure 6b, Supplementary Figure 6e, Supplementary Figure 6c, Supplementary Figure 5c, Supplementary Figure 6d, Figure 3, Supplementary Figure 7b, Supplementary Figure 7a.

   3_preterm_ARGs_dynamics_covariate_maaslin2.R
   This script is used to analyze the dynamics and covariates of ARGs in preterm infants, including Figure 4b, Figure 4a, Figure 4c, Supplementary Figure 9b, Supplementary Figure 9c, Supplementary Figure 10a, Supplementary Figure 10b, Figure 4e, Supplementary Figure 9d, Figure 4f.

   4_nec_diversity_signatures.R
   This script is used to analyze differences between cases and controls, including, Figure 4b, Figure 4a, Figure 4c, Supplementary Figure 9b, Supplementary Figure 9c, Supplementary Figure 10a, Supplementary Figure 10b, Figure 4e, Supplementary Figure 9d, Figure 4f, Supplementary Figure 11a, Supplementary Figure 11b, Supplementary Figure 11c, Figure 5a, Figure 5c, Figure 5b, Supplementary Figure 11d, Figure 5d, Supplementary Figure 11e.

   5_1_nec_rf_discovery.R
   This script is used to produce the AU-ROC from four public cohorts for discovery and validation based on features from subtype, species, subtype + species, respectively, with intra-cohort, combined cohorts of 10-times fivefold stratified cross-validation, and LOSO prediction. The output files are auc_cb_intra.csv, auc_loso_cb.csv, which are used for 5_2_nec_rf_validation.R.

   5_2_nec_rf_validation.R
   This script is used to analyze the predictive performance from public cohorts and validation cohort, including Supplementary Figure 12a, Figure 5e, Supplementary Figure 12b, Figure 6a, Figure 6b, Figure 6c.
