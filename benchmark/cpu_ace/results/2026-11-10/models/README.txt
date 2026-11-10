In references to: https://arxiv.org/abs/2501.06116

DATA:
combined_data_oct24_2023.extxyz - Multiphase Ti dataset.
ti64_combined_data_Apr08.extxyz - Multipahse Ti-6Al-4V dataset  (This corresponds to ACE models in pre-print).
ti64_combined_data_May30.extxyz - Multipahse Ti-6Al-4V dataset (This corresponds to ACE2 models in pre-print).

GAP models:
ti64_cd_Apr08_customSP_3700SP.xml* - This is a GAP model for the Ti64 dataset prior to targeted data for phonon benchmarks.
Ti_cd_oct23_2024_awfr_0.02_min_0.025_cwvr_0.05_min_0.05_SPbyConf_3300sp.xml* - This is a GAP model for the multiphase Ti dataset.

ACE models: (.json for ACE models as loaded in with Python ACEpotentials, .yace + .table for use in LAMMPS).
Ti_combined_data_oct24_2023_ace_order_3_degree_20_rcut_6.0_cw_moreV.* - ACE model for multiphase Ti
ti64_combined_data_Apr08_rcut_6.0_order_3_degree_15_custom_reg.* - ACE model for Ti-6Al-4V
ti64_combined_data_May30_rcut_6.0_order_3_degree_15_custom_reg.*  - ACE2 model for Ti-6Al-4V