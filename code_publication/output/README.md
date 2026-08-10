# output/

Analysis outputs are written here by the scripts. Expected files:

| Script | Output |
|---|---|
| 01_mouse_HSPC_processing.R | mouse_HSPC_annotated.rds, allmarkers_celltype_mouse.csv, gsva_res_mouse.csv, mouse_celltype_proportions.csv |
| 02_mouse_annotation_crossvalidation.R | mouse_label_transfer_predictions.rds |
| 03_mouse_HSC_subclustering.R | mouse_HSC.rds, HSC_module_correlation.csv |
| 04_mouse_Tcell_analysis.R | mouse_Tcell_CD4.rds, mouse_CD4Tcells.rds, Treg_DE.csv, Treg_GO_BP.csv, Treg_GSVA_matrix.xls, Treg_GSVA_diff.csv |
| 05_human_data_integration.R | human_bone_marrow.rds, gsva_res_human.csv, human_celltype_proportions.csv |
| 06_human_HSC_survival_analysis.R | human_HSC.rds, HSC_CD69high_vs_low_markers.csv, human_HSC_GSVA_diff.csv, HSC_CD69high_multivariate_Cox.csv |
| 07_spatial_transcriptomics_analysis.R | spatial/ST_WT/ and spatial/ST_Kras/ (RCTD.rds, Spatial_CellType.tsv, HSC_neighbor_cell_proportion.tsv) |
| 08_cell_cell_communication.R | reads cellphonedb/ST_WT and ST_Kras (CellPhoneDB output) |

CellPhoneDB must be run separately (Python) on the RCTD-annotated spots; place
its output (count_network.txt, pvalues.txt, means.txt, significant_means.txt)
in output/cellphonedb/ST_WT/ and output/cellphonedb/ST_Kras/ before running
script 08. See README.md for the CellPhoneDB command.