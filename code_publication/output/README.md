# output/

Analysis outputs are written here by the scripts. Expected files:

| Script | Output |
|---|---|
| 01_mouse_HSPC_processing.R | mouse_HSPC_annotated.rds, allmarkers_celltype_mouse.csv, gsva_res_mouse.csv, mouse_celltype_proportions.csv |
| 02_mouse_annotation_crossvalidation.R | mouse_label_transfer_predictions.rds |
| 03_mouse_HSC_subclustering.R | mouse_HSC.rds, HSC_module_gene_overlap.csv, HSC_TF_target_overlap.csv, HSC_module_correlation.csv |
| 04_mouse_Tcell_analysis.R | mouse_Tcell_CD4.rds, mouse_CD4Tcells.rds, Treg_DE.csv, Treg_GO_BP.csv, Treg_GSVA_matrix.tsv, Treg_GSVA_descriptive.csv |
| 05_human_data_integration.R | human_bone_marrow.rds, human_cluster_annotation_map.csv, gsva_res_human.csv, human_celltype_proportions.csv |
| 06_human_HSC_survival_analysis.R | human_HSC.rds, HSC_CD69high_vs_low_markers.csv, human_HSC_GSVA_descriptive.csv, HSC_CD69high_survival_table.csv, Table_S6_univariable_Cox.csv, Figure_8H_multivariable_Cox.csv, Figure_8H_model_audit.csv, Figure_8H_proportional_hazards_test.csv |
| 07_spatial_transcriptomics_analysis.R | spatial/ST_WT/ and spatial/ST_Kras/ (RCTD.rds, Spatial_CellType.tsv, HSC_neighbor_cell_proportion.tsv) |
| 08_prepare_cellphonedb_inputs.py | cellphonedb/ST_WT and ST_Kras (sparse h5ad, metadata, orthologue-mapping report) |
| 08_run_cellphonedb.py | standardised CellPhoneDB result tables and cellphonedb_run_metadata.json |
| 08_cell_cell_communication.R | CellPhoneDB network and HSC-centred dot-plot PDFs |
| 09_NB4_bulk_RNAseq.R | NB4_CD69OE_DESeq2_all_genes.csv, NB4_CD69OE_DESeq2_significant.csv, NB4_CD69OE_PCA.pdf, NB4_CD69OE_sessionInfo.txt |

CellPhoneDB must be run separately after `08_prepare_cellphonedb_inputs.py`;
place its output (`count_network.txt`, `pvalues.txt`, `means.txt` and
`significant_means.txt`) in the corresponding sample directory before running
the R visualisation script. See README.md for the versioned command.
