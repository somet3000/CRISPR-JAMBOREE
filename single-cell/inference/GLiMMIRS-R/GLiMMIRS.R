#' Negative Binomial Test
#'
#' @param mudata_input_fp Path to input MuData
#' @param mudata_output_fp Path to output MuData
#'
#' @return This function is called for its side effect of writing the output MuData.
# TODO: confirm that Low MOI and high MOI case are the same for negative binomial
perform_GLiMMIRS <- function(mudata_input_fp, mudata_output_fp) {
    
  # read in input mudata
  mudata <- MuData::readH5MU(mudata_input_fp)

  # Rename primary assay to 'counts' for both guide and gene
  if(is.null(SummarizedExperiment::assayNames(mudata[["gene"]]))){
    SummarizedExperiment::assayNames(mudata[['gene']]) <- 'counts'    
  } else{
    SummarizedExperiment::assayNames(mudata[['gene']])[[1]] <- 'counts'  
  }
  SummarizedExperiment::assayNames(mudata[['guide']])[[1]] <- 'counts'

  # get S phase and G2M phase genes from Seurat
  s_genes <- Seurat::cc.genes$s.genes
  g2m_genes <- Seurat::cc.genes$g2m.genes

  # separate expression matrix to get gene names
  expr_matrix <- SummarizedExperiment::assay(
    mudata[['gene']],
    'counts'
  )
  ensembl_gene_ids <- rownames(expr_matrix)

  # use BioMart to convert Ensembl gene IDs to HGNC symbols
  mart <- biomaRt::useMart('ensembl')
  mart <- biomaRt::useDataset('hsapiens_gene_ensembl', mart)
  gene_hgnc_symbols <- biomaRt::getBM(
    filters = 'ensembl_gene_id',
    attributes = c('ensembl_gene_id', 'hgnc_symbol'),
    values = ensembl_gene_ids,
    mart = mart
  )

  # merge with original gene symbols to retain gene ordering
  gene_ensembl_hgnc <- merge(
    data.frame(ensembl_gene_ids),
    gene_hgnc_symbols,
    all.x = TRUE,
    by.x = 'ensembl_gene_ids',
    by.y = 'ensembl_gene_id',
    sort = FALSE
  )

  # merge again to handle NA values (if applicable)
  gene_ensembl_hgnc <- merge(
    data.frame(ensembl_gene_ids),
    gene_ensembl_hgnc,
    by = 'ensembl_gene_ids',
    sort = FALSE
  )

  # input Ensembl gene IDs for genes without HGNC
  for (i in 1:nrow(gene_ensembl_hgnc)) {
    if (gene_ensembl_hgnc[i, 'hgnc_symbol'] == '') {
      gene_ensembl_hgnc[i, 'hgnc_symbol'] = gene_ensembl_hgnc[
        i, 
        'ensembl_gene_ids'
      ]
    }
  }

  # set rownames of expression matrix to HGNC symbols for Seurat
  rownames(expr_matrix) <- gene_ensembl_hgnc$hgnc_symbol

  # initialize Seurat object with count data
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = expr_matrix
  )

  # perform Seurat preprocessing
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = 'vst')
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = rownames(seurat_obj))
  
  # perform cell cycle scoring
  # seurat_obj <- Seurat::CellCycleScoring(
  #   seurat_obj,
  #   s.features = s_genes,
  #   g2m.features = g2m_genes,
  #   set.ident = TRUE
  # )

  # get MOI
  moi <- MultiAssayExperiment::metadata(mudata[["guide"]])$moi

  # In low-MOI case, extract control cells as those containing an NT gRNA
  if (moi == "low") {
    non_targeting_guides <- SummarizedExperiment::rowData(mudata[["guide"]]) |>
      as.data.frame() |>
      dplyr::filter(targeting == "FALSE") |>
      rownames()
    nt_grna_presence <- SummarizedExperiment::assay(
      mudata[["guide"]],
      "guide_assignment"
    )[non_targeting_guides, ] |>
      apply(MARGIN = 2, FUN = max)
    control_cells <- names(nt_grna_presence)[nt_grna_presence == 1]
  }

  # Extract pairs to test
  pairs_to_test <- MultiAssayExperiment::metadata(mudata)$pairs_to_test |> 
    as.data.frame()

  # Initialize test results data frame based on pairs_to_test
  test_results <- pairs_to_test |>
    dplyr::mutate(p_value = NA_real_, l2fc = NA_real_)
    
  # Carry out the negative binomial regression test for each pair
  for (pair_idx in 1:nrow(pairs_to_test)) {
        
    # Extract the gene and element to be tested
    gene_id <- pairs_to_test[pair_idx, "gene_id"]
    intended_target_name <- pairs_to_test[pair_idx, "intended_target_name"]

    # find vector with whether element is targeted or not
    grnas_targeting_element <- SummarizedExperiment::rowData(mudata[["guide"]]) |>
      as.data.frame() |>
    dplyr::filter(intended_target_name == !!intended_target_name) |>
      rownames()
    element_targeted <- SummarizedExperiment::assay(
      mudata[["guide"]],
      "guide_assignment"
    )[grnas_targeting_element, ] |>
      apply(MARGIN = 2, FUN = max)
    
    # define treatment cells for low MOI case
    treatment_cells <- names(element_targeted)[element_targeted == 1]

    # get gene expression for gene of interest
    gene_expression = SummarizedExperiment::assay(
      mudata[["gene"]],
      "counts"
    )[gene_id, ]

    # get total gene expression for each cell
    umis_per_cell <- SummarizedExperiment::colData(mudata[["gene"]])
    umis_per_cell <- umis_per_cell[, 'total_gene_umis']
        
    # create modeling data frame
    model_df <- data.frame(cbind(element_targeted, gene_expression, umis_per_cell))

    # subset cells in low MOI case - cells with guide and NT guide
    if (moi == "low") {
      model_df <- model_df[c(control_cells, treatment_cells), ]
    }
    
    # run negative binomial regression
    mdl <- MASS::glm.nb(gene_expression ~ element_targeted + offset(log(umis_per_cell)), model_df)

    # save L2FC and p-value
    mdl.coeffs <- summary(mdl)$coefficients
    test_results[pair_idx, "p_value"] <- mdl.coeffs['element_targeted', 'Pr(>|z|)']
    mdl.intercept <- mdl.coeffs['(Intercept)', 'Estimate']
    mdl.beta1 <- mdl.coeffs['element_targeted', 'Estimate']
    mdl.l2fc <- exp(mdl.intercept + mdl.beta1) / exp(mdl.intercept)
    test_results[pair_idx, 'l2fc'] <- mdl.l2fc
  }

  # Add output to MuData and write to disk
  mudata_output <- mudata
  MultiAssayExperiment::metadata(mudata_output)$test_results <- test_results
  MuData::writeH5MU(mudata_output, mudata_output_fp)
}
       