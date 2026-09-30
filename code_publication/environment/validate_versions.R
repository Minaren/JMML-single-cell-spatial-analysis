# Run from the code_publication root before executing the analysis scripts.

if (getRversion() < "4.5.0") {
  stop("R >= 4.5.0 is required for the GSVA 2.2.0/Bioconductor 3.21 environment.")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  stop("BiocManager is required to validate the Bioconductor release.")
}
if (as.character(BiocManager::version()) != "3.21") {
  stop("Bioconductor 3.21 is required; found ", BiocManager::version(), ".")
}

exact <- c(Seurat = "4.3.0", GSVA = "2.2.0")
missing <- names(exact)[!vapply(names(exact), requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required package(s): ", paste(missing, collapse = ", "))

observed <- vapply(names(exact), function(pkg) as.character(packageVersion(pkg)), character(1))
wrong <- names(exact)[observed != exact]
if (length(wrong)) {
  stop("Exact package-version mismatch: ",
       paste0(wrong, " expected ", exact[wrong], " found ", observed[wrong],
              collapse = "; "))
}

message("Core environment versions validated: R ", getRversion(),
        ", Bioconductor ", BiocManager::version(),
        ", Seurat ", observed[["Seurat"]], ", GSVA ", observed[["GSVA"]], ".")
