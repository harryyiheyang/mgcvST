#' MISO E13 spatial transcriptomics example
#'
#' MISO E13 covariates, three selected gene-expression vectors (Mapt, Map1b,
#' and Hist1h2ao), and the baseline and finer SPDE meshes. Rows of `expression`
#' align with rows of `covariates`. Use `spdePC_g999` with `pc_cutoff = 0.999`.
#'
#' @format A list with `covariates`, `expression`, and `meshes`.
#' @examples
#' data(MISO_E13)
#' head(MISO_E13$covariates)
#' head(MISO_E13$expression)
"MISO_E13"

#' Visium B spatial transcriptomics example
#'
#' Visium B covariates, three selected gene-expression vectors (mt-co3,
#' mt-co2, and BRAFhuman), and the baseline and finer SPDE meshes. Rows of
#' `expression` align with rows of `covariates`. Use `spdePC_g999` with
#' `pc_cutoff = 0.999`.
#'
#' @format A list with `covariates`, `expression`, and `meshes`.
#' @examples
#' data(Visium_B)
#' head(Visium_B$covariates)
#' head(Visium_B$expression)
"Visium_B"
