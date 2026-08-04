#' Visium human DLPFC data
#'
#' An example dataset
#'
#' @format ## `Visium_human_DLPFC`
#' A list with 3 elements
#' \describe{
#'   \item{counts}{Gene expression count matrix. Rows and columns indicate genes and spots, respectively.}
#'   \item{spatial_coords}{Spatial coordinates matrix. Rows indicate spots.}
#'   \item{cell_type_proportions}{Cell type proportion matrix estimated by RCTD. Rows and columns indicate spots and cell types, respectively.}
#' }
#' @source Visium human dorsolateral prefrontal cortex (DLPFC) sample 151673.
#'   The count matrix and spatial coordinates were obtained via the
#'   \pkg{spatialLIBD} resource (\url{https://research.libd.org/spatialLIBD/})
#'   and were originally generated and described in Maynard KR,
#'   Collado-Torres L, Weber LM, et al. (2021) Transcriptome-scale spatial gene
#'   expression in the human dorsolateral prefrontal cortex.
#'   \emph{Nature Neuroscience} \strong{24}, 425--436.
#'   \doi{10.1038/s41593-020-00787-0}. Cell type proportions were estimated by
#'   the authors of this package using RCTD with a reference dataset (\url{https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/human-frontal-cortex-ffpe-dataset/}) as the
#'   single-cell reference.
#'
#'   This dataset is redistributed here for illustration only. It is not
#'   covered by the package license and remains subject to the terms of the
#'   original data providers. Please cite Maynard et al. (2021) when using it.
"Visium_human_DLPFC"
