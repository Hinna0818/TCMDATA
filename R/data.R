#' Combined Human and Mouse TF-Target Interactions
#'
#' A comprehensive dataset containing high-confidence (levels A, B, and C) transcription factor (TF)
#' and target gene interactions for both Human and Mouse, derived from the DoRothEA database.
#'
#' @format A tibble with `r nrow(tf_targets)` rows and 5 columns:
#' \describe{
#'   \item{Species}{Character. The organism ("Human" or "Mouse").}
#'   \item{TF}{Character. Gene symbol of the transcription factor.}
#'   \item{Target}{Character. Gene symbol of the target gene.}
#'   \item{Confidence}{Character. Confidence level of the interaction source:
#'     \itemize{
#'       \item \strong{A}: High confidence (Literature-curated).
#'       \item \strong{B}: Moderate confidence (ChIP-seq evidence).
#'       \item \strong{C}: Medium confidence (TFBS predictions + ChIP-seq).
#'     }
#'   }
#'   \item{Mode_of_Regulation}{Numeric. Indicates if the TF activates or inhibits the target:
#'     \itemize{
#'       \item \strong{1}: Activation
#'       \item \strong{-1}: Inhibition
#'     }
#'   }
#' }
#' @source \url{https://saezlab.github.io/dorothea/}
#' @references
#' Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J.
#' Benchmark and integration of resources for the estimation of human transcription factor activities.
#' Genome Research. 2019. DOI: 10.1101/gr.240663.118.
"tf_targets"