#' Multilevel Compositional Data
#'
#' A simulated dataset containing multiple days of compositional data.
#'
#' @format A data table containing 9 variables.
#' \describe{
#'   \item{TST}{Total Sleep Time (minutes) --- repeated measure}
#'   \item{WAKE}{Wake time while in bed, trying to sleep (minutes) --- repeated measure}
#'   \item{MVPA}{Moderate to Vigorous Physical Activity (minutes) --- repeated measure}
#'   \item{LPA}{Light Physical Activity (minutes) --- repeated measure}
#'   \item{SB}{Sedentary Behavior (minutes) --- repeated measure}
#'   \item{ID}{A unique identifier for each individual}
#'   \item{Age}{Age in years --- baseline measure only}
#'   \item{Female}{Binary: whether participants identified as female (1) or not (0) --- baseline measure only}
#'   \item{STRESS}{Self report stress measures on a 0 to 10 scale --- repeated measure}
#' }
"mcompd"

#' Possible Pairwise Substitutions
#'
#' A dataset containing possible pairwise subsitutions.
#'
#' @format A data table containing 5 variables.
#' \describe{
#'   \item{V1}{first compositional variable}
#'   \item{V2}{second compositional variable}
#'   \item{V3}{third compositional variable}
#'   \item{V4}{fourth compositional variable}
#'   \item{V5}{fifth compositional variable}
#' }
"psub"

#' Sequential Binary Partition
#'
#' A matrix of sequential binary partition.
#'
#' @format A matrix with 5 columns and 4 rows.
#' \describe{
#'   \item{V1}{first compositional variable}
#'   \item{V2}{second compositional variable}
#'   \item{V3}{third compositional variable}
#'   \item{V4}{fourth compositional variable}
#'   \item{V5}{fifth compositional variable}
#' }
"sbp"