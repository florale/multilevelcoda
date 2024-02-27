#' Multilevel Compositional Data
#'
#' A simulated dataset containing multiple days of compositional data.
#'
#' @format A data table containing 10 variables.
#' \describe{
#'   \item{ID}{A unique identifier for each individual}
#'   \item{Time}{Recurrence time of repeated measures by individual}
#'   \item{Stress}{Self report stress measures on a 0 to 10 scale --- repeated measure}
#'   \item{TST}{Total Sleep Time (minutes) --- repeated measure}
#'   \item{WAKE}{Wake time while in bed, trying to sleep (minutes) --- repeated measure}
#'   \item{MVPA}{Moderate to Vigorous Physical Activity (minutes) --- repeated measure}
#'   \item{LPA}{Light Physical Activity (minutes) --- repeated measure}
#'   \item{SB}{Sedentary Behavior (minutes) --- repeated measure}
#'   \item{Age}{Age in years --- baseline measure only}
#'   \item{Female}{Binary: whether participants identified as female (1) or not (0) --- baseline measure only}
#' }
"mcompd"

#' Possible Pairwise Substitutions
#'
#' A dataset containing possible pairwise subsitutions.
#'
#' @format A data table containing 5 variables.
#' \describe{
#'   \item{TST}{first compositional variable}
#'   \item{WAKE}{second compositional variable}
#'   \item{MVPA}{third compositional variable}
#'   \item{LPA}{fourth compositional variable}
#'   \item{SB}{fifth compositional variable}
#' }
"psub"

#' Sequential Binary Partition
#'
#' A matrix of sequential binary partition.
#'
#' @format A matrix with 5 columns and 4 rows.
#' \describe{
#'   \item{TST}{first compositional variable}
#'   \item{WAKE}{second compositional variable}
#'   \item{MVPA}{third compositional variable}
#'   \item{LPA}{fourth compositional variable}
#'   \item{SB}{fifth compositional variable}
#' }
"sbp"

#' multilevelcoda Simulation Study results
#'
#' A list of 4 components
#'
#' @format A list with 5 columns and 4 rows.
#' \describe{
#'   \item{brmcoda_tab}{Simulation results for brmcoda() for tables}
#'   \item{sub_tab}{Simulation results for substitution() for tables}
#'   \item{brmcoda_plot}{Simulation results for brmcoda() for graphs}
#'   \item{sub_plot}{Simulation results for substitution() for graphs}
#' }
"sim"