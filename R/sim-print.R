#' Internal Printing and Summary Helpers
#'
#' Build compact summaries for `mlsim_data` objects and provide small utility
#' conversions used by print methods.
#'
#' @param x An `mlsim_data`, summary, metadata, or table object.
#' @param name Generator name.
#' @param metadata Generator metadata list.
#' @param ... Additional arguments passed to table printing.
#' @param row.names Logical; passed to [print()].
#'
#' @return Internal helper outputs used by print and summary methods.
#'
#' @examples
#' multilevelcoda:::.mlsim_na_character(NULL)
#' sim <- simulate_data(n = 2, generators = list(x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1)))
#' multilevelcoda:::.mlsim_generated_columns(sim)
#'
#' @keywords internal
#' @name multilevelcoda-internal-printing
NULL

#' @rdname multilevelcoda-internal-printing
.mlsim_na_character <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1L]])) {
    return(NA_character_)
  }
  as.character(x[[1L]])
}

#' @rdname multilevelcoda-internal-printing
.mlsim_na_integer <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1L]])) {
    return(NA_integer_)
  }
  as.integer(x[[1L]])
}

#' @rdname multilevelcoda-internal-printing
.mlsim_has_metadata <- function(x) {
  !is.null(x) && length(x) > 0L
}

#' @rdname multilevelcoda-internal-printing
.mlsim_generated_columns <- function(x) {
  design_cols <- c(x$metadata$group_id, x$metadata$obs_id, x$metadata$time_id)
  design_cols <- design_cols[!is.na(design_cols)]
  setdiff(names(x$data), design_cols)
}

#' @rdname multilevelcoda-internal-printing
.mlsim_design_summary <- function(x) {
  metadata <- x$metadata
  group_sizes <- metadata$group_sizes
  grouped <- !is.null(metadata$n_groups)

  data.table::data.table(
    n = nrow(x$data),
    n_cols = ncol(x$data),
    n_generated_cols = length(.mlsim_generated_columns(x)),
    n_groups = if (grouped) as.integer(metadata$n_groups) else NA_integer_,
    group_id = if (grouped) .mlsim_na_character(metadata$group_id) else NA_character_,
    group_size_min = if (grouped) as.integer(min(group_sizes)) else NA_integer_,
    group_size_median = if (grouped) stats::median(group_sizes) else NA_real_,
    group_size_max = if (grouped) as.integer(max(group_sizes)) else NA_integer_,
    obs_id = .mlsim_na_character(metadata$obs_id),
    time_id = .mlsim_na_character(metadata$time_id),
    seed = .mlsim_na_integer(metadata$seed),
    n_generators = length(metadata$generator_order)
  )
}

#' @rdname multilevelcoda-internal-printing
.mlsim_generator_summary_row <- function(name, metadata) {
  vars <- metadata$vars %||% character()

  data.table::data.table(
    generator = name,
    distribution = .mlsim_na_character(metadata$distribution),
    level = .mlsim_na_character(metadata$level),
    vars = paste(vars, collapse = ", "),
    n_vars = length(vars),
    parameter_level = .mlsim_na_character(metadata$parameter_level),
    parameter_count = .mlsim_na_integer(metadata$parameter_count),
    has_row_parameters = .mlsim_has_metadata(metadata$row_parameters),
    has_group_parameters = .mlsim_has_metadata(metadata$group_parameters),
    has_fixed_parameters = .mlsim_has_metadata(metadata$fixed_parameters),
    has_random_cov = .mlsim_has_metadata(metadata$random_cov),
    has_random_effects = .mlsim_has_metadata(metadata$random_effects) ||
      .mlsim_has_metadata(metadata$random_intercepts) ||
      .mlsim_has_metadata(metadata$scale_random_intercepts) ||
      .mlsim_has_metadata(metadata$random_cov),
    has_residuals = .mlsim_has_metadata(metadata$residual_sd) ||
      .mlsim_has_metadata(metadata$residual_cov) ||
      .mlsim_has_metadata(metadata$residual_cor),
    has_scale_model = .mlsim_has_metadata(metadata$scale_fixed_intercept) ||
      .mlsim_has_metadata(metadata$scale_random_intercepts),
    has_composition = isTRUE(metadata$compositional) ||
      .mlsim_has_metadata(metadata$parts) ||
      .mlsim_has_metadata(metadata$sbp),
    has_custom_output = .mlsim_has_metadata(metadata$custom_output_level)
  )
}

#' @rdname multilevelcoda-internal-printing
.mlsim_generator_summary <- function(x) {
  rows <- Map(.mlsim_generator_summary_row, names(x$generator_metadata), x$generator_metadata)
  data.table::rbindlist(rows, use.names = TRUE)
}

#' @rdname multilevelcoda-internal-printing
.mlsim_print_table <- function(x, ..., row.names = FALSE) {
  print(x, ..., row.names = row.names)
}

#' Summarize Simulated Data
#'
#' Build a compact summary of an `mlsim_data` object.
#'
#' @param object An `mlsim_data` object returned by [simulate_data()].
#' @param ... Ignored.
#'
#' @return A `summary.mlsim_data` object with `design` and `generators` tables.
#'
#' @examples
#' sim <- simulate_data(n = 3, generators = list(x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1)))
#' summary(sim)
#'
#' @method summary mlsim_data
#' @export
summary.mlsim_data <- function(object, ...) {
  structure(
    list(
      design = .mlsim_design_summary(object),
      generators = .mlsim_generator_summary(object)
    ),
    class = "summary.mlsim_data"
  )
}

#' Print Simulated Data
#'
#' Print a compact overview of an `mlsim_data` object.
#'
#' @param x An `mlsim_data` object returned by [simulate_data()].
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' sim <- simulate_data(n = 3, generators = list(x = gen_mvn("x", fixed_intercept = 0, residual_cov = 1)))
#' print(sim)
#'
#' @method print mlsim_data
#' @export
print.mlsim_data <- function(x, ...) {
  x_summary <- summary(x)
  design <- x_summary$design

  cat("<mlsim_data>\n")
  cat("  rows: ", design$n, "\n", sep = "")
  cat("  columns: ", design$n_cols, " (generated: ", design$n_generated_cols, ")\n", sep = "")
  cat("  seed: ", if (is.na(design$seed)) "<none>" else design$seed, "\n", sep = "")
  if (is.na(design$n_groups)) {
    cat("  grouping: none\n")
  } else {
    cat(
      "  grouping: ", design$group_id, " (", design$n_groups,
      " groups; size min/median/max: ", design$group_size_min, "/",
      design$group_size_median, "/", design$group_size_max, ")\n",
      sep = ""
    )
  }
  if (!is.na(design$time_id)) {
    cat("  time: ", design$time_id, "\n", sep = "")
  }
  cat("  generators: ", design$n_generators, "\n", sep = "")

  cat("\nGenerators:\n")
  generator_cols <- c("generator", "distribution", "level", "vars")
  generator_table <- as.data.frame(x_summary$generators[, generator_cols, with = FALSE])
  print(generator_table, row.names = FALSE)

  invisible(x)
}

#' Print a Simulated Data Summary
#'
#' @param x A `summary.mlsim_data` object.
#' @param ... Additional arguments passed to table printing.
#'
#' @return `x`, invisibly.
#'
#' @method print summary.mlsim_data
#' @export
print.summary.mlsim_data <- function(x, ...) {
  cat("<summary.mlsim_data>\n")
  cat("\nDesign:\n")
  .mlsim_print_table(x$design, ...)
  cat("\nGenerators:\n")
  .mlsim_print_table(x$generators, ...)
  invisible(x)
}
