## make Rcmd check happy
utils::globalVariables(c("i",  "..cols", ".", "To", ".SD", "t", "head",  "fitted",
                         "x", "object", "predict", "residuals", "tail", "vcov", "coef",
                         "Mean",  "CI_low", "CI_high", "From", "Delta", "pairs",
                         "spread", "value", "variable", "ID", "EffectType", "Level", "Reference",
                         "update"))


# expand grid data frame
expand.grid.df <- function(...) Reduce(function(...) merge.data.frame(..., by = NULL, all = TRUE), list(...))

# check sequence of number
is.sequential <- function(x) {
  all(length(x) > 2 & all(abs(diff(x)) == 1))
}
