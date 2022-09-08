## precompile vignettes as they are slow to run
## and depend on things like mc stan being available

## Must manually move image files from multilevelcoda/ to multilevelcoda/vignettes/ after knit
## note: do move them, don't copy, do not leave artifacts behind at top level pkg dir

library(knitr)
knit(
  "vignettes/substitution-model.Rmd.orig",
  "vignettes/substitution-model.Rmd")
knit(
  "vignettes/comp-outcome.Rmd.orig",
  "vignettes/comp-outcome.Rmd")
knit(
  "vignettes/comp-predictor.Rmd.orig",
  "vignettes/comp-predictor.Rmd")
