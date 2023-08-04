## precompile vignettes as they are slow to run
## and depend on things like mc stan being available

## Must manually move image files from multilevelcoda/ to multilevelcoda/vignettes/ after knit
## note: do move them, don't copy, do not leave artifacts behind at top level pkg dir

library(knitr)
knit(
  "vignettes/D-substitution-model.Rmd.orig",
  "vignettes/D-substitution-model.Rmd")
knit(
  "vignettes/C-composition-MMLM.Rmd.orig",
  "vignettes/C-composition-MMLM.Rmd")
knit(
  "vignettes/B-composition-MLM.Rmd.orig",
  "vignettes/B-composition-MLM.Rmd")
