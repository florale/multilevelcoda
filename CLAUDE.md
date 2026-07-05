# AGENTS.md

## Project type

This is an R package. Follow standard R package development workflows.

## Documentation rules

- Never edit files in `man/*.Rd` directly.
- Documentation source lives in Roxygen comments in `R/*.R`.
- To change documentation, edit the relevant `#'` Roxygen comments above the function definition.
- After documentation changes, run:

```r
  devtools::document()
```

## Vignette rules

CRAN rationale

This package may include precompiled vignettes because some examples are too slow or computationally expensive to run during CRAN checks.
The `.Rmd.orig` files contain the reproducible source workflow, while the generated `.Rmd` files are the CRAN-facing vignette inputs.

Before finishing work involving vignettes

When vignette source files are changed:

1. Edit the relevant `*.Rmd.orig` file.
2. Edit `precompile.R` only if the generation procedure needs to change.
3. Run the precompile script.
4. Confirm the generated `*.Rmd` changed only as a result of precompilation.
5. Do not hand-edit generated vignette outputs.

```md
## Generated vignette files are read-only

Treat generated vignette files as read-only.

Never edit `vignettes/*.Rmd` directly when there is a corresponding `vignettes/*.Rmd.orig`.

For example:

- `vignettes/foo.Rmd.orig` is editable source.
- `vignettes/foo.Rmd` is generated output and must not be edited manually.

All changes must be made in one of:

- the relevant `vignettes/*.Rmd.orig` file
- `precompile.R`

Then regenerate the `.Rmd` files using the precompilation workflow.

If asked to modify a generated `.Rmd` file, first check whether a matching `.Rmd.orig` exists. If it does, apply the change to the `.Rmd.orig` source and regenerate. Do not directly patch the `.Rmd`.