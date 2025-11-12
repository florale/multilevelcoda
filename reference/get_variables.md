# Extract variable names from an object

Generic function to extract variable names from a supported object.

## Usage

``` r
get_variables(object)

# S3 method for class 'complr'
get_variables(object)

# S3 method for class 'brmcoda'
get_variables(object)
```

## Arguments

- object:

  A `brmcoda` object

## Value

A list of variable names.

## Examples

``` r
# For a complr object:
# get_variables(complr_object)

# For a brmcoda object:
# get_variables(brmcoda_object)
```
