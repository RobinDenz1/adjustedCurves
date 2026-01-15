# Print Method for `curve_test` Objects

Prints some important parts of the output object.

## Usage

``` r
# S3 method for class 'curve_test'
print(x, digits=4, ...)
```

## Arguments

- x:

  An object of class `curve_test` created by the `adjusted_curve_test`
  function.

- digits:

  How many digits to use when rounding the results.

- ...:

  Currently not used.

## Value

Silently returns the `data.frame` which can be seen when calling the
function. `ABC` is an abbreviation for "area between the curves".

## Author

Robin Denz

## See also

[`adjusted_curve_test`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_curve_test.md),
[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md),
[`adjustedcif`](https://robindenz1.github.io/adjustedCurves/reference/adjustedcif.md)

## Examples

``` r
# See ?adjusted_curve_diff
```
