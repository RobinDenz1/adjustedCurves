# Plot Method for `curve_test` Objects

Produces either a spaghetti-plot of the bootstrapped difference curves
(`type="curves"`) or a kernel-density plot of the shifted bootstrap
distribution of the difference curve integrals (`type="integral"`).

## Usage

``` r
# S3 method for class 'curve_test'
plot(x, type="curves", xlab=NULL,
     ylab=NULL, title=NULL, ...)
```

## Arguments

- x:

  An object of class `curve_test` created by the `adjusted_curve_test`
  function.

- type:

  Either `"curves"` or `"integral"`, specifying what should be plotted.

- xlab:

  The label of the X-Axis. Set to `NULL` to use default label.

- ylab:

  The label of the Y-Axis. Set to `NULL` to use default label.

- title:

  The title of the plot. Set to `NULL` to use no title.

- ...:

  Currently not used.

## Details

When using `type="curves"` the black curve shows the observed curve of
the difference. When using `type="integral"` the red line shows the
observed integral of the curve of the difference.

Both graphics can be used to check if the assumptions of the test hold.
The bootstrap-shifted distribution of the integral of the difference
should approximately be normally distributed. If the kernel-density
estimate shown with `type="integral"` is clearly not normally
distributed, the estimated p-value might be wrong. Similarly, if the
curves of the differences do not vary randomly around the black line
when using `type="curves"`, the estimated p-value might be wrong. You
could also try to rerun the `adjustedsurv` or `adjustedcif` function
with a bigger number in `n_boot`.

## Value

Returns a `ggplot2` object.

## Author

Robin Denz

## See also

[`adjusted_curve_test`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_curve_test.md),
[`adjustedsurv`](https://robindenz1.github.io/adjustedCurves/reference/adjustedsurv.md),
[`adjusted_rmst`](https://robindenz1.github.io/adjustedCurves/reference/adjusted_rmst.md)

## Examples

``` r
# See ?adjusted_curve_test
```
