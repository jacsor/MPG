Variability decomposition across related mixture distributions
==============================================================

This package fits MPG, an algorithm to compare multiple related mixture distributions.

### Install
The package can be installed on Linux and Mac using `devtools`:

```S
install.packages('devtools')
library('devtools')
devtools::install_github('MPG', 'jacsor')
```

### Use
There are three functions in this package, and their descriptions are provided in the help files

```S
ans = fit(Y, C)
plotDiff(ans)
cal = calibrate(ans)
```

### Reference
To do.