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

### Example

n = c(250, 250)
p = 4
Y1 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 3, ncol = p))
Y2 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 4, ncol = p))
Y = rbind(Y1, Y2)
C = c( rep(1,sum(n)), rep(2,sum(n)))
ans = mpg(Y, C)  
plotDiff(ans, type = "weight")
plotDiff(ans, type = "shift")

### Reference
To do.