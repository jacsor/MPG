#' Plot differences across mixture models
#'
#'
#' @param object Object of class inheriting from "mpg".
#' @return A plot. 
#' @examples
#' 
#' n = c(250, 250)
#' p = 4
#' 
#' Y1 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 3, ncol = p))
#' Y2 = rbind( matrix( rnorm( n[1]*p), ncol = p), matrix( rnorm(n[2]*p) + 4, ncol = p))
#' Y = rbind(Y1, Y2)
#' C = c( rep(1,sum(n)), rep(2,sum(n)))
#' 
#' ans = mpg(Y, C)
#' plotDiff(ans, type = "weight")
#' plotDiff(ans, type = "shift")
plotDiff <- function( object, type = "weight", 
                      dim = c(1,2), group =  1:length(unique(object$data$C)),
                      xlim = range(object$data$Y[,dim[1]]), 
                      ylim = range(object$data$Y[,dim[2]]), ...  )
{
  if(class(object)!="MPG")
  {
    print("ERROR: Input's class should be MPG.")
    return(0)
  }
  colors = rainbow(150)[1:100]
  if(type == "weight")
    state = object$chain$Z > object$prior$K_0
  else if(type == "shift")
  {
    state = matrix(NA, nrow=nrow(object$chain$Z),  ncol=ncol(object$chain$Z))
    for( it in 1:nrow(state))
      state[it,] = ( object$chain$S[it, object$chain$Z[it,]] == 1)
  }
  else
  {
    print("ERROR: type is equal to 'weight' or 'shift'")
    return(0)
  }     
  
  
  size = ceiling(sqrt(length(group)))
  if( size * (size -1) == length(group) )
    mat = matrix( seq(1,  size*(size-1)) , ncol=size, byrow=TRUE  )
  else
    mat = matrix( seq(1,  size*size) , ncol=size, byrow=TRUE  ) 
  layout(mat = mat)
  
  if(is.null(colnames(object$data$Y)[dim[1]]))
    xlab = ""
  else
    xlab = colnames(object$data$Y)[dim[1]]
  
  if(is.null(colnames(object$data$Y)[dim[2]]))
    ylab = ""
  else
    ylab = colnames(object$data$Y)[dim[2]]
  
  for(j in group)
  {
    plot(object$data$Y[object$data$C==j,dim],
         col = colors[ round(99*colMeans(state[,object$data$C==j])) + 1 ], 
         pch = 16, xlim = xlim*c(1,1.1), ylim = ylim,
         xlab = xlab,
         ylab = ylab)
    par(new = TRUE)
    plot(rep(10,100), (1:100)/100, col = colors[1:100], bty='l', axes = FALSE,
         pch = 15, xlab = "", ylab = "", xaxt='n', cex = 2, 
         ylim = c(0,1), xlim = c(0,10))
    axis(side = 4, at=seq(0,1,by = .2), labels=seq(0,1,by = .2), las = 2)
    par(new = FALSE)
  }  
  
  layout(1)
  
  
}
