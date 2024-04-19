calc_suff_stat <- function(x, n, na.rm = TRUE){
  if(na.rm){
    xp    <- x[stats::complete.cases(x)]
    z     <- xp/n - 1/(2*n)
    odds  <- z/(1 - z)

    mean(log(odds))
  } else{
    z     <- x/n - 1/(2*n)
    odds  <- z/(1 - z)

    mean(log(odds))
  }
}

rank_array  <- function(arr, ...){
  rows  <- dim(arr)[1]
  cols  <- dim(arr)[2]
  n     <- dim(arr)[3]

  rankArray   <- array(NA, dim = c(rows, cols, n))
  for(r in 1:rows){
    for(l in 1:cols){
      rankArray[r,l,]   <- rank(arr[r,l,], ...)
    }
  }

  return(rankArray)

}
