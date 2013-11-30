check_data <- function(data){
  dat <- data
  cols <- names(dat)
  
  if (!("x" %in% cols)){
    stop("data must have a column 'x'.")
  }
  if (!("u" %in% cols)){
    stop("data must have a column 'u'.")
  }
  if (!("P" %in% cols)){
    stop("data must have a column 'P'.")
  }
  
  if (!("u_P" %in% cols)){
    dat$u_P <- (dat$u / dat$P)
  }
  
  if (sum(dat$u >= 0)){
    stop("'u' should be negative.")
  }
  
  return(dat)
}