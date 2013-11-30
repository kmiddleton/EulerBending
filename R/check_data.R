check_data <- function(data){
  dat <- data
  cols <- names(dat)
  
  flag <- TRUE
  
  if (!("x" %in% cols)){
    stop("data must have a column 'x'.")
    flag <- FALSE
  }
  if (!("u" %in% cols)){
    stop("data must have a column 'u'.")
    flag <- FALSE
  }
  if (!("P" %in% cols)){
    stop("data must have a column 'P'.")
    flag <- FALSE
  }
  return(flag)
}