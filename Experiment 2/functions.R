circleFun <- function(center = c(0,0), r = 2, npoints = 100){
  t <- seq(0,2*pi,length.out = npoints)
  x <- center[1] + r * cos(t)
  y <- center[2] + r * sin(t)
  return(data.frame(x = x, y = y))
}

readData <- function(filename) {
  dat <- read.csv(paste0('data/', filename)) 
  return(dat)
}