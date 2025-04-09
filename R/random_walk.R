adjCells <- function(r, xy){
  #r: vector - resolution of raster c(_,_)
  #xy: vector - xy location of center cell c(_,_)
  
  d <- c(rep(xy[1]-r[1], 3), rep(xy[1]+r[1],3), xy[1], xy[1],
         rep(c(xy[2]+r[2], xy[2], xy[2]-r[2]), 2), xy[2]+r[2], xy[2]-r[2])
  
  d <- matrix(d, ncol=2)

  
}

adjRowCol <- function(rowcol){
  
  x <- .subset2(rowcol, 1)
  y <- .subset2(rowcol, 2)
  
  d <- c(x - 1, y -1,
         x, y -1,
         x + 1, y -1,
         
         x - 1, y + 1,
         x, y + 1,
         x + 1, y  + 1,
         
         x - 1, y,
         x + 1, y)
  
  d <- matrix(d, ncol = 2, byrow = TRUE)
  
  return(d)
}



