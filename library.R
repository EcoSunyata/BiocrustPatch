library("tidyverse") # tidyverse tools
library("imagine") # for convolution2D 
library("EBImage") # image processing 
library("raster")  # for plotting

# define neighborhood kernel
make_kernel <- function(r) {
  n <- 2 * r + 1 # matrix size
  matrix <- matrix(0, n, n) # make zero matrix
  # fill the circle
  for (i in 1:n) {
    for (j in 1:n) {
      distance <- sqrt((i - r - 1)^2 + (j - r - 1)^2)
      if (distance <= r) {
        matrix[i, j] <- 1
      }
    }
  }
  return(matrix)
}


# define convolution function with different boundaries
myconvolution2D <- function(map,template,boundary="period"){
  row_ext <- ((dim(template)-1)/2)[1]
  col_ext <- ((dim(template)-1)/2)[2]
  row_x <- dim(map)[1]
  col_x <- dim(map)[2]
  if(boundary=="vanishing"){
    map_ext <- matrix(0,nrow=2*row_ext+row_x,ncol = 2*col_ext+col_x)
    map_ext[(1+row_ext):(row_x+row_ext),(1+col_ext):(col_x+col_ext)] <- map
    map_ext <- convolution2D(X = map_ext,kernel = template)
    row_z <- dim(map_ext)[1]
    col_z <- dim(map_ext)[2]
    return(map_ext[(1+row_ext):(row_z-row_ext),(1+col_ext):(col_z-col_ext)])
    }
  if(boundary=="period"){
    map <- map[c(((row_x-row_ext+1):row_x), (1:row_x), (1:row_ext)),]
    map <- map[,c(((col_x-col_ext+1):col_x), (1:col_x), (1:col_ext))]
    map_ext <- convolution2D(X = map,kernel = template)
    row_z <- dim(map_ext)[1]
    col_z <- dim(map_ext)[2]
    return(map_ext[(1+row_ext):(row_z-row_ext),(1+col_ext):(col_z-col_ext)])
  }
}


mydilate <- function(map,template,boundary="period"){
  row_ext <- ((dim(template)-1)/2)[1]
  col_ext <- ((dim(template)-1)/2)[2]
  row_x <- dim(map)[1]
  col_x <- dim(map)[2]
  if(boundary=="vanishing"){return(dilate(map, template))}
  if(boundary=="period"){
    map <- map[c(((row_x-row_ext+1):row_x), (1:row_x), (1:row_ext)),]
    map <- map[,c(((col_x-col_ext+1):col_x), (1:col_x), (1:col_ext))]
    map <- dilate(map, template)[(row_ext+1):(row_x+row_ext),
                                 (col_ext+1):(col_x+col_ext)]
    return(map)
  }
}


