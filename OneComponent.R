# Single component self-organization
# Author: Jingyao Sun
# Date: 2025-07-05
# Description: 
# Biocrust patches evolves within a 2-D lattice, with small patches gradually increasing in size. 
# Eventually, due to resource constraints, a regular pattern emerges in the later stages of development.
# --------------------------------------------------------------------- #
source("library.R")
# Define the iteration function
next_step <- function(map){
  e_moss <- map$moss # the place where moss live 
  p_moss <- myconvolution2D(e_moss,tml_nei) # the number of neighbor in local range
  k_moss <- ifelse(p_moss < K, T,F) # check if the number of neighbor is below the carrying capacity
  n_moss <- myconvolution2D(e_moss,tml_moss) # the number of neighbor in adjacent range
  
  expand <- dilate(e_moss, tml_moss) & k_moss & (!e_moss) # the place around the existing patches
  propa_p <- c(n_moss[expand]) ^2 * 0.0007
  expand[expand]<- sapply(propa_p, FUN = function(p){sample(c(1,0),1, prob=c(p,1-p),replace=T)})
  
  shrink <- !k_moss & e_moss
  shrink[shrink] <- sample(c(0,1),sum(shrink), prob=c(0.02,0.98),replace=T)
  
  map$moss <- ifelse(expand|shrink|(k_moss & e_moss), 1,0)
  disperse <- ifelse(sample(c(1,0), len^2, prob=c(0.00002,1-0.00002),replace=T),T,F)
  map$moss[disperse] <- 1
  
  return(map)
}

# Define the initial lattice
initialization <- function(){
  distri <- matrix(sample(c(1,0,-1),len*len,
                          prob=c(0.00005,1-2*0.00005,0.00005),replace=T),len,len)
  moss <- matrix(0,len,len)
  moss[distri==1] <- 1

  resource <- matrix(5,len,len)
  tml_pro <<- makeBrush(3, shape='box')
  tml_res <<- makeBrush(9, shape='disc')
  tml_res <<- tml_res/sum(tml_res)
  tml_moss <<- makeBrush(7, shape='disc')
  
  map <- list(moss=moss,resource=resource)
  return(map)
}
# --------------------------------------------------------------------- #
# Simulation
len <- 400
map <- initialization()
tml_nei <- makeBrush(101, shape='disc')
K <- 2300
map_list <- list()
map_list[[1]] <- map

#set.seed(2222)
for(i in 1:300){
  map_list[[i+1]] <- next_step( map_list[[i]])
  print(paste(i))
}

# --------------------------------------------------------------------- #
# Plot the change of density
density <- sapply(map_list,function(x){sum(x$moss)})/len^2
fig_data <- data.frame(step=0:300,density=density[1:301])
fig_data$mark <- "black"
fig_data$mark[fig_data$step %in% c(63,103,145,300)] <- "red"

ggplot(fig_data)+geom_point(aes(x=step,y=density,color=mark,size=mark))+
  scale_color_manual(values = c("grey30","#FD7446FF"))+
  scale_size_manual(values=c(1.5,3))+
  xlab("Step")+ylab("Density")+
  theme_bw()+
  theme(legend.position = "none")


# Plot the snapshot
temp <- map_list[[150]]$moss %>% as.matrix()
test <- temp %>% raster() %>% as.data.frame(xy=T) %>% mutate(layer=factor(layer))
ggplot(test)+geom_tile(aes(x=x,y=y,fill=layer))+
  coord_cartesian(expand = F)+theme_bw()+
  theme(legend.position = "none",axis.title = element_blank(),
        plot.margin = margin(t = 5, r = 10, b = 0, l = 0,unit = "pt"))+
  scale_fill_manual(values=c("white","#00BFC4"))




