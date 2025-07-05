# Self-organization of moss and lichen
# Author: Jingyao Sun
# Date: 2025-07-05
# Description: 
# Moss and lichen patches evolves within a 2-D lattice, with small patches gradually increasing in size. 
# Eventually, due to resource constraints, a regular pattern emerges in the later stages of development.
# --------------------------------------------------------------------- #
source("library.R")
# Define the iteration function
next_step <- function(map){
  e_moss <- map$moss
  e_lichen <- map$lichen
  
  # moss growth
  p_moss <- myconvolution2D(e_moss,tml_nei_m) # Number of neighbors in local range
  k_moss <- ifelse(p_moss < k_m, T,F) # Check if the local neighbors is below the carrying capacity
  n_moss <- myconvolution2D(e_moss,tml_moss,boundary="vanishing") # Number of neighbors in adjacent range
  
  expand_m <- mydilate(e_moss, tml_moss,boundary="vanishing") & 
    k_moss & (!e_moss) & (!e_lichen) # The place around the existing patches are potential to be occupied
  propa_p_m <- c(n_moss[expand_m]) *0.01 + 0.05  # Expanding probability
  expand_m[expand_m]<- sapply(propa_p_m, FUN = function(p){sample(c(T,F),1, prob=c(p,1-p),replace=T)})
  expand_m <- ifelse(expand_m>=1,T,F)
  
  shrink_m <- !k_moss & e_moss # The place where mosses are overcrowded
  shrink_m[shrink_m] <- sample(c(F,T),sum(shrink_m), prob=c(0.05,0.95),replace=T) 
  
  map$moss <- ifelse(expand_m|shrink_m|(k_moss & e_moss), 1,0)
  disperse_m <- ifelse(sample(c(T,F), len^2, prob=c(0.0000005,1-0.0000005),replace=T),T,F) # Patch dispersal
  map$moss[disperse_m] <- 1
  
  #lichen growth
  p_lichen <- myconvolution2D(e_lichen,tml_nei_l) # Number of neighbors in local range
  k_lichen <- ifelse(p_lichen < k_l, T,F) # # Check if the local neighbors is below the carrying capacity
  n_lichen <- myconvolution2D(e_lichen,tml_lichen,boundary="vanishing") # Number of neighbors in adjacent range
  
  expand_l <- mydilate(e_lichen, tml_lichen,boundary="vanishing") & 
    k_lichen  & (!e_lichen) & (!e_moss) # The place around the existing patches are potential to be occupied
  propa_p_l <- c(n_lichen[expand_l]) *0.01 # Expanding probability
  expand_l[expand_l]<- sapply(propa_p_l, FUN = function(p){sample(c(T,F),1, prob=c(p,1-p),replace=T)})
  expand_l <- ifelse(expand_l>=1,T,F)
  
  shrink_l <- !k_lichen & e_lichen
  shrink_l[shrink_l] <- sample(c(F,T),sum(shrink_l), prob=c(0.05,0.95),replace=T)
  
  map$lichen <- ifelse(expand_l|shrink_l|(k_lichen & e_lichen), 1,0)
  disperse_l <- ifelse(sample(c(T,F), len^2, prob=c(0.0000025,1-0.0000025),replace=T),T,F) # Patch dispersal
  map$lichen[disperse_l] <- 1
  
  #competition
  cross <- map$moss & map$lichen
  MoL   <- sample(c(1,0), sum(cross), prob=c(0.9,0.1),replace=T)
  map$moss[cross] <- MoL
  map$lichen[cross] <- (1-MoL)
  
  return(map)
}

# Define the initial lattice
initialization <- function(){
  distri <- matrix(sample(c(1,0,-1),len*len,
                          prob=c(0.0000005,1-0.0000025-0.0000005,0.0000025),replace=T),len,len)
  moss <- matrix(0,len,len)
  moss[distri==1] <- 1
  lichen <- matrix(0,len,len)
  lichen[distri==-1] <- 1
  resource <- matrix(5,len,len)
  map <- list(moss=moss,lichen=lichen,resource=resource)
  return(map)
}
# --------------------------------------------------------------------- #
# Simulation

len <- 400
map <- initialization()
tml_nei_m  <- makeBrush(101, shape='disc')
tml_nei_l  <- makeBrush(71, shape='disc')
tml_moss   <- makeBrush(7, shape='disc') # propagation
tml_lichen <- makeBrush(7, shape='disc') # propagation
map_list <- list()
map_list[[1]] <- map


set.seed(3234)
for(i in 1:300){
  k_m<- 800+1500/(1+exp(-0.1*(i-100)))
  if(i<100){
    k_l<- (500+ 500*exp(-(i-100)^2/(2*30^2)) )}
  if(i>=100){
    k_l<- (200+ 800*exp(-(i-100)^2/(2*30^2)) )}
  map_list[[i+1]] <- next_step( map_list[[i]])
  print(paste(i))
}

# The change of carrying capacity for moss (blue) and lichen (red)
i1 <- 1:100
i2 <- 101:300
y1 <- ((1/10)+(1/5))/(1+exp(-0.1*(c(i1,i2)-100))) 
y2 <- ((1/8)+ (1/8)*exp(-(i1-100)^2/(2*30^2)) ) 
y3 <- ((1/20)+ (1/5)*exp(-(i2-100)^2/(2*30^2)) ) 
plot(c(i1,i2),y1,type="l",lty=1,ylim=c(0,0.35),col="#00BFC4",xlab="Simulation steps",ylab="Density")
lines(c(i1,i2),c(y2,y3),lty=1,col="#F8766D")

sapply(map_list,function(x){sum(x$moss)}) %>%   `/`(len^2) %>% plot(.,col="#00BFC4",lwd=3,type="l")
sapply(map_list,function(x){sum(x$lichen)}) %>% `/`(len^2) %>% lines(.,col="#F8766D",lwd=3)



# Plot the snapshot
step <- 180
snapshot <- map_list[[step]]$moss*2 + map_list[[step]]$lichen %>% as.matrix() 
snapshot <- snapshot %>% raster() %>% as.data.frame(xy=T) %>% mutate(layer=factor(layer))

ggplot(snapshot) + geom_tile(aes(x=x*len,y=y*len,fill=layer))+
  coord_cartesian(expand = F)+theme_bw()+
  theme(legend.position = "none",axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_fill_manual(values=c("white","#F8766D","#00BFC4"))



