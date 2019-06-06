# En este script llevamos a cabo estudio de percolation in a lattice
# en estudiamos las activacion de los edges, la logica es la siguiente:
# Se selecciona un estado s cualquiera s=(0,1,0,1,1,....). Esto nos
# indica los vertex del lattice que quedaran activos y los que no.
# Para aquellos que quedaron activos, vemos si tienen vecinos activos.
# Para cada uno de los vecinos que esten activos, vemos si se conectan
# con probabilidad p.

# step 1:
# begin with an arbitrary configuration s. We associate to each pair of 
# neighbours a bond b_{i,k} being i and k part of the state s, a value 0 or 1. 
# If b_{i,k}=0, there is no link between the sites i and k; if b_{i,k}=1
# then the sites i and k are connected.
# These values are assigned according to the following conditional probabilities:
# P( b_{i,k}=0 | s_i != s_k ) = 1
# P( b_{i,k}=1 | s_i != s_k ) = 0
# P( b_{i,k}=0 | s_i = s_k ) = exp(-2*J_{ik}*\beta)
# P( b_{i,k}=1 | s_i = s_k ) = 1 - exp(-2*J_{ik}*\beta)
# Esto lo hacemos solo para los spins que estan en estado activo (s_i = 1),
# no para los que estan en activo 0 (a diferencia de lo que hace el algoritmo
# Swendsen-Wang)

# step 2:
# En el algoritmo de Swendsen-Wang, habria que identificar los spins que
# quedaron activos /e inactivos) y se les invierte el signo a los que 
# estan en un mismo clusters, se borran los edges, y comenzamos nuevamente al step 1.
# En nuestro caso, no hacemos eso, si no que repetimos el step 1 varias
# veces para el mismo estado y coleccionamos el tamano de los clusters <S>, 
# el numero de active bonds en los clusters, number of cluster per lattice, 
# the energies H of the clusters and
# los productos involucrados en el cluster mas grande.

# Creation date: 12-abr-19
# file name: percolation_bond_lattice_v3.R

# Notes:
# 12-abr-19: creation 
# 17-abr-19: repetimos el proceso muchas veces de estudiar los cluster que se forman
#           para un estado particular s.

# Recursos:
# https://www.youtube.com/watch?v=0SdzPJksV3Q
# https://rpubs.com/lgadar/generate-graphs






# # # # # # # # # # # # # # Carga de datos # # # # # # # # # # # # # #
rm(list = ls())
source("create_lattice_function.R") # ok
source("find_name_maxcluster_function.R")
source("energy_cluster_function.R")
source("correlation_length_function.R")
library(ggplot2)
library(igraph)
library(doMC)
library(purrr)

# 18-abr-19
# funcion para generar plots secuenciales
# function for creating file name with leading zeros
# makes it easier to process them sequentially
# https://davetang.org/muse/2015/02/12/animated-plots-using-r/
rename <- function(x){
  if (x < 10) {
    return(name <- paste('000',x,'plot.png',sep=''))
  }
  if (x < 100 && i >= 10) {
    return(name <- paste('00',x,'plot.png', sep=''))
  }
  if (x >= 100) {
    return(name <- paste('0', x,'plot.png', sep=''))
  }
}
# # # # # # # # # # # # # # Carga de datos # # # # # # # # # # # # # #


# # # # # # # # # ## # Parameters of the lattice # # # # # # # # # # #
# tamano del lattice:
size_lattice <- 20
# Para un lattice de N nodos, habran 2^N estados posibles y el lattice 
# tiene N^2 nodos.
number_of_edges <- 2*size_lattice*size_lattice-2*size_lattice 
# Distribucion ferromagnetica de acoples 
#J = rep(1, number_of_edges)
# Distribucion de los acoples: random [0,1]
#J <- runif(number_of_edges)
# Distribucion de los acoples: random [-1,1]
#J <- runif(number_of_edges)*2-1
#J[J < 0] <- 0
# custom discrete tistribution of J
J <- sample(x = c(0.1,2), number_of_edges, replace = T, prob = c(0.8, 0.2))
# Numero de veces de repeticion de la simulacion para cada estado
N = 50
# Random spin  generation 
# choosing a random spins to be active
set.seed(123)
#state <- as.numeric(runif(size_lattice*size_lattice)>0.5) # aprox la mitad de los vertices en on, la otra mitad en off
state <- rep(1, size_lattice*size_lattice) # todos los spins/vertices activados.
# # # # # # # # # ## # Parameters of the lattice # # # # # # # # # # #




# # # # # # # # # # # # Creation of the lattice # # # # # # # # # # #
# Lattice creation with L = size_lattice with weights
lat <- create_lattice(L=size_lattice)
xy <- lat[[1]]
edg <- lat[[2]]
g <- lat[[3]]

plot(g, layout=xy,
     #vertex.size=20,
     vertex.size=5,
     vertex.color = "red",
     #vertex.label.cex = 0.8,
     vertex.label="",
     edge.width=4)
# # # # # # # # # # # # Creation of the lattice # # # # # # # # # # #


# # # # # # #  # # # Lattice with active vertices # # # # # # # # # #
V(g)$state <- state
V(g)$color <- ifelse (V(g)$state, "red", "black") # black =off, red = on
E(g)$coupl <- J

plot(g, layout=xy,
     #vertex.size=20,
     vertex.size=5,
     vertex.color = V(g)$color,
     #vertex.label.cex = 0.8,
     vertex.label="",
     edge.width=4)
# # # # # # #  # # # Lattice with active vertices # # # # # # # # # #


# # # # # # # # # Lattice of connected active spins # # # # # # # # #
# Creation of the lattice with edges connected to all active spins
rw <- nrow(edg)
vec1 <- numeric(rw)
vec2 <- numeric(rw)
for (i in seq_along(vec1) ) {
  vec1[i] <- V(g)$state[which(V(g) == edg[i,1])]
  vec2[i] <- V(g)$state[which(V(g) == edg[i,2])]
}
link <- ifelse(vec1+vec2 == 2, 1, 0)
edg <- cbind(edg, link, J) # vertex1, vertex2, conextion 1=connected, 0 not connected
edg2 <- edg[edg[,"link"]==1,]
new_edg <- edg2[,c(1,2)]
new_g <- graph_from_edgelist(new_edg, directed = FALSE)
number_nodos_que_faltan <- size_lattice*size_lattice - max(new_edg)
if (number_nodos_que_faltan > 0) {
  new_g <- add_vertices(new_g, number_nodos_que_faltan )
} 
# solo si quiero colorear el strength del nodo
range1to100 <- function(x){1 + 99*(x-min(x))/(max(x)-min(x))}
palf <- colorRampPalette(c("lightgray", "dark red"))
colr <- palf(100);
# solo si quiero colorear el strength del nodo

E(new_g)$coupl <- edg2[,"J"]
V(new_g)$str <- strength(new_g, vids = V(new_g), weights = E(new_g)$coupl)
plot(new_g, layout=xy,
     #vertex.size=20,
     vertex.size=5,
     #vertex.color = V(g)$color,
     vertex.color = colr[round(range1to100(V(new_g)$str))], # solo si quiero colorear el strength del nodo
     #vertex.label.cex = 0.8,
     vertex.label="",
     edge.width=4)
# # # # # # # # # Lattice of connected active spins # # # # # # # # #


# NOTA: esto es lo que se REPITE N VECES.
# # # # # # # Lattice of connected active spins ramdom # # # # # # #
# Creation of the lattice with edges connected to  octive spins according
# probability p.
E(new_g)$p <- 1 - exp(-2*E(new_g)$coupl) # bond activation probability

# si p <= U, entonces hay conexion, si no, no hay conexion.
bond <- ifelse(E(new_g)$p >= runif(nrow(new_edg),0 ,1), 1, 0)
J <- E(new_g)$coupl
new_edg_ <- cbind(new_edg, bond, J) # vertex1, vertex2, conextion 1=connected, 0 not connected
edg3 <- new_edg_[new_edg_[,"bond"]==1,]
new_new_edg <- edg3[,c(1,2)]
new_new_g <- graph_from_edgelist(new_new_edg, directed = FALSE)
number_nodos_que_faltan <- size_lattice*size_lattice - max(new_new_edg)
if (number_nodos_que_faltan > 0) {
  new_new_g <- add_vertices(new_new_g, number_nodos_que_faltan )
} 

plot(new_new_g, layout=xy,
     #vertex.size=20,
     vertex.size=5,
     #vertex.color = V(g)$color,
     vertex.color = colr[round(range1to100(V(new_g)$str))], # solo si quiero colorear el strength del nodo
     #vertex.label.cex = 0.8,
     vertex.label="",
     edge.width=4)
E(new_new_g)$coupl <- edg3[,"J"]
# # # # # # # Lattice of connected active spins ramdom # # # # # # #





# # # # # # # # # # # # # # # Simulation # # # #  # # # # # # # # # #
# 17-abr-19
# repetimos el proces N veces, el proceso de formacion de clusters en 
# el lattice para un estado en particular.
almacen <- matrix(NA, ncol=6, nrow=N) #<S>, max(S), g(r), <k>, <E cluster>, <max E cluster>
vertex_names <-  vector("list", length = N)
cluster_sizes_list <- vector("list", length = N) # table with ns
for (it in 1:N) {
  cat("\r", "Processing data..Work done:", 100*it/N, "%")
  # Creation of the lattice with edges connected to  octive spins according
  # probability p.
  E(new_g)$p <- 1 - exp(-2*E(new_g)$coupl) # bond activation probability
  
  # si p <= U, entonces hay conexion, si no, no hay conexion.
  bond <- ifelse(E(new_g)$p >= runif(nrow(new_edg),0, 1), 1, 0)
  #print(bond)
  J <- E(new_g)$coupl
  new_edg_ <- cbind(new_edg, bond, J) # vertex1, vertex2, conextion 1=connected- 0 not connected, couplings
  edg3 <- new_edg_[new_edg_[,"bond"]==1, ]
  new_new_edg <- edg3[,c(1,2)]
  new_new_g <- graph_from_edgelist(new_new_edg, directed = FALSE)
  number_nodos_que_faltan <- size_lattice*size_lattice - max(new_new_edg)
  if (number_nodos_que_faltan > 0) {
    new_new_g <- add_vertices(new_new_g, number_nodos_que_faltan )
  } 
  E(new_new_g)$coupl <- edg3[,"J"]
  
  # almacenamiento de observables
  cluster_sizes <- clusters(new_new_g)
  #length correlation g(r=1)
  g_r <- correlation_length(new_new_g, L=size_lattice)
  # degree
  dg <- mean(degree(new_new_g))
  # para calcular la energia de todos los clusters
  clid <- 1:length(cluster_sizes$csize)
  all_energies <- map_dbl(clid, energy_cluster, members = cluster_sizes$membership)
  #drop energies equal to zero (cluster with size=1)
  all_energies <- all_energies[all_energies<0]
  
  almacen[it, ] <- c(mean(cluster_sizes$csize), max(cluster_sizes$csize), g_r, dg,
                     mean(all_energies), min(all_energies))
  vertex_names[[it]] <- find_name_maxcluster(cluster_sizes)
  cluster_sizes_list[[it]] <- table(cluster_sizes$csize)
  new_edg_ <- NULL
  
  #generacion de figuras para animar
  name <- rename(it)
  png(name)
  plot(new_new_g, layout=xy,
       #vertex.size=20,
       vertex.size=5,
       #vertex.color = V(g)$color,
       vertex.color = colr[round(range1to100(V(new_g)$str))], # solo si quiero colorear el strength del nodo
       #vertex.label.cex = 0.8,
       vertex.label="",
       edge.width=4,
       main = paste('Net lattice number= ', it, '<S>= ', round(mean(cluster_sizes$csize),2), '<E>= ', round(mean(all_energies),2)) )
  dev.off()
}

# analisis de los resultados:
almacen <- as.data.frame(almacen)
colnames(almacen) <- c("mean_s", "max_s", "gr", "k", "mean_E", "min_E")
colMeans(almacen)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
hist(almacen$mean_s,30)
hist(almacen$max_s,30)
hist(almacen$k,30)
hist(almacen$gr,30)
plot(almacen$mean_s, almacen$k)
plot(almacen$gr, almacen$k)
plot(almacen$mean_s, almacen$mean_E)
head(cluster_sizes_list)
head(vertex_names)
# en vertex_name habria que hacer analisis de cuantas veces aparece cada 
# producto en cada largest cluster.
# # # # # # # # # # # # # # # Simulation # # # #  # # # # # # # # # #
