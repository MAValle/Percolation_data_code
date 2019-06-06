# Es exactamente lo mismo que en percolation_bond_lattice_v4.R, pero simulo 
# con J positivos y negativos. Para aquellos bonds en los que J<0 hacemos que 
# J=0 (es decir, sin conexión) pero sólo para efectos de la simulación. No 
# obstante, la información de los bonds o edges que tienen J<0 se preserva 
# para calcular la energía de los clusters tomando en cuenta estos edges 
# con acoples negativos. Esto ocurrirá cuando queden dos nodos que son vecinos  
# y que sean parte de un mismo cluster (pero no conectaods entre sí), posean 
# acople negativo. 

# Procedure: 
# step 1: Create a lattice network of LXL with L^2 vertices.
# step 2: Assign to each edge of the lattice a distribution of J en [-1,1] (this is for energy computing only)
# step 3: Assign to each edge of the lattice the coupling J. If J<0, then J=0.
# step 4: identify the strength of each vertice on the lattice 
# step 5: start simulation connecting each neigbourhood vertex with probability p = 1 - exp(_2J)
#         with all vertex actives (s=1,1,1,,,,,,1)
# step 6: measure observables (for each cluster, their energy (with and without J<0), and size)
# step 7: drop or leave inactive one of the vertex of the lattice according to some rule and go to step 4.
# step 8: repeat step 4 until some proportion p of the vertex have been removed.
# step 9: repeat again from step 4 N times.


# Creation date: 29-abr-19
# file name: percolation_bond_lattice_v4a.R

# Notes:
# 29-abr-19: creation 
# 02-may-19: agregamos un script al final para simular conexion del lattice sola 1 vez, 
#           sin quitar nodos con mayor/menor strength , solo para analizar las energias.



# Recursos:
# https://www.youtube.com/watch?v=0SdzPJksV3Q
# https://rpubs.com/lgadar/generate-graphs






# # # # # # # # # # # # # # Carga de datos # # # # # # # # # # # # # #
rm(list = ls())
source("create_lattice_function.R") # ok
source("find_name_maxcluster_function.R")
source("energy_cluster_function.R")
source("correlation_length_function.R")
source("create_new_lattice_bonds_function.R")
library(ggplot2)
library(igraph)
library(doMC)
library(purrr)
library(latex2exp)

# 18-abr-19
# funcion para generar plots secuenciales
# function for creating file name with leading zeros
# makes it easier to process them sequentially
# https://davetang.org/muse/2015/02/12/animated-plots-using-r/
rename <- function(x){
  if (x < 10) {
    return(name <- paste('000',x,'plot.png',sep=''))
  }
  if (x < 100 && x >= 10) {
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
set.seed(123)
# Distribucion ferromagnetica de acoples 
#J0 <- Ji <- rep(1, number_of_edges)
# Distribucion de los acoples: random [0,1]
#J0 <- Ji <- runif(number_of_edges)
# Distribucion de los acoples: random [-1,1]
#J0 <- Ji <- runif(number_of_edges)*2-1   # J0 incluye los J < 0, mientras que Ji no los oncluye, es solo para la red.
#Ji[Ji < 0] <- 0
# custom discrete tistribution of J
#J0 <- Ji <- sample(x = c(0.1,2), number_of_edges, replace = T, prob = c(0.8, 0.2))
# Normal distribution of J truncated between 0 and 2
#library(truncnorm)
#J0 <- Ji <- rtruncnorm(number_of_edges, a=0, b=2, mean=0.5, sd=0.5)
# Normal distribution of J truncated between -1 and 2
library(truncnorm)
J0 <- Ji <- rtruncnorm(number_of_edges, a=-1, b=1, mean=0, sd=0.5)
# Numero de veces de repeticion de la simulacion de cada ciclo de simulaciones
N = 1 # si dejamos N=1 repetimos el ciclo solo 1 vez generando tantas graficas como nodos eliminemos
# number of nodes we remove one vertex at a time
fraction <- 0.9
K = round(size_lattice*size_lattice*fraction)
# Random spin  generation 
# choosing a random spins to be active
set.seed(123)
#state <- as.numeric(runif(size_lattice*size_lattice)>0.5) # aprox la mitad de los vertices en on, la otra mitad en off
state <- rep(1, size_lattice*size_lattice) # todos los spins/vertices activados.
# solo si quiero colorear el strength del nodo
range1to100 <- function(x){1 + 99*(x-min(x))/(max(x)-min(x))}
palf <- colorRampPalette(c("lightgray", "red"))
colr <- palf(100)
# # # # # # # # # ## # Parameters of the lattice # # # # # # # # # # #




# # # # # # # # # # # # Creation of the lattice # # # # # # # # # # #
# Lattice creation with L = size_lattice with weights
lat <- create_lattice(L = size_lattice)
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
# here we distinguish active vertex from inactive vertex
# Es importante aclarar que durante la simulacion, la red g no cambio,
# ni tampoco sus propiedades de nodos ni de edges.
# V(g)$color, V(g)$state cambian, E(g)$coupl, V(g)$str no cambian
V(g)$state <- state
E(g)$coupl <- Ji
E(g)$Je <- J0 # este es para calcular las energias de los clusters considerando J<0 (#28-abr-19)
V(g)$str <- strength(g, vids = V(g), weights = E(g)$coupl)
colores <- colr[round(range1to100(V(g)$str))]
V(g)$color <- ifelse(V(g)$state, colores, "black") # black =off, red = on
plot(g, layout=xy,
     vertex.size=5,
     #vertex.size=5,
     vertex.color = V(g)$color,
     vertex.label.cex = 0.8,
     vertex.label="",
     edge.width=4)
# # # # # # #  # # # Lattice with active vertices # # # # # # # # # #



# # # # # # # # # # # # ESTO YA NO SE EJECUTA # # # # # # # # # # # #
# Esto es lo que se va a repetir mcuhas veces
# # # # # # # # # Lattice of connected active spins # # # # # # # # #
# Creation of the lattice with edges connected to all active spins
# only adjacent vertexs which are both active can be connected between them, else, not.
# Here we distinguish the edges that will connect or not connect the vertex.
# 21-abr-19
new_lattice <- create_new_lattice_bonds(edg, g)
new_g <- new_lattice[[1]]
new_edg <- new_lattice[[2]]
plot(new_g, layout=xy,
     #vertex.size=20,
     vertex.size=5,
     #vertex.color = V(g)$color,
     #vertex.color = colr[round(range1to100(V(new_g)$str))], # solo si quiero colorear el strength del nodo
     vertex.color = ifelse(V(g)$state, colores, "black"), 
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
     #vertex.color = ifelse(V(g)$state, colr[round(range1to100(V(new_g)$str))], "black"),
     vertex.color = ifelse(V(g)$state, colores, "black"), 
     #vertex.label.cex = 0.8,
     vertex.label="",
     edge.width=4)
E(new_new_g)$coupl <- edg3[,"J"]
# # # # # # # Lattice of connected active spins ramdom # # # # # # #
# # # # # # # # # # # # ESTO YA NO SE EJECUTA # # # # # # # # # # # #





# # # # # # # # # # # # # # # Simulation # # # #  # # # # # # # # # #
# 17-abr-19
# repetimos el proces N veces, el proceso de formacion de clusters en 
# el lattice para un estado en particular.
almacen <- matrix(NA, ncol=7, nrow=K) #<S>, max(S), g(r), <k>, <E cluster>, <max E cluster> <E cluster con J<0 >
vertex_names <-  vector("list", length = K)
cluster_sizes_list <- vector("list", length = K) # table with ns
off_vertex_ids <- vector()
auxstr <- data.frame(str = V(g)$str, idvertex = 1:length(V(g)$str)) # generate the information
for (it in 1:K) {
  itk <- 100*it/K
  #format(itk, digits=1)
  cat("\r", "Processing data..Work done:", format(itk, digits=1), "%")
  #cat("\r", "Processing data..Work done:", 100*it/K, "%")
  
  #creation of the lattice with dropped vertexs (except for the first iteration)
  new_lattice <- create_new_lattice_bonds(edg, g)
  new_g <- new_lattice[[1]]
  new_edg <- new_lattice[[2]]
  
  # Creation of the lattice with edges connected to  octive spins according
  # probability p.
  E(new_g)$p <- 1 - exp(-2*E(new_g)$coupl) # bond activation probability
  # si p <= U, entonces hay conexion, si no, no hay conexion.
  bond <- ifelse(E(new_g)$p >= runif(nrow(new_edg),0, 1), 1, 0)
  bond <- ifelse(E(new_g)$p >= runif(nrow(new_edg),0, 1), 1, 0)
  #print(bond)
  J <- E(new_g)$coupl
  Je <- E(new_g)$Je #28-abr-19
  new_edg_ <- cbind(new_edg, bond, J, Je) # vertex1, vertex2, conextion 1=connected- 0 not connected, couplings (#28-abr-19)
  edg3 <- new_edg_[new_edg_[,"bond"]==1, ]
  if (is.vector(edg3)) {
    message("There are only one edge remaining. Iteration stopped!")
    break
  }
  new_new_edg <- edg3[,c(1,2)]
  new_new_g <- graph_from_edgelist(new_new_edg, directed = FALSE)
  number_nodos_que_faltan <- size_lattice*size_lattice - max(new_new_edg)
  if (number_nodos_que_faltan > 0) {
    new_new_g <- add_vertices(new_new_g, number_nodos_que_faltan )
  } 
  V(new_new_g)$str <- V(g)$str # como el lattice queda constante en numero de nodos, le asignamos los strength de la original.
  E(new_new_g)$coupl <- edg3[,"J"]
  E(new_new_g)$Je <- edg3[,"Je"] #28-abr-19
  
  # almacenamiento de observables
  cluster_sizes <- clusters(new_new_g)
  #length correlation g(r=1)
  g_r <- correlation_length(new_new_g, L=size_lattice)
  # degree
  dg <- mean(degree(new_new_g))
  # para calcular la energia de todos los clusters
  clid <- 1:length(cluster_sizes$csize)
  all_energies <- map_dbl(clid, energy_cluster, members = cluster_sizes$membership)
  all_energies_e <- map_dbl(clid, energy_cluster_e, members = cluster_sizes$membership) # energias considerando J < 0  # 28-abr-19
  #temp <- data.frame(e1=all_energies, e2=all_energies_e, s=cluster_sizes$csize) # analisis en caso que sea necesario
  #drop energies equal to zero (cluster with size=1)
  all_energies <- all_energies[all_energies != 0] #28-abr-19 antes era "<0" # 28-abr-19
  all_energies_e <- all_energies[all_energies_e != 0] #28-abr-19 antes era "<0"  # 28-abr-19
  almacen[it, ] <- c(mean(cluster_sizes$csize), max(cluster_sizes$csize), g_r, dg,
                     mean(all_energies), min(all_energies), mean(all_energies_e, na.rm = T)) # 28-abr-19
  vertex_names[[it]] <- find_name_maxcluster(cluster_sizes)
  cluster_sizes_list[[it]] <- table(cluster_sizes$csize)
  new_edg_ <- NULL
  
  #generacion de figuras para animar
  # name <- rename(it)
  # png(name)
  # plot(new_new_g, layout=xy,
  #      #vertex.size=20,
  #      vertex.size=5,
  #      #vertex.color = V(g)$color,
  #      vertex.color = ifelse(V(g)$state, colores, "black"), # solo si quiero colorear el strength del nodo
  #      #vertex.label.cex = 0.8,
  #      vertex.label="",
  #      edge.width=4,
  #      main = paste('Net lattice number=', it, '<S>=', round(mean(cluster_sizes$csize),2), '<E>=', round(mean(all_energies),2)) )
  # dev.off()
  
  # # # # # # # # # # # # # Dropping one vertex # # # ## # # # # # # # #
  #auxstr <- data.frame(str = V(new_g)$str, idvertex = 1:length(V(new_g)$str)) # generate the information
  #strmaxid <- which.max(auxstr[,"str"]) # find the rownumber with max strenght
  # strmaxid <- which.min(auxstr[,"str"]) # find the rownumber with min strenght
  strmaxid <- which(auxstr$str == min(auxstr$str, na.rm = TRUE)) #ind the rownumber with min strenght
  if (length(strmaxid) > 1) {
    strmaxid <- sample(strmaxid, 1)
  }
  node <- auxstr[strmaxid, "idvertex"] #find the id vertex with max strenght
  auxstr <- auxstr[!(auxstr$idvertex == node),] # drop the row with the max strength
  #d<-d[!(d$A=="B" & d$E==0),]
  off_vertex_ids <- c(off_vertex_ids, node) #aqui van quedando todos los nodos eliminados en orden
  rm(new_g, new_edg, new_new_g, new_new_edg)
  # generar red g -> no es necesario la dejamos igual solo modificamos lo que viene a continuacion:
  #V(g)$state <- state # vector de ones
  #V(g)$state[off_vertex_ids] <- 0
  V(g)$state[node] <- 0 # 22-abr-19
  #V(g)$str[off_vertex_ids] <- 0
  #V(g)$color <- ifelse(V(g)$state, "red", "black") # black =off, red = on
  # # # # # # # # # # # # # Dropping one vertex# # # ## # # # # # # # #
}

# analisis de los resultados:
almacen <- as.data.frame(almacen)
colnames(almacen) <- c("mean_s", "max_s", "gr", "k", "mean_E", "min_E", "mean_Ereal")
colMeans(almacen, na.rm=TRUE)
save.image("~/Dropbox/Research/Project_percolation/Percolation_data_code/results_from_percolation_bond_lattice_v4_220419.RData") 
save.image("~/Dropbox/Research/Project_percolation/Percolation_data_code/results_from_percolation_bond_lattice_v4_230419.RData") 
save.image("~/Dropbox/Research/Project_percolation/Percolation_data_code/results_from_percolation_bond_lattice_v4_230419b.RData") 
save.image("~/Dropbox/Research/Project_percolation/Percolation_data_code/results_from_percolation_bond_lattice_v4_230419c.RData") 
# size_lattice=20, J=random 0.1 y 2 con p=0.8 y 0.2 respectivamente, droping 200 vertex fron high to lowest strength.

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
plot(1:K, almacen$mean_s)
plot(1:K, almacen$k)
plot(1:K, almacen$mean_E)
plot(1:K, almacen$mean_Ereal)
plot(1:K, almacen$gr)
# check de mean_E con meanEreal
plot(almacen$mean_E, almacen$mean_Ereal)
# en vertex_name habria que hacer analisis de cuantas veces aparece cada 
# producto en cada largest cluster.
# # # # # # # # # # # # # # # Simulation # # # #  # # # # # # # # # #



# 02-may-19
# # # # # # # # # # # # # # Energy analysis # # # # # # # # # # # # # #
# Esto es solo la simulacion de 1 vez de percolacion del lattice para 
# analyzar la energia son considerar J<0, ves la energia considerando los
# acoples negativos.
# 1. ejecutar L42 a L147
# 2. creation of the lattice with dropped vertexs (except for the first iteration)
# in this case, all vertices are active. 
new_lattice <- create_new_lattice_bonds(edg, g)
new_g <- new_lattice[[1]]
new_edg <- new_lattice[[2]]
# 3. Creation of the lattice with edges connected to  octive spins according
# probability p.
E(new_g)$p <- 1 - exp(-2*E(new_g)$coupl) # bond activation probability
# si p <= U, entonces hay conexion, si no, no hay conexion.
bond <- ifelse(E(new_g)$p >= runif(nrow(new_edg),0, 1), 1, 0)
bond <- ifelse(E(new_g)$p >= runif(nrow(new_edg),0, 1), 1, 0)
J <- E(new_g)$coupl
Je <- E(new_g)$Je #28-abr-19
new_edg_ <- cbind(new_edg, bond, J, Je) # vertex1, vertex2, conextion 1=connected- 0 not connected, couplings (#28-abr-19)
edg3 <- new_edg_[new_edg_[,"bond"]==1, ]
new_new_edg <- edg3[,c(1,2)]
new_new_g <- graph_from_edgelist(new_new_edg, directed = FALSE)
number_nodos_que_faltan <- size_lattice*size_lattice - max(new_new_edg)
if (number_nodos_que_faltan > 0) {
  new_new_g <- add_vertices(new_new_g, number_nodos_que_faltan )
} 
V(new_new_g)$str <- V(g)$str # como el lattice queda constante en numero de nodos, le asignamos los strength de la original.
E(new_new_g)$coupl <- edg3[,"J"]
E(new_new_g)$Je <- edg3[,"Je"] #28-abr-19
# 4. almacenamiento de observables
cluster_sizes <- clusters(new_new_g)
# para calcular la energia de todos los clusters
clid <- 1:length(cluster_sizes$csize)
all_energies <- map_dbl(clid, energy_cluster, members = cluster_sizes$membership)
all_energies_e <- map_dbl(clid, energy_cluster_e, members = cluster_sizes$membership) # energias considerando J < 0  # 28-abr-19
temp <- data.frame(e1=all_energies, e2=all_energies_e, s=cluster_sizes$csize) # analisis en caso que sea necesario
# 5 analisis grafico 
plot(temp$e1, temp$e2, xlab="Energy - no negative Js", ylab="Total energy" )
plot(temp$s, temp$e1, ylab="Energy", xlab="Cluster size" )
points(temp$s, temp$e2, col="red", pch=19)
abline(coef = c(0,1))
text(7, -20, "red pts: total energy / white pts: Energy - no negative Js", cex = .8)
# # # # # # # # # # # # # # Energy analysis # # # # # # # # # # # # # #
