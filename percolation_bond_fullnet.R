# Hacemos lo mismo que en percolation_bond_lattice_v4a.R pero esta vez sobre una red 
# completa a la cual le asignamos J randoms entre -1 y 1. Luego podemos aplicar un 
# umbral para eliminar acoples muy pequeños. Luego conectamos nodos (todos activos) 
# solo para aquellos con J>0 con probabilidad p=1-exp(-2J) y calculamos la energia 
# de los clusters y otros observables.

# Procedure: 
# step 1: Create a complete network with make_full_graph(n) con n nodos o spins.
# step 1.1: Asignamos acoples a cada edge de acuerdo a distribucion
# step 2: conectamos los nodos con acoples positivos con probabilidad p = 1 - exp(-2J)
# step 3: computamos observables (size clusters, energies, length correlation, etc..)
# Lo podemos hacer a distintas temperaturas para ver a que temperatura percola.

# Creation date: 02-may-19
# file name: percolation_bond_fullnet.R

# Notes:
# 02-may-19: creation 
# 07-may-19: agregamos a la grafica un mst para ir viendo como se relaciona este con
#         los clusters que se generan a cierto nivel de probabilidad umbral. 
# 08-may-19: agregamos un loop para hacer la grafica de percolacion. P
# 09-may-19: agregamos la medicion de la energia de cada cluster.
# 11-may-19: agregamos probs <- rev(probs) para simular de mayor umbral a menor umbral p*
 



# # # # # # # # # # # # # # Carga de datos # # # # # # # # # # # # # #
rm(list = ls())
source("create_fullnet_function.R") 
source("create_fullnet_bonds_function.R")
source("energy_cluster_function.R")
#source("correlation_length_function.R")
#source("create_new_lattice_bonds_function.R")
library(ggplot2)
library(igraph)
library(doMC)
library(purrr)
library(latex2exp)


# 18-abr-19
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
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

# # # # # # # # # ## # Parameters of the simulation # # # # # # # # # # #
# tamano del lattice:
num_vertex <- 40
# Para un net completo de N nodos, habran N*(N-1)/2 edges
number_of_edges <- num_vertex*(num_vertex-1)/2
set.seed(123)
# umbral de probabilidad deltap
deltap = 0.857
# # # # # # # # distribution of couplings # # # # # # # # 
# Distribucion ferromagnetica de acoples 
#J0 <- Ji <- rep(1, number_of_edges)
# Distribucion de los acoples: random [0,1]
#J0 <- Ji <- runif(number_of_edges)
# Distribucion de los acoples: random [-1,1]
#J0 <- Ji <- runif(number_of_edges)*2-1   # J0 incluye los J < 0, mientras que Ji no los oncluye, es solo para la red.
#Ji[Ji < 0] <- 0
# Distribucion de exponencial de acoples positivos
#J0 <- Ji <- rexp(number_of_edges, 3) # van a dar entre 0 y app 2.5
# Distribucion de exponencial de acoples positivos y negativos centrados en 0
#J0 <- Ji <- c(rexp(number_of_edges*0.5, 3), -rexp(number_of_edges*0.5, 3))
#Ji[Ji < 0] <- 0
# custom discrete tistribution of J
#J0 <- Ji <- sample(x = c(0.1,2), number_of_edges, replace = T, prob = c(0.8, 0.2))
#Ji[Ji < 0] <- 0
# Normal distribution of J truncated between 0 and 2
#library(truncnorm)
#J0 <- Ji <- rtruncnorm(number_of_edges, a=0, b=2, mean=0.5, sd=0.5)
# Normal distribution of J truncated between -1 and 1
#library(truncnorm)
J0 <- Ji <- rtruncnorm(number_of_edges, a=-2, b=2, mean=0, sd=0.5)
Ji[Ji < 0] <- 0
# # # # # # # # distribution of couplings # # # # # # # # 
# Random spin  generation 
# choosing a random spins to be active
#state <- as.numeric(runif(num_vertex)>0.5) # aprox la mitad de los vertices en on, la otra mitad en off
state <- rep(1, num_vertex) # todos los spins/vertices activados.
# solo si quiero colorear el strength del nodo
range1to100 <- function(x){1 + 99*(x-min(x))/(max(x)-min(x))}
palf <- colorRampPalette(c("lightgray", "red"))
colr <- palf(100)
# # # # # # # # # ## # Parameters of the simulation # # # # # # # # # # #


# Step 1: crear complete network
lat <- create_fullnet(num_vertex)
xy <- lat[[1]]
edg <- lat[[2]]
g <- lat[[3]]
plot(g, layout=xy,
     vertex.size=5,
     vertex.frame.color = "white",
     vertex.color="orange",
     vertex.label="",
     vertex.arrow.mode=0,
     vertex.label.cex = 0.5,
     edge.width=1)


# step 1.1: we assign attributes and coloring of vertexs according to their strength (only J>0)
V(g)$state <- state
E(g)$coupl <- Ji 
E(g)$Je <- J0 # este es para calcular las energias de los clusters considerando J<0 (#28-abr-19)
V(g)$str <- strength(g, vids = V(g), weights = E(g)$coupl)
colores <- colr[round(range1to100(V(g)$str))]
V(g)$color <- ifelse(V(g)$state, colores, "black") # black =off, red = on
plot(g, layout=xy,
     vertex.size=5,
     vertex.frame.color = "white",
     vertex.color = V(g)$color,
     vertex.label.cex = 0.8,
     vertex.label="",
     vertex.arrow.mode=0,
     edge.width=1)

# step 1.2 building the MST
# Tranformar los acoples en una medida de distancia
library(gdata)
#M <- as_adjacency_matrix(g, type="both", names=TRUE, sparse=FALSE, attr="coupl") # matriz de adyacencia de acoples (en vez de coupl debiera ser Je)
M <- as_adjacency_matrix(g, type="both", names=TRUE, sparse=FALSE, attr="Je") # matriz de adyacencia de acoples (en vez de coupl debiera ser Je)
Ds <- sqrt(-M + 3) 
# Convertir la matriz de distancia en objeto igraph
netdis <- graph_from_adjacency_matrix(as.matrix(Ds), mode="upper", weighted=TRUE, diag=FALSE) # esta red 
E(netdis)$Je <- E(g)$Je # este es para calcular las energias de los clusters (incluye J < 0)
# es solo para rescatar los couplings de la red de acoples y pasarlas a g, y luego a mst_g.
# en mstdistance_and_energies.R hay codigos para hacer mst con ggnet.
mst_g <- minimum.spanning.tree(netdis, algorithm="prim")
plot(mst_g, 
     #layout=xy,
     vertex.size=5,
     vertex.frame.color = "white",
     vertex.color = V(g)$color,
     vertex.label.cex = 0.8,
     vertex.label="",
     vertex.arrow.mode=0,
     edge.color = ifelse(E(mst_g)$Je > 0, "darkgrey", "blue"),
     edge.width=2, 
     main="MST")

# step 1.3: we desconnect edges between then in case of having off vertexs
new_fullnet <- create_new_fullnet_bonds(edg, g)
new_g <- new_fullnet[[1]]
new_edg <- new_fullnet[[2]]
plot(new_g, layout=xy,
     vertex.size=5,
     vertex.frame.color = "white",
     vertex.color = V(g)$color,
     vertex.label.cex = 0.8,
     vertex.label="",
     vertex.arrow.mode=0,
     edge.width=1)

# step 2: conectamos los nodos con acoples positivos con probabilidad p = 1 - exp(-2J)
# Creation of the lattice with edges connected to  octive spins according
# probability p.
E(new_g)$p <- 1 - exp(-2*E(new_g)$coupl) # bond activation probability
# # # # # # # probability umbral deltap # # # # # # # 
# en vez de aplicar un umbral a los acoples, mejor aplicamos un umbral a las probabilidades
# de conexion: Para todas las probabilidades E(new_g)$p > deltap, se genera un 
# random number r entre 0 y 1. Si E(new_g)$p > r, entonces hay conexion. Si no, se
# deja sin edge. Para todas las probabilidades E(new_g)$p < deltap se deja sin conexion.
# Esto se asemeja a percolating net con distintas umbrales de probabilidad.
E(new_g)$p <- ifelse( E(new_g)$p >= deltap, E(new_g)$p, 0)
# si p >= r, entonces hay conexion, si no, no hay conexion.
bond <- ifelse(E(new_g)$p >= runif(nrow(new_edg),0, 1), 1, 0)
#print(bond)
J <- E(new_g)$coupl # solo J > 0
Je <- E(new_g)$Je #28-abr-19
new_edg_ <- cbind(new_edg, bond, J, Je) # vertex1, vertex2, conextion 1=connected- 0 not connected, couplings (#28-abr-19)
edg3 <- new_edg_[new_edg_[,"bond"]==1, ]
new_new_edg <- edg3[,c(1,2)]
new_new_g <- graph_from_edgelist(new_new_edg, directed = FALSE)
number_nodos_que_faltan <- num_vertex - max(new_new_edg)
if (number_nodos_que_faltan > 0) {
  new_new_g <- add_vertices(new_new_g, number_nodos_que_faltan )
} 
V(new_new_g)$str <- V(g)$str # como el lattice queda constante en numero de nodos, le asignamos los strength de la original.
E(new_new_g)$coupl <- edg3[,"J"]
E(new_new_g)$Je <- edg3[,"Je"] #28-abr-19
plot(new_new_g, layout=xy,
     vertex.size=5,
     vertex.frame.color = "white",
     vertex.color = V(g)$color,
     vertex.label.cex = 0.8,
     vertex.label="",
     vertex.arrow.mode=0,
     edge.width=1)


# step 3: computamos observables (size clusters, energies, length correlation, etc..)
# almacenamiento de observables
cluster_sizes <- clusters(new_new_g)
cs <- cluster_sizes$csize # size of the clusters
mean(cs[which(cs > 1)]) # media de los clusters sin tomar en cuenta los nodos apagados 
sum(cs > 1) # numero de cluster con size > 1

# degree
dg <- mean(degree(new_new_g))
dg2 <- degree(new_new_g)
mean(dg2[which(dg2 > 0)]) # degree real sin tomar en cuenta nodos apagados o con degree 0


# energy cluster
number_of_clusters <- length(cluster_sizes$csize)
# ejemplo: 
E <- energy_cluster(cl=10, members = cluster_sizes$membership)
# para calcular la energia de todos los clusters
clid <- 1:number_of_clusters
all_energies <- map_dbl(clid, energy_cluster, members = cluster_sizes$membership)
#drop energies equal to zero (cluster with size=1)
all_energies <- all_energies[all_energies != 0]
#mean(all_energies)
#max(all_energies)

# energy cluster taking J < 0.
# ejemplo: 
Er <- energy_cluster_f(cl=10, members = cluster_sizes$membership)
all_energies_f <- map_dbl(clid, energy_cluster_f, members = cluster_sizes$membership) # energias considerando J < 0  # 28-abr-19
all_energies_f <- all_energies_f[all_energies_f != 0] #28-abr-19 antes era "<0"  # 28-abr-19


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #




# 08-MAY-19
# # # # # # REPEAT THE PROCESS LOT OF TIMES FOR DIFFERENTS PROBABILITIES THRESHOLDS # # # # # # 
# Para ejecutar estas simulaciones, hay que cargar previamente los parametros de la simulacion al inicio.
# Numero de simulaciones para cada valor de threshold p
N <- 1
probs <- c(seq(0.1,0.7,0.05), seq(0.7,0.99, 0.001)) 
# 11-may-19
probs <- rev(probs) # en caso que quisiera simularlo desde un mayor umbral a uno menor (como en percolacion clasica)
#probs <- seq(0.1,0.95,0.005) # este es solo para probar
# Almacenamiento de datos:
# threshold probability pt, <k>, <S>, max(S), <E>, max(E) 
MyResults <- matrix(NA, ncol=7, nrow=length(probs))
# Nota: estamos dejando fijo los aciples y fijo el numero de spins (nodos) activos.
for (probab in seq_along(1:length(probs) ) ) {
  #cat(paste("iteration:", probab, "complete\n")) 
  deltap <- probs[probab]
  cat(paste("iteration:", probab, "Threshold p*:", deltap, "complete\n")) 
  #Almacenamiento de datos
  # aqui vamos guardando los resultados:
  result <- matrix(NA, ncol=6, nrow=N) 
  colnames(result) <- c("k", "S", "maxS", "E", "Ereal", "maxEreal")
  # Step 1: crear complete network
  lat <- create_fullnet(num_vertex)
  xy <- lat[[1]]
  edg <- lat[[2]]
  g <- lat[[3]]
  # step 1.1: we assign attributes and coloring of vertexs according to their strength (only J>0)
  V(g)$state <- state
  E(g)$coupl <- Ji
  E(g)$Je <- J0 # este es para calcular las energias de los clusters considerando J<0 (#28-abr-19)
  V(g)$str <- strength(g, vids = V(g), weights = E(g)$coupl)
  colores <- colr[round(range1to100(V(g)$str))]
  V(g)$color <- ifelse(V(g)$state, colores, "black") # black =off, red = on
  # step 1.2: build de the MST (no es necesario aqui)
  # step 1.3: we desconnect edges between then in case of having off vertexs
  new_fullnet <- create_new_fullnet_bonds(edg, g)
  new_g <- new_fullnet[[1]]
  new_edg <- new_fullnet[[2]]
  
  # step 2: conectamos los nodos con acoples positivos con probabilidad p = 1 - exp(-2J)
  # Creation of the lattice with edges connected to  octive spins according
  # probability p.
  E(new_g)$p <- 1 - exp(-2*E(new_g)$coupl) # bond activation probability
  # # # # # # # probability umbral deltap # # # # # # # 
  E(new_g)$p <- ifelse( E(new_g)$p >= deltap, E(new_g)$p, 0)
  for (it in 1:N) {
    # si p >= r, entonces hay conexion, si no, no hay conexion.
    bond <- ifelse(E(new_g)$p >= runif(nrow(new_edg),0, 1), 1, 0)
    #print(bond)
    J <- E(new_g)$coupl
    Je <- E(new_g)$Je #28-abr-19
    new_edg_ <- cbind(new_edg, bond, J, Je) # vertex1, vertex2, conextion 1=connected- 0 not connected, couplings (#28-abr-19)
    edg3 <- new_edg_[new_edg_[,"bond"]==1, , drop=F]
    if (nrow(edg3) == 0) {
      message("The network is now completely disconnected! Iteration has stopped.")
      break
    }
    new_new_edg <- edg3[,c(1,2), drop = F]
    new_new_g <- graph_from_edgelist(new_new_edg, directed = FALSE)
    number_nodos_que_faltan <- num_vertex - max(new_new_edg)
    if (number_nodos_que_faltan > 0) {
      new_new_g <- add_vertices(new_new_g, number_nodos_que_faltan )
    } 
    V(new_new_g)$str <- V(g)$str # como el lattice queda constante en numero de nodos, le asignamos los strength de la original.
    E(new_new_g)$coupl <- edg3[,"J"]
    E(new_new_g)$Je <- edg3[,"Je"] #28-abr-19
    
    # step 3: computamos observables (size clusters, energies, length correlation, etc..)
    # almacenamiento de observables
    cluster_sizes <- clusters(new_new_g)
    cs <- cluster_sizes$csize
    # degree
    dg <- mean(degree(new_new_g))
    dg2 <- degree(new_new_g)
    # energy of the clusters
    # para calcular la energia de todos los clusters
    clid <- 1:length(cluster_sizes$csize)
    all_energies <- map_dbl(clid, energy_cluster, members = cluster_sizes$membership)
    #drop energies equal to zero (cluster with size=1)
    all_energies <- all_energies[all_energies != 0]
    # energies taking J < 0:
    all_energies_f <- map_dbl(clid, energy_cluster_f, members = cluster_sizes$membership) # energias considerando J < 0  # 28-abr-19
    all_energies_f <- all_energies_f[all_energies_f != 0] #28-abr-19 antes era "<0"  # 28-abr-19
    
    result[it, ] <- c(dg, mean(cs, na.rm=T), max(cs, na.rm=T), mean(all_energies, na.rm=T), mean(all_energies_f, na.rm=T), max(all_energies_f, na.rm=T))
  }
  toputinside <- colMeans(result)
  MyResults[probab, ] <- c(deltap, toputinside)
}

# # # Analysis of the results:
#MyResults <- MyResults[complete.cases(MyResults), ]
colnames(MyResults) <- c("pThreshold", "k", "s", "maxs", "e", "ereal", "maxereal")
MyResults <- as.data.frame(MyResults)
MyResults$snorm <- MyResults$s/num_vertex
MyResults$qThreshold <- 1 - MyResults$pThreshold
#plot((MyResults$pThreshold), (MyResults$snorm), pch = "+", col="red", main="<S>")
plot((MyResults$qThreshold), (MyResults$snorm), pch = "+", col="red", main="<S>", xlab="1-p*", ylab="<s>/N")
#plot((MyResults$pThreshold), (MyResults$maxs), pch = "+", col="red", main="max(S)")
#plot((MyResults$qThreshold), (MyResults$maxs), pch = "+", col="red", main="max(S)", xlabel="1-p*")
#plot(MyResults$pThreshold, MyResults$k, pch = "+", col="red", main="<k>")
plot(MyResults$qThreshold, MyResults$k, pch = "+", col="red", main="<k>", xlab="1-p*", ylab="<k>")
plot(MyResults$qThreshold, MyResults$e, pch = "+", col="red", main="<E>", xlab="1-p*", ylab="<E>")
points(MyResults$qThreshold, MyResults$ereal, col="grey", pch="o")
plot(MyResults$snorm, MyResults$e, pch = "+", col="red", xlab="<S>", ylab="<E>")
plot(MyResults$snorm, MyResults$ereal, pch = "+", col="red", xlab="<S>", ylab="<Ereal>")
# # # # # # REPEAT THE PROCESS LOT OF TIMES FOR DIFFERENTS PROBABILITIES THRESHOLDS # # # # # # 