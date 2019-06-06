# Este script lo que hace es crear una red de Erdos-Renyi (ER) or random graph model
# con distintos valores de p (probabilidad de conexion entre nodos) la cual vamos
# cambiando y vamos viendo el tamano del cluster mayor.


# Creation date: 21-mar-19
# file name: percolation_ERnetwork.R

# Notes:
# 21-mar-19: creation 
# 21-mar-19: lo hacemos para varios valores de p
# 22-mar-19: hacemos una funcion para calcular el correlation lengh g(r=1)
#             g(r=1) es la probabilidad de encontrar a connected vertex at distance r of the reference
#              vertex. En un lattice tendremos, para r=1: 
# 27-mar-19: intentamos hacer un heatmap de la matriz de adyacencia de new_g



# Main process:
# Modelo de un grpah ER random graph model G(n,q) n es numero de nodos y q es prob de 
# coneccion entre un par de nodos elegidos al azar.
# Comenzamos con una red completa fully connected, de n nodos, por lo que habran
# n(n-1)/2 edges. This is a complete graph. Randomly  we delete edges from it.
# Sea p = 1-q la probabilidad de borrar un edge. Hacemos esto para distintos valores
# de q. 
# definimos lambda=np, y vamos graficando el promedio de los clusters de la red
# para varios valores de p, y graficamso entonces  lambda vs <s> y lambda vs correlation
# length.

# Recurso:
# https://www.youtube.com/watch?v=0SdzPJksV3Q
# https://rpubs.com/lgadar/generate-graphs


rm(list=ls())
library(ggplot2)
library(igraph)
load("results_from_percolation_ERnetwork_220319.Rda") # resultados con n = 100 nodos


# 21-mar-19
# 1. create a complete or full graph
g <- make_full_graph(n=300)
coords <- layout_(g, nicely())
plot(g, layout=coords,
     vertex.size=20,
     #vertex.size=4,
     vertex.label.cex = 0.8,
     #vertex.label="",
     edge.width=4)

E(g)
length(E(g)) # n(n-1)/2 edges
edg <- as_edgelist(g)

# 2. seleccionamos algunos edges de acuerdo al valor p
p <- 0.2
U <- runif(nrow(edg))
link <- ifelse(U < p, 1, 0) # 1 es conectado, 0 es desconectado
edg <- cbind(edg, U, link) 
# ahora removemos de edg, las filas con link = 0
edg2 <- edg[edg[,"link"]==1,]
new_edg <- edg2[,c(1,2)]

# 3. create the resulting graph
new_g <- graph_from_edgelist(new_edg, directed = FALSE)
plot(new_g, layout=coords,
     #vertex.size=20,
     vertex.size=2,
     #vertex.label.cex = 0.8,
     vertex.label="",
     edge.width=0.5)


# 3. get the clusters sizes
cluster_sizes <- clusters(new_g)
quantile(cluster_sizes$csize, probs = c(0.25, 0.5, 0.75))






# # # # # # # # # # # # ## # # # Heatmap of adjacency matrix # # # # # # # # # # # # # # # # # # #
# 27mar19
ady_M <- as_adjacency_matrix(new_g, type="both")
library(ggplot2)
ady_M <- as.matrix(ady_M)
library(lattice)
# https://www.rdocumentation.org/packages/dichromat/versions/1.1/topics/colorRampPalette
pal <- colorRampPalette(c("yellow", "black"), space = "rgb")
levelplot(ady_M, main="Adjacency matrix", xlab="", ylab="", col.regions=pal(4), cuts=3, at=seq(0,1,0.5))
# # # # # # # # # # # # ## # # # Heatmap of adjacency matrix # # # # # # # # # # # # # # # # # # #





# 23-mar-19
# # # # # Ahora lo hacemos para varios valores de probabilidad
n = 20
# Numero de simulaciones para cada valor de p
N <- 1000
probs <- seq(0.1,0.9,0.025)
# Almacenamiento de datos:
# mean of the degree <k>, size clusters: p25, p50 y p75 de la pdf del tamano de cluster, 
MyResults <- matrix(NA, ncol=5, nrow=length(probs))
for (probab in seq_along(1:length(probs) ) ) {
  #pb <- progress_bar$new(format = "  Progress [:bar] :percent eta: :eta", total = N*length(probs), clear = FALSE, width= 60) 
  cat(paste("iteration:", probab, "complete\n")) 
  p <- probs[probab]
  #Almacenamiento de datos
  # aqui vamos guardando los resultados:
  result <- matrix(NA, ncol=5, nrow=N) 
  # mean degree <k>, p25, p50 y p75 de la pdf del tamano de
  # de los clusters, correlation_length g(r=1)
  colnames(result) <- c("k", "p25", "p50", "p75", "g")
  # iteracion de simulaciones:
  for (n in seq_along(1:N) )  {
    # complete graph creation
    g <- make_full_graph(n=20)
    edg <- as_edgelist(g)
    
    
    # Simulation of random connection on the graph
    U <- runif(nrow(edg))
    link <- ifelse(U < p, 1, 0) # 1 es conectado, 0 es desconectado
    edg <- cbind(edg, U, link) 
    # ahora removemos de edg, las filas con link = 0
    edg2 <- edg[edg[,"link"]==1,]
    new_edg <- edg2[,c(1,2)]
    
    # creation of the ER graph
    new_g <- graph_from_edgelist(new_edg, directed = FALSE)
    
    # almacenamiento
    dg <- mean(degree(new_g))
    g_mean <- corr_length_net(new_g)
    cluster_sizes <- clusters(new_g)
    result[n, ] <- c(dg, quantile(cluster_sizes$csize, probs = c(0.25, 0.5, 0.75)), g_mean) 
    getmeans <- colMeans(result)
  }
  MyResults[probab, ] <- getmeans
}
colnames(MyResults) <- c("k", "p25", "p50", "p75", "g")
MyResults <- as.data.frame(MyResults)
MyResults$probability <- probs
#plot(MyResults$probability, MyResults$p50)
save(MyResults, file="results_from_percolation_ERnetwork_220319.Rda")



# Graph 1: probability of connection p versus mean of size of clusters
gr <- ggplot(MyResults, aes(x=probability, y=p50)) +
  geom_errorbar(aes(ymin=p25, ymax=p75), width=.007) + 
  geom_point() + 
  scale_shape_manual(values=c(3)) +
  xlab("Probability of connection") + 
  ylab("Cluster size <S>") +
  theme_bw()
gr
# Graph 2: probability of connection p versus mean of size of clusters
gr2 <- ggplot(MyResults, aes(x=k, y=p50)) +
  geom_point() + 
  scale_shape_manual(values=c(3)) +
  xlab("<k>") + 
  ylab("<S>") +
  theme_bw()
gr2
# Graph 4: mean degree <k> versus <s> mean of size of clusters
gr3 <- ggplot(MyResults, aes(x=k, y=p50) ) +
    geom_point() + 
    scale_shape_manual(values=c(3)) +
    xlab("k") + 
    ylab("<S>") +
    theme_bw()
gr3

# 22-mar-19
# # # # # # # correlation length
# Nota: Como no se trata de un lattice, y se trata de una red, es como una red con 
# dimension infinita, por lo que la correlation length aqui lo definimos como 
# la probabilidad de que un nodo tenga conexion a cualquier otro nodo de la red, 
# es decir, la degree del nodo dividido por el numero total de nodos de la red.
# inputs: new_g: objeto igraph
# outputs: g-mean = media de g(r=1) correlation length
corr_length_net <- function(new_g) {
  num_nodes <- length(V(new_g))
  G <-  vector(mode = "numeric", length = num_nodes)
  for (nodo in seq_along(1:num_nodes) ) {
    adyte <- as.numeric(neighbors(new_g, v=nodo))
    if ( is.null(adyte) ) {
      g <- 0
    } else {
      g <- degree(new_g, nodo)/num_nodes
    }
    G[nodo] <- g
  }
  g_mean <- mean(G)
  return(g_mean)
}
# ejemplo
# g_mean <- corr_length_net(new_g)




