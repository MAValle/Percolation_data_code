# Lo mismo que en percolation:ERnetwork.R, pero ahora para varios valores de 
# probabilidad p.


# Creation date: 23-mar-19
# file name: percolation_ERnetwork_v2.R

# Notes:
# 23-mar-19: creation 




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
source("correlation_length_net_function.R")
load("results_from_percolation_ERnetwork_220319.Rda") # resultados con n = 100 nodos




# 23-mar-19
# # # # # Ahora lo hacemos para varios valores de probabilidad
num_vertex = 100
# Numero de simulaciones para cada valor de p
N <- 1000
#probs <- seq(0.1,0.9,0.025) # este es el original
probs <- seq(0.1,0.9,0.5) # este es solo para probar
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
  # # # # iteracion de simulaciones:
  # complete graph creation
  g <- make_full_graph(n=num_vertex)
  edg <- as_edgelist(g)
  for (n in seq_along(1:N) )  {
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




# 25-mar-19 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# tratando de paralelizar de la L63 a la L79
internal_function <- function(g, edg, p) {
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
  result <- c(dg, quantile(cluster_sizes$csize, probs = c(0.25, 0.5, 0.75)), g_mean) 
  return(result)
  #result[n, ] <- c(dg, quantile(cluster_sizes$csize, probs = c(0.25, 0.5, 0.75)), g_mean) 
  #getmeans <- colMeans(result)
  #return(getmeans)
}
internal_function_for <- function(g, edg, p) {
  for (n in seq_along(1:N) ) {
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
  return(getmeans)
}

# utilizando internal_function_for (igual que el original) solo con dos probabiliades
system.time({ 
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
    # # # # iteracion de simulaciones:
    # complete graph creation
    g <- make_full_graph(n=num_vertex)
    edg <- as_edgelist(g)
    getmeans <- internal_function_for(g, edg, p)
    MyResults[probab, ] <- getmeans
  }
})
#user      system   elapsed 
#278.016  36.546   316.347  = 5.27 min


# Ahora en paralelo solo con dos probabilidades
# https://www.r-bloggers.com/the-wonders-of-foreach/
library(doMC)
registerDoMC(cores=4)
system.time({ 
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
  # # # # iteracion de simulaciones:
  # complete graph creation
  g <- make_full_graph(n=num_vertex)
  edg <- as_edgelist(g)
  result <- foreach(n = 1:N, .combine = rbind) %dopar% internal_function(g, edg, p)
  getmeans <- colMeans(result)
  MyResults[probab, ] <- getmeans
}
})
#user       system  elapsed 
#164.340   2.535    45.190  = 0.75 min.