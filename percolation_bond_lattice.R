# En este script llevamos a cabo estudio de percolation in a lattice
# en que la logica es la siguiente:
# seleccionamos o activamos aleatoriamente una fracción p de nodos del 
# lattice. Si dos nodos estan prendidos y son vecinos, entonces se
# forma un edge entre ellos. Luego estudiamos el degree, tamano de cluster y so
# percola para distintos valores de p.


# Creation date: 04-abr-19
# file name: percolation_bond_lattice.R

# Notes:
# 04-abr-19: creation 

# Recursos:
# https://www.youtube.com/watch?v=0SdzPJksV3Q
# https://rpubs.com/lgadar/generate-graphs


# Main process:
# Comenzamos lattice fully connected de tamano L con N = LXL nodos. 
# Luego comenzamos a seleccionar aleatoriamente una fraccion q=1-p de nodos
# 


# Luego comenzamos un proceso de dropeo o erase of edges en forma aleatoria 
# en que Randomly  we delete edges. En este caso el peso de cada link
# representa la probabilidad de que el link se borre. En este caso es como 
# tenerdistintas probabilidad p de borrar un edge. Por lo tanto tendremos 
# asociado a cada edge un p_ij que es el peso entre el nodo i y j.


# Objetivo: comparar el <k> degree y p, versis <S> para un lattice L=100
# con percolacion simple (from percolation_lattice_v3.R) y con percolacion
# cpn pesos (percolation_Wlattice.R)



# # # # # # # # # # # # # # Carga de datos # # # # # # # # # # # # # #
rm(list = ls())
source("create_lattice_function.R")
source("create_random_spins_activations_function.R")
source("percolation_detection_function.R")
source("correlation_length_function.R")
library(ggplot2)
library(igraph)
library(doMC)
load("results_from_percolation_bond_lattice_050419L50.Rda") # lattice L=20, con weights U(0,1)
# # # # # # # # # # # # # # Carga de datos # # # # # # # # # # # # # #




# Funcion que realiza computes necesarias para poder paralelizar
# 26-mar-19 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
internal_function <- function(g, edg, p, L) {
  # Simulation of random activation of spins (vertex) and removal of 
  # edges not connected to activated spins.
  new_g <- create_random_spins_activations(g, p=p)
  
  # Percolation detection
  detection <- percolation_detection(new_g, L=size_lattice)
  value <- detection[[1]] # FALSE no percolation, TRUE, percolation detected.
  cluster_sizes <- clusters(new_g)
  g_r <- correlation_length(new_g, L=size_lattice)
  dg <- mean(degree(new_g))
  
  result <- c(value, dg, min(cluster_sizes$csize), mean(cluster_sizes$csize), max(cluster_sizes$csize), g_r) 
  return(result)
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # Simulacion # # # # # # # # # # # # # # # #
# Simulating percolations for different values of connection probabilities p 
# Numero de simulaciones para cada valor de p
N <- 50
# tamano del lattice:
size_lattice <- 10
# valores de probabilidad p:
probs <- seq(0.1,0.9,0.025)
#probs <- c(seq(0.1,0.3, 0.02),seq(0.4,0.9,0.1))
# Almacenamiento de datos:
# percolation: false no percola, mean degree, true percola, p25, p50 y p75 de la pdf del tamano de cluster, <g(r)>
MyResults <- matrix(NA, ncol=6, nrow=length(probs))
# paralelizacion
for (probab in seq_along(1:length(probs) ) ) {
  #pb <- progress_bar$new(format = "  Progress [:bar] :percent eta: :eta", total = N*length(probs), clear = FALSE, width= 60) 
  cat(paste("iteration:", probab, "complete\n")) 
  p <- probs[probab]
  #Almacenamiento de datos
  # aqui vamos guardando los resultados:
  result <- matrix(NA, ncol=6, nrow=N) 
  # percolation: false (0) no percola (1), true percola, p25, p50 y p75 de la pdf del tamano de
  # de los clusters.
  # colnames(result) <- c("percolation", "k", "p25", "p50", "p75", "g_r")
  # iteracion de simulaciones:
  # Lattice creation with L = size_lattice with weights
  lat <- create_lattice(L=size_lattice)
  xy <- lat[[1]]
  edg <- lat[[2]]
  g <- lat[[3]]
  # # Paralelizacion
  result <- foreach(n = 1:N, .combine = rbind) %dopar% internal_function(g=g, edg=edg, p=p, L=size_lattice)
  getmeans <- colMeans(result)
  MyResults[probab, ] <- getmeans
}
# FIN Paralelizacion
colnames(MyResults) <- c("percolation", "k", "clustersize_min", "clustersize_mean", "clustersize_max", "gr")
MyResults <- as.data.frame(MyResults)
MyResults$probability <- probs
#plot(MyResults$probability, MyResults$p50)
MyResults$S <- MyResults$clustersize_max/size_lattice^2 #Normalized cluster size
save(MyResults, file="results_from_percolation_bond_lattice_050419L10.Rda")
# # # # # # # # # # # # # # # Simulacion # # # # # # # # # # # # # # # #




# # # # # # # # # # # # # # # Graficas# # # # # # # # # # # # # # # #
# Graph 1: probability of connection p versus mean of size of clusters
gr <- ggplot(MyResults, aes(x=probability, y=clustersize_max)) +
  #geom_errorbar(aes(ymin=p25, ymax=p75), width=.007) + 
  geom_point() + 
  scale_shape_manual(values=c(3)) +
  xlab("Probability of connection") + 
  ylab("Cluster size <S>") +
  theme_bw()
gr
# Graph 1b: probability of connection p versus mean of max size of clusters
grb <- ggplot(MyResults, aes(x=probability, y=S)) +
  #geom_errorbar(aes(ymin=p25, ymax=p75), width=.007) + 
  geom_point() + 
  scale_shape_manual(values=c(3)) +
  xlab("Probability of connection") + 
  ylab("max Cluster size <S>") +
  theme_bw()
grb
# Graph 2: probability of connection p versus mean of size of clusters
gr2 <- ggplot(MyResults, aes(x=k, y=S)) +
  geom_point() + 
  scale_shape_manual(values=c(3)) +
  xlab("<k>") + 
  ylab("<S/N>") +
  theme_bw()
gr2
# Graph 2: probability of connection p versus mean of size of clusters
gr3 <- ggplot(MyResults, aes(x=probability, y=percolation)) +
  geom_point() + 
  scale_shape_manual(values=c(3)) +
  xlab("p") + 
  ylab("% times that percolates") +
  theme_bw()
gr3
# Graph 2: probability of connection p versus mean correlation length
gr4 <- ggplot(MyResults, aes(x=probability, y=gr)) +
  geom_point() + 
  scale_shape_manual(values=c(3)) +
  xlab("p") + 
  ylab("Correlation length") +
  theme_bw()
gr4
# # # # # # # # # # # # # # # Graficas# # # # # # # # # # # # # # # #
