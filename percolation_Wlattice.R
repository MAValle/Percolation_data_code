# En este script llevamos a cabo estudio de percolation in a weigthed lattice
# Utilizamos la misma estructura algoritmica de percolation_lattice_v3.R

# Lo mismo que en percolation_lattice_v3.R, pero ahora para varios valores de 
# probabilidad p. La idea es que donde simulamos percolacion con pesos en la 
# red  con pik de distintas distribuciones, siempre entre 0 y 1.


# Creation date: 23-mar-19
# file name: percolation_Wlattice.R

# Notes:
# 29-mar-19: creation 

# Recursos:
# https://www.youtube.com/watch?v=0SdzPJksV3Q
# https://rpubs.com/lgadar/generate-graphs


# Main process:
# Comenzamos lattice fully connected de tamano L con n = LXL nodos. 
# A cada edge le asignaremos
# un peso entre 0 y 1 de acuerdo a alguna distribucion pre-establecido (por ejemplo
# una normal para simular un disordered glass).
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
#source("create_lattice_function.R")
source("create_weigthed_lattice_function.R")
#source("create_random_lattice_net_function.R")
source("create_random_Wlattice_net_function.R")
source("percolation_detection_function.R")
source("correlation_length_function.R")
library(ggplot2)
library(igraph)
library(doMC)
load("results_from_percolation_Wlattice_030419L20.Rda") # lattice L=20, con weights U(0,1)
load("results_from_percolation_Wlattice_030419L20b.Rda") # lattice L=20, con weights Normal truncado entre 0 y 1, con media 0.5

# # # # # # # # # # # # # # Carga de datos # # # # # # # # # # # # # #

# Funcion que realiza computes necesarias para poder paralelizar
# 26-mar-19 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
internal_function <- function(g, edg, p, L) {
  # Simulation of random connection on lattice
  S <- create_random_Wlattice_net(xy, edg, p=p, L=L)
  new_g <- S[[1]]
  
  # Percolation detection
  detection <- percolation_detection(new_g, L=size_lattice)
  value <- detection[[1]] # FALSE no percolation, TRUE, percolation detected.
  cluster_sizes <- clusters(new_g)
  g_r <- correlation_length(new_g, L=size_lattice)
  dg <- mean(degree(new_g))
  
  result <- c(value, dg, quantile(cluster_sizes$csize, probs = c(0.25, 0.5, 0.75)), g_r) 
  return(result)
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # Simulacion # # # # # # # # # # # # # # # #
# Simulating percolations for different values of connection probabilities p 
# Numero de simulaciones para cada valor de p
N <- 60
# tamano del lattice:
size_lattice <- 20
# valores de probabilidad p:
#probs <- seq(0.1,0.9,0.025)
probs <- c(seq(0.1,0.3, 0.02),seq(0.4,0.9,0.1))
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
  #lat <- create_lattice_weighted(L=size_lattice, dist='norm', moda=0.5, mi=0, ma=1)
  #lat <- create_lattice_weighted(L=size_lattice, dist='norm', moda=0.3, mi=0, ma=1)
  lat <- create_lattice_weighted(L=size_lattice, dist='norm', moda=0.8, mi=0, ma=1)
  xy <- lat[[1]]
  edg <- lat[[2]]
  g <- lat[[3]]
  # # Paralelizacion
  result <- foreach(n = 1:N, .combine = rbind) %dopar% internal_function(g=g, edg=edg, p=p, L=size_lattice)
  getmeans <- colMeans(result)
  MyResults[probab, ] <- getmeans
}
# FIN Paralelizacion
colnames(MyResults) <- c("percolation", "k",  "p25", "p50", "p75", "gr")
MyResults <- as.data.frame(MyResults)
MyResults$probability <- probs
#plot(MyResults$probability, MyResults$p50)
save(MyResults, file="results_from_percolation_Wlattice_030419L20b.Rda")
MyResults$S <- MyResults$p50/size_lattice^2 #Normalized cluster size
# # # # # # # # # # # # # # # Simulacion # # # # # # # # # # # # # # # #




# # # # # # # # # # # # # # # Graficas# # # # # # # # # # # # # # # #
# Graph 1: probability of connection p versus mean of size of clusters
gr <- ggplot(MyResults, aes(x=probability, y=p50)) +
  #geom_errorbar(aes(ymin=p25, ymax=p75), width=.007) + 
  geom_point() + 
  scale_shape_manual(values=c(3)) +
  xlab("Probability of connection") + 
  ylab("Cluster size <S>") +
  theme_bw()
gr
# Graph 2: probability of connection p versus mean of size of clusters
gr2 <- ggplot(MyResults, aes(x=k, y=S)) +
  geom_point() + 
  scale_shape_manual(values=c(3)) +
  xlab("<k>") + 
  ylab("<S/N>") +
  theme_bw()
gr2
# Graph 2: probability of connection p versus mean of size of clusters
gr3 <- ggplot(MyResults, aes(x=1-probability, y=percolation)) +
  geom_point() + 
  scale_shape_manual(values=c(3)) +
  xlab("p") + 
  ylab("% times that percolates") +
  theme_bw()
gr3
# # # # # # # # # # # # # # # Graficas# # # # # # # # # # # # # # # #
