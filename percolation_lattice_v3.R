# Este script lo que hace es simular percolation in a lattice of LxL = N
# Utilizamos las funciones originadas en percolation_lattice.R

# Creation date: 26-mar-19
# file name: percolation_lattice_v3.R

# Funcion clusters de igraph:
# La funcion clusters de igraph se tiene: $membership que nos da el valor numerico que
# identifica al cluser
# $scize: nos da el tamano de cada cluster
# $no: nos da el numero de cluster que hay.

# Que vamos a medir: Logic, p25, p50, p75 de la distribucion del tamano de los cluster.
# Logic indica TRUE or FALSE, si percola o no percola respectivamente.

# Notes:
# 13-mar-19: creation
# 17-mar-19: hamos aagregado al loop el calculo del correlation length con g=1 utilizando
#           la funcion correlation_lenght.
# 22-mar-19: le agregue a los resultados el degree <k> de los nodos del lattice.



rm(list = ls())
source("create_lattice_function.R")
source("create_random_lattice_net_function.R")
source("percolation_detection_function.R")
source("correlation_length_function.R")
library(ggplot2)
library(igraph)
library(doMC)
load("results_from_percolation_lattice_v2_160319.Rda") # resultados con L=10
load("results_from_percolation_lattice_v2_160319b.Rda") # resultados con L=8
load("results_from_percolation_lattice_v2_160319c.Rda") # resultados con L=12
load("results_from_percolation_lattice_v2_160319d.Rda") # resultados con L=20

# Funcion que realiza computes necesarias para poder paralelizar
# 26-mar-19 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
internal_function <- function(g, edg, p) {
  # Simulation of random connection on lattice
  S <- create_random_lattice_net(xy, edg, p=p, L=size_lattice)
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


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Simulating percolations for different values of connection probabilities p 
# Numero de simulaciones para cada valor de p
N <- 1000
# tamano del lattice:
size_lattice <- 20
# valores de probabilidad p:
probs <- seq(0.1,0.9,0.025)
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
  colnames(result) <- c("percolation", "k", "p25", "p50", "p75", "g(r)")
  # iteracion de simulaciones:
  # Lattice creation with L = size_lattice
  lat <- create_lattice(L=size_lattice)
  xy <- lat[[1]]
  edg <- lat[[2]]
  g <- lat[[3]]
  # # Paralelizacion
  result <- foreach(n = 1:N, .combine = rbind) %dopar% internal_function(g, edg, p)
  getmeans <- colMeans(result)
  MyResults[probab, ] <- getmeans
}
# FIN Paralelizacion
colnames(MyResults) <- c("percolation", "k",  "p25", "p50", "p75", "g(r)")
MyResults <- as.data.frame(MyResults)
MyResults$probability <- probs
#plot(MyResults$probability, MyResults$p50)
save(MyResults, file="results_from_percolation_lattice_v2_160319d.Rda")
MyResults$S <- MyResults$p50/size_lattice^2 #Normalized cluster size
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



# 16-mar-19
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Graph of mean of  p35, p50 and p75 (distribution of size of clusters)
# for different values of p.
library(reshape2)
df <- melt(MyResults, id='probability')
dfs <- subset(df, variable != "percolation")

# Graph 1: probability of connection p versus mean of size of clusters
p <- ggplot(MyResults, aes(x=probability, y=S)) +
  geom_errorbar(aes(ymin=p25/size_lattice^2, ymax=p75/size_lattice^2), width=.007) + 
  geom_point() + 
  scale_shape_manual(values=c(3)) +
  xlab("Probability of connection") + 
  ylab("Cluster size") +
  theme_bw()
p
# Graph 2: probability of connection p versus mean of times that percolates
ggplot(MyResults, aes(x=probability, y=percolation)) +
  geom_point() + 
  scale_shape_manual(values=c(3)) +
  xlab("Probability of connection") + 
  ylab("% of times that percolates") +
  theme_bw()
# Graph 3: probability of connection p versus correlation length g(r=1)
ggplot(MyResults, aes(x=MyResults$`g(r)`, y=percolation) ) +
  geom_point() + 
  scale_shape_manual(values=c(3)) +
  xlab("Correlation Length g(r=1)") + 
  ylab("% of times that percolates") +
  theme_bw()
# Graph 4: mean degree <k> versus <s> mean of size of clusters
ggplot(MyResults, aes(x=k, y=S) ) +
  geom_point() + 
  scale_shape_manual(values=c(3)) +
  xlab("<k>") + 
  ylab("<S>") +
  theme_bw()


# http://rpubs.com/MarkusLoew/226759  ---> para grafica con dos ejes.
# Con dos ejes al mismo tiempo:
p <- ggplot(MyResults, aes(x=probability, y=p50)) +
  geom_errorbar(aes(ymin=p25, ymax=p75), width=.007) + 
  geom_point() + 
  scale_shape_manual(values=c(3)) 
p <- p + geom_point( aes(x=probability, y=percolation*100, colour="red") )
p <- p + scale_y_continuous(sec.axis = sec_axis(~./100, name="% of times that percolates" ) )
#p <- p + scale_colour_manual(values=c("black", "red"))
p <- p + labs(y = "Cluster size", x = "Probability of connection")
#p <- p + theme(legend.position="none")
p <- p + theme_bw() + theme(legend.position="none")
p <- p + scale_x_continuous(breaks=seq(0.1,1,0.1))
p
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 




# LO QUE FALTA.
# * lo mismo pero para L=5, L=15
# lo mismo para redes de small world. como percolan?

# * falta calcular el exponente gamma_p (ver pag. 274 de netqork science de barabasi). En 
# este libro tambien calcular P_inf que es la probabilidad que un nodo pertenezca al 
# cluster mas grande.
# * Del paper Percolation on sparse networks de Kerri, me falta calcular el pi_i(s):
# la probabilidad de que un nodo i, pertenezca a un cluster de tamano s.
