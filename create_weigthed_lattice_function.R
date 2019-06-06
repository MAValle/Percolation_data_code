# Function to create a lattice of size L X L with random weigths 
# according to a specific distribution.

# Date: 300319

# recursos:
# https://www.howardrudd.net/how-tos/triangular-distribution/
# http://personality-project.org/r/distributions.html


# Para simular distintos valores de weigths (que serian los Jik 
# de un ising model) lo que hacemos es tomar valores de una uniforme
# entre -1 y 1; o  de una normal con media 0 y sigma=0.5 (asi 
# el 95.4% de los datos estaran entre -1 y 1).



# las critical probabilities asociadas a cada weigth, se calculan 
# considerando que P = K*exp(weight) (considerando modelo ising)
# y la constante K = 1/sum [exp(weight)]
# No obstante aqui haremos aqui una simulacion directa de probabilidades.
# La primera simplemente uniforme(0,1).
# https://www.rdocumentation.org/packages/truncnorm/versions/1.0-8/topics/truncnorm (normal truncada)
# La segunda de una normal truncada con parametros a=0, b=1, mode=0.5, sd=0.1
# La tercera de una normal truncada con parametros a=0, b=1, mode=0.2, sd=0.1
# La cuarta de una normal truncada con parametros a=0, b=1, mode=0.8, sd=0.1
# La probabilidad indica la probabilidad de que dos vertices se conecten.


# Inputs: 
# L es el tamano del lattice
# dist = indica si es 'unif' U(0,1) o 'tria' T(ma,mi,mo)
# Output:
# xy (coordenadas), edg (edge list), g (objeto igraph)

# para normal truncada
#library(truncnorm)
#y <- rtruncnorm(n=10000, a=0, b=1, mean = 0.5, sd = 0.2)



create_lattice_weighted <- function(L, dist, moda, mi, ma) {
  library(igraph)
  library(truncnorm)
  g <- make_lattice(length = L, dim = 2) # lattice de 5 X 5 nodos.
  xy <- layout_on_grid(g, dim = 2)
  # add weigths to the edges
  if (dist == "unif") {
    E(g)$weight <- runif(ecount(g), mi, ma)
  } else {
    E(g)$weight <- rtruncnorm(n=ecount(g), a=mi, b=ma, mean = moda, sd = 0.2)
    #E(g)$weight <- rnorm(ecount(g), mean=0, sd=0.5)
  }
  # get the critical probabilities
  #K <- 1/sum(exp(E(g)$weight))
  #E(g)$critical_probs <- K*exp(E(g)$weight)
  
  # get the adjacency list
  # https://www.r-graph-gallery.com/248-igraph-plotting-parameters/
  edg <- get.edgelist(g, names=TRUE)
  # the adjacency list with their weights and critial probabilities
  edg2 <- get.data.frame(g)
  edg2$critial_prob <- E(g)$critical_probs
  # matriz de incidencia o adyacencia con los pesos
  #get.adjacency(g, attr="weight")
  output <- list(xy, edg2, g) 
  return(output)
}
# example
#lat <- create_lattice_weighted(L=10, dist='norm', moda=0.5, mi=0, ma=1)
#xy <- lat[[1]]
#edg <- lat[[2]]
#g <- lat[[3]]
#plot(g, layout=xy,
#     vertex.size=4,
#     vertex.label.cex = 0.5,
#     vertex.label="",
#     edge.width=4)