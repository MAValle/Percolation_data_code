# Function to compute correlation length g(r) en which r=1 for each vertex, 
# for a lattice. 
# g(r=1) es la probabilidad de encontrar a connected vertex at distance r of the reference
# vertex. En un lattice tendremos, para r=1: 


# Inputs: L (size of the lattice) and new_g (net)
# Outputs: value indicating the mean of g(r).
correlation_length <- function(new_g, L) {
  # comenzamos por las esquinas
  G_corners <- vector(mode = "numeric", length = 4)
  corners <- c(1, L, L*L+(1-L), L*L)
  # vemos cuales son los nodos adyacentes de nodos corners
  for (n in seq_along(1:4) ) {
    nodo <- corners[n]
    adyte <- as.numeric(neighbors(new_g, v=nodo))
    if ( is.null(adyte) ) {
      g <- 0
    } else {
      g <- length(adyte)/2
    }
    G_corners[n] <- g
  }
  # seguimos por los nodos contorno
  contorno_in <- seq(1,L,1)
  contorno_in <- contorno_in[-1]
  G_contorno_in <- vector(mode = "numeric", length = length(contorno_in))
  for (n in seq_along(1:length(contorno_in) ) ) {
    nodo <- contorno_in[n]
    adyte <- neighbors(new_g, v=nodo)
    if ( is.null(adyte) ) {
      g <- 0
    } else {
      g <- length(adyte)/2
    }
    #g <- length(adyte)/2
    G_contorno_in[n] <- g
  }
  contorno_sup <- seq(L*L+(1-L),L*L,1)
  contorno_sup <- contorno_sup[-1]
  G_contorno_sup <- vector(mode = "numeric", length = length(contorno_sup))
  for (n in seq_along(1:length(contorno_sup)) ) {
    nodo <- contorno_sup[n]
    adyte <- neighbors(new_g, v=nodo)
    if ( is.null(adyte) ) {
      g <- 0
    } else {
      g <- length(adyte)/2
    }
    G_contorno_sup[n] <- g
  }
  contorno_iz <- seq(1,L*L+(1-L),L)
  contorno_iz <- contorno_iz[-1]
  G_contorno_iz <- vector(mode = "numeric", length = length(contorno_iz))
  for (n in seq_along(1:length(contorno_iz)) ) {
    nodo <- contorno_iz[n]
    adyte <- neighbors(new_g, v=nodo)
    if ( is.null(adyte) ) {
      g <- 0
    } else {
      g <- length(adyte)/2
    }
    G_contorno_iz[n] <- g
  }
  contorno_de <- seq(L,L*L,L)
  contorno_de <- contorno_de[-1]
  G_contorno_de <- vector(mode = "numeric", length = length(contorno_de))
  for (n in seq_along(1:length(contorno_de)) ) {
    nodo <- contorno_de[n]
    adyte <- neighbors(new_g, v=nodo)
    if ( is.null(adyte) ) {
      g <- 0
    } else {
      g <- length(adyte)/2
    }
    G_contorno_de[n] <- g
  }
  # Seguimos con los nodos del interior
  all_nodes <- seq(1,L*L,1)
  remover <- c(corners, contorno_de, contorno_in, contorno_iz, contorno_sup)
  insiders <- all_nodes[! all_nodes %in% remover]
  G_insiders <- vector(mode = "numeric", length = length(insiders))
  for (n in seq_along(1:length(insiders)) ) {
    nodo <- insiders[n]
    adyte <- neighbors(new_g, v=nodo)
    if ( is.null(adyte) ) {
      g <- 0
    } else {
      g <- length(adyte)/2
    }
    G_insiders[n] <- g
  }
  G <- c(G_corners, G_contorno_in, G_contorno_sup, G_contorno_iz, G_contorno_de, G_insiders )
  final_output <- mean(G)
  return(final_output)
}
# ejemplo:
# correlation_length(new_g, L)
