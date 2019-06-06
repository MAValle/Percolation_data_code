# Function to compute correlation length g(r) in a Network


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

