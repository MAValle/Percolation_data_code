# Function to create a lattice with random connections
# Date: 040319


# Note: 
# 18-mar-19: he agregado el codigo:
#             new_g <- add_vertices(new_g, length(seq(max(new_edg)+1,L*L,1)) )
#             que incorpora al objeto new_g, los nodos que se pierden del objeto g
#             

# Inputs: xy = layout of the lattice
#         edg = layout of the lattice
#         p = probability of NO-connction  beetween  edges (for this case is probability of no-connection)
#         L = size of the lattice

# outputs: new_g = igraph object

create_random_lattice_net <- function(xy, edg, p, L) {
  #colnames(edg) <- c("V1", "V2")
  U <- runif(nrow(edg)) 
  link <- ifelse(p > runif(nrow(edg),0,1), 1, 0) 
  edg <- cbind(edg, link) 
  edg2 <- edg[edg[,"link"]==1,]
  new_edg <- edg2[,c(1,2)]
  new_g <- graph_from_edgelist(new_edg, directed = FALSE)
  # Ahora hay que agregar los nodos que faltan. Si por ejemplo la ultima edge de new_edg
  # es 58 con 59, el objeto new_g tendra hasta el nodo 59.
  number_nodos_que_faltan <- L*L - max(new_edg)
  if (number_nodos_que_faltan == 0) {
    output <- list(new_g) 
  } else {
    new_g <- add_vertices(new_g, number_nodos_que_faltan )
    output <- list(new_g) 
  }
  return(output)
}
#Example
#S <- create_random_lattice_net(xy, edg, p=0.5)
# new_g <- S[[1]]
#plot(new_g, layout=xy,
#     vertex.size=4,
#     vertex.label.cex = 0.5,
#     vertex.label="",
#     edge.width=4)