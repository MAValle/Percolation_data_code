# Function to create a lattice with random connections from a lattice with
# weighted conenctions
# Date: 030419

# name of teh file: create_random_Wlattice_net_function.R

# Note: 
# 18-mar-19: 
#             

# Inputs: xy = layout of the lattice
#         edg = edgelist with the weights
#         p = probability of conection of edges
#         L = size of the lattice

# outputs: new_g = igraph object

create_random_Wlattice_net <- function(xy, edg, p, L) {
  # edg$weight es la prob de que dos nodos se conecten. Si p<edg$weight, entonce se conectan.
  link <- ifelse(p < edg$weight, 1, 0) # 1 es conectado, 0 es desconectado
  edg <- cbind(edg, link) 
  edg2 <- edg[edg[,"link"]==1,]
  new_edg <- edg2[,c(1,2)]
  new_g <- graph_from_edgelist(as.matrix(new_edg), directed = FALSE)
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
#S <- create_random_Wlattice_net(xy, edg, p=0.5, L=L)
# new_g <- S[[1]]
#plot(new_g, layout=xy,
#     vertex.size=4,
#     vertex.label.cex = 0.5,
#     vertex.label="",
#     edge.width=4)