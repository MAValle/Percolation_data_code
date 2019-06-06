
# Esta funcion recibe edg y g de un complete net (g) y lo que hace es
# dejar conectado los edge entre nodos que solo son vecinos y que esten activados (spin on).

# 05-may-19

# Funcion
# inputs:
#	  edg: edge list matrix the comes from function create_fullnet 
# 	g: igraph object that was created with create_fullnet 
# output:
#   new_edg: edge list matrix with two columns v1 and v2
#   new_g: igraph object 
create_new_fullnet_bonds <- function(edg, g) {
  rw <- nrow(edg)
  vec1 <- numeric(rw)
  vec2 <- numeric(rw)
  for (i in seq_along(vec1) ) {
    vec1[i] <- V(g)$state[which(V(g) == edg[i,1])]
    vec2[i] <- V(g)$state[which(V(g) == edg[i,2])]
  }
  link <- ifelse(vec1+vec2 == 2, 1, 0)
  edg <- cbind(edg, link, Ji, J0) # vertex1, vertex2, conextion 1=connected, 0 not connected # 28-abr-19
  #edg <- cbind(edg, link) # vertex1, vertex2, conextion 1=connected, 0 not connected
  edg2 <- edg[edg[,"link"]==1,]
  new_edg <- edg2[,c(1,2)]
  new_g <- graph_from_edgelist(new_edg, directed = FALSE)
  number_nodos_que_faltan <- num_vertex - max(V(new_g))
  if (number_nodos_que_faltan > 0) {
    new_g <- add_vertices(new_g, number_nodos_que_faltan )
  } 
  V(new_g)$str <- V(g)$str # como el lattice queda constante en numero de nodos, le asignamos los strength de la original.
  E(new_g)$coupl <- edg2[,"Ji"] # 28-abr-19 original: E(new_g)$coupl <- edg2[,"J0"]
  E(new_g)$Je <- edg2[,"J0"] # 28-abr-19 
  #V(new_g)$str <- strength(new_g, vids = V(new_g), weights = E(new_g)$coupl)
  return(list(new_g, new_edg))
}
# ejemplo:
# new_fullnet <- create_new_fullnet_bonds(edg, g)
# new_g <- new_fullnet[[1]]
# new_edg <- new_fullnet[[2]]
# plot(new_g, layout=xy,
#      vertex.size=5,
#      vertex.frame.color = "white",
#      vertex.color = V(g)$color,
#      vertex.label.cex = 0.8,
#      vertex.label="",
#      vertex.arrow.mode=0,
#      edge.width=1)

