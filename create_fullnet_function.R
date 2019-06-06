

# Esta funcion creat un red completa. 

# 05-may-19

# step 1: complete graph creation
# Inputs: num_vertex
# outputs: xy = layout of the lattice
#          edg <- edgelist
#          g = igraph object
create_fullnet <- function(num_vertex) {
  # https://rstudio-pubs-static.s3.amazonaws.com/337696_c6b008e0766e46bebf1401bea67f7b10.html
  # https://kateto.net/netscix2016.html
  library(igraph)
  g <- make_full_graph(n=num_vertex) # complete net with num_vertex nodes
  xy <- layout_in_circle(g)
  # https://www.r-graph-gallery.com/248-igraph-plotting-parameters/
  edg <- get.edgelist(g, names=TRUE)
  output <- list(xy, edg, g) 
  return(output)
}
# example
# lat <- create_fullnet(num_vertex)
# xy <- lat[[1]]
# edg <- lat[[2]]
# g <- lat[[3]]
# plot(g, layout=xy,
#      vertex.size=5,
#      vertex.frame.color = "white",
#      vertex.color="orange",
#      vertex.label="",
#      vertex.arrow.mode=0,
#      vertex.label.cex = 0.5,
#      edge.width=1)

