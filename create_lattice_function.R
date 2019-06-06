# Function to create a lattice of size L X L
# Date: 040319

# Inputs: L
# outputs: xy = layout of the lattice
#          edg <- edgelist
#          g = igraph object

create_lattice <- function(L) {
  library(igraph)
  g <- make_lattice(length = L, dim = 2) # lattice de 5 X 5 nodos.
  xy <- layout_on_grid(g, dim = 2)
  # https://www.r-graph-gallery.com/248-igraph-plotting-parameters/
  edg <- get.edgelist(g, names=TRUE)
  output <- list(xy, edg, g) 
  return(output)
}
# example
#lat <- create_lattice(L=10)
#xy <- lat[[1]]
#edg <- lat[[2]]
#g <- lat[[3]]
#plot(g, layout=xy,
#     vertex.size=4,
#     vertex.label.cex = 0.5,
#     vertex.label="",
#     edge.width=4)