# Function  companion for percolation_bond_lattice.R.
# Given a lattice with on-spins that are connected 
# between them, we get the names of the largest cluster .

# date: 18-abr-19

# notes:


# FUNCTION to find names of the largest cluster
# find the membersho id of the largest cluster
# input:
# number_of_clusters: numero de cluster en el lattice
find_name_maxcluster <- function(cluster_sizes) {
  vtx <- vector("list", length = cluster_sizes$no)
  for (n in 1:cluster_sizes$no) {
    #m[n] <- length( V(new_new_g)[cluster_sizes$membership == n] )
    vtx[[n]] <- V(new_new_g)[cluster_sizes$membership == n]
  }
  library(purrr)
  idlist <- which.max(map_dbl(vtx, length))
  vertices_clustermax <- vtx[[idlist]]
  return( vertices_clustermax )
}
# ejemplo
# cluster_sizes <- clusters(g)
# vertices_names <- find_name_maxcluster(cluster_sizes)