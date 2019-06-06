# Function to compute the energy of a cluster in a lattice


# 19-abr-19
# Funcion
# calculo de la energia de un cluster
# input:
# cl: numero de cluster
# members: vector con el id identificatorio del cluster al cual pertenece el vertice
# output:
# energy: energia del cluster = E = -2*sum(todos los acoples de los vertices conectados)
# nota que en cada pasada por nodo, tomo doble cada edge, por lo que al dividir por dos, 
# solo nos queda que E = sum(acoples de los vertices conectados.)
energy_cluster <- function(cl, members) {
  vert_ids <- which(members == cl)
  energy <- 0
  for (v in seq_along(vert_ids)) {
    vertice <- vert_ids[v]
    edgess <- E(new_new_g)[from(vertice)]
    acoples <- E(new_new_g)$coupl[edgess]
    partial_energy <- -sum(acoples)
    energy <- energy + partial_energy 
  }
  return(energy)
}
# number_of_clusters <- length(cluster_sizes$csize)
# ejemplo: E <- energy_cluster(cl=15, members = cluster_sizes$membership)
# para calcular la energia de todos los clusters
#clid <- 1:length(cluster_sizes$csize)
#all_energies <- map_dbl(clid, energy_cluster, members = cluster_sizes$membership)
#drop energies equal to zero (cluster with size=1)
#all_energies <- all_energies[all_energies>0]
#mean(all_energies)
#max(all_energies)


# 28-abr-19
# Designed for percolation_bond_lattice_v4a.R en lattice
# This function computes the energy of the energy of the cluster same as energy_cluster function,
# but it includes the J < 0 attributed of the edges which is located in E(new_g)$Je
energy_cluster_e <- function(cl, members) {
  vert_ids <- which(members == cl)
  if (length(vert_ids) < 4) { # ver pag. 208 libro de apuntes.
    energy <- NA
  } else {
    ng2 <- induced_subgraph(new_g, vert_ids)
    energy <- -1*sum(E(ng2)$Je)
  }
  return(energy)
}



# 12-amy-19
# Designed for percolation_bond_fullnet.R
# This function computes the energy of the energy of the cluster same as energy_cluster function,
# but it includes the J < 0 attributed of the edges which is located in E(new_g)$Je
energy_cluster_f <- function(cl, members) {
  vert_ids <- which(members == cl)
  if (length(vert_ids) == 1) { # si se trata de cluster de size = 1, la energia es nula.
    energy <- 0
  } else {
    ng2 <- induced_subgraph(new_g, vert_ids)
    energy <- -1*sum(E(ng2)$Je)
  }
  return(energy)
}