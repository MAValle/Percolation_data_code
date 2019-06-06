# Function  companion for percolation_bond_lattice.R.
# Given a fully connected lattice, the function drop edges
# of the lattice if they are incident to off-spins. The function
# randomly activate a fraction p of the spins(vertices) of the
# lattice. The edges incident to these spins remains.

# date: 05-abr-19

# notes:


# Input:
#   g: objeto igraph lattice fully connected of size L
#   p = proportion of spins activated 
# output
#   new_g: objeto igraph lattice of size L, with edges removed according
#   to the spins offs.
create_random_spins_activations <- function(g, p) {
  on_spins <- sort(sample(V(g), p*length(V(g))))
  off_spins <- setdiff(V(g), on_spins)
  # # # Borro los edges de g que son adyacentes a los nodos que estan en off_spins
  drop_edges <- vector(mode="numeric", length=0)
  for ( i in seq_along(1:length(off_spins))) {
    #drop_edges <-c(drop_edges, E(g)[inc(V(g)[i])]) 
    drop_edges <-c(drop_edges, E(g)[inc(V(g)==off_spins[i])]) 
  }
  drop_edges <- unique(drop_edges)
  new_g <- delete.edges(g, drop_edges )
  return(new_g)
}
# ejemplo:
#new_g <- create_random_spins_activations(g, p=p)