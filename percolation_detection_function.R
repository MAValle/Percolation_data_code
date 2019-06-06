# Function to detect percolation in random lattice network
# The percolation is when there is (are) some path from any of the
# vertexs of the first row (lower lattice) to any vertex of the upper lattice.


# Inputs: L and new_g
# outputs: state_of_perc : TRUE or FALSE: TRUE a percolation detected.
#          msg: message  
#          stack : list os visited vertex
percolation_detection <- function(new_g, L) {
  net <- new_g
  V(new_g)$name <-V(new_g)
  V(net)$name <- V(new_g)
  intento = 0
  vertex_sink <- seq(L*L-(L-1),L*L, by=1) #estos son los nodos sink (salida de la percolacion)
  get_out_of_the_loop <- FALSE
  msg <- "No percolation has been detected"
  state_of_perc <- FALSE
  while ( (L %in% V(net)$name) && (get_out_of_the_loop == FALSE) ) { #meintras el nodo L=10 este presente en la red, seguimos
    #while ( intento < 30 ) {
    v <- V(net)$name[1] # Initial node of the modified net
    stack <- vector(mode="numeric", length=0) # aqui vamos poniendo y quitando los nodos por los que viajamos
    stack <- c(stack, v)
    #intento <- intento + 1
    #m1 <- sprintf("Search number %i", intento)
    #print(stack)
    #cat(m1)
    #while ( length(stack) != 0 && intento < 30 ) {
    while ( (length(stack) != 0) && (L %in% V(net)$name) ) {  
      intento <- intento + 1
      #print(intento)
      #print(stack)
      v <- stack[length(stack)] #asignamos a v el ultimo vertice agregado a stack
      #ady <- adjacent_vertices(net, v)
      ady <- neighbors(net, which(V(net)$name==v)) #nos de el nombre de los nodos adyacentes
      ady <- ady$name #como net tiene atributo llamado "name", el resultado de neigbors hereda este atributo
      n <- length(ady)
      if ( n == 0 ) { #hay que sacar el nodo v de stack y borrarlo como nodo de net
        id <- which(stack == v)
        stack <- stack[-id]
        net <- delete_vertices(net, which(V(net)$name==v))
      } else if ( n == 1 ) { # # # # # # # # # # # # # # # # # # # # # # # #
        #id <- which(stack == V(net)$name[ady]) #vemos si el nodo adyacente esta en stack
        id <- which(stack == ady[1] ) 
        #id <- which(stack == ady )
        #Si elemento de ady no esta en stack, ponerlo en el stack
        if ( length(id)==0 ) {
          stack <- c(stack, ady[1] )
          #stack <- c(stack, ady)
        } else {
          #Si elemento de ady esta en stack entonces. v hay sacarlo del stack  y borrarlo de la red net)
          # En este caso, v corresponde a un nodo terminal y hay que quitarlo de la red
          stack <- stack[-c(which(stack==v))]
          net <- delete_vertices(net, which(V(net)$name==v))
        }
      } else if ( n > 1 ) { # # # # # # # # # # # # # # # # # # # # # # # # 
        #temp <- as.numeric(ady)
        num_nodos_ady_in_stack <- 0
        for (j in 1:n) {
          #buscamos en que posicion de stack esta el vertice de ady
          #id <- which(stack == V(net)$name[temp[j]]) #original
          #id <- which(stack == V(net)[temp[j]]) #original
          #id <- which(stack == ady[j])
          id <- which(stack == ady[j] ) 
          #Si elemento de ady no esta en stack, ponerlo en el stack
          if ( length(id)==0 ) {
            #stack <- c(stack, V(net)$name[temp[j]]) #original
            stack <- c(stack, ady[j] )
          } else {
            # Si elemento de ady esta en stack. dejarlo quieto, (esto ya no: sacarlo y delete ese elemento)
            num_nodos_ady_in_stack <- num_nodos_ady_in_stack + 1
            #stack <- stack[-id]
            #net <- delete_vertices(net, v)
          }
        }#fin del loop
        # Si todos los elementos de ady estan en stack, entonces encontramos un loop ->
        # en ese caso, el nodo v hay que sacarlo de stack y eliminarlo de la red!!
        if ( num_nodos_ady_in_stack == n  ) { #recuerda que n es el numero de nodos adyacentes de v
          stack <- stack[-c(which(stack==v))]
          net <- delete_vertices(net, which(V(net)$name==v))
        } # fin de chequeo de loop
      } #fin else if 
      
      #Tenemos que detectar si el sistema percola o no:
      #Si cualquier nodo del 91 al 100 (del LXL-1 al LXL) aparece en stack, entonces percola
      if ( sum(stack %in% vertex_sink) > 0 ) {
        get_out_of_the_loop <- TRUE
        msg <- "A percolation has been detected!"
        state_of_perc <- TRUE
        break
      }
    }#fin while
    
  }# fin codigo 
  sprintf("%s", msg)
  output <- list(state_of_perc, msg, stack)
  return(output)
}
# example
#S <- percolation_detection(new_g, L)
