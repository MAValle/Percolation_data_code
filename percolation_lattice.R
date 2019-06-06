# Este script lo que hace es crear una red lattice de L X L a la cual simulamos
# bond percolation, es decir con cierta aprobabilidad p, hacemos que 
# los espacios se unan con edges con probabilidad p. Luego analizamos el tamano
# de los clusters que se forman.

# El segundo paso, es hacer lo mismo, con la diferencia es que los edges que
# hay en cada espacio, tienen grosores distintos. Mayor grosor indica que 
# que la probabilidad de union de esos espacio es mayor.


# Creation date: 25-feb-19
# file name: percolation_lattice.R

# Notes:
# 28-feb-19: el primer intento para ver si el sistema percola fracasa. Lo deje
#           inclompleto debido a que el codigo se complica.
# 01-mar-19: comienzo 2do intento para programar si el sistema percola, esta vez
#           utilizando nocion de nodos adyacentes.
# 05-mar-19: finalmente logre programa la deteccion de percolation. Todas
#           las funciones requeridas se pasaron a un archivo R cuyo nombre
#           es el nombre de la respectiva funcion: por ejemplo:
#           create_lattice_function.R, create_random_lattice_net_function.R, 
#           percolation_detection_function.R.
# 17-mar-19: comenzamos a probar un codigo para calcular el g(r) o correlation
#           length.
# 27-mar-19: intentamos hacer un heatmap de la matriz de adyacencia de new_g

# Primero tenemos que crear un lattice de L X L
# Luego se unen los espacios de lattice con probabilidad p
# luego contamos el tamano y numero de clusters
# Un cluster se forma cuando dos espacioes contiguos quedan unidos

#https://www.r-bloggers.com/percolation-threshold-on-a-square-lattice/
rm(list=ls())
library(ggplot2)
library(igraph)

# # # # # # # # # # # # # # # # # # Creation of the lattice net # # # # # # # # # # # # # # # # # #
# create_lattice_function.R
L <- 100
g <- make_lattice(length = L, dim = 2) # lattice de 5 X 5 nodos.
xy <- layout_on_grid(g, dim = 2)
# https://www.r-graph-gallery.com/248-igraph-plotting-parameters/
plot(g, layout=xy,
     #vertex.size=20,
     vertex.size=0.5,
     #vertex.label.cex = 0.8,
     vertex.label="",
     edge.width=4)
####adj <- get.adjacency(g, type="both", sparse=F)
#A solo los elementos de adj que son 1, podemos hacer que permanezcan en 1 (conectado)
# o desconectarlos (igual a cero) Esto lo podemos hacer con un random number u. 
# if u >0.5, entonces dejamos el 1, si no, lo dejamos en 0 (desconectado). 
# En este caso estamso usando p=0.5
edg <- get.edgelist(g, names=TRUE)
V(g)
E(g)
# # # # # # # # # # # # # # # # # # Creation of the lattice net # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # ## # # # Creation of random lattice net # # # # # # # # # # # # # # # # # #
# create_random_lattice_net_function.R
# nos podemos dar cuenta en en edg aparecen los nodos que solo estan conectados.
# Con estos podemos definir si se conectan o no y a partir de aqui, 
colnames(edg) <- c("V1", "V2")
p <- 0.25
U <- runif(nrow(edg))
link <- ifelse(U < p, 1, 0) # 1 es conectado, 0 es desconectado
edg <- cbind(edg, U, link) 
# ahora removemos de edg, las filas con link = 0
edg2 <- edg[edg[,"link"]==1,]
new_edg <- edg2[,c(1,2)]
# ahora a partir de new_edg, actualizamos la red g
new_g <- graph_from_edgelist(new_edg, directed = FALSE)
plot(new_g, layout=xy,
     #vertex.size=20,
     vertex.size=0.5,
     #vertex.label.cex = 0.8,
     vertex.label="",
     edge.color = "red",
     edge.width=4)
# # # # # # # # # # # # ## # # # Creation of random lattice net # # # # # # # # # # # # # # # # # #



# Ahora contar numero de cluster y tamano de cada cluster
https://igraph.org/r/doc/components.html
https://igraph.org/r/doc/decompose.html
http://cneurocvs.rmki.kfki.hu/igraph/doc/R/clusters.html
https://sites.fas.harvard.edu/~airoldi/pub/books/BookDraft-CsardiNepuszAiroldi2016.pdf

clusters(new_g)

# procedimeinto
# * observables: S = mean cluster size, 
# * PI(L,p)=probability that a lattice of
# linear dimension L percolates at conceptration p (percolate is when at least one cluster
# connects the top line or plane with bottom line or plane). Otra forma de definirlo:
# probability that a random n-by-n system with site vacancy probability p percolates. 
# We refer to this quantity as the percolation probability. 
# To estimate its value, we simply run a number of experiments.

# * p_av =average concentration at which for the first time a percolating cluster
# connects top and bottom of the cluster. Esto se puede hacer utilizando a
# recursive scheme known as depth-first search.

# * length correlation = https://en.wikipedia.org/wiki/Correlation_function_(statistical_mechanics)

# * Adaptive plot:  plots the percolation probability as a function of the site 
# vacancy probability p for an n-by-n system. It uses a recursive approach that 
# produces a good-looking curve at relatively low cost: https://introcs.cs.princeton.edu/java/24percolation/

# Recursos:
# http://kateto.net/networks-r-igraph
# https://en.wikipedia.org/wiki/Depth-first_search
# https://introcs.cs.princeton.edu/java/24percolation/


# Procedimiento:
# primero tomamo p=1/2 y vemos if liquid percolate. 
# si lo hace, decrease p un 1/4, de otra forma increase p a un 1/4.
# repeat this division until con suficiente accuracy we get the concentration
# p = p_av e que aparece the first spanning cluster (para L grande p_c -> p_av)
# para cada una de las repeticiones, calcular la media de S (tamano del cluster).


save.image(file='data_testing_040319_funcionando.RData')
library(prodlim)
load('data_testing_040319.RData')


# # # # # # # # # # # # # # # # # # percolation detection # # # # # # # # # # # # # # # # # #
# percolation_detection_function.R
#Mar2019
# Get node from a graph new_g (debe ser del 1 al 10 o del 1 a L)
net <- new_g
V(new_g)$name <-V(new_g)
V(net)$name <- V(new_g)
intento = 0
vertex_sink <- seq(L*L-(L-1),L*L, by=1) #estos son los nodos sink (salida de la percolacion)
get_out_of_the_loop <- FALSE
msg <- "No percolation has been detected"
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
    print(intento)
    print(stack)
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
      break
    }
  }#fin while
  
}# fin codigo 
sprintf("%s", msg)
# # # # # # # # # # # # # # # # # # percolation detection # # # # # # # # # # # # # # # # # #


# 17-mar-19
# # # # # # # # # # # # ## # # Creation of correlation length # # # # # # # # # # # # # # # #
# Creamos un codigo para calular g(r) o correlation length r, en donde r=1, es decir
# calculamos la probabilidad de encontrar a connected vertex at distance r of the reference
# vertex. En un lattice tendremos, para r=1: 
# 1) los vertex de las esquinas con dos posibles conexiones: 1, L, L*L+(1-L) y L*L
# 2) los vertex del contorno del lattice con 3 posibles conexiones: 
#   base inferior: seq(1,L,1)
#   base superior: seq(L*L+(1-L),L*L,1)
#   contorno izquierdo: seq(1,L*L+(1-L),L)
#   contorno derecho: seq(L,L*L,L)
# 3) los vertex del interior del lattice con 4 posibles conexiones
# Procedimiento:
# Nota: se asume que la red new_g (cpercolada o no) ya esta cargada con los nombres de nodos
# 1. para cada nodo, verificar si es nodo del lattice tipo 1, 2 o 3 y luego encontrar sus vecinos 
# 2. de los nodos que son vecinos, encontrar los que estan conectados
# 3. calcular g(r) para ese nodo como la suma de los nodos conectados dividido por el numero total de vecinos
# Funcion:
# Input: L, tamano del lattice
#        new_g: objeto igraph con la red percolada o no percolada
# Output: media de g(r)  <g(r=1)>
correlation_length <- function(new_g, L) {
  # comenzamos por las esquinas
  G_corners <- vector(mode = "numeric", length = 4)
  corners <- c(1, L, L*L+(1-L), L*L)
  # vemos cuales son los nodos adyacentes de nodos corners
  for (n in seq_along(1:4) ) {
    nodo <- corners[n]
    adyte <- neighbors(new_g, v=nodo)
    adyte <- adyte$name
    if ( is.null(adyte) ) {
      g <- 0
    } else {
      g <- length(adyte)/2
    }
    G_corners[n] <- g
  }
  # seguimos por los nodos contorno
  contorno_in <- seq(1,L,1)
  contorno_in <- contorno_in[-1]
  G_contorno_in <- vector(mode = "numeric", length = length(contorno_in))
  for (n in seq_along(1:length(contorno_in) ) ) {
    nodo <- contorno_in[n]
    adyte <- neighbors(new_g, v=nodo)
    adyte <- adyte$name
    if ( is.null(adyte) ) {
      g <- 0
    } else {
      g <- length(adyte)/2
    }
    #g <- length(adyte)/2
    G_contorno_in[n] <- g
  }
  contorno_sup <- seq(L*L+(1-L),L*L,1)
  contorno_sup <- contorno_sup[-1]
  G_contorno_sup <- vector(mode = "numeric", length = length(contorno_sup))
  for (n in seq_along(1:length(contorno_sup)) ) {
    nodo <- contorno_sup[n]
    adyte <- neighbors(new_g, v=nodo)
    adyte <- adyte$name
    if ( is.null(adyte) ) {
      g <- 0
    } else {
      g <- length(adyte)/2
    }
    #g <- length(adyte)/2
    G_contorno_sup[n] <- g
  }
  contorno_iz <- seq(1,L*L+(1-L),L)
  contorno_iz <- contorno_iz[-1]
  G_contorno_iz <- vector(mode = "numeric", length = length(contorno_iz))
  for (n in seq_along(1:length(contorno_iz)) ) {
    nodo <- contorno_iz[n]
    adyte <- neighbors(new_g, v=nodo)
    adyte <- adyte$name
    if ( is.null(adyte) ) {
      g <- 0
    } else {
      g <- length(adyte)/2
    }
    #g <- length(adyte)/2
    G_contorno_iz[n] <- g
  }
  contorno_de <- seq(L,L*L,L)
  contorno_de <- contorno_de[-1]
  G_contorno_de <- vector(mode = "numeric", length = length(contorno_de))
  for (n in seq_along(1:length(contorno_de)) ) {
    nodo <- contorno_de[n]
    adyte <- neighbors(new_g, v=nodo)
    adyte <- adyte$name
    if ( is.null(adyte) ) {
      g <- 0
    } else {
      g <- length(adyte)/2
    }
    #g <- length(adyte)/2
    G_contorno_de[n] <- g
  }
  # Seguimos con los nodos del interior
  all_nodes <- seq(1,L*L,1)
  remover <- c(corners, contorno_de, contorno_in, contorno_iz, contorno_sup)
  insiders <- all_nodes[! all_nodes %in% remover]
  G_insiders <- vector(mode = "numeric", length = length(insiders))
  for (n in seq_along(1:length(insiders)) ) {
    nodo <- insiders[n]
    adyte <- neighbors(new_g, v=nodo)
    adyte <- adyte$name
    if ( is.null(adyte) ) {
      g <- 0
    } else {
      g <- length(adyte)/2
    }
    #g <- length(adyte)/2
    G_insiders[n] <- g
  }
  G <- c(G_corners, G_contorno_in, G_contorno_sup, G_contorno_iz, G_contorno_de, G_insiders )
  final_output <- mean(G)
  return(final_output)
}
# ejemplo:
correlation_length(new_g, L)
# # # # # # # # # # # # ## # # Creation of correlation length # # # # # # # # # # # # # # # #




# # # # # # # # # # # # # # # #  # percolation path graph # # # # # # # # # # # # # # # # # #
#040319
# Tratando de graficar el percolation path
# # utilizando all_simple_paths:
# OJO: all_simple_paths toma mucho tiempo cuando L=20, de hecho se queda sin memoria.
#percolation_path <- all_simple_paths(new_g, stack[1], stack[length(stack)]) #me da el path con los nombres de los vertices
#percolation_path_names <- percolation_path[[4]] #nombres de los vertex del path de percolation
# # Utilizando get.shortest.paths: (este es mucho mas rapido!)
percolation_path_names <- get.shortest.paths(new_g, from=stack[1], to=stack[length(stack)], mode = "out")
percolation_path_names <- percolation_path_names$vpath[[1]]
#percolation_path_names <- stack

edg <- get.edgelist(g, names=TRUE) ####### Importante: edgelist original del lattice
#necesito dejar solo los edges que conectan cada vertice de percolation_path_names, que se llamara edg3
g_ring <- graph.ring(length(percolation_path_names))
V(g_ring)$name <- percolation_path_names
g_ring <- delete_edges(g_ring, E(g_ring)[length(E(g_ring))]) #sacamos el ultimo edge
edg_of_ring <- get.edgelist(g_ring, names=TRUE) ###### Importante: edgelist of percolation path

# En el vector link vamos poniendo si el edge de edg se retira o se mantiene
link <- rep(0, nrow(edg)) # 1 es conectado, 0 es desconectado
for ( i in 1:nrow(edg)) {
  edge_is <- c(edg[i,1], edg[i,2]) # sacamos un edge de edg (el lattice)
  for (j in 1:nrow(edg_of_ring)) { # vamos sacando edges del edg_of_ring
    edge_in <- c(edg_of_ring[j,1], edg_of_ring[j,2])
    if ( sum(edge_is %in% edge_in) == 2 ) {
      link[i] = 1
      break
    }
  }
}

edg <- cbind(edg, link) 
# ahora removemos de edg, las filas con link = 0
edg3 <- edg[edg[,"link"]==1,]
edg4 <- edg3[,c(1,2)]

percol_g <- graph_from_edgelist(edg4, directed = FALSE)
plot(percol_g, layout=xy,
     #vertex.size=12,
     vertex.size=4,
     vertex.label.cex = 0.5,
     vertex.label="",
     edge.width=4,
     edge.color="red")
# # # # # # # # # # # # # # # #  # percolation path graph # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # ## # # # Heatmap of adjacency matrix # # # # # # # # # # # # # # # # # # #
# 27mar19
ady_M <- as_adjacency_matrix(new_g, type="both")
library(ggplot2)
ady_M <- as.matrix(ady_M)
# alternativa 1
image(ady_M, xlab = 'Matrix rows', ylab = 'Matrix columns') # demora mucho tiempo
# alternativa 2
#Build the palette and plot it
library(lattice)
# https://www.rdocumentation.org/packages/dichromat/versions/1.1/topics/colorRampPalette
pal <- colorRampPalette(c("yellow", "black"), space = "rgb")
system.time(
  levelplot(ady_M, main="Adjacency matrix", xlab="", ylab="", col.regions=pal(4), cuts=3, at=seq(0,1,0.5))
) # 40.3 sec.
# alternativa 3
# https://cran.r-project.org/web/packages/heatmap3/vignettes/vignette.pdf
# https://stackoverflow.com/questions/5667107/how-can-i-make-a-heatmap-with-a-large-matrix
# https://rdrr.io/cran/heatmap3/man/heatmap3.html
library(heatmap3)
system.time(
heatmap3(ady_M, useRaster=TRUE, Rowv=NA,Colv=NA, labRow="", labCol="")
)
#con una matriz de 5000X5000 se demora 7 s.
# # # # # # # # # # # # ## # # # Heatmap of adjacency matrix # # # # # # # # # # # # # # # # # # #




































# Esto no funciono..
# Funcion para buscar e identificar una conexion al nodo v
# Input: edgelist de la red con columna V1 y V2
# Input: v que es el nombre del nodo a buscar
# Output: * FALSE en caso que nodo v no tenga conexion
#         * si v tiene conexion retorna vector con nombre de los nodos a los cuales esta conectado
search_connection <- function(edgelist, v){
  # vector de salida:
  salida <- vector(mode="numeric", length=0)
  # identificamos filas en que aparece el vactor v
  #id <-  which(apply(edgelist, 1, function(r) any(r %in% v)))
  id <- which(edgelist[,1] == v, arr.ind=TRUE)
  if ( is.integer(id) && length(id) == 0L ) {
    salida <- integer(0)
  } else {
    salida <- id
  }
  return(salida)
}
# Ejemplo_
# search_connection(edgelist=new_edg, v=2)

# # # # # # # # # # # # Intento para ver si el sistema percola # # # # # # # # # # # #
# Feb2019
######### Initial settings
# chequeo de percolation
perc_edg <- matrix(NA, nrow=nrow(new_edg), ncol=2) # aqui vamos cuardando los edges de la percolacion
# con dupl_edg vamos a ir trabajando. Solo nos interesa la percolacion
# cuando parte de un nodo de la primera fila.
dupl_edg <- new_edg
v <- dupl_edg[2,1] #comenzamos con primer nodo conectado en new_edg


#############---falta otro while  cccc que vuelve de L186
test <- vector(mode="numeric", length=nrow(new_edg)) #1=hay que borrar, 0 no hay que borrar (va de la mano con dupl_edg)
#first_row_edg <- new_edg[new_edg[,1] <= 10,] #10 corresponde a L (tamano del lattice)
first_row_id <- which(new_edg[,1] <= 10)
fila <- 1
n <- 0
# vamos chequeando si un nodo tiene mas de un edge, va de la mano id con id con perc_edg
mas_de_un_edge <- vector(mode="numeric", length=nrow(new_edg)) 
idx <- 1 #cualquier numero distinto de cero
n <- length(perc_edg[!is.na(perc_edg)]) + 1
while ( length(idx) != 0) {
  #n <- n + 1
  idx <- search_connection(edgelist=dupl_edg, v=v) # nos entrega el id de dupl_edg de los edges que tiene el nodo v
  mas_de_un_edge[n] <- length(idx) # para testear en caso que el nodo v tenga mas de un edge
  if ( length(idx) != 0 ) {
    perc_edg[n, ] <- dupl_edg[idx[1], ] # guardamos el edge seleccionado que supuestamente percola
    v <- dupl_edg[idx[1], 2]
    print(n)
    print(perc_edg[n, ])
  }
  n <- n + 1
}
#nos salimos el while porque ya no hay mas edges
perc_edg <- perc_edg[complete.cases(perc_edg), ,drop=FALSE] # Dejamos todos los edges encontrados a partir del nodo v
mas_de_un_edge[mas_de_un_edge < 2] <- NA
mas_de_un_edge <- mas_de_un_edge[complete.cases(mas_de_un_edge )]
if ( length(mas_de_un_edge) == 0 ) {
  #encontrar el edge en  perc_edg en dupl_edg para marcarlo como fallido con library(prodlm)
  if ( nrow(perc_edg) == 1 ) {
    fnd <- row.match(perc_edg[1,], dupl_edg[,1:2])
    test[fnd] <- 1 
  } else {
    fnd <- vector(mode="numeric", length=0)
    for (i in 1:nrow(perc_edg)) {
      fnd <- c(fnd,row.match(perc_edg[i,], dupl_edg[,1:2]))
      test[fnd] <- 1 
    }
  }
  perc_edg <- matrix(NA, nrow=nrow(new_edg), ncol=2) # parte perc_edg desde cero
} else {
  elimino_desde <- length(mas_de_un_edge) + 1
  drop_edg <- perc_edg[c(elimino_desde:nrow(perc_edg)), ]
  perc_edg <- perc_edg[-c(elimino_desde:nrow(perc_edg)), ]
  #encontrar estos edges de drop_edg en dupl_edg para marcarlo como fallido con library(prodlm)
  fnd <- vector(mode="numeric", length=0)
  for (i in 1:nrow(drop_edg)) {
    fnd <- c(fnd,row.match(drop_edg[i,], dupl_edg[,1:2]))
  }
  test[fnd] <- 1 
}
#veamos que dupl_edg hay que eliminar edges que no percolan, utilizando el indicador de test
quita <- which( test == 1 )
dupl_edg <- dupl_edg[-c(quita), ]
perc_edg <- matrix(NA, nrow=nrow(new_edg), ncol=2) #partimos con per_edg de cero.
# # # # # # # # # # # # Intento para ver si el sistema percola # # # # # # # # # # # #
