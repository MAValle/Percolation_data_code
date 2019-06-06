# Hacemos lo mismo que en percolation_bond_fullnet.R pero esta vez consideramos
# los antibonds en la formacion de la percolacion. Es decir, consideramos los bonds
# con J<0.
# De igial manera que en percolation_bond_fullnet.R  asignamos J randoms entre -1 y 1. 
# a una red completa. Luego podemos aplicar un  umbral para eliminar acoples muy pequeños. 
# (esto ultimo es opcional) Luego conectamos nodos (todos activos) con probabilidad 
# p=1-exp(-2J) para J>0 y con  y p=1-exp(-2|J|)  para J<0. 
# Luego calculamos la calculamos la energia de los clusters y otros observables.

# Procedure: 
# step 1: 
# step 1.1: 
# step 2: 
# step 3: 



# Creation date: 05-jun-19
# file name: percolation_bond_fullnet_V2.R

# Notes:
# 05-jun-19: creation 

 



# # # # # # # # # # # # # # Carga de datos # # # # # # # # # # # # # #
rm(list = ls())
source("create_fullnet_function.R") 
source("create_fullnet_bonds_function.R")
source("energy_cluster_function.R")
#source("correlation_length_function.R")
#source("create_new_lattice_bonds_function.R")
library(ggplot2)
library(igraph)
library(doMC)
library(purrr)
library(latex2exp)

