## Script list for percolation project 

###Date: apr 18, 2019



* **percolation_lattice_v3.R**: percolation on lattice.

  Associated functions:

  * create_lattice in **create_lattice_function.R**

  * create_random_lattice_net in **create_random_lattice_net_function.R**
  * percolation_detection in **percolation_detection_function.R**
  * correlation_length in **correlation_length_function.R**



* **percolation_ERnetwork_v3.R**: percolation on ER net.

  Associated functions:

  * correlation_length in **correlation_length_function.R**



* **percolation_Wlattice.R**: percolation on a weighted lattice.

  Associated functions:

  * create_lattice_weighted in **create_weigthed_lattice_function.R**

  * create_random_Wlattice_net in **create_random_Wlattice_net_function.R**

  * percolation_detection in **percolation_detection_function.R**

  * correlation_length in **correlation_length_function.R**

* **percolation_Wnetwork.R**: percolation on a weighted net (not worked on yet). I think that a percolation on a weighted ER net following the same idea of percolation_Wlattice.R will be the same, with a change in critical percolation probability.



* **percolation_bond_lattice.R**: percolation on a lattice activating the sites or vertexs.

  Associated functions:

  * create_lattice in **create_lattice_function.R**

  * create_random_spins_activations in **create_random_spins_activations_function.R**
  * percolation_detection in **percolation_detection_function.R**
  * correlation_length in **correlation_length_function.R**



* **percolation_bond_lattice_v3.R**: percolation in a lattice activating the spins or vertex and each edge has activation probability to the other active spin.

  Associated functions:

  * create_lattice in **create_lattice_function.R**
  * correlation_length in **correlation_length_function.R**
  * find_name_maxcluster in **find_name_maxcluster_function.R**
  * energy_cluster in **energy_cluster_function.R**



* **percolation_bond_lattice_v4.R**: Similar a percolation_bond_lattice_v3.R, pero se simula con una distribucion de "J" Normal o custom donde algunos nodos tendran alto strength y otros con bajos  strength ($J>0$), y luego simulo la percolacion quitando nodos (dejandolos inactivos en el lattice) desde  el mas alto strength hasta al mas bajo, desde el mas bajo al mas alto, y luego en forma aleatoria. Se conecta con el MST (esto último aún pendiente).

  Associated functions:

  *  create_lattice in **create_lattice_function.R**
  * correlation_length in **correlation_length_function.R**
  * energy_cluster in **energy_cluster_function.R**
  * create_new_lattice_bonds in **create_new_lattice_bonds_function.R**



* **percolation_bond_lattice_v4a.R**: Es exactamente lo mismo que en percolation_bond_lattice_v4.R, pero simulo con J positivos y negativos. Para aquellos bonds en los que $J<0$ hacemos que $J=0$ (es decir, sin conexión) pero sólo para efectos de la simulación. No obstante, la información de los bonds o edges que tienen $J<0$ se preserva para calcular la energía de los clusters tomando en cuenta estos edges con acoples negativos. Esto ocurrirá cuando queden dos nodos que son vecinos  y que sean parte de un mismo cluster (pero no conectaods entre sí), posean acople negativo. 

  Associated functions: Las mismas que en percolation_bond_lattice_v4.R.



+ **percolation_bond_fullnet.R**: Hacemos lo mismo que en percolation_bond_lattice_v4a.R pero esta vez sobre una red completa a la cual le asignamos J randoms entre -1 y 1. Luego podemos aplicar un umbral para eliminar acoples muy pequeños. Luego conectamos nodos (todos activos) solo para aquellos con $J>0$ con probabilidad $p=1-\exp(-2J)$ y calculamos la energia de los clusters y otros observables.

  Associated functions:

  * create_fullnet in **create_fullnet_function.R**
  * create_fullnet_bonds in **create_fullnet_bonds_function.R**
  * energy:cluster_f in **energy_cluster_function.R**  



* **percolation_Wlattice_v2.R**: percolation on a weighted network (network products) under different temperatures. (aun no lo hago pero es la idea 3, pag.190)



* **percolation_bond_fullnet_v2.R**: necesito uno que sea como percolation_bond_fullnet pero que vaya creando los clusters con los antibonds. Vemos que esto es complicado por las razones expuestas en pag. 229-230.



