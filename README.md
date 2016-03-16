# multiband_GWLAN

Ciao!

Questa repository contiene la mia versione del simulatore per multiband GWLAN, rispetto al draft del 4 febbraio. 

Il tool viene eseguito chiamando lo script 'launcher.sh' , che per ogni istanza di ogni configurazione (queste ultime si possono consultare in 'gwlan-configs.dat') invoca "instance_maker.cc". 
Il generatore di istanze genera un file, il quale viene letto da "algo_main.cc" , che inizializza e chiama le procedure del solver (definite in "gwlan_lib.cc") per tale istanza, per poi salvare le soluzioni in un apposito file.
Infine viene chiamato "fetch_statistics.sh" dalla cartella "scripts" per generare la tabella riassuntiva.


###
TUTTI I FILE NECESSARI A ESEGUIRE IL PROGETTO SONO INCLUSI NELLA REPOSITORY. Compilare chiamando il comando "make" sulla shell bash per generare i files binari (occhio ai path di 'CPLEX' e 'concert' specificati nel Makefile).
###

