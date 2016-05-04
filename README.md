# mgsolve
Multi Grid solver for elleptic PDEs

In order to run this programm use the following command:
```sh
$ ./mgsolve <number levels> <number cycles>
```
This command executes n V(2,1)-cycles
The program writes the solution of the PDE to file "solution.txt"

#Bonus task
In order to execute the probelm given in the bonus task, use the define macro GHOST in the file mgsolve.cpp and set it to 1 (=true), then recompile and execute. 

This file includes:
  - Makefile
  - mgsolve.cpp
  - Auswertung.pdf




