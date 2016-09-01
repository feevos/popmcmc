# popmcmc++
Population mcmc c++ library.

This is a population mcmc c++ library. It incorporates basically 4 algorithms, 3 reliable and tested. 
DE,DES based on differential evolution- Braak et al. (2006)
and GWS - Goodman and Weare 2010, the stretch move.  There is also a new experimental one (PCX), 
inspired from the PCX crossover (Deb et al. 2002). Parallel tempering versions are also provided 
(Algos: DE_PT, DES_PT, GWS_PT, PCX_PT). 

DEPENDS: intel TBB. 

The library is fully parallelized using intel TBB library. There is the option of selecting number of threads, 
or running in single thread mode. 

Jumping in the examples directory should get you started easily. 

The only file that you need to check to understand how the library works is mcmc.hpp, it should be kind (-ish) self explanatory. 

Promise I'll add proper documentation quite soon. 

compile: g++ -std=c++14 main.cpp -o main -ltbb 
