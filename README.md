# popmcmc++
Population mcmc c++ library.

This is a population mcmc c++ library. It incorporates 4 algorithms, 3 reliable and tested and one experimental. 
Tested algorithms: DE,DES based on differential evolution- Braak et al. (2006), Streens 2002 -  
and GWS (Goodman and Weare 2010), the stretch move.  

Experimental one (PCX), inspired from the PCX crossover (Deb et al. 2002). Parallel tempering versions are also provided 
(Algos: DE,DES,GW_Stretch,PCX,DE_PT, DES_PT, GW_Stretch_PT, PCX_PT). 

Dependencies: intel TBB. 

The library is fully parallelized using intel TBB library. There is the option of selecting number of threads, 
or running in single thread mode, without the use of TBB.  

Jumping in the examples directory should get you started easily. The algorithm outputs variable chains and last column has loglikelihood values. Start with example Rosenbrock, it's better documented than others. Example multiv_gaussian has also an experimental implementation of combo stepper, that *needs debugging* (on my TODO list), *avoid using* it. 


compile: g++ -std=c++14 main.cpp -o main -ltbb -O3 

Software tested with gcc 5.4.0  
