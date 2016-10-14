# popmcmc++
Population mcmc c++ library.

**Beta testing mode - developement branch**

This is a population mcmc c++ library. It incorporates 4 algorithms, 3 reliable and tested and one experimental **(PCX - not finalized, cannot use it)**. 
Tested algorithms: DE,DES based on differential evolution- Braak et al. (2006), Streens 2002 -  
and GWS (Goodman and Weare 2010), the stretch move.  

 Parallel tempering versions are also provided 
(Algos: DE, DES, GW_Stretch, DE_PT, DES_PT, GW_Stretch_PT). 

Dependencies: intel TBB. 

The library is fully parallelized using intel TBB library. There is the option of selecting maximum number of threads (TBB actually decides how many will use), or running in single thread mode, without the use of TBB.  

Jumping in the examples directory should get you started easily. There is a walkthrough example in the wiki. The algorithm outputs variable chains and last column has loglikelihood values. Start with example Rosenbrock, it's better documented than others. Example multiv_gaussian has also an experimental implementation of combo stepper, that *needs debugging* (on my TODO list), *avoid using* it. 


compile: g++ -std=c++14 main.cpp -o main -ltbb -O3 

Software tested with gcc 5.4.0  
