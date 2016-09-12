/*
    This file is part of popmcmc++.

    popmcmc++ is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    popmcmc++ is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with popmcmc++. If not, see <http://www.gnu.org/licenses/>.

*/




#ifndef _PCX_stepper_
#define _PCX_stepper_

#include <vector>
#include <set>

#include "../globals.h"
#include "../individual.hpp"
#include "../tools/vec_ops.hpp"

namespace mcmc{
namespace steppers{


/**
	
	PCX - inspired stepper 

*/


template <class T> // This will always be std::vector<mcmc::individual>  - remove template at later stage. 
class PCX_stepper  {

	private:
 	
		std::vector<std::vector<double> > PriorRange; // PriorRange for scaling / rescaling variables

		int Nvars; // Total number of variables  
		int Npop;  // Population 
		T propvecx; 



		/* Random integer, for selection of individuals from population */
		std::uniform_int_distribution<> dis_int; // (0, Npop-1);	// Integer in the range 0,...,Nwalkers-1.
        	std::normal_distribution<> gauss;// (-epsilon, epsilon);// from a=1 to b=2. Call with: "dis(gen)" without " " characters.
	

		double sigma; // Variance of gaussian normal 	
		/* Total number of parents contributing to crossover, 2 < mu < Npop */
		int mu; 
		

		string name; 

		// Necessary for normalization (scaling) of all variables to [0,1] 
        	double scale_vars (double &x, int &idx)
			{
                	return (x - PriorRange[idx][0] ) / (  PriorRange[idx][1] -  PriorRange[idx][0] ) ;           
                	};

        	double restore_vars (double &x, int &idx)
				{
                		return  PriorRange[idx][0] +   x * (  PriorRange[idx][1] -  PriorRange[idx][0] ) ;                                         

        		};




	public: 
		PCX_stepper(const int &_Nvars, const int &_Npop);

		void set_mu(int &_mu){mu = _mu;}
		void set_PriorRange(std::vector<std::vector<double> > &_PriorRange){PriorRange = _PriorRange; } /**< Sets PriorRange. If this is not set inside some algo, then it doesn't have default values, sampler will give erroneous results */ 

		/*
		For a weighted version of mean estimate, I need information of the likelihood function as well. 
		So vector<vector<double> > X is not appropriate choise. Still, I do not know if a weighted version satisfied detailed balance. I cannot think of a reason why it doesn't? 
		*/
		T propose( std::vector<T> &X, int &j) ;
		
		string get_name(){return name;}; 	

};










template <class T>
PCX_stepper<T>::PCX_stepper(const int &_Nvars, const int &_Npop ): Nvars (_Nvars), Npop(_Npop) {
	
	mu = 5; // Default choice. 

	name = "PCX_stepper"; 

	
	// Set up distributions. 
	uniform_int_distribution<>::param_type newParams2(0, Npop-1);
	dis_int.param(newParams2);


	// Variance of normal distribution
	sigma=0.1;
	normal_distribution<>::param_type newParams1(0.0,sigma);
	gauss.param(newParams1);
		
	// Set up proper size of proposed vector, and gamma step 	
	propvecx.Vars.resize(Nvars); 



}



/*
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	---------------- DEBUGGING MODE  ---------------

// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*/

/**
	Algorithm: 
		a. Select a pool of random mu (mate_pool) individuals from the population. 
		b. Select random individual X. 
		c. Calculate the mean vector of all variables. 
		d. Calculate distance d = pool_mean - X 
		e. Calculate average perpendicular distance Di of all individuals i in pool from d and find <D> (average distance). 
		f. Calculate the normal (to d) basis e(i) of the mate_pool individuals. 
		h. Proposed offspring: 
		
				Y = P[j] + gauss_1 * d * gauss_2 * \sum_i <D> e(i); 

*/

template <class T>
T  PCX_stepper<T>::propose(std::vector<T> &P, int &idx_p) 
	{

	// Assign X[idx_p] to be the PARENT: 
	vector<double>  Xparent = P[idx_p].Vars;
	for (auto i=0; i < Nvars; ++i)
			Xparent[i] = scale_vars(Xparent[i],i);	


	// a. Select a pool of random mu+1 (mate_pool+random) individuals from the population, different from the parent. 
	
	std::set<int> pool_indices; 
	int R1; 
	do{
		do{
			R1=dis_int(gen);
		}while (R1==idx_p);
		pool_indices.insert(R1);

	}while (pool_indices.size() < (mu+1)); 


	std::vector<vector<double> > mate_pool;
	mate_pool.reserve(mu+1);
	for (auto &it: pool_indices)
		mate_pool.push_back(P[it].Vars);

	// Scale vars of mate_pool: 
	for (auto i=0; i < mu+1; ++i)
		for (auto j=0; j < Nvars; ++j)
			mate_pool[i][j] = scale_vars( 	mate_pool[i][j], j); 


	// b. Select random individual X (last entry from the pool). 
	std::vector<double> X = mate_pool.back(); 


// @@@@@@@@@@@@ SAFE GUARD: remove last element from mate_pool. @@@@@@@@@@@@@@@@@ -- remove after debugging. 
	mate_pool.pop_back(); 
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
	

	// c. Calculate the mean vector of individuals in the mate pool. 
        // Mean vector g
        std::vector<double> g(Nvars,0.0);  
        for (auto kk=0; kk < Nvars; ++kk)
                {
                for (auto i=0; i <mu; ++i)
                        g[kk]+= mate_pool[i][kk];
                g[kk] /= double (mu);
                }

	// d. difference vector d, from X to mean. 
        std::vector<double> d(Nvars,0.0);
        for (auto kk=0; kk < Nvars; ++kk)
                d[kk] = X[kk] - g[kk];




	// e. Calculate perpendicular distances Di for each of the rest mate_pool (except parent) to the line d 
 

      
	mcmc::tools::vec_ops myOps; 


        // Construct Parallel  + perpendicular vectors (necessary for decomposition of each individual to vertical and parallel)
 
        std::vector<std::vector<double> > Xp(mu,std::vector<double> (Nvars,0.0)); // Parallel to d 
        std::vector<std::vector<double> > Xv(mu,std::vector<double> (Nvars,0.0)); // Vertical to d 



                // Construct parallel first:    
                {
        	for (auto kk=0; kk < mu; ++kk)
                        Xp[kk] = myOps.projection(mate_pool[kk],d);

                // Now construct Verticals (vertical components of mate_pool vectors, are the perpendicular distances). 
        	for (auto kk=0; kk < mu; ++kk)
                        for (auto j=0; j < Nvars; ++j)
                                Xv[kk][j] = mate_pool[kk][j] - Xp[kk][j];

                }



        	// Find average distance D:
                // Perpendicular distances of each of the rest of the mu-1 individuals (length of Xv vector). 
                std::vector<double> D(mu);
                for (auto i=0; i < mu; ++i)
                        D[i] = myOps.norm(Xv[i]);


                double Davg=0.0;
                for (auto &it: D)
                        Davg += it;
                Davg /= double (mu);



		// Calculate orthonormal basis ev from the set of the rest parents. Use Eigenvector decomposition 
        myOps.Ortho_GS_base(Xv); // This operation makes Xv orthonormal




        // Evaluate new vector y from parent xp According to PCX proposed value 
        std::vector<double> Vertical(Nvars,0.0);
        for (auto i=0; i < Nvars; ++i)
                {
                // Term with sum of normal vectors 
                for (auto j=0; j < mu; ++j)
                        Vertical[i] += Davg*Xv[j][i];

                }



        double w1,w2;

        w1=gauss(gen);
        w2=gauss(gen);

        for (auto i=0; i < Nvars; ++i){
                        propvecx.Vars[i] = Xparent[i] + w1 * d[i] + w2 * Vertical[i];
			propvecx.Vars[i] = restore_vars(propvecx.Vars[i],i); // Scaling back to physical dimensions. 
	}

        return propvecx;



}






}//end steppers namespace 
}// end mcmc namespace 


#endif
