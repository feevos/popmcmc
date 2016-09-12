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


#ifndef _pop_
#define _pop_

#include <iostream>
#include <fstream>
#include <vector>



#include "../individual.hpp"
#include "../parallel" // Parallelization routines TBB
#include "../macros.h"
#include "../globals.h"

namespace mcmc{
namespace population{


/*
This class does not participate in expensive operations, I can use virtual base class. 
*/
template < class LogLikelihood> 
class pop  { 
	private: 
		std::vector<mcmc::individual> pop_vec; /**< Population -once initialized is all we care about */
	protected:	
		bool pop_init; /**< Knots if population is initialized */
	public:	
		LogLikelihood logP;  /**< LogLikelihood function  */
	private:	
		std::vector<std::vector<double> > PriorRange;  /**< Range of real variables */

	protected:	
		int Nvars; /**< ConstLength_vars.size() */
		int Npop; /**< Total number of individuals. */

	private:	
		void populate_random(int &idx_start, int &idx_end);  
		void read_pop(std::string &_flname_in);

		void talk_to_me(){
			std::cout<< "Population initiated, proceeding to sampling..." << std::endl; 

			}

	
	public: 
		pop(LogLikelihood &_logP, std::vector<std::vector<double> > &_PriorRange): logP(_logP), PriorRange(_PriorRange), Nvars(_PriorRange.size()){ pop_init = false;}  

		virtual ~pop(){}

		// TODO: Add threaded versions in initialization. 
		void init(int &_Npop); /**< Random initializer */
		void init(std::vector<mcmc::individual> &_pop_vec);  /**< Start from known population */ 
		void init(int &_Npop , mcmc::individual &_indiv);  /**< Start from a known good solution by perturbing around the _indiv */ 
		void init(int &_Npop, std::vector<mcmc::individual> &_pop_vec);  // Initiate using all of _pop vector, and complete difference Npop - _pop.size() with random initiates Npop > _pop.size()
		void init(std::string &flname_read); /**< Initiate population directly from file */ 
		//void init(int &_Npop, std::string &flname_read); /**< Initiate population directly from file,  complete difference Npop - _pop.size() with random initiates  */


		/*
			Get instantatiated vector of individuals 
		*/

		std::vector<mcmc::individual> get_pop() 
			{
			if(pop_init == false){
				DEBUG(pop_init);
				std::cerr<< " Request for non initialized population, aborting ..." << std::endl;
				throw 0; 
				}
			
			return pop_vec;
			}

		
		void set_pop(std::vector<mcmc::individual> &_pop_vec ) // Set population values, WITHOUT calculating logP. 
			{
			pop_vec = _pop_vec; 
			}


		int get_Npop(){return Npop;}
		int get_Nvars(){return Nvars;}
		std::vector<std::vector<double> > get_PriorRange() 
			{
			return 
				PriorRange; 
			}

		pop<LogLikelihood> * get_full_pop(){return this;} 

};





template <class LogLikelihood> 
void pop<LogLikelihood>::populate_random(int &idx_start, int &idx_end){


	if (idx_start >= idx_end){
		std::cerr << idx_start << " >= " << idx_end << " can't populate, aborting ... " << std::endl; 
		DEBUG(idx_start);
		DEBUG(idx_end);
		throw 0;
	}else if ( (Npop-1) <idx_end  || idx_start < 0 ) {
		std::cerr << " Error in range of indices,  can't populate, aborting ... " << std::endl; 
		DEBUG(idx_start);
		DEBUG(idx_end);
		DEBUG(Npop); 
		throw 0;

	} else {




/*	
        // Fill in  
        std::uniform_real_distribution<> unif_vars(0.0,1.0);

	// ******************* TODO: PARALLELIZE TBB *********************************
	// Initiate random population  
	for (auto i=idx_start; i <= idx_end; ++i)
		{
		pop_vec[i].Vars.resize(Nvars);

                 // Insert some random values for Real Vars :
                 for (auto j=0; j < Nvars; ++j)
                 	pop_vec[i].Vars[j] =  unif_vars(gen) * (PriorRange[j][1] - PriorRange[j][0])   + PriorRange[j][0]  ;


                 // Now evaluate and store the logP value for each chromosome - LENGTHY operation 
                 pop_vec[i].loglkhood=logP(pop_vec[i]);
		}
	// ******************* TODO: PARALLELIZE TBB *********************************
	
*/


	auto populate = [&] (tbb::blocked_range<int> &r){
		
	
        	// Fill in  
	        std::uniform_real_distribution<> unif_vars(0.0,1.0);

		for (auto i=r.begin(); i != r.end(); ++i)
			{
			pop_vec[i].Vars.resize(Nvars);

        	         // Insert some random values for Real Vars :
	                 for (auto j=0; j < Nvars; ++j)
        	         	pop_vec[i].Vars[j] =  unif_vars(gen) * (PriorRange[j][1] - PriorRange[j][0])   + PriorRange[j][0]  ;


                	 // Now evaluate and store the logP value for each chromosome - LENGTHY operation 
	                 pop_vec[i].loglkhood=logP(pop_vec[i]);
			}



		};



	tbb::affinity_partitioner ap; 
	tbb::parallel_for(tbb::blocked_range<int>(idx_start,idx_end+1), populate , ap);

	}






}





// This function reads directly 
template <class LogLikelihood> 
void pop<LogLikelihood>::read_pop ( std::string &_flname_in){

	
        std::string tmpString;
        std::ifstream txtFile(_flname_in);

        std::vector<std::vector<double> > tMatrix;
        tMatrix.reserve(1000000); // Reserve some large number of rows for reading data. Won't need in general more than 10^6 input read. 

        std::vector<double> temp_vals;
        if(txtFile.is_open())
                {
                while(getline(txtFile, tmpString))
                        {
			if (tmpString.front() != '#')	// Do not read lines that start with # 
				{
                        	temp_vals.clear();  // Empty vector<T> that holds line. 
	                        std::stringstream lineStream (tmpString); // Pass string to stringstream for easier processing. 
	                        copy(std::istream_iterator<mcmc::individual>(lineStream), std::istream_iterator<mcmc::individual>(), std::back_inserter(temp_vals)); // Copy values to vector<T> 
                        	// Pass values to temporary Matrix. 
        	                tMatrix.push_back(temp_vals);
				} // End if 
                        }

                txtFile.close();
                }

	// Pass tmatrix -->  matrix 
//      Matrix = std::move(tMatrix);
	tMatrix.shrink_to_fit();  // release unnecessary memory 

	if (tMatrix.size() != Npop) {
		std::cerr<< "File size different than declared population size, aborting ..." << std::endl; 
		DEBUG(tMatrix.size());
		DEBUG(Npop);
		throw(1); 

	}	


	pop_vec.resize(tMatrix.size()); 


	// ******************* TODO: PARALLELIZE TBB *********************************
	for (auto i=0; i < Npop; ++i)
		{
		pop_vec[i].Vars.resize(Nvars); 
		for (auto j=0; j < Nvars; ++j)
			pop_vec[i].Vars[j] = tMatrix[i][j]; 

		pop_vec[i].loglkhood = logP(pop_vec[i]);  // This is a slow operation. 
		}	
	// ******************* TODO: PARALLELIZE TBB *********************************


	// Parallel routine for evaluation of new likelihood. 
	auto populate = [&] (tbb::blocked_range<int> &r){

		for (auto i=r.begin(); i != r.end(); ++i)
			{		
			pop_vec[i].Vars.resize(Nvars); 
			for (auto j=0; j < Nvars; ++j)
				pop_vec[i].Vars[j] = tMatrix[i][j]; 

	                pop_vec[i].loglkhood=logP( pop_vec[i] ); // Expensive operation. 
			}
		}; 
	

	tbb::affinity_partitioner ap; 
	tbb::parallel_for(tbb::blocked_range<int>(0,Npop), populate , ap);



}







template <class LogLikelihood> 
void pop<LogLikelihood>::init (int &_Npop){

	Npop = _Npop; 

	pop_vec.resize(Npop);

	int idx_start = 0;
	int idx_end = pop_vec.size()-1;
	populate_random(idx_start,idx_end);  // Parallelized routine. 
	
	pop_init = true; 
	talk_to_me(); 

}



template <class LogLikelihood> 
void pop<LogLikelihood>::init (std::vector< mcmc::individual > &_pop_vec ){
	
	Npop =  static_cast<int>( _pop_vec.size()) ; 

	pop_vec = _pop_vec; 

/*
	// ******************* PARALLELIZE TBB *********************************
		// Copy given population to part of pop
                for (auto i=0; i <  Npop; ++i)
                        {
                        pop_vec[i] = _pop_vec[i];
                        pop_vec[i].loglkhood = logP( pop_vec[i]);
                        }
	// ******************* PARALLELIZE TBB *********************************
*/
	
	// Parallel routine for evaluation of new likelihood. 
	auto populate = [&] (tbb::blocked_range<int> &r){

		for (auto i=r.begin(); i != r.end(); ++i)
			{		
                	 // Evaluate and store the logP value for each individual 
	                 pop_vec[i].loglkhood=logP( pop_vec[i] ); // Expensive operation. 
			}
	}; 
	

	tbb::affinity_partitioner ap; 
	tbb::parallel_for(tbb::blocked_range<int>(0,Npop), populate , ap);


	pop_init = true; 
	talk_to_me(); 
}


/*
	Parallelized, create random centered at individual. 
*/
template <class LogLikelihood> 
void pop<LogLikelihood>::init (int &_Npop,   mcmc::individual &_indiv ){
	
	Npop  = _Npop;
	if (Nvars != static_cast<int>(_indiv.Vars.size() )){
		std::cerr<< "Incompatible variables size of PriorRange and _indiv, aborting ..." << std::endl;
		DEBUG(Nvars);
		DEBUG(_indiv.Vars.size() );
		throw 0; 

	}



	pop_vec.resize(Npop); 
	auto populate = [&] (tbb::blocked_range<int> &r){


		// Makes sure that each call is an individual random number 
		// Construct random gaussians with variance at the same magnitude as each variable. 
		std::vector<std::normal_distribution<> > gauss_vec(Nvars);   
		
		double exponent; 
		double sigma; 
	        for (auto j=0; j < Nvars; ++j)
			{
			exponent = std::round(log10( std::abs(_indiv.Vars[j]) ));
			sigma	 = (1.e-3) * pow(10.0,exponent); 

			std::normal_distribution<>::param_type newParams1(0.0,sigma);
			gauss_vec[j].param(newParams1); 
			}


		for (auto i=r.begin(); i != r.end(); ++i)
			{		
			pop_vec[i].Vars.resize(Nvars);
        	         // Insert some random values for Real Vars :
	                 for (auto j=0; j < Nvars; ++j)
				pop_vec[i].Vars[j] = _indiv.Vars[j] * ( gauss_vec[j](mcmc::gen) +1.);
                	 // Now evaluate and store the logP value for each chromosome - LENGTHY operation 
	                 pop_vec[i].loglkhood=logP(pop_vec[i]); // Expensive operation. 
			}
		};


	tbb::affinity_partitioner ap; 
	tbb::parallel_for(tbb::blocked_range<int>(0,Npop), populate , ap);

	pop_init = true; 
	talk_to_me(); 
}





template <class LogLikelihood> 
void pop<LogLikelihood>::init(int &_Npop, std::vector< mcmc::individual > &_pop_vec ){

	Npop=_Npop; 
        int Nkeep = _pop_vec.size();  // Keep these from the restart file. 
        int Ndiff = Npop - Nkeep;  // Difference between initial population and requested population. 



        if ( Ndiff < 0  ){
                std::cout << "Vector of individuals exceeds size of population, need a different initialization function, aborting... " << std::endl;
                DEBUG(Ndiff);
                throw 0;


        } else if (Ndiff ==0) {

		init(_pop_vec); 

        } else if (Ndiff > 0){


                pop_vec.resize(Npop);

/*
	// *******************  PARALLELIZE TBB *********************************
		// Copy given population to part of pop
                for (auto i=0; i <  Nkeep; ++i)
                        {
                        pop_vec[i] = _pop_vec[i];
                        pop_vec[i].loglkhood = logP( pop_vec[i] );
                        }
	// *******************  PARALLELIZE TBB *********************************
*/


		// Parallel routine for evaluation of new likelihood. 
		auto populate = [&] (tbb::blocked_range<int> &r){

			for (auto i=r.begin(); i != r.end(); ++i)
				{		
                		 // Evaluate and store the logP value for each individual 
	                        pop_vec[i] = _pop_vec[i];
		                pop_vec[i].loglkhood=logP( pop_vec[i] ); // Expensive operation. 
				}
		}; 
	

		tbb::affinity_partitioner ap; 
		tbb::parallel_for(tbb::blocked_range<int>(0,Nkeep), populate , ap);


		// Initiate rest as random 
		int idx_end = pop_vec.size()-1;
		populate_random(Nkeep,idx_end); 
		}

	
	pop_init = true; 
	talk_to_me(); 
}





template <class LogLikelihood> 
void pop<LogLikelihood>::init( std::string &flname_read){

	read_pop(flname_read); 

	pop_init = true; 
	talk_to_me(); 
}








} // end population
} // End mcmc 
#endif
