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




#ifndef _popmcmc_
#define _popmcmc_
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

// Available steppers 
#include "steppers/GW_Stretch_stepper.hpp"
#include "steppers/DES_stepper.hpp"
#include "steppers/DE_stepper.hpp"
#include "steppers/PCX_stepper.hpp"

// Available MH ratios 
#include "mh_ratios/GW_Stretch_mh.hpp"
#include "mh_ratios/std_mh.hpp"


// Population and basic algo information
#include "population/pop.hpp"
#include "population/pop_pt.hpp"
#include "algos/mcmc_algo.hpp"
#include "algos/mcmc_algo_pt.hpp"
// ---------- PCX experimental ----------
#include "algos/mcmc_algo_pcx.hpp"
#include "algos/mcmc_algo_pcx_pt.hpp"


namespace mcmc{



template<class logLkhood>
class popmcmc final { 
	private: 
		bool algo_set; 
		bool algo_alloc; 
		/**< Does the evolution for all  generations */
		base::mcmc_algo_base * myalgo; 


// --------------------------------------------------- 	+
		string ALGO_NAME; //			|
// --------------------------------------------------- 	+


		logLkhood logP; 
		vector<vector<double> > PriorRange; 

		void initiate_algo ( mcmc::population::pop<logLkhood>  &mypop);
		void initiate_algo_pt ( mcmc::population::pop_pt<logLkhood>  &mypop,int &Nswap);

	public: 	
		popmcmc(logLkhood _logP,  vector<vector<double> > &_PriorRange,int Nthread= -1);
		~popmcmc();


/**
		I don't like the idea that I have a lot of init files in here. I'd like to find some way to make it more modular. Transfer them to a friend - templated class ? 
*/
// ----------------------------------------------------------------------------------------------------------- 

		/**< Initializer of population */
		void init(std::string &flname_read); /**< Initiate population directly from file */ 
		void init(int &_Npop); /**< Random initializer */
		void init(int &_Npop, mcmc::individual &_myIndiv);
		void init(std::vector<mcmc::individual> &_pop_vec);  /**< Start from known population */ 
		void init(int &_Npop, std::vector<mcmc::individual> &_pop_vec);  // Initiate using all of _pop vector, and complete difference Npop - _pop.size() with random initiates Npop > _pop.size()



		/**< Parallel tempering constructors */
		void init(std::string &flname_read,int dim_T, int Nswap); 
		void init(int &_Npop,int dim_T, int Nswap); /**< Random initializer - PT  */
		void init(std::vector<std::vector<mcmc::individual> > &_pop_pt_vec, int dim_T, int Nswap);  /**< Start from known population - PT*/ 
		void init(int &_Npop,mcmc::individual&_myIndiv, int dim_T, int Nswap);  /**< Start from known population - PT*/ 
		void init(int &_Npop, std::vector<std::vector<mcmc::individual> > &_pop_pt_vec,int dim_T, int Nswap);  // Initiate using all of _pop vector, and complete difference Npop - _pop.size() with random initiates Npop > _pop.size()

// ----------------------------------------------------------------------------------------------------------- 




		void set_algo(string _ALGO_NAME); 
		void set_algo(base::mcmc_algo_base * _myalgo); /**< constructor for arbitrary mcmc algo. */		

		void set_flname_out(string &_flname_out, int _write_binary=0); /* Set output filename */
		
		// All sampling will be done in parallel 
		void sample_single(int &Nsteps); /**< Sample till death */ 
		void sample(int &Nsteps); /**< Sample till death in parallel */ 


};



template <class logLkhood>
void popmcmc<logLkhood>::initiate_algo ( mcmc::population::pop<logLkhood> &mypop){

	if (ALGO_NAME == GW_Stretch ){

	typedef steppers::template GW_Stretch_stepper< mcmc::individual>  GWS_stepper; 
	typedef mh_ratios::GW_Stretch_mh GWS_mh_ratio; 

	myalgo = new  mcmc::algos::mcmc_algo<logLkhood,GWS_stepper,GWS_mh_ratio>( mypop ); // This refers to the population inheritance 
	algo_alloc = true; 

	}else if (ALGO_NAME == DES ){

	typedef steppers::template DES_stepper< mcmc::individual>  tDES_stepper; 
	typedef mh_ratios::std_mh std_mh_ratio; 

	myalgo = new  mcmc::algos::mcmc_algo<logLkhood,tDES_stepper,std_mh_ratio>( mypop );

	algo_alloc = true; 
	} else if ( ALGO_NAME == DE ){


	typedef steppers::template DE_stepper< mcmc::individual>  tDE_stepper; 
	typedef mh_ratios::std_mh std_mh_ratio; 

	myalgo = new  mcmc::algos::mcmc_algo<logLkhood,tDE_stepper,std_mh_ratio>( mypop );

	algo_alloc = true;  
	
	} else if ( ALGO_NAME == PCX ){


	typedef steppers::template PCX_stepper< mcmc::individual> tPCX_stepper; 
	typedef mh_ratios::std_mh std_mh_ratio; 

// @@@@@@@@@@@@@@@@@@ SPECIALIZATION: NEEDS to be integrated with other steppers with more careful design. @@@@@@@@@@@@@@@@@@@@@@@@
	myalgo = new  mcmc::algos::mcmc_algo_pcx<logLkhood,tPCX_stepper,std_mh_ratio>( mypop ); 
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


	algo_alloc = true;  



	} else {
		std::cerr<< "Either algo and PT are incompatible, or you haven't set algo BEFORE initiating population, aborting ... " << std::endl; 
		DEBUG (ALGO_NAME);
		throw 0;
	}
	


	algo_set = true; 



}


// ************************* PARALLEL TEMPERING ALGOS ********************** 

template <class logLkhood>
void popmcmc<logLkhood>::initiate_algo_pt ( mcmc::population::pop_pt<logLkhood> &mypop, int &Nswap){

	
	

	if (ALGO_NAME == GW_Stretch_PT ){

	typedef steppers::template GW_Stretch_stepper< mcmc::individual>  GWS_stepper; 
	typedef mh_ratios::GW_Stretch_mh GWS_mh_ratio; 

	myalgo = new  mcmc::algos::mcmc_algo_pt<logLkhood,GWS_stepper,GWS_mh_ratio>( mypop , Nswap); // This refers to the population inheritance 
	algo_alloc = true; 

	}else if (ALGO_NAME == DES_PT ){

	typedef steppers::template DES_stepper< mcmc::individual>  tDES_stepper; 
	typedef mh_ratios::std_mh std_mh_ratio; 

	myalgo = new  mcmc::algos::mcmc_algo_pt<logLkhood,tDES_stepper,std_mh_ratio>( mypop , Nswap);

	algo_alloc = true; 
	} else if ( ALGO_NAME == DE_PT ){


	typedef steppers::template DE_stepper< mcmc::individual>  tDE_stepper; 
	typedef mh_ratios::std_mh std_mh_ratio; 


	myalgo = new  mcmc::algos::mcmc_algo_pt<logLkhood,tDE_stepper,std_mh_ratio>( mypop , Nswap);

	algo_alloc = true; 

	} else if ( ALGO_NAME == PCX_PT ){


	typedef steppers::template PCX_stepper< mcmc::individual>  tPCX_stepper; 
	typedef mh_ratios::std_mh std_mh_ratio; 


	myalgo = new  mcmc::algos::mcmc_algo_pcx_pt<logLkhood,tPCX_stepper,std_mh_ratio>( mypop , Nswap);

	algo_alloc = true; 

	} else {
		std::cerr<< "Either algo and PT are incompatible, or you haven't set algo BEFORE initiating population, aborting ... " << std::endl; 
		DEBUG (ALGO_NAME);
		throw 0;
	}





	algo_set = true; 


}









template <class logLkhood> 
popmcmc<logLkhood>::popmcmc( logLkhood _logP, vector<vector<double> > &_PriorRange, int Nthread ):logP(_logP)
	{

	if (Nthread == -1){
		tbb::task_scheduler_init init;

	}else if (Nthread >=1){

		tbb::task_scheduler_init init(Nthread);
		
	};


	PriorRange = _PriorRange; 

	algo_set = false; 
	algo_alloc = false; 

	// Random seed 
	gen.seed(rd());

	
	uniform_real_distribution<>::param_type newParams1(0.0,1.0);
	unif_real.param(newParams1);

}


template <class logLkhood> 
popmcmc<logLkhood>::~popmcmc(){
	if (algo_alloc){
		delete myalgo; 
	}
		
}


/**< Initializer of population */


template <class logLkhood> 
void popmcmc<logLkhood>::init(int &_Npop){
	mcmc::population::pop<logLkhood> mypop(logP,PriorRange); 
	mypop.init(_Npop); 
	initiate_algo(mypop); 
}


template <class logLkhood> 
void popmcmc<logLkhood>::init(int &_Npop, mcmc::individual &_myIndiv){
	mcmc::population::pop<logLkhood> mypop(logP,PriorRange); 
	mypop.init(_Npop,_myIndiv); 
	initiate_algo(mypop); 

}


template <class logLkhood> 
void popmcmc<logLkhood>::init(std::vector<mcmc::individual> &_pop){
	mcmc::population::pop<logLkhood> mypop(logP,PriorRange); 
	mypop.init(_pop); 
	initiate_algo(mypop); 

}

template <class logLkhood> 
void popmcmc<logLkhood>::init(int &_Npop, std::vector<mcmc::individual> &_pop){

	mcmc::population::pop<logLkhood> mypop(logP,PriorRange); 
	mypop.init(_Npop,_pop); 
	initiate_algo(mypop); 

}


template <class logLkhood> 
void popmcmc<logLkhood>::init(std::string &flname_read){ /**< Initiate population directly from file */ 
	mcmc::population::pop<logLkhood> mypop(logP,PriorRange); 
	mypop.init(flname_read); 
	initiate_algo(mypop); 
}


// ***********  Parallel tempering constructors  ************

template <class logLkhood> 
void popmcmc<logLkhood>::init(std::string &flname_read,int dim_T, int Nswap){ /**< Initiate population directly from file */ 

	mcmc::population::pop_pt<logLkhood> mypop_pt(logP,PriorRange,dim_T); 
	mypop_pt.init(flname_read); 
	initiate_algo_pt(mypop_pt, Nswap); 

}

template <class logLkhood> 
void popmcmc<logLkhood>::init(int &_Npop,int  dim_T,int Nswap ){ 
	mcmc::population::pop_pt<logLkhood> mypop_pt(logP,PriorRange,dim_T); 
	mypop_pt.init(_Npop); 
	initiate_algo_pt(mypop_pt,Nswap); 
}



/**< Start from known population - PT*/ 
template <class logLkhood> 
void popmcmc<logLkhood>::init(std::vector<std::vector<mcmc::individual> > &_pop_pt_vec,int dim_T,  int Nswap){

	mcmc::population::pop_pt<logLkhood> mypop_pt(logP,PriorRange,dim_T); 
	mypop_pt.init(_pop_pt_vec); 
	initiate_algo_pt(mypop_pt,Nswap); 
}


/**< Start from known population - PT*/ 
template <class logLkhood> 
void popmcmc<logLkhood>::init(int &_Npop, mcmc::individual &_myIndiv,int dim_T,  int Nswap){

	mcmc::population::pop_pt<logLkhood> mypop_pt(logP,PriorRange,dim_T); 
	mypop_pt.init(_Npop, _myIndiv); 
	initiate_algo_pt(mypop_pt,Nswap); 
}




/**< Start from known population - PT*/ 
template <class logLkhood> 
void popmcmc<logLkhood>::init(int &_Npop, std::vector<std::vector<mcmc::individual> > &_pop_pt_vec,int dim_T , int Nswap){
	mcmc::population::pop_pt<logLkhood> mypop_pt(logP,PriorRange,dim_T); 
	mypop_pt.init(_Npop,_pop_pt_vec); 
	initiate_algo_pt(mypop_pt, Nswap); 
}





/**
	This is the constructor of various algorithms 
*/
template <class logLkhood> 
void popmcmc<logLkhood>::set_algo( string _ALGO_NAME)
	{	ALGO_NAME = _ALGO_NAME;

}


template <class logLkhood> 
void popmcmc<logLkhood>::set_algo(base::mcmc_algo_base * _myalgo) /**< constructor for arbitrary mcmc algo. */		{

	myalgo=_myalgo;
	algo_set = true; 
	algo_alloc = false; 
}



template <class logLkhood> 
void popmcmc<logLkhood>:: set_flname_out(string &_flname_out, int _write_binary ){ /* Set output filename */
	myalgo->set_flname_out(_flname_out,_write_binary);
}




template <class logLkhood> 
void popmcmc<logLkhood>::sample_single(int &Nsample) /**< Sample till death */ 
	{
	
	myalgo->sample_single(Nsample);
	
}

	
template <class logLkhood> 
void popmcmc<logLkhood>::sample(int &Nsample) /**< Sample till death */ 
	{
	myalgo->sample(Nsample);
}

} // end of namespace 
#endif
