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

#ifndef _mcmc_algo_pcx_
#define _mcmc_algo_pcx_

#include "../population/pop.hpp"
#include "../individual.hpp"


#include "../base/mcmc_algo_base.hpp"
#include "../parallel"

namespace mcmc{
namespace algos{



template<class LogLkhood,class stepper, class mh_ratio >
class mcmc_algo_pcx: public  mcmc::base::mcmc_algo_base { /**< This inheritance is use to declare type of CRTP classes as members */

	typedef  mcmc::population::pop<LogLkhood> mcmc_pop; 
	
	private: 

		long int counter; 

		LogLkhood logP;  	/**< LogLikelihood function  */

		stepper * my_step;  	/**< Memory will be allocated */
		mh_ratio *  my_mh; 	/**< Memory will be allocated */

		bool step_alloc;
		bool mh_alloc;

		vector<vector<double> > PriorRange; 
		
		std::vector <  mcmc::individual >  tP_1; /**< Temporary instance of individuals: it is needed for detailed balance, before and after modification */
		std::vector <  mcmc::individual >  tP_2; /**< Temporary instance of individuals: it is needed for detailed balance, before and after modification  */
	

		string flname_out; 
		ofstream mcmc_out; 
		bool flname_set;
		int write_binary_flag;  // Binary flag

		int Nvars;
		int Npop; 
		
		void update()	  /**< Updates population tP_1 from tP_2 */
			{
			tP_1 = tP_2; 
			}

		
		void write_txt();
		void write_binary(); 

	public:
		

		mcmc_algo_pcx(mcmc_pop &_P ); 
		mcmc_algo_pcx(mcmc_pop &_P,stepper * _my_step, mh_ratio * _my_mh ); 

		~mcmc_algo_pcx()
			{
			if (step_alloc)
				delete my_step; 
			if (mh_alloc)
				delete my_mh; 

			mcmc_out.close(); 
			}

		void set_flname_out(string & _flname_out, int _write_binary=0); 
		void write();


		// ************************ TBB parallelized routines **********************************************
		void evolve_single(); /**< Evolves the population P for a single generation */ 
		void evolve(); /**< Evolves the population P for a single generation using Nthreads */ 
		//void evolve(std::vector<double> &beta, int Nthreads); // Evolves using parallel tempering  for a single generation  // This will probably go in a different algorithm: I need different population.  
		// *************************************************************************************************


		
		std::vector<mcmc::individual> get_pop()
			{
			return tP_1;
			} 



		// All sampling will be done in parallel 
		void sample_single(int &Nsteps); /**< Sample till death */ 
		void sample(int &Nsteps); /**< Sample till death in parallel */ 

};

template<class LogLkhood,class stepper, class mh_ratio >
void mcmc_algo_pcx<LogLkhood, stepper,  mh_ratio >::write_txt() {/* Writes variables, loglkhood, accept/ratio in flname_out */


	for (auto i=0; i < Npop; ++i)
		{
		for (auto j=0; j < Nvars; ++j)
			mcmc_out<< tP_1[i].Vars[j]<< " ";
		mcmc_out<< tP_1[i].loglkhood<<"\n"; 

		// TODO: Add accept/reject ratio 

		}


}


template<class LogLkhood,class stepper, class mh_ratio >
void mcmc_algo_pcx<LogLkhood, stepper,  mh_ratio >::write_binary(){/* Writes variables, loglkhood, accept/ratio in flname_out */


	double val; 
	for (int i=0; i <Npop; ++i)
		{
		for (int j=0; j < Nvars; ++j)
			{
			val=tP_1[i].Vars[j];
        		mcmc_out.write((char *)&val,sizeof(val));
			}
			val=tP_1[i].loglkhood;
        		mcmc_out.write((char *)&val,sizeof(val));

			
		// TODO: Add accept/reject ratio 

		}





}




template<class LogLkhood,class stepper, class mh_ratio >
mcmc_algo_pcx<LogLkhood, stepper,  mh_ratio >::mcmc_algo_pcx(mcmc_pop &_P ): logP(_P.logP)
	{	

	counter			=	1;

	Nvars 			=  	_P.get_Nvars() ; // <--------------------------------
	Npop 			=   	_P.get_Npop() ;  // <--------------------------------

	write_binary_flag	=	0; // By default, writes in ascii file. 
	flname_set		=	false;

	PriorRange		=	_P.get_PriorRange(); // <---------------------------

	tP_1 			=	_P.get_pop();  // <---------------------------------
	tP_2			=	tP_1; 

	my_step 		= 	new   stepper(PriorRange, Nvars,Npop);

	my_mh 			= 	new   mh_ratio ;

	
	step_alloc		= 	true; 
	mh_alloc		=	true;


}


template<class LogLkhood,class stepper, class mh_ratio >
mcmc_algo_pcx<LogLkhood,stepper,mh_ratio>::mcmc_algo_pcx(mcmc_pop &_P,stepper * _my_step, mh_ratio * _my_mh ): logP(_P.logP), my_step(_my_step), my_mh(_my_mh){

	counter			=	1;

	Nvars 			=  	_P.get_Nvars() ; // <--------------------------------
	Npop 			=   	_P.get_Npop() ;  // <--------------------------------

	write_binary_flag	=	0; // By default, writes in ascii file. 
	flname_set		=	false;

	PriorRange		=	_P.get_PriorRange(); // <---------------------------

	tP_1 			=	_P.get_pop();  // <---------------------------------
	tP_2			=	tP_1; 

	step_alloc		= 	false; 
	mh_alloc		=	false;


}

template<class LogLkhood,class stepper, class mh_ratio >
void mcmc_algo_pcx<LogLkhood,stepper,mh_ratio>:: set_flname_out(string &_flname_out, int _write_binary){ 


	/* Set output filename */
	flname_out = _flname_out; 
	 if ( _write_binary == 1  ) {
		write_binary_flag=1;
		mcmc_out.open(flname_out.c_str(),std::ios::out );
		flname_set = true; 
	} else {
		mcmc_out.open(flname_out.c_str(),std::ios::out | std::ios::binary);
		flname_set = true; 
	}

}



template<class LogLkhood,class stepper, class mh_ratio >
void mcmc_algo_pcx<LogLkhood, stepper,  mh_ratio >::write(){
	switch (write_binary_flag){
	case 1:
		write_binary();
		break;
	default:
		write_txt();
		break; 
	}


}

template<class LogLkhood,class stepper, class mh_ratio >
void mcmc_algo_pcx<LogLkhood, stepper,  mh_ratio >::evolve_single(){



	// ********************************* PARALLELIZE TBB ********************************
	mcmc::individual  tY;

	bool within_PriorRange; 
	for (auto j=0; j< Npop; ++j)
		{

		within_PriorRange = true;  

		// Propose new individual. 
		tY = my_step->propose( tP_1 , j);

		// Accept / Reject 
		// Check if its values are within PriorRange 
		for (auto i=0; i < Nvars; ++i)
                        if( (tY.Vars[i] > PriorRange[i][1]) || (tY.Vars[i]< PriorRange[i][0])  )
				{
				within_PriorRange = false; 
				}

		
		if (within_PriorRange == false ) {// Reject proposed value  
			tP_2[j] =  tP_1[j] ; // Keep old value 
		} else {
			tY.loglkhood= logP(tY); 
			if ( my_mh->accept( tY ,tP_1[j])==true ){ 
				tP_2[j] = tY; // Accept proposed value 
			}else {
				tP_2[j] =  tP_1[j]; // Keep old value 
			}
		}	
		}
	
	// *******************************************************************************************************


		// Update population -- SLOW operation, would love to avoid it.  
		update();  // Copies values of tP_1 <-- tP_2

	
}





/**
This is the most important function of the algorithm. It is the SAME independent of algorithmic structure. ***** So it should NOT repeated.***** 

*/

template<class LogLkhood,class stepper, class mh_ratio >
void mcmc_algo_pcx<LogLkhood, stepper,  mh_ratio >::evolve(){

	

	auto  evolve_par = [&](tbb::blocked_range<int> &r)->void{

		mcmc::individual  tY;
		bool within_PriorRange; 

		for (auto j=r.begin(); j!= r.end(); ++j)
			{
			within_PriorRange = true;  

			// Propose new individual. 
			tY = my_step->propose( tP_1 , j);

			// Accept / Reject 
			// Check if its values are within PriorRange 
			for (auto i=0; i < Nvars; ++i)
        	                if( (tY.Vars[i] > PriorRange[i][1]) || (tY.Vars[i]< PriorRange[i][0])  )
					{
					within_PriorRange = false; 
					}

		
			if (within_PriorRange == false ) {// Reject proposed value  
				tP_2[j] =  tP_1[j] ; // Keep old value 
			} else {
				tY.loglkhood= logP(tY); 
				if ( my_mh->accept( tY ,tP_1[j])==true ){ 
					tP_2[j] = tY; // Accept proposed value 
				}else {
					tP_2[j] =  tP_1[j]; // Keep old value 
				}
			}	
			}
	
		}; // END lambda function 


	static tbb::affinity_partitioner ap; 
	tbb::parallel_for(tbb::blocked_range<int>(0,Npop), evolve_par , ap);


	// Update population -- SLOW operation, would love to avoid it.  
	update();  // Copies values of tP_1 <-- tP_2

}




template<class LogLkhood,class stepper, class mh_ratio >
void mcmc_algo_pcx<LogLkhood,stepper,mh_ratio>::sample_single(int &Nsteps) /**< Sample till death */ 
	{
	if (flname_set == true ){
	for (auto i=0; i < Nsteps; ++i)
		{

		// Vocal evolution
                if(i%int(0.05*Nsteps)==0) 
                        cout<<"Percentage completed  --> "<<100.0*double(i)/double(Nsteps)<<"%" << endl; 
		evolve_single(); 
		write ();
		counter++;
		}
	}else {
		std::cerr<< "Either flname_out or algo not set, yet you request sampling. WTF?  Aborting ..." << std::endl; 
		DEBUG (flname_set);
		throw 0;

	}
	
}


template<class LogLkhood,class stepper, class mh_ratio >
void mcmc_algo_pcx<LogLkhood,stepper,mh_ratio>::sample(int &Nsteps) /**< Sample till death */ 
	{
	if (flname_set == true ){
	for (auto i=0; i < Nsteps; ++i)
		{
		// Vocal evolution
                if(i%int(0.05*Nsteps)==0) 
                        cout<<"Percentage completed  --> "<<100.0*double(i)/double(Nsteps)<<"%" << endl; 
		evolve(); 
		write ();
		counter++;
		}
	}else {
		std::cerr<< "flname_out or algo not set, aborting ..." << std::endl; 
		DEBUG (flname_set);
		throw 0;

	}
	
}









}
}
#endif
