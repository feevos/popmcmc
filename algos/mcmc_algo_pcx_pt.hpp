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


#ifndef _mcmc_algo_pcx_pt_
#define _mcmc_algo_pcx_pt_

#include "../population/pop_pt.hpp"
#include "../individual.hpp"


#include "../base/mcmc_algo_base.hpp"
#include "../parallel"

namespace mcmc{
namespace algos{



template<class LogLkhood,class stepper, class mh_ratio >
class mcmc_algo_pcx_pt: public  mcmc::base::mcmc_algo_base { /**< This inheritance is use to declare type of CRTP classes as members */

	typedef  mcmc::population::pop_pt<LogLkhood> mcmc_pop; 
	

	private: 
		double Tref; 

		LogLkhood logP;  	/**< LogLikelihood function  */

		stepper * my_step;  	/**< Memory will be allocated */
		mh_ratio *  my_mh; 	/**< Memory will be allocated */

		vector<vector<double> > PriorRange; 
		
		vector<double> beta; 
		int Nswap; 

		std::vector < std::vector <  mcmc::individual > >  tP_1; /**< Temporary instance of individuals: it is needed for detailed balance, before and after modification */
		std::vector < std::vector <  mcmc::individual > > tP_2; /**< Temporary instance of individuals: it is needed for detailed balance, before and after modification  */
	


		
		// ########### SWAP proposal specifics ###############################
		
		uniform_int_distribution<> dist_int_beta;  /**< Gives a random integer value for selection of temperature */
		uniform_int_distribution<> dist_int; 	/**< Selects random population individuals */
		
		void swap( std::vector<std::vector<mcmc::individual> > &P,  int &i, vector<double> &_beta);

		// ##################### END SWAP ####################################



		string flname_out; 
		ofstream mcmc_out; 
		bool flname_set;
		int write_binary_flag;  // Binary flag


		int Nvars;
		int Npop; 
		int dim_T;
	
		void update()	  /**< Updates population tP_1 from tP_2 */
			{
			tP_1 = tP_2; 
			}


		void write_txt();
		void write_binary(); 

	public:
		

		mcmc_algo_pcx_pt(mcmc_pop &_P, int Nswap); 
		mcmc_algo_pcx_pt(mcmc_pop &_P, stepper * _my_step, mh_ratio * _my_mh, int Nswap ); 

		~mcmc_algo_pcx_pt()
			{
			delete my_step; 
			delete my_mh; 
			mcmc_out.close(); 
			}



		void set_flname_out(string & _flname_out, int _write_binary=0); /**< write ascii by default */
		void write();

		// ************************ TBB parallelized routines **********************************************
		void evolve_single(); /**< Evolves the population P for a single generation */ 
		void evolve(); /**< Evolves the population P for a single generation using Nthreads */ 
		//void evolve(std::vector<double> &beta, int Nthreads); // Evolves using parallel tempering  for a single generation  // This will probably go in a different algorithm: I need different population.  
		// *************************************************************************************************

		std::vector<std::vector<mcmc::individual> >get_pop()
			{
			return tP_1;
			} 


		// All sampling will be done in parallel 
		void sample_single(int &Nsteps); /**< Sample till death */ 
		void sample(int &Nsteps); /**< Sample till death in parallel */ 

};



template<class LogLkhood,class stepper, class mh_ratio >
void mcmc_algo_pcx_pt<LogLkhood, stepper,  mh_ratio >::write_txt() {/* Writes variables, loglkhood, accept/ratio in flname_out */


	// a. Find temperature: beta[T]==1; 
	//auto pos = std::find(beta.begin(),beta.end(),Tref);
	for (auto i=0; i < Npop; ++i)
		{
		for (auto j=0; j < Nvars; ++j)
			mcmc_out<< tP_1.back()[i].Vars[j]<< " ";
		mcmc_out<< tP_1.back()[i].loglkhood<<"\n"; 

		// TODO: Add accept/reject ratio 

		}


}


template<class LogLkhood,class stepper, class mh_ratio >
void mcmc_algo_pcx_pt<LogLkhood, stepper,  mh_ratio >::swap(std::vector<std::vector<mcmc::individual> > &P,  int &i, vector<double> &_beta){

	if(i%Nswap==0){ // Try it every Nswap steps. 

	// conjugate temperatures to swap 
	int TT= dist_int_beta(gen); 
	int TT2=TT+1;

	// Choose two chains at random to swap  
	int jj = dist_int(gen);
	int jj2; 

	do {		
	jj2 = dist_int(gen);
	}while (jj2==jj);



	double swaplogratio=log(unif_real(mcmc::gen) );
	double logM =  (_beta[TT] - _beta[TT2]) * ( P[TT2][jj].loglkhood - P[TT][jj2].loglkhood );		


	if(swaplogratio<=logM){ // swap temperatures. 
	std::swap(P[TT2][jj],P[TT][jj2]);  
	}
	};// End of IF Swap proposal. 


}




template<class LogLkhood,class stepper, class mh_ratio >
void mcmc_algo_pcx_pt<LogLkhood, stepper,  mh_ratio >::write_binary(){/* Writes variables, loglkhood, accept/ratio in flname_out */

/*
	// a. Find temperature: beta[T]==1; 
	//auto pos = std::find(beta.begin(),beta.end(),Tref);

	double val; 
	for (int i=0; i <Npop; ++i)
		{
		for (int j=0; j < Nvars; ++j)
			{
			val=tP_1[*pos][i].Vars[j];
        		mcmc_out.write((char *)&val,sizeof(val));
			}
			val=tP_1[*pos][i].loglkhood;
        		mcmc_out.write((char *)&val,sizeof(val));

			
		// TODO: Add accept/reject ratio 

		}
*/


	// a. Find temperature: beta[T]==1; 

	double val; 
	for (int i=0; i <Npop; ++i)
		{
		for (int j=0; j < Nvars; ++j)
			{
			val=tP_1.back()[i].Vars[j];
        		mcmc_out.write((char *)&val,sizeof(val));
			}
			val=tP_1.back()[i].loglkhood;
        		mcmc_out.write((char *)&val,sizeof(val));

			
		// TODO: Add accept/reject ratio 

		}




}




template<class LogLkhood,class stepper, class mh_ratio >
mcmc_algo_pcx_pt<LogLkhood, stepper,  mh_ratio >::mcmc_algo_pcx_pt(mcmc_pop &_P, int _Nswap ): logP(_P.logP)
	{

	Tref=1.0;

	Nvars 		=  	_P.get_Nvars();  // <--------------------------------
	Npop 		=   	_P.get_Npop() ;  // <--------------------------------
	dim_T		=	_P.get_dim_T() ; // <--------------------------------

	Nswap		=	_Nswap; 

// ################################# Am not using this in here, not necessary  #########################################
	// Set up temperatures. These can be in future versions different than linear. 
	beta.resize(dim_T); 
	for (auto i=0; i < dim_T; ++i)
		beta[i] = double (i+1) / double (dim_T); 

	
	std::uniform_int_distribution<>::param_type newParams2(0,dim_T-2);
	dist_int_beta.param(newParams2); 

	std::uniform_int_distribution<>::param_type newParams3(0,Npop-1);
	dist_int.param(newParams3); 


	PriorRange	=	_P.get_PriorRange();  // <--------------------------------

	tP_1 		=	_P.get_pop(); 	      // <--------------------------------
	tP_2		=	tP_1; 

	my_step 	= 	new   stepper(Nvars,Npop);
	
	// Set PriorRange 
	my_step->set_PriorRange(PriorRange);

	my_mh 		= 	new   mh_ratio; 


}


template<class LogLkhood,class stepper, class mh_ratio >
mcmc_algo_pcx_pt<LogLkhood, stepper,  mh_ratio >::mcmc_algo_pcx_pt(mcmc_pop &_P, stepper * _my_step, mh_ratio * _my_mh, int Nswap ): logP(_P.logP), my_step(_my_step), my_mh(_my_mh){


	Tref=1.0;

	Nvars 		=  	_P.get_Nvars();  // <--------------------------------
	Npop 		=   	_P.get_Npop() ;  // <--------------------------------
	dim_T		=	_P.get_dim_T() ; // <--------------------------------

// ################################# Am not using this in here, not necessary  #########################################
	// Set up temperatures. These can be in future versions different than linear. 
	beta.resize(dim_T); 
	for (auto i=0; i < dim_T; ++i)
		beta[i] = double (i+1) / double (dim_T); 


	PriorRange	=	_P.get_PriorRange();  // <--------------------------------

	tP_1 		=	_P.get_pop(); 	      // <--------------------------------
	tP_2		=	tP_1; 



}




template<class LogLkhood,class stepper, class mh_ratio >
void mcmc_algo_pcx_pt<LogLkhood,stepper,mh_ratio>:: set_flname_out(string &_flname_out, int _write_binary){ 


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
void mcmc_algo_pcx_pt<LogLkhood, stepper,  mh_ratio >::write(){
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
void mcmc_algo_pcx_pt<LogLkhood, stepper,  mh_ratio >::evolve_single(){


	// ********************************* PARALLELIZE TBB ********************************
	for (auto T=0; T < dim_T; ++T){
		mcmc::individual  tY;

		bool within_PriorRange; 
		for (auto j=0; j< Npop; ++j)
			{
			within_PriorRange = true;  

			// Propose new individual. 
			tY = my_step->propose( tP_1[T] , j);

			// Accept / Reject 
			// Check if its values are within PriorRange 
			for (auto i=0; i < Nvars; ++i)
	                        if( (tY.Vars[i] > PriorRange[i][1]) || (tY.Vars[i]< PriorRange[i][0])  )
					{
					within_PriorRange = false; 
					}

		
			if (within_PriorRange == false ) {// Reject proposed value  
				tP_2[T][j] =  tP_1[T][j] ; // Keep old value 
			} else {
				tY.loglkhood= logP(tY); 
				if ( my_mh->accept( tY ,tP_1[T][j],beta[T]) ){ 
					tP_2[T][j] = tY; // Accept proposed value 
				}else {
					tP_2[T][j] =  tP_1[T][j]; // Keep old value 
				}
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
void mcmc_algo_pcx_pt<LogLkhood, stepper,  mh_ratio >::evolve(){

	

	auto  evolve_par = [&](tbb::blocked_range2d<int> &r)->void{

	for (auto T=r.rows().begin(); T != r.rows().end(); ++T){
		mcmc::individual  tY;
		bool within_PriorRange; 

		for (auto j=r.cols().begin(); j!= r.cols().end(); ++j)
			{
			within_PriorRange = true;  

			// Propose new individual. 
			tY = my_step->propose( tP_1[T] , j);

			// Accept / Reject 
			// Check if its values are within PriorRange 
			for (auto i=0; i < Nvars; ++i)
        	                if( (tY.Vars[i] > PriorRange[i][1]) || (tY.Vars[i]< PriorRange[i][0])  )
					{
					within_PriorRange = false; 
					}

		
			if (within_PriorRange == false ) {// Reject proposed value  
				tP_2[T][j] =  tP_1[T][j] ; // Keep old value 
			} else {
				tY.loglkhood= logP(tY); 
				if ( my_mh->accept( tY ,tP_1[T][j], beta[T]) ){ 
					tP_2[T][j] = tY; // Accept proposed value 
				}else {
					tP_2[T][j] =  tP_1[T][j]; // Keep old value 
				}
			}	
			} // end of for j loop 
	

		} // End of for T loop 
		}; // END lambda function 


	static tbb::affinity_partitioner ap; 
	tbb::parallel_for(tbb::blocked_range2d<int>(0,dim_T,0,Npop), evolve_par , ap);


	// Update population -- SLOW operation, would love to avoid it.  
	update();  // Copies values of tP_1 <-- tP_2





	
}





template<class LogLkhood,class stepper, class mh_ratio >
void mcmc_algo_pcx_pt<LogLkhood,stepper,mh_ratio>::sample_single(int &Nsteps) /**< Sample till death */ 
	{
	if (flname_set == true ){
	for (auto i=0; i < Nsteps; ++i)
		{
		// Vocal evolution
                if(i%int(0.05*Nsteps)==0) 
                        cout<<"Percentage completed  --> "<<100.0*double(i)/double(Nsteps)<<"%" << endl; 
		evolve_single(); 
		swap(tP_1, i,beta);
		write ();
		}
	}else {
		std::cerr<< "Either flname_out or algo not set, yet you request sampling. WTF?  Aborting ..." << std::endl; 
		DEBUG (flname_set);
		throw 0;

	}
	
}


template<class LogLkhood,class stepper, class mh_ratio >
void mcmc_algo_pcx_pt<LogLkhood,stepper,mh_ratio>::sample(int &Nsteps) /**< Sample till death */ 
	{
	if (flname_set == true ){
	for (auto i=0; i < Nsteps; ++i)
		{
		// Vocal evolution
                if(i%int(0.05*Nsteps)==0) 
                        cout<<"Percentage completed  --> "<<100.0*double(i)/double(Nsteps)<<"%" << endl; 
		evolve(); 
		swap(tP_1, i,beta);
		write();
		}
	}else {
		std::cerr<< "flname_out not set, aborting ..." << std::endl; 
		DEBUG (flname_set);
		throw 0;

	}
	
}





}
}
#endif
