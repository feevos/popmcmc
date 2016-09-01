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

#ifndef _GW_Stretch_algo_
#define _GW_Stretch_algo_

#include "../steppers/GW_Stretch_stepper.hpp"
#include "../mh_ratios/GW_Stretch_mh.hpp"
#include "../population/pop.hpp"
#include "../individual.hpp"




namespace mcmc{
namespace algos{

/**
	Necessary typedefs in order to have class polymorhpism for CRTP 
*/
typedef  mcmc::steppers::GW_Stretch_stepper<mcmc::individual>  GWS_stepper; 
typedef  mcmc::mh_ratios::GW_Stretch_mh GWS_mh_ratio; 

typedef  mcmc::base::stepper_base< mcmc::individual, GWS_stepper> GWS_stepper_base;  /**< stepper_base */
typedef  mcmc::base::mh_base < mcmc::individual, GWS_mh_ratio > GWS_mh_ratio_base;  	/**< mh_base */




template<class LogLkhood>
class GW_Stretch_algo:public  Any< GWS_stepper_base>, Any<GWS_mh_ratio_base>{ /**< This inheritance is use to declare type of CRTP classes as members */
	using  Any< GWS_stepper_base>::Any;  // What's the point? I could just use the standard templates ... :(
	using  Any<GWS_mh_ratio_base>::Any;

	typedef  mcmc::population::population< LogLkhood > GW_pop; 
	

	private: 

		LogLkhood logP;  		/**< LogLikelihood function  */

		GWS_stepper_base * my_step;  	/**< Memory will be allocated */
		GWS_mh_ratio_base *  my_mh; 	/**< Memory will be allocated */


		vector<vector<double> > PriorRange; 
		
		std::vector <  mcmc::individual >  tP_1; /**< Temporary instance of individuals: it is needed for detailed balance, before and after modification */
		std::vector <  mcmc::individual >  tP_2; /**< Temporary instance of individuals: it is needed for detailed balance, before and after modification  */
	
		int Nvars;
		int Npop; 
		
		void update()	  /**< Updates population tP_1 from tP_2 */
			{
			tP_1 = tP_2; 
			}

	public:
		
		std::vector<mcmc::individual> get_pop(){return tP_1;} 

		GW_Stretch_algo(GW_pop &_P ); 

		~GW_Stretch_algo()
			{
			delete my_step; 
			delete my_mh; 
			}


		// ************************ TBB parallelized routines **********************************************
		void evolve(); /**< Evolves the population P for a single generation */ 
//		void evolve(std::vector<double> &beta); // Evolves using parallel tempering  for a single generation  // This will probably go in a different algorithm: I need different population.  
		// *************************************************************************************************

};



template<class LogLkhood>
GW_Stretch_algo<LogLkhood>::GW_Stretch_algo(GW_pop &_P ): logP(_P.logP)
	{
	Nvars 		=  	_P.get_Nvars() ; 
	Npop 		=   	_P.get_Npop() ; 

	PriorRange	=	_P.get_PriorRange(); 

	tP_1 		=	_P.get_pop(); 
	tP_2		=	tP_1; 

	my_step 	= 	new   GWS_stepper(Nvars,Npop);
	my_mh 		= 	new   GWS_mh_ratio ();


}


/**
This is the most important function of the algorithm. It is the SAME independent of algorithmic structure. ***** So it should NOT repeated.***** 

*/

template<class LogLkhood>
void GW_Stretch_algo<LogLkhood>::evolve(){

	mcmc::individual  tY;

	bool within_PriorRange; 
	for (auto j=0; j< Npop; ++j)
		{
		within_PriorRange = true;  

		// Propose new individual. 
		tY = my_step->propose( tP_1 , j);

		// Accept / Reject 
		// ********************************* Wrap this step inside a function? *********************************
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
			if ( my_mh->accept( tY ,tP_1[j]) ){ 
				tP_2[j] = tY; // Accept proposed value 
			}else {
				tP_2[j] =  tP_1[j]; // Keep old value 
			}
		}	
		// *******************************************************************************************************
		}

		// Update population 
		update();  // Copies values of tP_1 <-- tP_2





		


}









}
}
#endif
