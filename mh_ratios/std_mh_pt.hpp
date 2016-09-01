
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

#ifndef _std_mh_pt_
#define _std_mh_pt_


#include "../globals.h"
#include "../individual.hpp"


namespace mcmc{
namespace mh_ratios{

class std_mh_pt{ 
	private: 
		double logmh_ratio(mcmc::individual  &Ytp1, mcmc::individual &Xt, double &beta_T); 
		int Npop; 
		int Nswap; 
		int Nbeta; 

		
		int TT,TT2,jj,jj2;  /**< helping variables for efficiency */
		double logM, swaplogratio; 

		uniform_int_distribution<> dist_int_beta;  /**< Gives a random integer value for selection of temperature */
		uniform_int_distribution<> dist_int; 	/**< Selects random population individuals */
		
	public:

		std_mh_pt(vector<double> &_beta, int &_Npop, int _Nswap):  Npop(_Npop), Nswap(_Nswap)
			{
			Nbeta = static_cast<int> (_beta.size()); 
			std::uniform_real_distribution<>::param_type newParams1(0.0,1.);
			mcmc::unif_real.param(newParams1);

			std::uniform_int_distribution<>::param_type newParams2(0,Nbeta-2);
			dist_int_beta.param(newParams2); 

			std::uniform_int_distribution<>::param_type newParams3(0,Npop-1);
			dist_int.param(newParams3); 
			} 

		
		/**< MH accept_reject for detailed balance */
		bool accept(mcmc::individual  &Ytp1, mcmc::individual &Xt, double & beta_T);
		void swap(std::vector<std::vector<mcmc::individual> > &P,  int &i,vector<double> &beta);

};


double std_mh_pt::logmh_ratio ( mcmc::individual &Ytp1, mcmc::individual &Xt, double &beta_T)
	{
	// This is standard mh ratio 
	return  
		beta_T * (Ytp1.loglkhood - Xt.loglkhood);
}




bool std_mh_pt::accept(mcmc::individual  &Ytp1, mcmc::individual &Xt, double &beta_T)
	{                    
	double trialrandom=std::log(mcmc::unif_real(gen)); // Proposed random value. 
	if (trialrandom <= logmh_ratio( Ytp1, Xt, beta_T) ){
		return true; 
	}else {
		return false; 
	}

}



void std_mh_pt::swap(std::vector<std::vector<mcmc::individual> > &P,  int &i, vector<double> &beta){

	if(i%Nswap==0){ // Try it every Nswap steps. 

	// conjugate temperatures to swap 
	TT= dist_int_beta(gen); 
	TT2=TT+1;

	// Choose two chains at random to swap  
	jj = dist_int(gen);
	
	do {		
	jj2 = dist_int(gen);
	}while (jj2==jj);



	swaplogratio=log(unif_real(mcmc::gen) );
	logM =  (beta[TT] - beta[TT2]) * ( P[TT2][jj].loglkhood - P[TT][jj2].loglkhood );		


	if(swaplogratio<=logM){ // swap temperatures. 
	//std::swap(beta[TT2],beta[TT]);  
	std::swap(P[TT2][jj],P[TT][jj2]);  
	}
	};// End of IF Swap proposal. 


}




}//end of mh_ratios
}// End of namespace mcmc 
#endif 
