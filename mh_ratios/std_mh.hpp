
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

#ifndef _std_mh_
#define _std_mh_


#include "../globals.h"
#include "../individual.hpp"


namespace mcmc{
namespace mh_ratios{

class std_mh { 
	private: 
		double logmh_ratio(mcmc::individual  &Ytp1, mcmc::individual &Xt); 
		double logmh_ratio(mcmc::individual  &Ytp1, mcmc::individual &Xt,double &beta_T); 

	public:

		std_mh()
			{
			std::uniform_real_distribution<>::param_type newParams1(0.0,1.);
			mcmc::unif_real.param(newParams1);
			} 

		
		/**< MH accept_reject for detailed balance */
		bool accept(mcmc::individual  &Ytp1, mcmc::individual &Xt);
		bool accept(mcmc::individual  &Ytp1, mcmc::individual &Xt, double &beta_T);


};

double std_mh::logmh_ratio ( mcmc::individual &Ytp1, mcmc::individual &Xt)
	{
	// This is standard mh ratio 
	return  
		Ytp1.loglkhood - Xt.loglkhood;
}



double std_mh::logmh_ratio ( mcmc::individual &Ytp1, mcmc::individual &Xt, double &beta_T)
	{
	// This is standard mh ratio 
	return  
		beta_T * (Ytp1.loglkhood - Xt.loglkhood);
}



bool std_mh::accept(mcmc::individual  &Ytp1, mcmc::individual &Xt)
	{
                      
	double trialrandom=std::log(mcmc::unif_real(gen)); // Proposed random value. 
	if (trialrandom <= logmh_ratio( Ytp1, Xt) ){
		return true; 
	}else {
		return false; 
	}


}


bool std_mh::accept(mcmc::individual  &Ytp1, mcmc::individual &Xt, double &beta_T)
	{                    
	double trialrandom=std::log(mcmc::unif_real(gen)); // Proposed random value. 
	if (trialrandom <= logmh_ratio( Ytp1, Xt, beta_T) ){
		return true; 
	}else {
		return false; 
	}

}




}
}// End of namespace mcmc 
#endif 
