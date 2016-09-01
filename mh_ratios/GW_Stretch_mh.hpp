
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

#ifndef _GW_Stretch_mh_
#define _GW_Stretch_mh_

#include <cmath>

#include "../globals.h"
#include "../individual.hpp"
#include "../base/mh.hpp"


namespace mcmc{
namespace mh_ratios{

class GW_Stretch_mh: public mcmc::base::mh_base< mcmc::individual, GW_Stretch_mh>{
	private:
		double tNvars; 
		double logmh_ratio(mcmc::individual  &Ytp1, mcmc::individual &Xt); /**< MH ratio */


	public:
		GW_Stretch_mh()
			{
			std::uniform_real_distribution<>::param_type newParams1(0.0,1.);
			mcmc::unif_real.param(newParams1);
			} 


		/**< MH accept_reject for detailed balance */
		bool accept(mcmc::individual  &Ytp1, mcmc::individual &Xt);


};


double GW_Stretch_mh::logmh_ratio ( individual &Ytp1, individual &Xt)
	{
	tNvars = static_cast<double> (Ytp1.Vars.size());
	return  
		(tNvars - 1.)*std::log(Ytp1.Z) + Ytp1.loglkhood - Xt.loglkhood;
}




bool GW_Stretch_mh::accept(mcmc::individual  &Ytp1, mcmc::individual &Xt)
	{		
	double trialrandom=std::log(mcmc::unif_real(gen)); // Proposed random value. 
	if (trialrandom <= logmh_ratio( Ytp1, Xt) ){
		return true; 
	}else {
		return false; 
	}


}






}// end mh_ratios 
}// End of namespace mcmc 
#endif 
