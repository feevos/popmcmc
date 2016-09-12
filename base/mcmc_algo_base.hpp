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

#ifndef _mcmc_algo_base_
#define _mcmc_algo_base_

#include "../individual.hpp"

namespace mcmc{
namespace base{

class mcmc_algo_base{
	public:
		virtual ~mcmc_algo_base(){}

		virtual void set_flname_out(string & _flname_out, int _write_binary)=0; /**< Set output filename */
		
		virtual void sample_single(int &Nsteps)=0;	/**< Sample using one thread, without TBB */ 
		virtual void sample(int &Nsteps)=0; 		/**< Sample till Nsteps in parallel */ 

};




}
}
#endif
