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


#ifndef _pop_base_
#define _pop_base_

#include <vector>
#include "../individual.hpp"

namespace mcmc{
namespace base{

template <class T>
class pop_base{ 

	
	public: 
		virtual ~pop_base(){}

	
		
		virtual int get_PriorRange()=0;
		virtual int get_Nvars()=0;
		virtual int get_dim_T()=0;

		/*
			Get instantatiated vector of individuals 
		*/
		virtual T get_pop() =0; 
		

};



} // end base
} // End mcmc 
#endif
