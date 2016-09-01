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


#ifndef _stepper_
#define _stepper_


#include <vector> 

namespace mcmc{
namespace base{

/**
	Abstract strategy - static polymorphism 
*/
template <typename T, class spec_stepper> // I DONT NEED the T template: individual will be fixed for all 
class stepper_base{
	public:
		T propose( std::vector<T> &X, int &j)
			{
			return static_cast<spec_stepper *> (this) -> propose(X,j); 
			} 

};



}
} // end namespace 

#endif
