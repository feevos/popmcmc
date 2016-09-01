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

#ifndef _mh_
#define _mh_

namespace mcmc{
namespace base{


// Abstract strategy interface  - static polymorphism 
template <class T, class spec_mh>
class mh_base{
	public:
		// Metropolis hastis ratio 
		bool accept(T &Ytp1, T &Xt){
			return 
				static_cast< spec_mh * > (this) -> accept(Ytp1,Xt); 
		} 

};

}
}

#endif
