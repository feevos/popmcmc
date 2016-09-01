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
#include <cstdlib>


namespace mcmc{

using namespace std; 



/* Abstract strategy interface. */
template <class T>
class stepper_base{
	public:
		virtual T propose( vector<T> &X, int &j) = 0; 

};





// I don't think this is necessarY! 
/* Implimentation of a specific stepper  */
template <class T>
class stepper{
	private: 
	/* Some specific stepper */
	stepper_base<T> * spec_stepper; 	

	public: 
		explicit stepper(stepper_base<T> * _stepper): spec_stepper(_stepper){};
		T propose(vector<T> &X, int &j) 
			{
			return 
				spec_stepper->propose(X,j);
			}



};









} // end namespace 

#endif
