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


#ifndef _DES_stepper_
#define _DES_stepper_

#include <string>


#include "../globals.h"
#include "../base/stepper.hpp"


namespace mcmc{
namespace steppers{


/**
	DES - Streens, 2002 - MCMC Sampling usind direct search optimization 	

	// Main difference with DE - MCMC is the coefficient gamma in front of the difference: here it is gaussian. 
*/

template <class T>
class DES_stepper: public mcmc::base::stepper_base<T, DES_stepper<T> >  {

	private:
 
		T propvecx; 
		int Nvars; // Total number of variables  
		int Npop;  // Population 

        	std::normal_distribution<> gauss;// (0, log(2.));  
		std::uniform_int_distribution<> dis_int; // (0, Npop-1);	// Integer in the range 0,...,Nwalkers-1.
		

		std::string name; 

	public: 
		DES_stepper(const int &_Nvars, const int &_Npop);
		T propose(std::vector<T> &X, int &j) ;
	
	std::string get_name(){return name;}; 	

};



template <class T>
DES_stepper<T>::DES_stepper(const int &_Nvars, const int &_Npop):Nvars (_Nvars), Npop(_Npop) {
	
	name = "DES_stepper";

	// Set up distributions. 
	gen.seed(rd());

	std::normal_distribution<>::param_type newParams1(0.0,log(2.));
	gauss.param(newParams1);

	std::uniform_int_distribution<>::param_type newParams2(0, Npop-1);
	dis_int.param(newParams2);
			
	// Set up proper size of proposed vector, and gamma step 	
	propvecx.Vars.resize(Nvars); 

}


/* Proposes new value for individual j */
template <class T>
T DES_stepper<T>::propose(std::vector<T> &X, int &j) 
	{
	int R1,R2; 
	do{
		R1=dis_int(gen);
	}while (R1==j);
	do{
		R2=dis_int(gen);
	}while (R2 == R1 || R2 == j);

        for(auto i=0;i<Nvars;i++)
                propvecx.Vars[i]=X[j].Vars[i]+gauss(gen)*(X[R2].Vars[i]-X[R1].Vars[i]);

	// loglkhood has not been assigned at this stage. 
	return propvecx;

}














}
}



#endif
