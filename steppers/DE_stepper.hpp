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


#ifndef _DE_stepper_
#define _DE_stepper_

#include <string>

#include "../globals.h"


namespace mcmc{
namespace steppers{


/**
	Differential Evolution stepper - Braak, 2006, 	DOI 10.1007/s11222-006-8769-1

*/

template <class T>
class DE_stepper  {

	private:
 
		double epsilon;// (1.e-4);	 
		double gamma; //= 2.38 *pow(2.0* static_cast<double>(Nparams),-0.5 );
		T propvecx; 
		int Nvars; // Total number of variables  
		int Npop;  // Population 

		// These should go in globals section 
        	std::uniform_real_distribution<> Uprop;// (-epsilon, epsilon);// from a=1 to b=2. Call with: "dis(gen)" without " " characters.
		std::uniform_int_distribution<> dis_int; // (0, Npop-1);	// Integer in the range 0,...,Nwalkers-1.
		
		
		std::string name; 

	public: 
		DE_stepper(const int &_Nvars, const int &_Npop);
		T propose(std::vector<T> &X, int &j) ;

		std::string get_name(){return name;}
};



template <class T>
DE_stepper<T>::DE_stepper(const int &_Nvars, const int &_Npop):Nvars (_Nvars), Npop(_Npop) {	

	name = "DE_stepper"; 

	gen.seed(rd());
	epsilon = 1.e-4;

	// Set up distributions. 

	std::uniform_real_distribution<>::param_type newParams1(-epsilon,epsilon);
	Uprop.param(newParams1);

	std::uniform_int_distribution<>::param_type newParams2(0, Npop-1);
	dis_int.param(newParams2);
			
	// Set up proper size of proposed vector, and gamma step 	
	propvecx.Vars.resize(Nvars); 
	gamma = 2.38 *pow(2.0* static_cast<double>(Nvars),-0.5 );

}


template <class T>
T DE_stepper<T>::propose(std::vector<T> &X, int &j) 
	{
	
	int R1, R2; 
	do{
		R1=dis_int(gen);
	}while (R1==j);
	do{
		R2=dis_int(gen);
	}while (R2 == R1 || R2 == j);

        for(auto i=0;i<Nvars;i++)
                propvecx.Vars[i]=X[j].Vars[i] + gamma*(X[R2].Vars[i]-X[R1].Vars[i]) + Uprop(gen);

	// loglkhood has not been assigned at this stage. 
	return propvecx;

}














}
}



#endif
