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




/*
	// ################################### MORE #####################################
	Short description: this class allows creation of a specific algorithm, without the user to have to define each individually / manually. 

*/


#ifndef _mcmc_algo_
#define _mcmc_algo_


// Put all these somewhere else. 
#include "stepper.hpp"
#include "DE_stepper.hpp"
#include "DES_stepper.hpp"
#include "GW_Stretch_stepper.hpp"

/*
// NON working algorithms as yet 
#include "GW_Walk_stepper.hpp"
#include "PCX_stepper.hpp"
*/

namespace mcmc{



enum stepper_type{
	DE,
	DES,
	GW_Stretch,
	GW_Walk_mcmc,	// To be implemented.
	PCX		// To be implemented.
//	, HMCMC	 	// To be implemented. 
};



/*
	This is implementation of Strategy design pattern in order to select between available steppers. Works better than Builder that implements a single strategy. 

*/
template <class T> 
class mcmc_algo{
	private: 
		/* Specific stepper to be implemented */
		stepper_base<T> * mystepper; 
		bool initiated; 
		
		// Nvariables and individuals in population 
		int Nvars; 
		int Npop; 


	public: 
		// Constructor 
		mcmc_algo(int &_Nvars, int &_Npop); 
		/* No copy / assign, protection for pointer memory assignment */
		mcmc_algo(const mcmc_algo<T> &G) = delete; // No copy
		mcmc_algo& operator=(mcmc_algo<T> &G) = delete; // No assign

		~mcmc_algo(); 
		void set_stepper(stepper_type &Type); 

		T propose(vector<T> &X,int &j);


		int get_Nvars(){return Nvars;};
		int get_Npop(){return Npop;};

};



template<class T>
mcmc_algo<T>::mcmc_algo( int &_Nvars, int &_Npop): Nvars(_Nvars), Npop(_Npop)
	{
	mystepper = nullptr;
	initiated = false; 
}



template<class T>
mcmc_algo<T>::~mcmc_algo()
	{
	if (initiated == true)
		delete mystepper; 
}

template<class T>
void mcmc_algo<T>::set_stepper(stepper_type &Type)
	{

	if (Type == DE){
		mystepper = new DE_stepper<T>(Nvars,Npop);
		initiated = true;	

	} else if (Type == DES){	
		mystepper = new DES_stepper<T>(Nvars,Npop);
		initiated = true;	

	}else if (Type == GW_Stretch){
		mystepper = new GW_Stretch_stepper<T>(Nvars,Npop);
		initiated = true;	
	}
/*
	}else if (Type == GW_Walk){
		mystepper = new GW_Walk_stepper<T>(Nvars,Npop);
		initiated = true;	
	
	}else if (Type == PCX){
		mystepper = new PCX_stepper<T>(Nvars,Npop);
		initiated = true;	

	}
*/
}

template<class T>
T mcmc_algo<T>::propose (vector<T> &X, int &j)
	{
	return 
		mystepper -> propose (X,j);

}




} // END of namespace mcmc 

#endif 
