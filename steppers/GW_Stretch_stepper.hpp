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


#ifndef _GW_Stretch_stepper_
#define _GW_Stretch_stepper_

#include <string>

#include "../globals.h"
#include "../base/stepper.hpp"


namespace mcmc{
namespace steppers{

/** 
	Ensemble samplers with affine invariance, Goodman&Weare stepper, DOI:  10.2140/camcos.2010.5.65, 2010.  

*/
template <class T> 
class GW_Stretch_stepper: public mcmc::base::stepper_base<T, GW_Stretch_stepper<T> >  {

	private:
 
		T propvecx; 
		int Nvars; // Total number of variables  
		int Npop;  // Population 

		std::uniform_int_distribution<> dis_int; // (0, Npop-1);	// Integer in the range 0,...,Nwalkers-1.
		std::uniform_real_distribution<> Zunif_real; // (0., 1.);


		// -------------- TODO: Wrap inside a structure outside of this class-------- 
		// Auxiliary variables for g(z) random distribution. 
		double a;			
		double inva;				
		double _2inva; // Optimized for faster evaluation. 
						
		/* Random number created from g(z) */  	
		double Zzrandpdf();				
		// --------------------------------------------------------------------------


		std::string name; 

	public: 
		GW_Stretch_stepper(int &_Nvars, int &_Npop);
		/** Custom set of alpha variable */
		void set_alpha(const double &_a){a=_a;}; 

		/* Stepper function */
		T propose(std::vector<T> &X, int &j) ;
		

//		double get_Z(){return ZZ;}


		std::string get_name(){return name;}; 	

};



template <class T> 
double GW_Stretch_stepper<T>::Zzrandpdf() {
	double xZran=Zunif_real(gen);
	double xZran_sq = xZran*xZran;
	// Reference formula: 
	//(1. + x*(2.0*(a-1.) + (1. + (-2. + a)*a)*x))/a;

        return
		inva + 2.*xZran -_2inva*xZran- 2.*xZran_sq +xZran_sq*a + inva*xZran_sq;
}



template <class T>
GW_Stretch_stepper<T>::GW_Stretch_stepper(int &_Nvars, int &_Npop):Nvars (_Nvars), Npop(_Npop) {
	
	name = "GW_Stretch_stepper";

	gen.seed(rd());

	// ===========>      <==============================================================
        a=2.0;//related to random distribution that creates the affine step. 
	// ===========>      <============================================================== 

	// Helping variables for optimized evaluation: 
	inva  = 1./a;
	_2inva = 2.*inva;


	// Set up distributions. 
	std::uniform_real_distribution<>::param_type newParams1(0.0,1.0);
	Zunif_real.param(newParams1);

	std::uniform_int_distribution<>::param_type newParams2(0, Npop-1);
	dis_int.param(newParams2);
			
	// Set up proper size of proposed vector	
	propvecx.Vars.resize(Nvars); 

}



template <class T>
T GW_Stretch_stepper<T>::propose(std::vector<T> &X, int &j) 
	{
	
	// Random walker kk != j 
	int kk; 
	double ZZ;

         do{
            kk=dis_int(gen);// -1 to include the zzero index value (unif(n)>0)
	} while(kk==j); /** Defining a random integer different than j */

        ZZ=Zzrandpdf();// Random value. 


	// Attention: same ZZ for all variables (vector coordinates).
        for(auto i=0;i<Nvars;i++)
		{
                propvecx.Vars[i]=X[kk].Vars[i]+ZZ*(X[j].Vars[i]-X[kk].Vars[i]);
		}
	// The class T should also store the value of the random Z to be used in MH ratio. 
	propvecx.Z=ZZ; 
	
	// loglkhood not assigned yet. 
	return propvecx;

}














}// end steppers namespace 
}// end mcmc namespace 



#endif
