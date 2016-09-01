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



#ifndef _DES_GWS_PCX_stepper_
#define _DES_GWS_PCX_stepper_

#include <vector>
#include <set>

#include "../globals.h"
#include "../individual.hpp"
#include "DES_stepper.hpp"
#include "GW_Stretch_stepper.hpp"
#include "PCX_stepper.hpp"


namespace mcmc{
namespace steppers{



template <class T>
class DES_GWS_PCX_stepper{
	private:
		mcmc::steppers::DES_stepper<T>		* 	stepper1; 
		mcmc::steppers::GW_Stretch_stepper<T> 	*	stepper2;
		mcmc::steppers::PCX_stepper<T>		*	stepper3;

		std::vector<double> weights; 
		int Nweights; 
		vector<double> cumprob; // Vector of cumulative probability

	public: 
		
		DES_GWS_PCX_stepper(vector<double> &_weights, std::vector<std::vector<double> > &_PriorRange, const int &_Nvars, const int &_Npop, int _mu=5);
		~DES_GWS_PCX_stepper(); 

		T propose( std::vector<T> &X, int &j) ;



};


template <class T>
DES_GWS_PCX_stepper<T>::DES_GWS_PCX_stepper(vector<double> &_weights, std::vector<std::vector<double> > &_PriorRange, const int &_Nvars, const int &_Npop, int _mu):weights(_weights){
	stepper1 = new DES_stepper<T>(_Nvars,_Npop);
	stepper2 = new GW_Stretch_stepper<T>(_Nvars,_Npop); 
	stepper3 = new PCX_stepper<T>(_PriorRange,_Nvars,_Npop,_mu);
	
	Nweights = static_cast<int> (weights.size());

	// Normalize weights: 
	double sum=0.0; 
	for (auto &it: weights)
		sum += it; 
	
	for (auto &it: weights)
		it /= sum; 

	cumprob.resize(Nweights);
	for (auto i=0; i < Nweights; ++i)
		{
		cumprob[i]=0.0;
		for (auto j=0; j <= i; ++j)
			cumprob[i] += weights[j]; 
		}		

	uniform_real_distribution<>::param_type newParams(0.0, 1.0);
	mcmc::unif_real.param(newParams);

}




template <class T>
DES_GWS_PCX_stepper<T>::~DES_GWS_PCX_stepper(){
	delete stepper1;
	delete stepper2;
	delete stepper3;
}



template <class T>
T DES_GWS_PCX_stepper<T>::propose( std::vector<T> &X, int &j) {


	double u = mcmc::unif_real(mcmc::gen);
	
	if (  u < cumprob[0]){
		return stepper1->propose(X,j);
	} else if ( u < cumprob[1] ){
		return stepper2->propose(X,j);
	}  else if ( u < cumprob[2] ){
		return stepper3->propose(X,j);
	};
		
}







}
}
#endif
