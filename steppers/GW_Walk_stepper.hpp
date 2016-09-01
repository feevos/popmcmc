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





// ***********************************************************************************
// ***********************************************************************************
// ***********************************************************************************
// ##########################  NEEDS DEBUGGING - INCOMPLETE ##########################
// ***********************************************************************************
// ***********************************************************************************
// ***********************************************************************************





#ifndef _GW_Walk_stepper_
#define _GW_Walk_stepper_

#include "../base/stepper.hpp"


namespace mcmc{


/** 
	Ensemble samplers with affine invariance, Goodman&Weare stepper, DOI:  10.2140/camcos.2010.5.65, 2010.  

*/
class GW_Walk_stepper: public stepper_base  {

	private:
 
		vector<double> propvecx; 
		int Nvars; // Total number of variables  
		int Npop;  // Population 

		std::random_device rd;
		std::mt19937 gen;

		std::uniform_int_distribution<> dis_int; // (0, Npop-1);	// Integer in the range 0,...,Nwalkers-1.
		std::normal_distribution gauss;

		/* Number of parents for the evaluation of the Walk proposal */
		int mu; 

		double tempSum;
		vector<double> W; 

	public: 
		GW_Walk_stepper(const int &_Nvars, const int &_Npop);
		/** Total number of parents  */
		void set_mu(const int &_mu){mu=_mu;}; 

		/* Stepper function */
		vector<double> propose(vector<vector<double> > &X, int &j) ;
		

};




GW_Walk_stepper::GW_Walk_stepper(const int &_Nvars, const int &_Npop):Nvars (_Nvars), Npop(_Npop) , gen(rd()){
	

	// ===========>      <==============================================================
        mu=5;//related to random distribution that creates the affine step. 
	// ===========>      <============================================================== 

	// Set up distributions. 
	normal_distribution<>::param_type newParams1(0.0,0.1);
	gauss.param(newParams1);

	uniform_int_distribution<>::param_type newParams2(0, Npop-1);
	dis_int.param(newParams2);
			
	// Set up proper size of proposed vector	
	propvecx.resize(Nvars); 
	W.resize(Nvars);

}




vector<double>  GW_Walk_stepper::propose(vector<vector<double> > &X, int &j) 
	{
	// Evaluate a set of parent vectors 
	set<int> idxs;
	int kk; 
	do {
		kk =   dis_int(gen);	
	} while (kk == j);

	idxs.insert(kk);

	for (auto i=1; i<mu; ++i)
		{
		do{
			kk = dis_int(gen);	
		}while (kk == j ||  (std::find(idxs.begin(),idxs.end(),kk)) );

		idxs.insert(kk);

		}

	// Calculate mean vector 
	vector<double> Xs(Nvars);
	double denum = static_cast<double> (mu);
	tempSum;
	for (auto i=0; i<Nvars; ++i)
		{	
		tempSum=0.0;
		for (auto &it:idxs)
			tempsum += X[it][i];	
		Xs[i] = tempsum/denum;
		}

	// Calculate proposal vector W 
	vector<double> ZZrands(mu);
	for (auto i=0; i< mu; ++i)
		ZZrands[i]=gauss(gen);


// ***************************** INCOMPLETE  *********************************************










	// Proposed vector: 
        for(auto i=0;i<Nvars;i++)
               propvecx[i]=X[j][i]+W[i];


	return null_ptr;
	//return propvecx;

}















}



#endif
