#include "../../mcmc.hpp"
//#include "../../steppers/DES_GWS_PCX_combo_stepper.hpp"
#include "../../steppers/DES_PCX_combo_stepper.hpp"
#include "../../mh_ratios/std_mh_pt.hpp"

#include <cmath>
#include <chrono>
#include <vector>

using namespace std;
using namespace mcmc; 



class gaussian{
	private: 
		double mu;
		double sigma; 

	public: 
		gaussian(){}
		void set_params(double &_mu, double &_sigma)
			{
			mu = _mu;
			sigma = _sigma; 
			}
		gaussian(double &_mu, double &_sigma):mu(_mu),sigma(_sigma){}
		double operator()(double &x){
			
			return 
				-0.5*log(2.0*M_PI*sigma*sigma)
				-0.5*pow((x-mu)/sigma,2);
		}
};

class mult_gaussian{
	private: 

		int Ngauss;
		vector<double> weights;
		vector<double> sigmas;
		vector<double> means; 

		vector<gaussian> mgauss;  


	public: 
		mult_gaussian(vector<double> &_weights, vector<double> &_sigmas, vector<double> &_means)
			{

			Ngauss = _weights.size(); 

			weights = _weights;
			sigmas  = _sigmas;
			means	= _means; 

			// Normalize weights: 
			double sum=0.0; 
			for (auto &it: weights)
				sum += it; 
			for (auto &it: weights)
				it /= sum; 

			mgauss.resize(weights.size());
			for (auto i=0; i < weights.size(); ++i)
				mgauss[i].set_params(means[i],sigmas[i]);

			}

		double operator() (mcmc::individual &X)
			{
			double sum = 0.0; 
			for (auto i=0; i < Ngauss; ++i)
				sum += mgauss[i](X.Vars[i]);

			return sum; 
			}



};




int main(){


	
	mcmc::gen.seed(mcmc::rd());	

	int Npop=1000;
	int Nvars=4;
	int Nsample=1000; 

	vector<double> means {-2.0, 1.0, 1.5, 1.5};
	vector<double> sigmas {.2, 0.1, .15, 0.25};
	vector<double> weights {1.,2.,2.,1.};

	
	std::vector<vector<double> > PriorRange(Nvars); 
	for (auto i=0; i < Nvars; ++i)
		{
		PriorRange[i].resize(2); 
		PriorRange[i][0]=-5.; 
		PriorRange[i][1]=5.; 
	
		}

	mult_gaussian logP(weights,sigmas,means);
	
	
	mcmc::popmcmc<mult_gaussian> mysampler(logP,PriorRange);

	//mysampler.set_algo(mcmc::PCX);
	//mysampler.set_algo(mcmc::GW_Stretch);
	//mysampler.set_algo(mcmc::DES);

	
	//mysampler.init(Npop);  



/*
	// ********************************** NEEDS DEBUGGING DOES NOT WORK ***************************************
	// -------------------------------------------------------------------------------------------------------+
	//													//|
	//	Constructing custom stepper(DES_PCX_combo) 							//|
	//													//|
														//|
	// a. Set up population: 										//|
	population::pop<mult_gaussian> mypop (logP,PriorRange);							//|
	mypop.init(Npop);											//|
														//|
														//|
	vector<double> sweights {0.5,0.5};									//|
	typedef  mcmc::steppers::DES_PCX_stepper<mcmc::individual> stepper; 					//|
	stepper mystepper(sweights,PriorRange,Nvars,Npop);							//|
	mcmc::mh_ratios::std_mh my_mh; 										//|
	// Custom Algo												//|
	algos::mcmc_algo_pcx<mult_gaussian,stepper,mh_ratios::std_mh> myalgo(mypop,&mystepper,&my_mh); 		//|
	mysampler.set_algo(&myalgo);										//|
														//|
	// -------------------------------------------------------------------------------------------------------+
*/




	
	// -------------------------------------------------------------------------------------------------------+
	//													//|
	//	Constructing custom stepper: 									//|
	//													//|
														//|
	// a. Set up population: 										//|
	int dim_T=10;												//|
	population::pop_pt<mult_gaussian> mypop (logP,PriorRange,dim_T);					//|
	mypop.init(Npop); // (10,30)										//|
														//|
														//|
	vector<double> sweights {0.5,0.5};									//|
	typedef  mcmc::steppers::DES_PCX_stepper<mcmc::individual> stepper; 					//|
	stepper mystepper(sweights,PriorRange,Nvars,Npop);							//|

	
	int Nbeta = 10; 
	int Nswaps = 30; 
	std::vector<double> beta(Nbeta);  
	for (auto i=0; i < Nbeta; ++i)
		beta[i] = double (i) / double (Nbeta-1); 

	mcmc::mh_ratios::std_mh my_mh; 								//|
	// Custom Algo													//|
	algos::mcmc_algo_pcx_pt<mult_gaussian,stepper,mh_ratios::std_mh> myalgo(mypop,&mystepper,&my_mh,Nswaps);	//|
	mysampler.set_algo(&myalgo);											//|
															//|
	// -------------------------------------------------------------------------------------------------------+



	




	
	//string flname_out = "mcmc_DES.chains"; 
	//string flname_out = "mcmc_GWS.chains"; 
	//string flname_out = "mcmc_DES_PCX_combo.chains"; 
	string flname_out = "mcmc_DES_PCX_combo_PT.chains"; 
	//string flname_out = "mcmc_PCX_lambda_10.chains"; 
	mysampler.set_flname_out(flname_out);

	
	auto t1=std::chrono::high_resolution_clock::now();
	//mysampler.sample_single(Nsample); 
	mysampler.sample(Nsample);  // parallel version 
	auto t2=std::chrono::high_resolution_clock::now();
        auto Dt = t2-t1;

        std::cout <<    std::chrono::duration_cast<std::chrono::milliseconds>(Dt).count()<<  endl;



}

