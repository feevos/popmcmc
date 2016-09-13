#include "../../mcmc.hpp"
#include <cmath>
#include <chrono>

using namespace std;



// This is how you define a Functor for the logLikelihood. ATTENTION: it *must* be copy constructable. 
class ros_logP{

	public:
		double operator()(mcmc::individual &X){
			vector<double> theta = X.Vars; 
			return
        	        // In this example Rosenbrock density. 
                	-(100.0*pow(theta[1]- pow(theta[0],2) ,2) + pow(1.- theta[0],2) )/20.0;
	}
};



int main(){


	int Npop=1000;
	int Nvars=2;
	int Nsample=2000; 


	vector<vector<double> > PriorRange(Nvars); 
	for (auto i=0; i < Nvars; ++i)
		{
		PriorRange[i].resize(2); 
		PriorRange[i][0]=-30.0; 
		PriorRange[i][1]= 30.0; 
		}

	ros_logP logP;
	
	// Constructor of population mcmc sampler. 
	mcmc::popmcmc<ros_logP> mysampler(logP,PriorRange);



	// Available methods, 

// ------------- Standard samplers -- must be combined with standard init functions (i.e. non - tempered). 

	//mysampler.set_algo(mcmc::GW_Stretch);
	//mysampler.set_algo(mcmc::DE);
	mysampler.set_algo(mcmc::DES);
	//mysampler.set_algo(mcmc::PCX);

	// Characteristic init functions. See file mcmc.hpp for all constructors.  
	mysampler.init(Npop);   
// -------------------------------------------------------------------------------------------
	
// ---------------- Parallel Tempering algorithms --------------------------------------------
	//mysampler.set_algo(mcmc::GW_Stretch_PT);
	//mysampler.set_algo(mcmc::DE_PT);
	//mysampler.set_algo(mcmc::DES_PT);
	//mysampler.set_algo(mcmc::PCX_PT);
	
	// Characteristic init functions. 10: dimension of temperatures, 30: Nswap, propose swap every Nswap iterations. 
	//mysampler.init(Npop,20,30);  
// ------------------------------------------------------------------------------------------	


	// Set filename for output. 
	//string flname_out = "caramba_PCX"; 
	//string flname_out = "caramba_PCX_PT"; 
	//string flname_out = "caramba_GWS_1"; 
	//string flname_out = "caramba_GWS_8_Pt"; 
	//string flname_out = "caramba_GWS"; 
	//string flname_out = "caramba_DES_8_shuffle"; 
	//string flname_out = "caramba_GWS_8_shuffle"; 
	string flname_out = "caramba_DES_8"; 
	//string flname_out = "caramba_DES_PT_1"; 
	//string flname_out = "caramba_DE"; 
	//string flname_out = "caramba_DE_PT"; 
	//string flname_out = "caramba_DES_PT"; 
	mysampler.set_flname_out(flname_out);


	
	auto t1=std::chrono::high_resolution_clock::now(); // Time the evaluation. 
	// ********* SAMPLE ************
	//mysampler.sample_single(Nsample); 
	mysampler.sample(Nsample);  // parallel version - the number corersponds to number of threads. 
	// ****************************
	auto t2=std::chrono::high_resolution_clock::now();
        auto Dt = t2-t1;

	// Output on screen time taken to sample. 
	std::cout << "Time elapsed (msec): ";
        std::cout <<    std::chrono::duration_cast<std::chrono::milliseconds>(Dt).count()<<  endl;
	std::cout<< "mcmc chains written in file: " << flname_out << std::endl; 

}
