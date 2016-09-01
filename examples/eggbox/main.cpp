#include "../../mcmc.hpp"
#include <cmath>
#include <chrono>

using namespace std;


class egg_logP{

	public:
		double operator()(mcmc::individual &X){
			vector<double> theta = X.Vars; 

			double chi=1.;
			for (auto i=0; i < X.Vars.size(); ++i)
				chi *= cos(0.5*X.Vars[i]); 

			return
				pow(chi+2.0,5); 
	}

};


int main(){
	

	
	int Npop=1000;
	int Nvars=10;
	int Nsample=1000; 


	
	vector<vector<double> > PriorRange(Nvars); 
	for (auto i=0; i < Nvars; ++i)
		{
		PriorRange[i].resize(2); 
		PriorRange[i][0]=0.0; 
		PriorRange[i][1]= 10.0*M_PI; 
		}

	egg_logP logP;
	
	mcmc::popmcmc<egg_logP> mysampler(logP,PriorRange);

	//mysampler.set_algo(mcmc::PCX);
	//mysampler.set_algo(mcmc::PCX_PT);
	//mysampler.set_algo(mcmc::DES_PT);
	//mysampler.set_algo(mcmc::DES);


	mysampler.init(Npop,10,30);  
	//mysampler.init(Npop);   


	string flname_out = "egg_DES_dim10_PT"; 
	//string flname_out = "egg_DES_dim2"; 
	//string flname_out = "egg_DES_dim2_PT"; 
	mysampler.set_flname_out(flname_out);


	
	auto t1=std::chrono::high_resolution_clock::now();
	//mysampler.sample(Nsample); 
	mysampler.sample(Nsample,8);  // parallel version 
	auto t2=std::chrono::high_resolution_clock::now();
        auto Dt = t2-t1;

        std::cout <<    std::chrono::duration_cast<std::chrono::milliseconds>(Dt).count()<<  endl;






















}
