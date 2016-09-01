
#ifndef _vec_ops_
#define _vec_ops_

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>


namespace mcmc{
namespace tools{



class vec_ops{

	private: 

		std::vector<std::vector<double> > XGS; 

	public: 
		
		/* Necessary functions for Gram-Schmidt orthonormalization */      // <---- Add in external file.  
		double inner_prod(std::vector<double> &u, std::vector<double> &v);
		std::vector<double> cross_prod (std::vector<double> &u, std::vector<double> &v); 
		double norm  (std::vector<double> &u);
		double norm_sq (std::vector<double> &u);
		void normalize (std::vector<double> &u);

		std::vector<double> projection (std::vector<double> &u, std::vector<double> &v);

		/* Constructs orthonormal basis Gram-Schmidt process */
		void Ortho_GS_base (std::vector<std::vector<double> > &X); 
		std::vector<std::vector<double> > get_GS() {return XGS;}; 

		

};




double vec_ops::inner_prod(std::vector<double> &u, std::vector<double> &v) {
			
		auto N = static_cast<int> (u.size());
		if (v.size() !=  u.size()){
			std::cerr << "Incompatible vectors dimensionality for projection, aborting ..." << std::endl; 
			throw(0); 
			}
	
		double result = 0.0; 
		for (auto i=0; i < N; ++i)
			result += u[i] * v[i];

		return 
			result; 

}




std::vector<double> vec_ops::cross_prod(std::vector<double> &u, std::vector<double> &v) {
			
		auto N = u.size();
		if (v.size() != N){
			std::cerr << "Incompatible vectors dimensionality for projection, aborting ..." << std::endl; 
			throw(0); 
			}
	

		std::vector<double> temp { 
			u[1]*v[2] - u[2]*v[1] ,  
			u[2]*v[0] - u[0]*v[2] ,  
			u[0]*v[1] - u[1]*v[0] 
		};

		return std::move(temp); 


}






double vec_ops::norm  (std::vector<double> &u)
	{
		double norm=0.0;
		for (auto &it: u)
			norm += it*it;

		norm = sqrt(norm);

		return norm;
}


// Squared norm. 
double vec_ops::norm_sq (std::vector<double> &u)
	{
	double norm=0.0;
	for (auto &it: u)
		norm += it*it;

	return norm;
}


// Normalizes a vector 
void vec_ops::normalize (std::vector<double> &u){

		double norm_val = norm(u);
		if (norm_val != 0.0){
			for (auto &it: u)
				it /= norm_val; 	

		}

}



// Euclidean projection operator: projection of V (second) along u (first)
std::vector<double> vec_ops::projection (std::vector<double> &u, std::vector<double> &v)
		{
		auto N = static_cast<int> (u.size());

		double norm_u = norm_sq(u); 
		std::vector<double> temp(N,0.0);
		if (norm_u !=0.0){
		double factor =  inner_prod(u,v) /  norm_u ; 
			for (auto j=0; j <N; ++j)
				temp[j] =  factor * u[j]; 
		}
		// Now this vector is the projected along u vector 
		return temp; 
}




// Constructs orthonormal basis Gram-Schmidt process:  
void vec_ops::Ortho_GS_base (std::vector<std::vector<double> > &X)  {

		auto Nvectors = static_cast<int> ( X.size() );
		auto Ndim = static_cast<int> (X[0].size()); 

		// Numerically stable basis. 
		for (auto i=0; i < Nvectors; ++i)
			{
			normalize(X[i]);	
			for (auto j=i+1; j <Nvectors; ++j)
				{
				std::vector<double> proj = projection(X[i],X[j]);
				for (auto kk =0; kk < Ndim; ++kk)
					X[j][kk]  = X[j][kk] - proj[kk];
				}
			}




}





}// End namespace tools
} // End namespace EA 

#endif
