
#ifndef _vec_ops_
#define _vec_ops_




class vec_ops{

	private: 

		vector<vector<double> > XGS; 

	public: 
		
		/* Necessary functions for Gram-Schmidt orthonormalization */      // <---- Add in external file.  
		double inner_prod(vector<double> &u, vector<double> &v);
		double norm  (vector<double> &u);
		double norm_sq (vector<double> &u);
		void normalize (vector<double> &u);

		vector<double> projection (vector<double> &u, vector<double> &v);

		/* Constructs orthonormal basis Gram-Schmidt process */
		void Ortho_GS_base (vector<vector<double> > &X); 
		vector<vector<double> > get_GS() {return XGS;}; 

		

};




double vec_ops<T>::inner_prod(vector<double> &u, vector<double> &v) {
			
		auto N = u.size();
		if (v.size() != N){
			cerr << "Incompatible vectors dimensionality for projection, aborting ..." << endl; 
			throw(0); 
			}
	
		double result = 0.0; 
		for (auto i=0; i < N; ++i)
			result += u[i] * v[i];

		return 
			result; 

}


double vec_ops<T>::norm  (vector<double> &u)
	{
		double norm=0.0;
		for (auto &it: u)
			norm += it*it;

		norm = sqrt(norm);

		return norm;
}


// Squared norm. 
double vec_ops::norm_sq (vector<double> &u)
	{
	double norm=0.0;
	for (auto &it: u)
		norm += it*it;

	return norm;
}


// Normalizes a vector 
void vec_ops::normalize (vector<double> &u){
		vector<double> temp(u.size(),0.0);

		double norm_val = norm(u);
		if (norm_val != 0.0){
			for (auto &it: u)
				it /= norm_val; 	

		}

}



// Euclidean projection operator. 
template <class T>
vector<double> vec_ops<T>::projection (vector<double> &u, vector<double> &v)
		{
		auto N = u.size();

		double norm_u = norm_sq(u); 
		vector<double> temp(N,0.0);
		if (norm_u !=0.0){
		double factor =  inner_prod(u,v) /  norm_u ; 
			for (auto j=0; j <N; ++j)
				temp[j] =  factor * u[j]; 
		}
		// Now this vector is the projected along u vector 
		return temp; 
}




// Constructs orthonormal basis Gram-Schmidt process:  
template <class T>
void vec_ops<T>::Ortho_GS_base (vector<vector<double> > &X)  {

		auto Nvectors = X.size();
		auto Ndim = X[0].size(); 

		// Numerically stable basis. 
		for (auto i=0; i < Nvectors; ++i)
			{
			normalize(X[i]);	
			for (auto j=i+1; j <Nvectors; ++j)
				{
				vector<double> proj = projection(X[i],X[j]);
				for (auto kk =0; kk < Ndim; ++kk)
					X[j][kk]  = X[j][kk] - proj[kk];
				}
			}


		XGS = std::move(X);


}






#endif
