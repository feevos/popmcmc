#ifndef _composite_stepper_
#define _composite_stepper_



namespace mcmc {


template <class spec_stepper> 
struct combo_unit{
	spec_stepper unit_stepper; 	/* Specific stepper  */
	double weight; 			/* Weight factor, 0.  <= weight <= 1, else it is re-normalized  */

};






/*
enum stepper_type{
	DE,
	GW_Stretch,
	GW_Walk_mcmc,
	PCX,
	HMCMC	 // To be implemented. 
};
*/


struct weighted_stepper{
	stepper_type Type;	// Stepper type 	
	double weight; 		// Relative probability of realization 
};


template < class spec_stepper,   class T> 
class combo_stepper:public stepper<T> {
	private: 

		/** Vector of weighted steppers to be used */ 
		vector<combo_unit< spec_stepper > > vec_of_steppers; 



	public: 

	//explicit combo_stepper(); 
	explicit combo_stepper( vector<combo_unit< spec_stepper > >   &_vec_of_steppers );


	// This results, according to some probability, a unique proposed value. 
	T propose (vector<T> &X , int &j); 



};


template < class spec_stepper,   class T> 
combo_stepper< spec_stepper,T>::combo_stepper( vector<combo_unit< spec_stepper > >   &_vec_of_steppers ): vec_of_steppers(& _vec_of_steppers)
	{












}










}
#endif
