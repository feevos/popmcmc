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

#ifndef _individual_
#define _individual_

#include <vector>
using namespace std;

namespace mcmc{

struct individual{
	vector<double> Vars; // Vector of parameters. 
	double loglkhood; // Value of loglikelihood
	double Z; // Value of random in case GW_Stretch is used. 
};

} // end of namespace

#endif
