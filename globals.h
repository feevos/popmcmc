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


#ifndef _globals_
#define _globals_

#include <random>
#include <string> 
namespace mcmc{

// Implemented algorithm names 
const std::string  GW_Stretch = "GW_Stretch"; 
const std::string  DE = "DE"; 
const std::string  DES = "DES"; 
const std::string  PCX = "PCX"; 

const std::string  GW_Stretch_PT = "GW_Stretch_PT"; 
const std::string  DE_PT = "DE_PT"; 
const std::string  DES_PT = "DES_PT"; 
const std::string  PCX_PT = "PCX_PT"; 



std::random_device rd;
std::mt19937 gen;

std::uniform_real_distribution<> unif_real; /**< It is used inside steppers */




}
#endif 
