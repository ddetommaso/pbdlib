/**
Copyright (C) 2014, Davide De Tommaso

This file is part of PbDLib.

    PbDLib is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PbDLib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with PbDLib.  If not, see <http://www.gnu.org/licenses/>.
*/

/*! \file demonstration.h
\brief Demonstration class
The class Demonstration model a single demonstration in PbD framework

\author Davide De Tommaso
\bug No known bugs.
*/


#ifndef DEMONSTRATION_H
#define DEMONSTRATION_H

#include "armadillo"
#include <fstream>

#include "pbdlib/datapoints.h"
#include "pbdlib/parameters.h"

using namespace arma;

namespace pbdlib
{

class Demonstration
{
    private:
        uint nVARS, nPOINTS, nPARAMS;        
        Parameters params;
        std::vector<std::string> params_names;
        Datapoints data;

    public:
        Demonstration(std::string path);
        Demonstration(uint _nVARS, uint _nPOINTS, uint _nPARAMS);

        Datapoints&                 getDatapoints();
        Parameters&                 getParameters();
        std::vector<Parameters>&    getParams();
        std::vector<std::string>&   getParamNames();        

        void                        setParamNames(const std::vector<std::string> &paramsnames);        
        void                        saveInFile(std::string path);
        void                        setParams(const Parameters &parameters);

};

} //end of pbdlib namespace

#endif
