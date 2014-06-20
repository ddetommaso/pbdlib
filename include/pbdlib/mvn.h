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

/*! \file mvn.h
\brief GaussianDistribution class
The class GaussianDistribution model a multivariate Gaussian distribution.

\author Davide De Tommaso, Milad Malekzadeh
\bug No known bugs.
*/


#ifndef MVN_H
#define MVN_H

#define THRESHOLD_MIN std::numeric_limits<float>::epsilon()
#define THRESHOLD_MAX 1.0e60

#define PI 3.14159265359

#include "armadillo"
#include "pbdlib/datapoints.h"

using namespace arma;

namespace pbdlib
{


class GaussianDistribution
{
    private:
        uint nVARS;
        mat SIGMA;
        colvec MU;

    public:
        GaussianDistribution(uint _nVARS);
        GaussianDistribution(colvec& _MU, mat& _SIGMA);

        uint                    getNumVARS();
        colvec&                 getMU();
        mat&                    getSIGMA();
        colvec                  getPDFValue(const mat& SAMPLES);        

        void                    setMU(const colvec& _MU);
        void                    setSIGMA(const mat& _SIGMA);
        void                    setNumVARS(uint numvars);

};

} //end of namespace pbdlib

#endif
