/**
Copyright (C) 2014, Davide De Tommaso, Milad Malekzadeh

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

#include "pbdlib/mvn.h"

namespace pbdlib
{


GaussianDistribution::GaussianDistribution(colvec& _MU, mat& _SIGMA)
{
    nVARS = _MU.n_rows;
    MU = _MU;
    SIGMA = _SIGMA;
}


GaussianDistribution::GaussianDistribution(uint _nVARS)
{
    nVARS = _nVARS;
}


uint GaussianDistribution::getNumVARS()
{
    return nVARS;
}


colvec &GaussianDistribution::getMU()
{
    return MU;
}


mat& GaussianDistribution::getSIGMA()
{
    return SIGMA;
}


void GaussianDistribution::setMU(const colvec &_MU)
{
    MU = _MU;
}


void GaussianDistribution::setSIGMA(const mat &_SIGMA)
{
    SIGMA = _SIGMA;
}


void GaussianDistribution::setNumVARS(uint numvars)
{
    nVARS = numvars;
}


colvec GaussianDistribution::getPDFValue(const mat &SAMPLES)
{
        colvec Probs(SAMPLES.n_cols);   
        mat D_Tmp = trans(SAMPLES) - repmat(trans(MU),SAMPLES.n_cols,1);        
        mat invTmp;

        invTmp = inv(SIGMA);
        Probs = sum((D_Tmp * invTmp) % D_Tmp, 1);
        Probs = exp(-0.5*arma::abs(Probs)) / sqrt(pow((2*PI),SIGMA.n_cols) * (fabs(det(SIGMA)) + THRESHOLD_MIN));

        return Probs;
}

} //end of pbdlib namespace
