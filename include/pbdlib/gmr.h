/**
Copyright (C) 2014, Tohid Alizadeh, Milad Malekzadeh, Leonel Rozo, Davide De Tommaso

This file is part of PbDLib (Programming-by-demonstration C++ Library).

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

/*! \file gmr.h
\brief gmr class
    This class is the implementation of the Gaussian mixture regression (GMR).

\author Tohid Alizadeh, Milad Malekzadeh, Leonel Rozo, Davide De Tommaso
\bug No known bugs.
*/

#ifndef GMR_H
#define GMR_H

#include "pbdlib/gmm.h"
#include "armadillo"


using namespace arma;

namespace pbdlib
{

class GMR
{
    private:
        GMM_Model *task_model;
        Datapoints *data_in;
        /*!
          COMPONENTS: A vector of Gaussian components, each component has PRIORS, MU and SIGMA
        */
        std::vector<GaussianDistribution> COMPONENTS;
        /*!
          Priors: Vector of priors (mixing coefficients)
        */
        rowvec Priors;
        /*!
          Mu: Centers of Gaussian components
        */
        mat Mu;
        /*!
          Sigma: Covariance matrices of Gaussian components in a cube format
        */
        cube Sigma;
        /*!
          Sigma_y: Covariance matrices of the predicted output of GMR
        */
        cube Sigma_y;
        /*!
          y: Predicted output of GMR (predicted centers of output Gaussian)
        */
        mat y;

    public:
        GMR(){}
        // The GMM_Model should be already learnt. (SIGMA, MU, PRIORS)
        GMR(GMM_Model* model);

        /*!
          get_Sigma_y(): Gives Sigma_y
        */
        cube&                               get_Sigma_y();
        /*!
          get_y(): Gives y
        */
        mat&                                get_y();
        /*!
          getCOMPONENTS(): Gives the vector of COMPONENTS
        */
        std::vector<GaussianDistribution>&  getCOMPONENTS();

        void                                setGMMModel(GMM_Model* gmmmodel);
        /*!
          regression(): Implementation of GMR
        */
        void                                regression(Datapoints* data_in);
        /*!
          saveGMRInFiles(): Saves the output in two .txt files (GMR_Y.txt and GMR_SigmaY.txt)
        */
        void                                saveGMRInFiles();
};

} //end pbdlib namespace

#endif


