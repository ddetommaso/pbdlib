/*
Copyright (C) 2014, Tohid Alizadeh, Milad Malekzadeh, Leonel Rozo, Davide De Tommaso

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

/*! \file pgmm.h
\brief pgmm class
    This class is the implementation of the task parameterized Gaussian mixture models (TP-GMM) ([Calinon, 2012], [Calinon, ICRA14]).

\author Tohid Alizadeh, Milad Malekzadeh, Leonel Rozo, Davide De Tommaso
\bug No known bugs.
*/
#ifndef PGMM_H
#define PGMM_H

#include "armadillo"
#include "pbdlib/gmm.h"

using namespace arma;

namespace pbdlib
{


#define REALMIN 2.2251e-200
#define REALMAX 1.7977e200

class PGMM_Model
{
    /*!
      nVARS: Number of variables
      nSATAES: Number of Gaussian components of the model (This sould be set by the user)
      nPARAMS: Number of task parameters (This sould be set by the user)
      GMMS: A vector of GMM models containing for each task parameter, the center and covariance of the PGMM
      PRIORS: A vector containing the mixing coefficients of the Gaussian components

    */
    uint nVARS, nSTATES, nPARAMS;
    std::vector<GMM_Model> GMMS;              //Added to have seperate GMMS for each frame (This might be better than previous solution)
    rowvec PRIORS;
    std::vector<std::string> vars_names;
    std::vector<std::string> frames_names;
    /*!
      DEMONSTRATIONS: A vector of demonstrations, each demonstration contains the trajectory and the task parameters

    */
    std::vector<Demonstration> DEMONSTRATIONS;//All the demonstrations will be loaded in this variable


    /*!
      initTensorGMMTimeBased(): It is used to initialize the centers and the covariances of GMMS using the equally spaced time domain

    */
    void    initTensorGMMTimeBased();
    /*!
      EM_tensorGMM(): The implementation of the Expectation-Maximization algorithm for pgmm

    */
    void    EM_tensorGMM();

    public:
        /*!
          EM_TensorLearn(): Calls initTensorGMMTimeBased() &
                            EM_tensorGMM()
        */

        uint EM_TensorLearn();  //EM algorithm with the last formulation (Tensor)

        PGMM_Model(std::vector<Demonstration> &demos, uint nSTATES, uint nPARAMS);
        PGMM_Model(uint nVARS, uint nSTATES, uint nPARAMS);
        ~PGMM_Model(){}

        /*!
          getNumVARS(): Gives nVARS
        */
        uint                            getNumVARS();
        /*!
          getNumSTATES(): Gives nSTATES
        */
        uint                            getNumSTATES();
        /*!
          getNumFRAMES(): Gives nPARAMS
        */
        uint                            getNumFRAMES();
        /*!
          getPRIORS(): Gives the PRIORS
        */
        rowvec&                         getPRIORS();
        void                            setVARSNames(const std::vector<std::string>& vars);
        std::vector<std::string>&       getVARSNames();
        void                            setFRAMESNames(const std::vector<std::string>& frames);
        std::vector<std::string>&       getFRAMESNames();
        /*!
          loadPGMMfromMATLAB(): From a set of given .txt files reads the PRIORS, MU and SIGMA of the GMMS
          in the MuFileName.txt there is a matrix containing centers like this: [MU_1 MU_2 MU_3 ...]
          in the SigmaFileName.txt there is matrix containing the covariance matrices like this: [SIGMA_1_1 SIGMA_1_2 ... SIGMA_i_j] ; i=1:nPARAMS and j=1:nSTATES
        */
        void                            loadPGMMfromMATLAB(std::string priorsFileName, std::string MuFileName, std::string SigmaFileName);
        /*!
          prodGMM(): By providing task parameters (A & b) using MU and SIGMA of the GMMS gives us product of linearly transformed Gaussians
                    This results are used for reproduction using GMR.
                    [Calinon, 2012] Eq. (3) - Product of linearly transformed Gaussians
        */
        GMM_Model*                      prodGMM(std::vector<Parameters> params);		// [Calinon, 2012] Eq. (3) - Prod. of linearly transformed Gaussians
        /*!
          getGMMS(): Gives GMMS
        */
        std::vector<GMM_Model>&         getGMMS();

        GMM_Model*                      getGMM(std::vector<mat> A, std::vector<colvec> b);
};

} // end of namespace pbdlib
#endif
