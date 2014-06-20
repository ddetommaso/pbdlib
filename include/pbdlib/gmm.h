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

/*! \file gmr.h
\brief GMM_model class
    The class GMM_model allows to use a Gaussian Mixture Model and to learn the parameters from task demonstrations

\author Davide De Tommaso, Milad Malekzadeh
\bug No known bugs.
*/

#ifndef GMM_H
#define GMM_H

#include "pbdlib/datapoints.h"
#include "pbdlib/demonstration.h"
#include "pbdlib/mvn.h"
#include "armadillo"

using namespace arma;


namespace pbdlib
{

class GMM_Model
{
    private:
        uint nVARS, nSTATES;
        std::vector<GaussianDistribution> COMPONENTS;
        rowvec PRIORS;
        std::vector<Demonstration> DEMONSTRATIONS;
        std::vector<std::string> vars_names;

        void    learnKMEANS();      // [MacQueen, 1967]
        bool    EM_isfinished(double l_old, double l_new); // [Calinon, 2009] p.48
        uint    EM(double likelihood); // [Calinon, 2009] p.48



    public:
        GMM_Model(std::vector<Demonstration> &demos, uint _nSTATES);
        GMM_Model(uint _nSTATES, uint _nVARS);
        GMM_Model(const std::string &priors_path, const std::string &mu_path, const std::string &sigma_path,const std::string &vars_path);
        ~GMM_Model(){}

        uint                                EM_learn();
        std::vector<std::string>&           getVARSNames();
        uint                                getIndexOfVARName(const std::string& varname);
        uint                                getNumVARS();
        uint                                getNumSTATES();
        rowvec&                             getPRIORS();
        std::vector<GaussianDistribution>&  getCOMPONENTS();
        double                              getProbability(const colvec& sample); // [Calinon, 2009] p.34 (2.1)
        double                              getLikelihood(const mat& SAMPLES);  // [Calinon, 2009] p.35 (2.2)        
        bool                                addDemo(Demonstration &demo);


        void                                setPRIORS(const rowvec& priors);
        void                                setCOMPONENTS(const std::vector<GaussianDistribution>& components);
        void                                setVARSNames(const std::vector<std::string>& vars);
        void                                saveInFiles();
};

} //end of pbdlib namespace

#endif



/*
 |
 |  >> REFERENCES <<
 |
 |    [Calinon, 2009]
 |            @book{Calinon09book,
 |              author="S. Calinon",
 |              year="2009",
 |              title="Robot Programming by Demonstration: A Probabilistic Approach",
 |              publisher="EPFL/CRC Press"
 |            }
 |
 |
 |    [MacQueen, 1967]
 |           @inproceedings{macqueen1967some,
 |              title={Some methods for classification and analysis of multivariate observations},
 |              author={MacQueen, James and others},
 |              booktitle={Proceedings of the fifth Berkeley symposium on mathematical statistics and probability},
 |              volume={1},
 |              number={281-297},
 |              pages={14},
 |              year={1967},
 |              organization={California, USA}
 |           }
 |
 |
*/
