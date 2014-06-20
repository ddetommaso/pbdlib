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

/*! \file datapoints.h
\brief Datapoints class
The class Datapoints allows to access samples of n-points and n-variables
stored in a Armadillo matrix.

\author Davide De Tommaso
\bug No known bugs.
*/


#ifndef DATAPOINTS_H
#define DATAPOINTS_H

#include "armadillo"
#include <fstream>

using namespace arma;

namespace pbdlib
{

class Datapoints
{
    private:
        uint nVARS, nPOINTS;
        mat data; // [nVARS]x[nPOINTS]
        mat dataAll; // Reserved for PGMM [nVARSxnFRAMES]x[nPOINTS]
        std::vector<std::string> vars_names;

    public:
        Datapoints(){}
        Datapoints(uint _nVARS, uint _nPOINTS);


        /*!
          getData() provides in output an Armadillo mat of size [nVARS]x[nPOINTS]

        */
        mat&                        getData();
        uint                        getNumVARS();
        uint                        getNumPOINTS();
        std::string&                getVarName(uint varIndex);
        uint                        getIndexVarname(std::string varname);
        std::vector<std::string>&   getVarNames();
        mat&                        getDataAll();

        void                        setData(mat &_data);
        void                        setDataAll(mat &dataAll);
        void                        setVarNames(std::vector<std::string> &names);                
        void                        saveInFile(std::string path);
        void                        loadFromFile(std::string path, bool is_transpose = false);
};

} //end of pbdlib namespace

#endif
