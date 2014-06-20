/**
Copyright (C) 2014, Tohid Alizadeh, Milad Malekzadeh, Davide De Tommaso

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

/*! \file parameters.h
\brief Parameters class

The class Parameters allows the user to define and use the task parameters encoded as the projection matrix A (frame of reference) and the position vector b.
During demonstration the user should record the data together with these task parameters (A and b) for each time-step. These parameters may be constant during
each demonstration.

b: Position vector
A: Projection matrix
Frame: Frame of reference (task parameters)
nPARAMS: Number of task parameters - according to the task this should be determined before
nVARS: Number of variables
nPOINTS: Number of data-points
params: The set of task parameters


\author Tohid Alizadeh, Milad Malekzadeh, Davide De Tommaso
\bug No known bugs.
*/


#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <fstream>
#include "armadillo"
#include "pbdlib/datapoints.h"


using namespace arma;

namespace pbdlib
{

struct Frame
{
    vec b;
    mat A;
};

class Parameters
{
    private:
        std::vector< std::vector<Frame> > params;
        uint nPARAMS;
        uint nVARS;
        uint nPOINTS;

    public:
        Parameters(){}
        ~Parameters(){}
        Parameters(uint nVARS, uint nPARAMS);
        Parameters(uint nVARS, uint nPARAMS, uint nPOINTS);

        /*!
          getParams() Returns the task parameters for a given time-index (point_index)

        */
        std::vector<Frame>&                 getParams(uint point_index);
        /*!
          loadFromFile() Loads the paramerts from the given path
          The file is a .txt file containing b in the first row and the corresponding A in the next rows and the same for the other frames:
          [b_1;A_1;b_2;A_2;...]
          The .txt file should contain a matrix with the following number of rows and coloumns:
          n_cols = nVARS
          n_rows = nPARAMS * (nVARS + 1)

        */
        void                                loadFromFile(std::string path);

};

} // end of pbdlib namespace

#endif // PARAMETERS_H
