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

/*! \file test_datapoint.cpp
\brief testing datapoint class
Testing basic features of the Datapoint class

\author Davide De Tommaso and Milad Malekzadeh
\bug No known bugs.
*/

#include "pbdlib/datapoints.h"
#include <sstream>


using namespace pbdlib;
using namespace arma;

int main(int argc, char **argv)
{
    mat A = randu<mat>(3,100);
    Datapoints datapoint(3,100);
    datapoint.setData(A);

    cout << "\n Data \n " << datapoint.getData();
    datapoint.saveInFile("data01.txt");

    datapoint.loadFromFile("data01.txt");
    cout << "\n Data \n " << datapoint.getData();

    return 0;
}
