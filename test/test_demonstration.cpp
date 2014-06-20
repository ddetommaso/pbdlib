/**
Copyright (C) 2014, Milad Malekzadeh, Davide De Tommaso

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

/*! \file test_demonstration.cpp
\brief testing Demonstration class
Testing basic features of the Demonstration class

\author Milad Malekzadeh and Davide De Tommaso
\bug No known bugs.
*/

#include "pbdlib/demonstration.h"
#include <sstream>

using namespace pbdlib;
using namespace arma;

int main(int argc, char **argv)
{
    mat A = randu<mat>(3,100);
    Datapoints datapoint(3,100);
    datapoint.setData(A);

    cout << "\n Data \n " << datapoint.getData();

    datapoint.saveInFile("data02.txt");

    Demonstration demo =  Demonstration(4,200,0);
    std::vector<Demonstration> demos;
    std::vector<std::string> vars;

    vars.push_back("t");
    vars.push_back("x");
    vars.push_back("y");
    vars.push_back("z");

    demo.getDatapoints().setVarNames(vars);

    demo.getDatapoints().loadFromFile("data02.txt");

    cout << endl << "Data is: " << demo.getDatapoints().getData() << endl;
    cout << endl << "nbVar is: " << demo.getDatapoints().getNumVARS() << endl;
    cout << endl << "nbDatapoints is: " << demo.getDatapoints().getNumPOINTS() << endl;

    demos.push_back(demo);

    demo.saveInFile("demo01.txt");

    return 0;
}
