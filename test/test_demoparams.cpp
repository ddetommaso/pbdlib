/**
Copyright (C) 2014, Tohid Alizadeh, Davide De Tommaso

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
\brief testing Demonstration and Parameters classes
Testing basic features of the Demonstration and Parameters classes

\author Milad Malekzadeh and Davide De Tommaso
\bug No known bugs.
*/

#include "pbdlib/demonstration.h"
#include "pbdlib/pgmm.h"
#include "pbdlib/gmr.h"
#include <sstream>

using namespace pbdlib;
using namespace arma;

int main(int argc, char **argv)
{

    uint nPARAMS = 2;
    uint nStates = 3;
    uint nVars = 3;
    uint nDemos = 4;

    std::vector<std::string> varNames;
    varNames.push_back("t");    // time is input (t) and the attractor position (2D vector [x y]) is the output.
    varNames.push_back("x");
    varNames.push_back("y");

    //......... loading the demonstrations and task parameters from the txt files ......
    std::vector<Demonstration> demos;
    Demonstration demo =  Demonstration(nStates,200,nPARAMS);  //(3 vars, 200 datapoints, 2 frames of references)


    char filename[256];

    cout << "Loading the demonstrations and the task parameters ...";

    for (int m=0; m<nDemos; m++){   //Loading demos in the loop
          Parameters param = Parameters(3,2,1);  //nVAR, nPAR, nPOINTS (set to 1 to show that the task parameters are fixed all the time)
          sprintf(filename, "../../data/pgmm/Data0%d.txt",m+1);

          demo.getDatapoints().loadFromFile(filename, true);
          demo.getDatapoints().setVarNames(varNames);

          cout << endl << "Data is: " << demo.getDatapoints().getData() << endl;
          sprintf(filename, "../../data/pgmm/Param0%d.txt",m+1);
          param.loadFromFile(filename);
          demo.setParams(param);
          cout << endl << "Parameter is: " << param.getParams(0).at(0).A << endl;
          demos.push_back(demo);
      }

      cout<<endl<<"Demonstrations and task parameters are loaded successfully !!!!"<<endl;

}
