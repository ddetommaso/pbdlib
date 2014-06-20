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

/*! \file repro_gmr.cpp
\brief Learning GMM model and testing reproductions with GMR
Learning a GMM model from a demonstration saved in the file data_txyz.txt and after

\author Davide De Tommaso and Milad Malekzadeh
\bug No known bugs.
*/

#include "pbdlib/gmm.h"
#include "pbdlib/gmr.h"
#include <sstream>

using namespace pbdlib;

int main(int argc, char **argv)
{
    std::vector<Demonstration> demos;
    Demonstration demo =  Demonstration(4,200,0);
    std::vector<std::string> vars,vars2;
    vars.push_back("t");
    vars.push_back("x");
    vars.push_back("y");
    vars.push_back("z");

    demo.getDatapoints().loadFromFile("../../data/data_txyz.txt");

    demo.getDatapoints().setVarNames(vars);


    cout << endl << "Data is: " << demo.getDatapoints().getData() << endl;
    cout << endl << "nbVar is: " << demo.getDatapoints().getNumVARS() << endl;
    cout << endl << "nbDatapoints is: " << demo.getDatapoints().getNumPOINTS() << endl;

    getchar();
    demos.push_back(demo);

    GMM_Model *gmm;
    gmm = new GMM_Model(demos, 3);
    gmm->setVARSNames(vars);


    cout << "\n EM rounds = " << gmm->EM_learn();

    cout << "\n MU = " << endl << gmm->getCOMPONENTS().at(0).getMU();
    cout << "\n MU = " << endl << gmm->getCOMPONENTS().at(1).getMU();
    cout << "\n MU = " << endl << gmm->getCOMPONENTS().at(2).getMU();

    cout << "\n SIGMA = " << endl << gmm->getCOMPONENTS().at(0).getSIGMA();
    cout << "\n SIGMA = " << endl << gmm->getCOMPONENTS().at(1).getSIGMA();
    cout << "\n SIGMA = " << endl << gmm->getCOMPONENTS().at(2).getSIGMA();



    Datapoints *Repros;
        Repros = new Datapoints(1,200);


        vars2.push_back("t");

        Repros->setVarNames(vars2);


        mat data = demo.getDatapoints().getData().submat(0,0,0,199);
        Repros->setData(data);
        GMR* reg = new GMR(gmm);
        reg->regression(Repros);


        cout << "\n SIGMA = "  << reg->getCOMPONENTS().at(199).getSIGMA();
        cout << "\n MU = "  << reg->getCOMPONENTS().at(199).getMU();

  return 0;
}
