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


#include "pbdlib/parameters.h"

namespace pbdlib
{

Parameters::Parameters(uint _nVARS, uint _nPARAMS, uint _nPOINTS)
{
    nVARS = _nVARS;
    nPARAMS = _nPARAMS;
    nPOINTS = _nPOINTS;
    params.reserve(_nPOINTS);
}

std::vector<Frame>& Parameters::getParams(uint point_index)
{
    return params.at(point_index);
}

void Parameters::loadFromFile(std::string path)
{
    mat paramTmp;
    if( !paramTmp.load(path, raw_ascii) )
    {
        std::cout << "\n [ERROR]::Parameters::readParamsFromTxtFile(std::string path) if( paramTmp.load(path, raw_ascii) )   ... else .";
        return;
    }

    uint i,j;
    std::vector<Frame> params_tmp;
    Frame frame_tmp;


    for(i=0; i<nPOINTS; i++)
    {
        for(j=0; j<nPARAMS; j++)
        {
            frame_tmp.b = trans(paramTmp.row( j*(nVARS+1)+i*nPARAMS*(nVARS+1) ));

            frame_tmp.A = paramTmp.rows( j*(nVARS+1)+1+i*nPARAMS*(nVARS+1) , j*(nVARS+1)+nVARS+i*nPARAMS*(nVARS+1) );
            params_tmp.push_back(frame_tmp);
        }
        params.push_back(params_tmp);
        params_tmp.clear();
    }

}

} //end of pbdlib namespace

