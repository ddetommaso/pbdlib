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

#include "pbdlib/datapoints.h"

namespace pbdlib
{

Datapoints::Datapoints(uint _nVARS, uint _nPOINTS)
{
    nVARS = _nVARS;
    nPOINTS = _nPOINTS;
    data = zeros(_nVARS,_nPOINTS);
}


void Datapoints::setData(mat& _data)
{
    data = _data;
}

mat& Datapoints::getData()
{
    return data;
}

void Datapoints::setDataAll(mat& dataAll)
{
    this->dataAll = dataAll;
}

mat& Datapoints::getDataAll() //Reserved for the PGMM version. dataAll is the data projected in all the frames of references
{
    return this->dataAll;
}



uint  Datapoints::getNumVARS()
{
    return nVARS;
}

uint  Datapoints::getNumPOINTS()
{
    return nPOINTS;
}

std::string& Datapoints::getVarName(uint varIndex)
{
    if( varIndex < nVARS )
        return vars_names.at(varIndex);
    else
        std::cout << "\n [ERROR]::Datapoints::getVarName  if( varIndex < nVARS ) ... else .";

    return vars_names.at(0);
}

void Datapoints::setVarNames(std::vector<std::string>& names)
{
    if(names.size() == nVARS)
        vars_names = names;
    else
        std::cout << "\n [ERROR]::Datapoints::setVarNames  if(names.size() == nVARS) ... else .";
}


std::vector<std::string>& Datapoints::getVarNames()
{
    return vars_names;
}

uint Datapoints::getIndexVarname(std::string varname)
{
    std::vector<std::string>::iterator it;
    it = std::find (vars_names.begin(), vars_names.end(), varname);
    return std::distance(vars_names.begin(), it);
}


void Datapoints::loadFromFile(std::string path, bool is_transpose)
{
    mat _data;
    _data.load(path, raw_ascii);
    if(is_transpose)
        _data = _data.t();
    setData(_data);
}


void Datapoints::saveInFile(std::string path)
{
    data.save(path, raw_ascii);
}


} // end of pbdlib namespace
