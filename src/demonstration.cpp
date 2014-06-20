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

#include "pbdlib/demonstration.h"

namespace pbdlib
{

Demonstration::Demonstration(uint _nVARS, uint _nPOINTS, uint _nPARAMS)
{
    nVARS = _nVARS;
    nPOINTS = _nPOINTS;
    nPARAMS = _nPARAMS;

    data = Datapoints(nVARS,nPOINTS);
}

Datapoints& Demonstration::getDatapoints()
{
    return data;
}


void Demonstration::setParamNames(const std::vector<std::string>& paramsnames)
{
    if(paramsnames.size() == nPARAMS)
        params_names = paramsnames;
    else
        std::cout << "\n [ERROR]::Demonstration::setParamNames if(paramsnames.size() == nPARAMS) ... else .";
}

std::vector<std::string>& Demonstration::getParamNames()
{
    return params_names;
}


Parameters& Demonstration::getParameters()
{
    return params;
}

void Demonstration::setParams(const Parameters &parameters)
{
    params = parameters;
}


void Demonstration::saveInFile(std::string path)
{
    data.getData().save(path, raw_ascii);
}

} // end of pbdlib namespace
