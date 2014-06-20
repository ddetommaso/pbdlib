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

#include "pbdlib/gmr.h"


namespace pbdlib
{

GMR::GMR(GMM_Model* model)
{
    task_model = NULL;
    setGMMModel(model);
}

void GMR::regression(Datapoints* data_in)
{ 
    uint i,j,k;

    COMPONENTS.clear();

    if(task_model == NULL)
        return;

    if(data_in->getVarNames().size() == 0)
    {
        std::cout << "\n [ERROR]::GMR::regression() if(data_in->getVarNames().size() == 0) ... else .";
        return;
    }

    urowvec in,out;
    bool found = false;

    in.zeros(1, data_in->getNumVARS());
    out.zeros(1, task_model->getNumVARS() - data_in->getNumVARS());

    for(j=0; j<data_in->getNumVARS(); j++)
        for(i=0; i<task_model->getNumVARS(); i++)
            if ( data_in->getVarName(j).compare( task_model->getVARSNames().at(i) ) == 0)
                in(j) = i;


    for(i=0, k=0; i<task_model->getNumVARS(); i++)
    {
        for(j=0; j<data_in->getNumVARS(); j++)
            if ( data_in->getVarName(j).compare( task_model->getVARSNames().at(i) ) == 0)
                found = true;

        if(!found)
            out(k++) = i;
        found = false;

    }

    //std::cout << "\n IN = " << in;
    //std::cout << "\n OUT = " << out;

    mat x = data_in->getData();
    // P+Q = D ?
    if( (in.n_elem + out.n_elem) == task_model->getNumVARS())
    {
        if(x.n_rows == in.n_elem)
        {
            uint i,j;
            ucolvec icol;
            colvec Mu_tmp;
            mat Sigma_tmp, beta, yj_tmp, Sigmaj_y_tmp;
            mat Pxi = zeros(data_in->getNumPOINTS(),task_model->getNumSTATES());
            y = zeros(out.n_elem, data_in->getNumPOINTS());
            Sigma_y = zeros(out.n_elem, out.n_elem, data_in->getNumPOINTS());

            // Compute the influence of each GMM component, given input x


            GaussianDistribution* Gtmp;
            for(i=0; i<task_model->getNumSTATES(); i++)
            {
                icol=i;
                Mu_tmp = Mu.submat(in,icol);
                Sigma_tmp = Sigma.slice(i).submat(in,in);
                //Mu_tmp.print("Mu(i) = ");			// DEBUGGING
                //Sigma_tmp.print("Sigma(i) = ");	// DEBUGGING
                Gtmp = new GaussianDistribution(Mu_tmp, Sigma_tmp);
                Pxi.col(i) = Priors(i) * (Gtmp->getPDFValue(data_in->getData()));

            }
            //Priors.print("Priors = ");	// DEBUGGING
            //data_in->getData().print("t = ");	// DEBUGGING
            //Pxi.print("Pxi = ");	// DEBUGGING
            beta = (Pxi / repmat(sum(Pxi,1),1,task_model->getNumSTATES())).t();
            //beta.print("beta = ");	// DEBUGGING

            colvec mu;
            for(i=0; i<data_in->getNumPOINTS(); i++)
            {
                for(j=0; j<task_model->getNumSTATES(); j++)
                {
                    icol=j;
                    yj_tmp = Mu(out,icol) + Sigma.slice(j).submat(out,in)*inv(Sigma.slice(j).submat(in,in)) * (x.col(i)-Mu(in,icol));
                    //yj_tmp.print("y_tmp = "); 	// DEBUGGING
                    y.col(i) = y.col(i) + (beta(j,i) * yj_tmp);

                    Sigmaj_y_tmp = Sigma.slice(j).submat(out,out) - (Sigma.slice(j).submat(out,in)*inv(Sigma.slice(j).submat(in,in))*Sigma.slice(j).submat(in,out));
                    Sigma_y.slice(i) = Sigma_y.slice(i) + ((beta(j,i)*beta(j,i)) * Sigmaj_y_tmp);
                }

                mu = y.col(i);
                COMPONENTS.push_back( GaussianDistribution(mu, Sigma_y.slice(i)) );
            }
        }
    }


}

cube&  GMR::get_Sigma_y()
{
        return Sigma_y;
}

mat&  GMR::get_y()
{
        return y;
}

void GMR::setGMMModel(GMM_Model* gmmmodel)
{
    task_model = gmmmodel;

    Priors = task_model->getPRIORS();

    Mu = zeros( task_model->getNumVARS(), task_model->getNumSTATES() );
    Sigma = zeros(task_model->getNumVARS(), task_model->getNumVARS(), task_model->getNumSTATES());

    for(int i=0; i<task_model->getNumSTATES(); i++)
    {
        Mu.col(i) = task_model->getCOMPONENTS().at(i).getMU();
        Sigma.slice(i) = task_model->getCOMPONENTS().at(i).getSIGMA();
    }

}

std::vector<GaussianDistribution>&  GMR::getCOMPONENTS()
{
    if(COMPONENTS.size() != 0)
        return COMPONENTS;
}


void GMR::saveGMRInFiles()
{
    y.save("GMR_Y.txt", raw_ascii);
    Sigma_y.save("GMR_SigmaY.txt", raw_ascii);
}

} //end of pbdlib namespace
