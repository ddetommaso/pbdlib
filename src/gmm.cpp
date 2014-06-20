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

#include "pbdlib/gmm.h"

namespace pbdlib
{

GMM_Model::GMM_Model(std::vector<Demonstration> &demos, uint _nSTATES)
{    
    DEMONSTRATIONS = demos;
    nVARS = demos.at(0).getDatapoints().getNumVARS();
    nSTATES = _nSTATES;
    PRIORS = rowvec(_nSTATES);
    COMPONENTS.reserve(_nSTATES);
    setVARSNames( demos.at(0).getDatapoints().getVarNames());
}

GMM_Model::GMM_Model(uint _nSTATES, uint _nVARS)
{
    COMPONENTS.reserve(_nSTATES);
    nVARS = _nVARS;
    nSTATES = _nSTATES;
    PRIORS = rowvec(_nSTATES);
}


GMM_Model::GMM_Model(const std::string &priors_path, const std::string &mu_path, const std::string &sigma_path,const std::string &vars_path)
{
        mat priors, mu, sigma;

        priors.load(priors_path, raw_ascii);
        mu.load(mu_path, raw_ascii);
        sigma.load(sigma_path, raw_ascii);

        nVARS = mu.n_rows;
        nSTATES = priors.n_elem;

        std::ifstream varsfile(vars_path.c_str());
        std::string varsnames;
        std::getline (varsfile,varsnames);

        std::vector<std::string> vars;

        std::stringstream ss(varsnames); // Insert the string into a stream
        std::string buf;

        while (ss >> buf)                      
              vars.push_back(buf);        

        setVARSNames(vars);

        mat _SIGMA = zeros(nVARS, nVARS);
        colvec _MU = zeros(nVARS,1);
        rowvec _PRIORS(1, nSTATES);
        std::vector<GaussianDistribution> components;      


        for(uint i=0; i<nSTATES; i++)
        {
            _PRIORS(0,i) = priors(0,i);

            for(uint j=0; j<nVARS; j++)
            {
                _MU(j,0) = mu(j,i);
                for(uint k=0; k<nVARS; k++)
                    _SIGMA(j,k) = sigma(j,k + i*(nSTATES+1));
            }
            components.push_back(GaussianDistribution(_MU, _SIGMA));

        }


        setPRIORS(_PRIORS);
        setCOMPONENTS(components);

}


void GMM_Model::saveInFiles()
{
    mat priors(1, nSTATES);
    mat mu(nVARS, nSTATES);
    mat sigma(nVARS, nVARS*nSTATES);
    std::ofstream varsfile ("GMM_vars.txt");

    for(uint i=0; i<nSTATES; i++)
    {
            priors(0,i) = getPRIORS()(i);
            for(uint j=0; j<nVARS; j++)
            {
                mu(j,i) = getCOMPONENTS().at(i).getMU()(j,0);
                for(uint k=0; k<nVARS; k++)
                    sigma(j,k + i*(nSTATES+1)) = getCOMPONENTS().at(i).getSIGMA()(j,k);
            }
    }

    for(uint i=0; i<nVARS; i++)
    {
        varsfile << vars_names.at(i);
        if(i<nVARS-1)
            varsfile << " ";
    }

    priors.save("GMM_priors.txt", raw_ascii);
    mu.save("GMM_mu.txt", raw_ascii);
    sigma.save("GMM_sigma.txt", raw_ascii);

    varsfile.close();
}

void GMM_Model::setPRIORS(const rowvec& priors)
{
    if(priors.n_elem == nSTATES)
        PRIORS = priors;
    else
        std::cout << "\n [ERROR]::GMM_Model::setPRIORS if(priors.n_elem == nSTATES) ... else .";
}


void GMM_Model::setCOMPONENTS(const std::vector<GaussianDistribution>& components)
{
    if(components.size() == nSTATES)
        COMPONENTS = components;
    else
    	std::cout << "\n [ERROR]::GMM_Model::setCOMPONENTS if(components.size() == nSTATES) ... else .";
}

uint GMM_Model::getIndexOfVARName(const std::string& varname)
{
        std::vector<std::string>::iterator it;
        it = std::find (vars_names.begin(), vars_names.end(), varname);
        return std::distance(vars_names.begin(), it);
}

uint GMM_Model::getNumVARS()
{
    return nVARS;
}

uint GMM_Model::getNumSTATES()
{
    return nSTATES;
}

rowvec&  GMM_Model::getPRIORS()
{
    return PRIORS;
}


std::vector<GaussianDistribution>&  GMM_Model::getCOMPONENTS()
{
    return COMPONENTS;
}

double GMM_Model::getProbability(const colvec& sample)
{
    double P = 0.0;
    for(uint k=0; k<nSTATES; k++)
        P += PRIORS.at(k)*as_scalar( COMPONENTS.at(k).getPDFValue(sample) );
    return P;
}



void GMM_Model::learnKMEANS()
{

    mat DemosTmp = DEMONSTRATIONS.at(0).getDatapoints().getData();

    for(int i=1; i<DEMONSTRATIONS.size(); i++)
        DemosTmp.insert_cols(DEMONSTRATIONS.at(0).getDatapoints().getNumPOINTS(), DEMONSTRATIONS.at(i).getDatapoints().getData());

    //Criterion to stop the EM iterative update
    double cumdist_threshold = 1e-10;
    uint maxIter = 100;

    //Initialization of the parameters
    uint nbVar = DemosTmp.n_rows;
    uint nbData = DemosTmp.n_cols;
    double cumdist_old = -std::numeric_limits<double>::max();
    uint nbStep = 0;

    //srand (time(NULL));

    urowvec idTmp = sort_index(randn<urowvec>(nbData));

    uvec allrows = linspace<uvec>(0, nbVar-1, nbVar);

    mat Mu = DemosTmp.submat(allrows, idTmp.cols(0,nSTATES-1));
    cube Sigma = zeros(nbVar,nbVar,nSTATES);


    //k-means iterations
    while(true)
    {
        mat distTmp = zeros(nSTATES,nbData);

        //E-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for (uint i=0; i<nSTATES; i++)
        {
            //Compute distances
            distTmp.row(i) = sum((DemosTmp - repmat(Mu.col(i), 1, nbData)) % (DemosTmp - repmat(Mu.col(i), 1, nbData)));

        }

        vec vTmp = zeros<vec>(nbData,1);
        vec idList = zeros<vec>(nbData,1);
        uword index;
        vec tempVec;
        double tempId;
        for (uint i=0; i<nbData; i++)
        {
            tempVec = distTmp.col(i);
            //cout << endl << "tempVec " << i << endl << tempVec << endl;

            vTmp.row(i) = min(tempVec);
            //cout << endl << "vTmp " << i << endl << vTmp.row(i) << endl;

            tempId = tempVec.min(index);
            //cout << endl << "tempId " << i << endl << index << endl;

            idList.row(i) = index;
        }
        //cout << endl << "vTmp" << endl << vTmp << endl;

        //cout << "idList" << idList << endl;

        double cumdist = sum(vTmp);
        //cout << endl << "cumdist: " << cumdist << endl;


        //M-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        uvec idTmp1;
        rowvec priors = zeros(1,nSTATES);
        mat sigTmp = zeros(nbVar,nbVar);
        for (uint i=0; i<nSTATES; i++)
        {
            idTmp1 = find(idList == i);
            //cout << endl << "idTmp1    " << idTmp1.n_elem << endl;
            priors(i) = idTmp1.n_elem;
            //cout << endl << "Priors_tpm :" << priors(i) << endl;
            Mu.col(i) = mean(DemosTmp.submat(allrows, idTmp1), 1);
            //cout << endl << "Mu_tpm :" << Mu.col(i) << endl;
            sigTmp = cov( trans(join_rows(DemosTmp.submat(allrows, idTmp1),DemosTmp.submat(allrows, idTmp1))) ,0 );
            //cout << endl << "sigTmp" << sigTmp << endl;
            Sigma.slice(i) = sigTmp;
            //cout << endl << "Sigma.slice(i)" << Sigma.slice(i) << endl;
        }

        priors = priors/nbData;

        //cout << endl << "Priors_tpm :" << priors << endl;
        //cout << endl << "Mu_tmp :" << Mu << endl;
        //cout << endl << "Sigma_tmp :" << Sigma << endl;

        //cout << endl << "cumdist_old   " << cumdist_old << endl;

        //Stopping criterion %%%%%%%%%%%%%%%%%%%%
        if (fabs(cumdist-cumdist_old) < cumdist_threshold)
        {
            cout << "%%%%%%%%%%% KMEANS %%%%%%%%%%%%%%%%%%" << endl;
            cout << endl << "Priors_kmeans :" << endl << priors << endl;
            cout << endl << "Mu_kmeans :" << endl << Mu << endl;
            cout << endl << "Sigma_kmeans :" << endl << Sigma << endl;
            cout << "%%%%%%% KMEANS FINISHED %%%%%%%%%%%%%" << endl;

            GaussianDistribution *GD;
            colvec mu;

            for(int i=0; i<nSTATES; i++)
            {
                PRIORS(i) = priors(i);
                mu = Mu.col(i);                
                COMPONENTS.push_back(GaussianDistribution(mu, Sigma.slice(i)));
            }

            return;
        }
        cumdist_old = cumdist;
        nbStep = nbStep + 1;

        //cout << endl << "nbStep   " << nbStep << endl;



    }

}

uint GMM_Model::EM_learn()
{
    learnKMEANS();
    return EM(std::numeric_limits<double>::min());
}



double GMM_Model::getLikelihood(const mat &SAMPLES)
{
    double L=0.0;
    colvec s(SAMPLES.n_rows);

    for(uint j=0; j<SAMPLES.n_cols; j++)
    {
        s = SAMPLES.col(j);
        L += log( getProbability( s ) );
    }
    return L;
}

bool GMM_Model::EM_isfinished(double l_old, double l_new)
{
    //cout << "\n error=" <<  fabs((l_new/l_old) - 1.0);

    if ( fabs((l_new/l_old) - 1.0)  <= THRESHOLD_MIN)
        return true;
    return false;
}


uint GMM_Model::EM(double likelihood)
{    
    double likelihood_new;

    uint k,i;
    rowvec E;
    colvec mu_tmp;
    mat Pxi, Pix, Pix_tmp, DataTmp1, sigma_tmp;

    mat demos = DEMONSTRATIONS.at(0).getDatapoints().getData();

    for(int i=1; i<DEMONSTRATIONS.size(); i++)
        demos.insert_cols(DEMONSTRATIONS.at(0).getDatapoints().getNumPOINTS(), DEMONSTRATIONS.at(i).getDatapoints().getData());


    // [BEGIN] E_STEP

        rowvec prior = PRIORS;

        Pxi = zeros(demos.n_cols, nSTATES);
        for(k=0; k<nSTATES; k++)
        {

            Pxi.col(k) = COMPONENTS.at(k).getPDFValue( demos );
            //cout << "Mu " << endl << COMPONENTS.at(k)->getMU();
            //cout << "Sigma " << endl << COMPONENTS.at(k)->getSIGMA();

            //cout << "Pxi " << endl << Pxi.col(k);
        }


        // compute posterior probability p(i|x)
        Pix_tmp = repmat(prior, demos.n_cols, 1) % Pxi;
        Pix = Pix_tmp / repmat(sum(Pix_tmp, 1), 1, nSTATES);
        // compute cumulated posterior probability
        E = sum(Pix);

        //cout << endl << "//////////////// EM ////////////////////" << endl;
        //cout << endl << "EM_E: " << endl << E << endl;


    // [END] E_STEP

    // [BEGIN] M_STEP

        for (i=0; i<nSTATES; i++)
        {
            // update the priors
            PRIORS(i) = E.col(i) / demos.n_cols;
            // update the centers
            mu_tmp = demos * Pix.col(i) / E.col(i);
            COMPONENTS.at(i).setMU(mu_tmp);
            // update the covariance matrices
            DataTmp1 = demos - repmat( COMPONENTS.at(i).getMU(), 1,demos.n_cols);

            //cout << "DataTmp1: " << endl << DataTmp1 << endl;

            sigma_tmp = (repmat(trans(Pix.col(i)), nVARS, 1) % DataTmp1 * trans(DataTmp1)) / E.col(i);
            sigma_tmp = sigma_tmp  + 1E-5 * eye(nVARS,nVARS);

            //cout << "sigma_tmp: " << endl << sigma_tmp << endl;

            COMPONENTS.at(i).setSIGMA( sigma_tmp );
        }

    // [END] M_STEP
        /*
        cout << "EM_Priors: " << endl << PRIORS << endl;
        cout << "Sum_EM_Priors: " << endl << sum(PRIORS) << endl;

        cout << "EM_MU: " << endl << COMPONENTS.at(0)->getMU() << endl;
        cout << "EM_MU: " << endl << COMPONENTS.at(1)->getMU() << endl;
        cout << "EM_MU: " << endl << COMPONENTS.at(2)->getMU() << endl;

        cout << "EM_SIGMA: " << endl << COMPONENTS.at(0)->getSIGMA() << endl;
        cout << "EM_SIGMA: " << endl << COMPONENTS.at(1)->getSIGMA() << endl;
        cout << "EM_SIGMA: " << endl << COMPONENTS.at(2)->getSIGMA() << endl;
        */

    likelihood_new = getLikelihood( demos );

    if( EM_isfinished(likelihood, likelihood_new) )
        return 1;

    return  EM(likelihood_new)+1;
}

void GMM_Model::setVARSNames(const std::vector<std::string>& vars)
{
    if(vars.size() == nVARS)
        vars_names = vars;
    else
        std::cout << "\n [ERROR]::GMM_Model::setVARSNames if(vars.size() == nVARS) ... else .";
}

std::vector<std::string>& GMM_Model::getVARSNames()
{
    return vars_names;
}


bool GMM_Model::addDemo(Demonstration& demo)
{

    if(demo.getDatapoints().getNumVARS() == nVARS)
        DEMONSTRATIONS.push_back(demo);
    else
        return false;

    return true;
}

} //end of pbdlib namespace
