/**
Copyright (C) 2014, Tohid Alizadeh, Milad Malekzadeh, Leonel Rozo

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

#include "pbdlib/pgmm.h"


namespace pbdlib
{


PGMM_Model::PGMM_Model(uint nVARS, uint nSTATES, uint nPARAMS){
	this->nVARS = nVARS;
	this->nSTATES = nSTATES;
    	this->nPARAMS = nPARAMS;
	this->PRIORS.resize(this->nSTATES);
  	this->GMMS.reserve(this->nPARAMS);
}


PGMM_Model::PGMM_Model(std::vector<Demonstration> &demos, uint nSTATES, uint nPARAMS)
{
    this->DEMONSTRATIONS = demos;
    this->nVARS = demos.at(0).getDatapoints().getNumVARS();
    this->nSTATES = nSTATES;
    this->PRIORS = rowvec(nSTATES);
    this->setVARSNames( demos.at(0).getDatapoints().getVarNames());
}


void PGMM_Model::loadPGMMfromMATLAB(std::string priorsFileName, std::string MuFileName, std::string SigmaFileName){
	cout << "Loading PGMM..." << endl;

	vec     priorsFile;     // To load the priors vector of the PGMM
	mat     MuFile;        // To load all the Mu of the PGMM
	mat     SigmaFile;     // To load all the Sigma of the PGMM
	uint	  ind1, ind2;     // Auxiliar indexes

	// Loading file containing the priors of the model
	priorsFile.load(priorsFileName, raw_ascii);
		// priorsFile.print();	// DEBUGGING

	for(uint i = 0 ; i < nSTATES ; i++)
		PRIORS(i) = priorsFile(i);
	//PRIORS.print("Priors = "); // DEBUGGING

	// Loading files containing the Zmu and Zsigma
	MuFile.load(MuFileName,raw_ascii);
		// MuFile.print(); // DEBUGGING
		// cout << "Mu data cols = " << MuFile.n_cols << ", rows = " << MuFile.n_rows << endl; // DEBUGGING
	SigmaFile.load( SigmaFileName,raw_ascii);
		// SigmaFile.print(); // DEBUGGING
		// cout << "Sigma data cols =" << SigmaFile.n_cols << ", rows = " << SigmaFile.n_rows << endl; // DEBUGGING

    for (uint m = 0; m < nPARAMS ; m++){
		// Resizing the ZSigma vector and initializing Zmu

		//ZMU.push_back(zeros(nVARS, nSTATES));
		//ZSIGMA.push_back(zeros(nVARS, nVARS, nSTATES));

		// Storing loaded data
	    for (uint i = 0; i < nSTATES; i++){
	    	// Mu
	      GMMS.at(m).getCOMPONENTS().at(i).setMU(MuFile.col((m*nSTATES) + i)); // Getting the i-th col (Zmu_i) for m-th frame
//	    	ZMU.at(m).col(i) = ZmuFile.col((m*nSTATES) + i);	// Getting the i-th col (Zmu_i) for m-th frame
	    	//ZMU.at(m).print("ZMu= "); //DEBUGGING
	    	// Sigma
	        ind1 = (m * nSTATES * nVARS) + (i * nVARS);        		// Initial index
	        ind2 = (m * nSTATES * nVARS) + ((i+1) * nVARS) - 1; 	// Last index
          GMMS.at(m).getCOMPONENTS().at(i).setSIGMA(SigmaFile.cols(ind1, ind2));
//	        ZSIGMA.at(m).slice(i) = ZsigmaFile.cols(ind1, ind2);   // Getting the i-th matrix (ZSigma_i) for the m-th frame
	        //ZSIGMA.at(m).print("ZSigma = "); // DEBUGGING
	    }
	}
	// cout << "ZMU size = " << ZMU.size() << endl; // DEBUGGING
	// cout << "ZSIGMA size = " << ZSIGMA.size() << endl; // DEBUGGING
	cout << "PGMM loaded from text files!" << endl;
}


GMM_Model* PGMM_Model::getGMM(std::vector<mat> A, std::vector<colvec> b){
    GMM_Model*	gmm;	// Resulting GMM
    gmm = new GMM_Model(nSTATES, nVARS);
    std::vector<GaussianDistribution> componentsTmp;	// Auxiliar Gaussian dist. vector
    componentsTmp.reserve(nSTATES);
    GaussianDistribution *singleStateTmp;				// Auxuliar Gaussian distribution
    colvec localMu = zeros(nVARS, 1);
    colvec MuTmp;
    mat    localSigma = zeros(nVARS, nVARS);
    mat    SigmaTmp;

    // cout << "A's size = " << A.size() << endl; // DEBUGGING
    // cout << "b's size = " << b.size() << endl; // DEBUGGING

    //Projections in the different candidate frames (see Eq. (3) Calinon et al, Humanoids'2012 paper)
    for (uint i = 0 ; i < nSTATES ; i++){
        MuTmp = zeros(nVARS, 1);
        SigmaTmp = zeros(nVARS, nVARS);

        for (uint m = 0; m < nPARAMS ; m++){
            //A.at(m).print("A = "); // DEBUGGING
            //b.at(m).print("b = "); // DEBUGGING
            //localMu = A.at(m) * ZMU.at(m).col(i) + b.at(m);		// Projecting Zmu of state "i" at frame "m" by Milad
            localMu = A.at(m) * GMMS.at(m).getCOMPONENTS().at(i).getMU() + b.at(m);		// Projecting Zmu of state "i" at frame "m"
            //localSigma = A.at(m) * ZSIGMA.at(m).slice(i) * trans(A.at(m));	// Projecting Zsigma of state "i" at frame "m" by Milad
            localSigma = A.at(m) * GMMS.at(m).getCOMPONENTS().at(i).getSIGMA() * trans(A.at(m));	// Projecting Zsigma of state "i" at frame "m"
            SigmaTmp += inv(localSigma);						// Accumulative product of Gaussians (Covariance matrix)
            MuTmp += inv(localSigma) * localMu;					// Accumulative product of Gaussians (mean)
        }
        // Saving the Gaussian distribution for state "i"
        singleStateTmp = new GaussianDistribution(nVARS);
        singleStateTmp->setSIGMA(inv(SigmaTmp));
        singleStateTmp->setMU(inv(SigmaTmp) * MuTmp);
        //singleStateTmp->getMU().print("MU = "); // DEBUGGING
        //singleStateTmp->getSIGMA().print("SIGMA = "); // DEBUGGING
        componentsTmp.push_back(*singleStateTmp);
    }
    gmm->setPRIORS(this->PRIORS);
    gmm->setCOMPONENTS(componentsTmp);
    return 		gmm;
}
/*----------------------------------------------------------------------------*/
GMM_Model* PGMM_Model::prodGMM(std::vector<Parameters> params){
    GMM_Model*	gmm;	// Resulting GMM
    gmm = new GMM_Model(nSTATES, nVARS);
    std::vector<GaussianDistribution> componentsTmp;	// Auxiliar Gaussian dist. vector
    componentsTmp.reserve(nSTATES);
    GaussianDistribution *singleStateTmp;				// Auxuliar Gaussian distribution
    colvec localMu = zeros(nVARS, 1);
    colvec MuTmp;
    mat    localSigma = zeros(nVARS, nVARS);
    mat    SigmaTmp;

    std::vector<mat> A;
    std::vector<colvec> b;

    A.resize(nPARAMS);
    b.resize(nPARAMS);
    for (uint m=0; m<nPARAMS; m++)
        {
          A[m] = params.at(0).getParams(0).at(m).A;
          b[m]= params.at(0).getParams(0).at(m).b;
        }

    // cout << "A's size = " << A.size() << endl; // DEBUGGING
    // cout << "b's size = " << b.size() << endl; // DEBUGGING

    //Projections in the different candidate frames (see Eq. (3) Calinon et al, Humanoids'2012 paper)
    for (uint i = 0 ; i < nSTATES ; i++){
        MuTmp = zeros(nVARS, 1);
        SigmaTmp = zeros(nVARS, nVARS);

        for (uint m = 0; m < nPARAMS ; m++){
            //A.at(m).print("A = "); // DEBUGGING
            //b.at(m).print("b = "); // DEBUGGING
            //localMu = A.at(m) * ZMU.at(m).col(i) + b.at(m);		// Projecting Zmu of state "i" at frame "m" by Milad
            localMu = A.at(m) * GMMS.at(m).getCOMPONENTS().at(i).getMU() + b.at(m);		// Projecting Zmu of state "i" at frame "m"
            //localSigma = A.at(m) * ZSIGMA.at(m).slice(i) * trans(A.at(m));	// Projecting Zsigma of state "i" at frame "m" by Milad
            localSigma = A.at(m) * GMMS.at(m).getCOMPONENTS().at(i).getSIGMA() * trans(A.at(m));	// Projecting Zsigma of state "i" at frame "m"
            SigmaTmp += inv(localSigma);						// Accumulative product of Gaussians (Covariance matrix)
            MuTmp += inv(localSigma) * localMu;					// Accumulative product of Gaussians (mean)
        }
        // Saving the Gaussian distribution for state "i"
        singleStateTmp = new GaussianDistribution(nVARS);
        singleStateTmp->setSIGMA(inv(SigmaTmp));
        singleStateTmp->setMU(inv(SigmaTmp) * MuTmp);
        //singleStateTmp->getMU().print("MU = "); // DEBUGGING
        //singleStateTmp->getSIGMA().print("SIGMA = "); // DEBUGGING
        componentsTmp.push_back(*singleStateTmp);
    }
    gmm->setPRIORS(this->PRIORS);
    gmm->setCOMPONENTS(componentsTmp);
    return 		gmm;
}
/*----------------------------------------------------------------------------*/
uint PGMM_Model::EM_TensorLearn()    //Tohid: I should implement the Tensor EM here.
{
    initTensorGMMTimeBased();          //Tohid: equally spaced time based initialization for the tensor GMM
    EM_tensorGMM();
    return true;  //Tohid: I just set this as a temporary solution, to test other things.
}
/*----------------------------------------------------------------------------*/
void PGMM_Model::initTensorGMMTimeBased()     //Tohid: I should implement the time based initialization of the Tensor GMM here.
{
    mat DataAll;
    mat Atmp, Mu;
    mat btmp;
    cube Sigma;
    vec TimingSep, idTmp2, Priors;
    double diagRegMat = 1e-4;
    mat DemosTmp = DEMONSTRATIONS.at(0).getDatapoints().getData();


    cout << "\n size= " << DEMONSTRATIONS.size();
    cout << "\n npoints= " << DEMONSTRATIONS.at(0).getDatapoints().getNumPOINTS();
    cout << "\n data= " << DEMONSTRATIONS.at(0).getDatapoints().getData();
    for(int i=0; i<DEMONSTRATIONS.size(); i++)
        DemosTmp.insert_cols(DEMONSTRATIONS.at(0).getDatapoints().getNumPOINTS(), DEMONSTRATIONS.at(i).getDatapoints().getData());



    cout<<"\n Demonstrations are loaded now!!!"<<endl;
    cout<<"\n DemosTmp.n_cols = "<<DemosTmp.n_cols<<endl;

    //cout<<"DemosTmp.n_rows = "<<DemosTmp.n_rows<<endl;
    //cout<<DEMONSTRATIONS.at(0).getParameters().getParams(0).at(1).b<<endl;
    //cout<<DEMONSTRATIONS.at(1).getParameters().getParams(0).at(1).b<<endl;
//    cout<<"End of the initTensorPGMMtimeBased!!!"<<endl;
//    cout<<DEMONSTRATIONS.at(0).getParameters().getParams(1).at(1).b<<endl;
//    cout<<DEMONSTRATIONS->getParameters().getParams(1).at(0).b<<endl;
//    Demo->getParameters().getParams(199).at(1).A
    //cout<<"Parameters.getParams(0).at(0).b =  "<<Parameters.params.at(0).getParams(1).at(1).b<<endl;
    uint nDemos = DEMONSTRATIONS.size();  //Tohid: I should determine the number of demonstrations, given the demos
    uint nData = DEMONSTRATIONS.at(0).getDatapoints().getNumPOINTS();  //tohid: I should determine the number of Data here.
    uint nPARAMS = DEMONSTRATIONS.at(0).getParameters().getParams(0).size();
    uint nVars = DemosTmp.n_rows;
//    cout<<"nDemos =    "<<nDemos<<endl;     cout<<"nData =     "<<nData<<endl;
//   cout<<"nPARAMS =   "<<nPARAMS<<endl;    cout<<"nVars =     "<<nVars<<endl;
//    cout<<"nVARS =     "<<nVARS<<endl;
    DataAll = zeros(nPARAMS*nVars, nData*nDemos);
    mat DataAllTmp1 = DataAll;
    mat DataAllTmp2 = zeros(nPARAMS*nVars, nData);
    for (uint n=0; n<nDemos; n++){
        DemosTmp = DEMONSTRATIONS.at(n).getDatapoints().getData();
        for(int m=0; m<nPARAMS; m++){
          //if(DEMONSTRATIONS.at(n).getFrames().at(m)==1){  //This is to check if the frame "m" was available during the demonstration.
            Atmp = DEMONSTRATIONS.at(n).getParameters().getParams(0).at(m).A;
            btmp = DEMONSTRATIONS.at(n).getParameters().getParams(0).at(m).b;
            DataAllTmp2.rows(m*nVars,(m+1)*nVars-1) = Atmp*DemosTmp+repmat(btmp,1, nData);
            //DataAll.submat(m*nVars,i*nData, (m+1)*nVars-1,(i+1)*nData-1) = Atmp*DemosTmp+repmat(btmp,1, nData);
          //}
        }
        DataAll.cols(n*nData, (n+1)*nData-1) = DataAllTmp2;
        DEMONSTRATIONS.at(n).getDatapoints().setDataAll(DataAllTmp2);
        //cout<<"DEMONSTRATIONS.at(i).getDatapoints().getDataAll()"<<endl;
        //cout<<DEMONSTRATIONS.at(i).getDatapoints().getDataAll()<<endl;
//        mat DataTMP1 = DEMONSTRATIONS.at(i).getDatapoints().getData();
//        mat DataTMP2 = DEMONSTRATIONS.at(i).getDatapoints().getDataAll();
//        cout<<"DataTMP1.size:  "<<DataTMP1.n_rows<<" x "<<DataTMP1.n_cols<<endl;
//        cout<<"DataTMP2.size:  "<<DataTMP2.n_rows<<" x "<<DataTMP2.n_cols<<endl;
    }
    TimingSep = linspace(min(DataAll.row(0)), max(DataAll.row(0)),nSTATES+1);
//    cout<<"TimingSep:  "<<endl<<TimingSep<<endl;
    Mu = zeros(nVars*nPARAMS,nSTATES);
    Sigma = zeros(nVars*nPARAMS,nVars*nPARAMS,nSTATES);
    Priors = zeros(nSTATES);
    for (int i=0; i<nSTATES; i++){
        idTmp2=zeros(nData*nDemos);
        uvec idtmp = find(DataAll.row(0)>=TimingSep(i));
        int m=0;
        for (int n=0; n<idtmp.n_elem; n++){
            if (DataAll(0,idtmp(n))<=TimingSep(i+1)){
                idTmp2(m) = idtmp(n);
                DataAllTmp1.col(m) = DataAll.col(idtmp(n));
                m+=1;
            }
        }
        mat DataAllTmp = DataAllTmp1.cols(0,m-1);
        Mu.col(i) = sum(DataAllTmp,1)/m;
        Sigma.slice(i) = cov(trans(DataAllTmp))+eye(nVars*nPARAMS,nVars*nPARAMS)*diagRegMat;
        Priors(i) = m;
    }
    Priors = Priors/sum(Priors);
    // In order to shape the GMMs of the PGMM
    for (int m=0; m<nPARAMS; m++){
        //cout<<"m =   : "<<m<<endl;
        GMM_Model*	gmm;	// Resulting GMM
        gmm = new GMM_Model(nSTATES, nVARS);
        std::vector<GaussianDistribution> componentsTmp;	// Auxiliar Gaussian dist. vector
        componentsTmp.reserve(nSTATES);
        GaussianDistribution *singleStateTmp;				// Auxuliar Gaussian distribution
        //colvec localMu = zeros(nVARS, 1);
        colvec MuTmp;
        mat    localSigma = zeros(nVARS, nVARS);
        mat    SigmaTmp;
        for (int i=0; i<nSTATES; i++){
            //cout<<"i =   : "<<i<<endl;
            MuTmp = Mu.submat(m*nVARS,i,(m+1)*nVARS-1,i);
            SigmaTmp = Sigma.slice(i);
            localSigma = SigmaTmp.submat(m*nVARS,m*nVARS, (m+1)*nVARS-1,(m+1)*nVARS-1);
            //cout<<"SigmaTmp =   "<<endl<<SigmaTmp<<endl;
            //cout<<"localSigma = "<<endl<<localSigma<<endl;
            //getchar();
            // Saving the Gaussian distribution for state "i"
            singleStateTmp = new GaussianDistribution(nVARS);
            singleStateTmp->setSIGMA((localSigma));
            singleStateTmp->setMU(MuTmp);
//            cout<<"singleStateTmp->setMU(MuTmp):    "<<MuTmp<<endl;
//            getchar();
            //singleStateTmp->getMU().print("MU = "); // DEBUGGING
            //singleStateTmp->getSIGMA().print("SIGMA = "); // DEBUGGING
            componentsTmp.push_back(*singleStateTmp);
        }
        //cout<<"...... "<<this->PRIORS<<endl;
        gmm->setPRIORS(trans(Priors));
        //gmm->setPRIORS(this->PRIORS);
        gmm->setCOMPONENTS(componentsTmp);
        GMMS.push_back(*gmm);
    }
//    cout<<"GMMS.size():    "<<GMMS.size()<<endl;
}
/*----------------------------------------------------------------------------*/
void PGMM_Model::EM_tensorGMM()
{
    //Parameters of the EM algorithm
    uint nbMinSteps = 5; //Minimum number of iterations allowed
    uint nbMaxSteps = 100; //Maximum number of iterations allowed
    double maxDiffLL = 1E-4; //Likelihood increase threshold to stop the algorithm
    double diagRegularizationFactor = 1E-4;
    uint nDATA = 0;     //nDATA is the total number of data points including all the demonstrations
    uint nDemos = DEMONSTRATIONS.size();
    uint nPoints = DEMONSTRATIONS.at(0).getDatapoints().getNumPOINTS();
    nPARAMS = DEMONSTRATIONS.at(0).getParameters().getParams(0).size();
    mat DataAll= zeros(nVARS*nPARAMS, nDemos*nPoints);
    mat DataTmp = zeros(nVARS, nDemos*nPoints);
    mat DataTmp2 = zeros(nVARS, nDemos*nPoints);
    for (uint n=0; n<DEMONSTRATIONS.size(); n++){
        nDATA+=DEMONSTRATIONS.at(n).getDatapoints().getNumPOINTS();
        DataAll.cols(n*nPoints, (n+1)*nPoints-1) = DEMONSTRATIONS.at(n).getDatapoints().getDataAll();
    }
    std::vector<mat> GAMMA0;
    mat GAMMA = zeros(nVARS, nDATA);
    mat GAMMA2 = zeros(nVARS, nDATA);
    vec LL = zeros(nbMaxSteps);
    GAMMA0.resize(nPARAMS);
    for (uint m=0;m<nPARAMS; m++){
        GAMMA0.at(m) = zeros(nVARS, nDATA);
    }
    cout<<"Start to learn PGMM using tensor EM"<<endl;
    for (uint nbIter=0; nbIter<nbMaxSteps; nbIter++){
        // Expectation-Step:
        mat L = ones(nVARS,nDATA);
        for (uint m=0; m<nPARAMS; m++){
          //if(DEMONSTRATIONS.at(n).getFrames().at(m)==1){  //This is to check if the frame "m" was available during the demonstration.
            DataTmp = DataAll.rows(m*nVARS, (m+1)*nVARS-1);
            for(uint i=0; i<nSTATES; i++){
                GAMMA0.at(m).row(i) = GMMS.at(m).getPRIORS().at(i)*trans(GMMS.at(m).getCOMPONENTS().at(i).getPDFValue(DataTmp));
                L.row(i) = L.row(i)%GAMMA0.at(m).row(i);
            }
          //}
        }
        GAMMA = L/repmat(sum(L,0)+REALMIN,L.n_rows,1);
//        GAMMA2 = GAMMA/repmat(sum(GAMMA,1)+REALMIN,1,GAMMA.n_cols);
        GAMMA2 = GAMMA/repmat(sum(GAMMA,1),1,GAMMA.n_cols);//I removed the realmin from here to make it the same as the Matlab code
        // Maximization-Step:
        vec MuTmp = zeros(nVARS);
        mat SigmaTmp = zeros(nVARS, nVARS);
        for (uint m=0; m<nPARAMS; m++){
            DataTmp = DataAll.rows(m*nVARS, (m+1)*nVARS-1);
            for(uint i=0; i<nSTATES; i++){
                //update PRIORS
                PRIORS(i) = sum(sum(GAMMA.row(i)))/GAMMA.n_cols;
                //update Mu
                MuTmp = DataTmp*trans(GAMMA2.row(i));
                GMMS.at(m).getCOMPONENTS().at(i).setMU(MuTmp);
                //update Sigma
                DataTmp2 = DataTmp - repmat(MuTmp, 1, DataTmp.n_cols);
                SigmaTmp = DataTmp2*diagmat(GAMMA2.row(i))*trans(DataTmp2)+eye(nVARS,nVARS)*diagRegularizationFactor;
                GMMS.at(m).getCOMPONENTS().at(i).setSIGMA(SigmaTmp);
            }
            GMMS.at(m).setPRIORS(PRIORS);
        }
        //  Compute average log-likelihood
        mat LL_Tmp;
        LL(nbIter) = 0;
        for(int n=0; n<nDATA; n++){
          LL_Tmp = log(sum(L.col(n),0));
          LL(nbIter) += as_scalar(sum(LL_Tmp,1));
        }
        LL(nbIter)/= nDATA;
        //  Stop the algorithm if EM converged (small change of LL)
        if (nbIter>nbMinSteps){
            if (LL(nbIter)-LL(nbIter-1)<maxDiffLL || nbIter==nbMaxSteps-1){
                cout<<"EM converged after "<<nbIter<<" iterations"<<endl;
                cout<<"The results of the model after the learning by EM:"<<endl;
                cout<<"PRIORS"<<endl;
                cout<<PRIORS<<endl;
                for (int i=0; i<nSTATES; i++){
                  for(int m=0; m<nPARAMS; m++){
                    cout<<"GMMS.at("<<m<<").getCOMPONENTS().at("<<i<<").getMU():     "<<endl;
                    cout<<GMMS.at(m).getCOMPONENTS().at(i).getMU()<<endl;
                    cout<<"GMMS.at("<<m<<").getCOMPONENTS().at("<<i<<").getSIGMA():  "<<endl;
                    cout<<GMMS.at(m).getCOMPONENTS().at(i).getSIGMA()<<endl;
                  }
                }
                break;
            }
        }
        if (nbIter==nbMaxSteps){
          cout<<"EM converged after "<<nbMaxSteps<<" iterations"<<endl;
        }
    }
}

uint PGMM_Model::getNumVARS()
{
    return this->nVARS;
}

std::vector<GMM_Model>& PGMM_Model::getGMMS(){
    return this->GMMS;
}


uint PGMM_Model::getNumSTATES()
{
    return this->nSTATES;
}

uint PGMM_Model::getNumFRAMES()
{
    return this->nPARAMS;
}

rowvec&  PGMM_Model::getPRIORS()
{
    return this->PRIORS;
}

std::vector<std::string>& PGMM_Model::getVARSNames()
{
    return this->vars_names;
}

void PGMM_Model::setVARSNames(const std::vector<std::string>& vars)
{
    if(vars.size() == nVARS)
        this->vars_names = vars;
    else
        std::cout << "\n [ERROR]::PGMM_Model::setVARSNames if(vars.size() == nVARS) ... else .";
}

void PGMM_Model::setFRAMESNames(const std::vector<std::string>& frames)
{
    if(frames.size() == nPARAMS)
        this->frames_names = frames;
    else
        std::cout << "\n [ERROR]::PGMM_Model::setFRAMESNames if(frames.size() == nPARAMS) ... else .";
}

std::vector<std::string>& PGMM_Model::getFRAMESNames()
{
    return this->frames_names;
}



} //end of namespace pbdlib
