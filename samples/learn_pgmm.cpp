/*
* Copyright (c) 2013
* - Davide De Tommaso @ dtmdvd[at]gmail[dot]com
* - Milad Malekzadeh @ milad[dot]malekzadeh[at]gmail[dot]com
* - Leonel Rozo @ leonel[dot]rozo[at]iit[dot]it
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of the <organization> nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY <copyright holder> ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "pbdlib/pgmm.h"
#include "pbdlib/gmr.h"
#include "pbdlib/datapoints.h"
#include "pbdlib/parameters.h"
#include <sstream>

using namespace pbdlib;

int main(int argc, char **argv)
{
  // Model variables
  uint nFrames = 2;
  uint nStates = 3;
  uint nVars = 3; //Tohid: t, x, y
  uint nDemos = 4;//Tohid: number of Demonstrations, for easier loading from the data files
  std::vector<std::string> varNames;
  varNames.push_back("t");    // time is input (t) and the attractor position (2D vector [x y]) is the output.
  varNames.push_back("x");
  varNames.push_back("y");
  //..................................................................................
  //......... loading the demonstrations and task parameters from the txt files ......
  std::vector<Demonstration> demos;
  Demonstration demo =  Demonstration(3,200,2);  //(3 vars, 200 datapoints, 2 frames of references)

  std::vector<Parameters> Params;


  char filename[256];
  cout<<"Loading the demonstrations and the task parameters ...";
  for (int m=0; m<nDemos; m++){   //Loading demos in the loop
      Parameters param = Parameters(3,2,1);  //nVAR, nPAR, nPOINTS (set to 1 to show that the task parameters are fixed all the time)
      sprintf(filename, "../../data/pgmm/Data0%d.txt",m+1);
      demo.getDatapoints().loadFromFile(filename, true);
      demo.getDatapoints().setVarNames(varNames);
      sprintf(filename, "../../data/pgmm/Param0%d.txt",m+1);
      param.loadFromFile(filename);
      Params.push_back(param);
      demo.setParams(param);
      demos.push_back(demo);
  }
  cout<<endl<<"Demonstrations and task parameters are loaded successfully !!!!"<<endl;
  cout<<"Press any key to continue (learning the PGMM model for the loaded demonstrations)"<<endl;
  getchar();

  //..................................................................................
  //......... Learning the PGMM model parameters using the EM algorithm ......
  cout<<"Learning model parameters using EM ...";
  PGMM_Model *pgmm;
  //pgmm = new PGMM_Model(nVars, nStates, nFrames);
  pgmm = new PGMM_Model(demos, nStates, nFrames);
  //    pgmm->loadPGMMfromMATLAB(priorsName, ZmuName, ZsigmaName);   //Tohid: commented, because we do not have Zmu and ZSigma anymore
  //cout<<"demos.size() =  "<<demos.size()<<endl;
  cout << "\n EM rounds = " << pgmm->EM_TensorLearn()<<endl;
  cout<<"Learning is done and the model is learned for the provided demonstrations!!!"<<endl;
  // saving the learned model parameters
  mat MuTmp = zeros(nVars, nStates*nFrames);
  mat SigmaTmp = zeros(nVars, nVars*nStates*nFrames);
  mat Priors = zeros(1,nStates);
  std::vector<GMM_Model> GMMSresult;
  GMMSresult = pgmm->getGMMS();
  for(uint m=0; m<nFrames; m++){
    for (uint i=0;i<nStates; i++){
      MuTmp.col(m*nStates+i) = GMMSresult.at(m).getCOMPONENTS().at(i).getMU();
      SigmaTmp.cols(m*nVars*nStates+i*nVars, m*nVars*nStates+i*nVars+nVars-1) = GMMSresult.at(m).getCOMPONENTS().at(i).getSIGMA();
    }
  }
  sprintf(filename, "../../data/pgmm/Mu0%d.txt",2);
  MuTmp.save(filename, raw_ascii);
  sprintf(filename, "../../data/pgmm/Sigma0%d.txt",2);
  SigmaTmp.save(filename, raw_ascii);
  sprintf(filename, "../../data/pgmm/Priors0%d.txt",2);
  Priors = pgmm->getPRIORS();
  Priors.save(filename, raw_ascii);
  cout<<"the learned model parameters are written into text files succesfully."<<endl;
  cout<<"Press any key to continue (performing the GMR for the new set of task parameters)"<<endl;
  getchar();
  //..................................................................................

  // Path for model files
  //std::string priorsName = "../../data/pgmm/Priors02.txt";
  //std::string MuName = "../../data/pgmm/Mu02.txt";
  //std::string SigmaName = "../../data/pgmm/Sigma02.txt";

  // Reproduction variables (GMR)
  GMR			*gmr = new GMR();
  Datapoints	*repro = new Datapoints(1,1);				// Only one datapoint passed at each time step
  mat			reproData, auxYgmr;
  std::vector<std::string> reproVarNames;
  reproVarNames.push_back("t");							      // time is the query point
  repro->setVarNames(reproVarNames);
  // loading new parameters for the reproduction from the text file
  // In this case we have considered that the parameters are the same all along the reproduction. They can be varying.
  std::vector<mat> A_tmp;
  mat auxA;
  std::vector<colvec> b_tmp;
  colvec auxb;
  sprintf(filename, "../../data/pgmm/ParamRepro0%d.txt",1);
  mat ParamTmp;
  ParamTmp.load(filename, raw_ascii);
  b_tmp.push_back(trans(ParamTmp.row(0)));
  A_tmp.push_back(ParamTmp.rows(1,nVars));
  b_tmp.push_back(trans(ParamTmp.row(nVars+1)));
  A_tmp.push_back(ParamTmp.rows(nVars+2,2*nVars+1));
  for(uint tn = 1 ; tn <= 200 ; tn++){
      // Computing the resulting gmm given the set of parameters {A,b}
      GMM_Model* gmm;
      gmm = pgmm->getGMM(A_tmp, b_tmp);
      gmm->setVARSNames(varNames);
      // Printing model components
      cout << "Resulting GMM given the set of parameters 'A' and 'b'" << endl;
      for(uint i = 0 ; i < nStates ; i++){
              cout << "State #" << i << ":"<< endl;
              gmm->getCOMPONENTS().at(i).getMU().print("Mu = ");
              gmm->getCOMPONENTS().at(i).getSIGMA().print("Sigma = ");
      }
      // Computing GMR
      gmr->setGMMModel(gmm);
      mat repDataPt;
      repDataPt << (tn * 0.01);	// time step as input
      repDataPt.print("tn = ");
      repro->setData(repDataPt);
      gmr->regression(repro);
      gmr->get_y().print("y = ");

      cout << "Please press [ENTER] to continue or [CTRL+C] to finish" << endl;
      char key;
      std::cin.ignore(1);
  }
  return 0;
}




