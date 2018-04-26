// ################################################################################
// # Program to perform the full angular analysis of the decay B0 --> K*0 mu+ mu- #
// ################################################################################

#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH3.h>
#include <TF2.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TCutG.h>
#include <Math/Functor.h>

#include <TMinuit.h>
#include <RooMinimizer.h>
#include <RooNLLVar.h>
#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooAbsPdf.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooArgSet.h>
#include <RooFitResult.h>
#include <RooPolynomial.h>
#include <RooMCStudy.h>
#include <RooMinuit.h>
#include <RooNLLVar.h>
#include <RooWorkspace.h>
#include <RooConstVar.h>
#include <RooRandom.h>
#include <RooDataHist.h>
#include <RooFunctorBinding.h>
#include <RooStats/RooStatsUtils.h>

#include <ctime>
#include <iostream>
#include <utility>
#include <sstream>

#include "ReadParameters.h"
#include "Utils.h"
#include "B0KstMuMuTreeContent.h"
using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;
using std::ios_base;
using std::pair;
using std::make_pair;
using namespace RooFit;


// ####################
// # Global constants #
// ####################
#define NBINS        20
#define MULTYIELD     1. // Multiplication factor to the number of entry in toy-MC
#define NCOEFFPOLYBKG 5  // Maximum number of coefficients (= degree) of the polynomial describing the combinatorial bkg

#define nJPSIS 230000.0
#define nJPSIB   2500.0
#define nPSIPS  15000.0
#define nPSIPB   1500.0
 
#define _USE_MATH_DEFINES

// ##########################################
// # Internal flags to control the workflow #
// ##########################################
#define MAKEmumuPLOTS false
#define SETBATCH      false
#define PLOTBKG       false    //plot sideband bkg
#define PROFILENLL    false   //[true= plot the profile likelihood]
#define PLOT          true  //[true= plot the results]
#define SAVEPOLY      false // ["true" = save bkg polynomial coefficients in new parameter file; "false" = save original values]
#define SAVEPLOT      true   //2015-12-17
#define RESETsigANG   false // Reset signal angular parameters before starting the fit
#define RESETcomANG   false // Reset combinatorial bkg angular parameters before starting the fit
#define FULLTOYS      false // Run generation-and-fit toys
#define FUNCERRBAND   false // Show the p.d.f. error band
#define MINIMIZER     "Minuit" // Minimizer type for 3D MODEL actual fit ["Minuit"; "Minuit2"]
#define GENPARAMS     "All" // Option to generate parameters for parameter file: "All" "misTagFrac" "FlP5pFsAs" "combBkgAng"
#define SCAN          true    //control initial value when scan
#define doBKG         false     // true = do bkg uncertainty by control initial value of P5p P1 As5
// ##################
// # External files #
// ##################
#define PARAMETERFILEIN  "/python/ParameterFile.txt"
#define PARAMETERFILEOUT "ParameterFileOut.txt"


// ############################################
// # Global variables from configuration file #
// ############################################
double PsiYieldGoodTag, PsiYieldGoodTagErr;
double PsiYieldMisTag,  PsiYieldMisTagErr;
double LUMI;

string CTRLfitWRKflow;
string ParameterFILE;

vector<double> q2Bins;
vector<double> cosThetaKBins;
vector<double> cosThetaLBins;
vector<double> phiBins;

vector<vector<string>*>             fitParam;    // Vector containing the pointers to the vectors containing the starting values for the fit
vector<vector<unsigned int>*>       configParam; // Vector containing the pointers to the vectors containing the configuration parameters for the fit

// ####################
// # Global variables #
// ####################
TTree* theTree;
Utils* Utility;
B0KstMuMuTreeContent* NTuple;

ofstream fileFitResults;
ofstream fileFitSystematics;

double* q2BinsHisto;


// ####################################
// # Useful variables from the NTuple #
// ####################################
RooDataSet* toy;
RooDataSet* toy2;
RooDataSet* SingleCandNTuple_JPsi;
RooDataSet* SingleCandNTuple_PsiP;
RooDataSet* SingleCandNTuple_RejectPsi;
RooDataSet* SingleCandNTuple;
RooRealVar* B0mass;
RooRealVar* mumuMass;
RooRealVar* mumuMassE;
RooRealVar* ctK;                //CosThetaKArb
RooRealVar* ctL;                //CosThetaMuArb
RooRealVar* phi;                //PhiKstMuMuPlaneArb
RooRealVar* truthMatchSignal;
RooRealVar* rightFlavorTag;


// #################################
// # Variables and pdf for the fit #
// #################################

// ##################
// # Signal B0 mass #
// ##################
RooRealVar* meanS;

RooRealVar* sigmaS1;
RooAbsPdf*  MassS1;

RooRealVar* sigmaS2;
RooAbsPdf*  MassS2;

RooRealVar* fracMassS;
RooAbsPdf*  MassSignal;

// #################
// # Signal angles #
// #################
RooRealVar* FlS;
RooRealVar* P4pS;
RooRealVar* P5pS;
RooRealVar* P6pS;
RooRealVar* P8pS;
RooRealVar* FsS;
RooRealVar* AsS;
RooRealVar* P1S;
RooRealVar* P2S;
RooRealVar* P3S;              
RooRealVar* As5S;             
RooAbsPdf*  AngleS;


// ####################
// # Total signal pdf #
// ####################
RooAbsPdf* Signal;
RooAbsPdf* SignalT;

// ####################################
// # Combinatorial background B0 mass #
// ####################################
RooRealVar* var1;
RooRealVar* var2;
RooAbsPdf*  BkgMassExp1;
RooAbsPdf*  BkgMassExp2;

RooRealVar* fracMassBExp;
RooAbsPdf*  BkgMassComb;

// #########################
// # Mistag signal B0 mass #
// #########################
RooRealVar* sigmaMisTag1;
RooAbsPdf*  MassMisTag1;
RooRealVar* sigmaMisTag2;
RooAbsPdf*  MassMisTag2;
RooRealVar* fracMisTag;
RooAbsPdf*  MassMisTag;
RooAbsPdf*  AngleMisTag;

// ##############################
// # Peaking background B0 mass #
// ##############################
RooRealVar* meanR1;
RooRealVar* sigmaR1;
RooRealVar* meanR2;
RooRealVar* sigmaR2;
RooAbsPdf*  BkgMassRPeak1;
RooAbsPdf*  BkgMassRPeak2;

RooRealVar* fracMassBRPeak;
RooAbsPdf*  BkgMassRPeak;

RooRealVar* meanL1;
RooRealVar* sigmaL1;
RooRealVar* meanL2;
RooRealVar* sigmaL2;
RooAbsPdf*  BkgMassLPeak1;
RooAbsPdf*  BkgMassLPeak2;

RooRealVar* fracMassBLPeak;
RooAbsPdf*  BkgMassLPeak;

RooRealVar* fracMassBPeak;
RooAbsPdf*  BkgMassPeak;

// #####################
// # Background angles #
// #####################
RooRealVar* p1Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleP1;
RooRealVar* c1Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleC1;

RooRealVar* p2Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleP2;
RooRealVar* c2Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleC2;

RooRealVar* p3Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleP3;
RooRealVar* c3Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleC3;

RooAbsPdf* BkgAnglesC;
RooAbsPdf* BkgAnglesP;

// ########################
// # Total background pdf #
// ########################
RooAbsPdf* BkgMassAngleComb;
RooAbsPdf* MassAngleMisTag;
RooAbsPdf* BkgMassAnglePeak;

RooRealVar* nSig;
RooRealVar* nBkgComb;  //07-26
RooRealVar* nMisTagFrac;
RooRealVar* nBkgPeak;

// ##################################
// # Total pdf for B0 --> K*0 mu mu #
// ##################################
RooAbsPdf* TotalPDFRejectPsi;

// ##############################################
// # Vector containing the Gaussian constraints #
// ##############################################
RooArgSet vecConstr;


// #######################
// # Function Definition #
// #######################

bool CheckGoodFit               (RooFitResult* fitResult, TPaveText* paveText = NULL);
RooRealVar* GetVar              (RooAbsPdf* pdf, string varName);
bool GetValueAndErrors          (RooAbsPdf* pdf, string varName, stringstream* myString, double& val, double& errLo, double& errHi);
void SetValueAndErrors          (RooAbsPdf* pdf, string varName, double multi, stringstream* myString, double* val, double* errLo, double* errHi);
void PrintVariables             (RooArgSet* setVar, string type);
void ClearVars                  (RooArgSet* vecVars);
void CloseAllAndQuit            (TApplication* theApp, TFile* NtplFile);

void SetStyle                   ();
string MakeName                 (const RooDataSet* data, unsigned int ID);
void DrawString                 (double Lumi, RooPlot* myFrame = NULL);
void AddGaussConstraint         (RooArgSet* vecConstr, RooAbsPdf* pdf, string varName);
RooAbsPdf* MakeAngWithEffPDF (unsigned int q2BinIndx, RooRealVar* y, RooRealVar* z,RooRealVar* p, unsigned int FitType, bool useEffPDF, RooArgSet* VarsAng);
unsigned int CopyFitResults     (RooAbsPdf* pdf, unsigned int q2BinIndx, vector<vector<string>*>* fitParam);

void MakeDatasets               (B0KstMuMuTreeContent* NTuple, unsigned int FitType);

//###############
// 3-D models   #
//###############

void InstantiateGen3AnglesFit     (RooAbsPdf** TotalPDF,
				bool useEffPDF,
				RooRealVar* y, RooRealVar* z,RooRealVar* p,
				string fitName, unsigned int FitType,
				vector<vector<unsigned int>*>* configParam,
				vector<vector<string>*>* fitParam,
				unsigned int q2BinIndx);
RooFitResult* Make3AnglesFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* y, RooRealVar* z,RooRealVar* p,  unsigned int FitType, RooArgSet* vecConstr, TCanvas* Canv, unsigned int ID);

void Iterative3AnglesFitq2Bins (RooDataSet* dataSet,
                                    bool useEffPDF,
                                    RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                    int specBin,
                                    unsigned int FitType,
                                    vector<double>* q2Bins,
                                    vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
                                    RooArgSet* vecConstr,
                                    unsigned int ID = 0);


// ###########################
// # Function Implementation #
// ###########################
bool CheckGoodFit (RooFitResult* fitResult, TPaveText* paveText)
// ####################################################
// # Covariance matrix quality:                       #
// # -1 : "Unknown, matrix was externally provided"   #
// #  0 : "Not calculated at all"                     #
// #  1 : "Approximation only, not accurate"          #
// #  2 : "Full matrix, but forced positive-definite" #
// #  3 : "Full, accurate covariance matrix"          #
// ####################################################
{
  if (fitResult != NULL)
    {
      if ((fitResult->covQual() == 3) && (fitResult->status() == 0))
	{
	  if (paveText != NULL) paveText->AddText("Fit status: GOOD");
	  return true;
	}
      else
	{
	  if (paveText != NULL) paveText->AddText("Fit status: BAD");
	  return false;
	}
    }


  return false;
}


RooRealVar* GetVar (RooAbsPdf* pdf, string varName)
{
  return (RooRealVar*)(pdf->getVariables()->find(varName.c_str()));
}


bool GetValueAndErrors (RooAbsPdf* pdf, string varName, stringstream* myString, double& val, double& errLo, double& errHi)
{
  if (GetVar(pdf,varName.c_str()) != NULL)
    {
      val   = GetVar(pdf,varName.c_str())->getVal();
      errLo = GetVar(pdf,varName.c_str())->getErrorLo();
      errHi = GetVar(pdf,varName.c_str())->getErrorHi();

      (*myString) << val << "   " << errHi << "   " << errLo << "   ";
      return true;
    }
  else (*myString) << "0.0   0.0   0.0";


  return false;
}


void SetValueAndErrors (RooAbsPdf* pdf, string varName, double multi, stringstream* myString, double* val, double* errLo, double* errHi)
// #############################################################################
// # If the error is an empty string --> setAsymError = -1/+1 and setError = 1 #
// # If both errLo and errHi are 0.0 --> setAsymError = -1/+1 and setError = 1 #
// #############################################################################
{
  string tmpStr;


  if (myString->str().empty() == false)
    {
      tmpStr.clear();
      (*myString) >> tmpStr;
      *val = atof(tmpStr.c_str()) * multi;

      tmpStr.clear();
      (*myString) >> tmpStr;
      if (tmpStr.empty() == true) *errLo = -1.0;
      else *errLo = atof(tmpStr.c_str()) * multi;

      tmpStr.clear();
      (*myString) >> tmpStr;
      if (tmpStr.empty() == true) *errHi = 1.0;
      else *errHi = atof(tmpStr.c_str()) * multi;
    }

  if ((pdf != NULL) && (GetVar(pdf,varName) != NULL))
    {
      pdf->getVariables()->setRealValue(varName.c_str(),*val);
      if ((*errLo == 0.0) && (*errHi == 0.0))
      	{
	  GetVar(pdf,varName)->setAsymError(0.0,0.0);
	  GetVar(pdf,varName)->setError(0.0);
      	}
      else
      	{
      	  GetVar(pdf,varName)->setAsymError(*errLo,*errHi);
	  GetVar(pdf,varName)->setError((*errHi - *errLo) / 2.);
      	}
    }
}



void PrintVariables (RooArgSet* setVar, string type)
{
  RooRealVar* tmpVar;
  int nEleSet = setVar->getSize();


  if (type == "vars")
    {
      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout <<   "@@@ Printing variables @@@" << endl;
      cout <<   "@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      
      if (setVar != NULL)
	{
	  TIterator* it = setVar->createIterator();
	  for (int i = 0; i < nEleSet; i++)
	    {
	      tmpVar = (RooRealVar*)it->Next();
	      cout << "Variable: " << i;
	      cout << "\tname: "   << tmpVar->GetName();
	      cout << "\tvalue: "  << tmpVar->getVal();
	      cout << "\terr: "    << tmpVar->getError();
	      cout << "\tErrLo: "  << tmpVar->getErrorLo();
	      cout << "\tErrHi: "  << tmpVar->getErrorHi() << endl;
	    }
	}
    }
  else if (type == "cons")
    {
      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout <<   "@@@@@@@@@ Printing constraints @@@@@@@@@" << endl;
      cout <<   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

      if (setVar != NULL)
	{
	  TIterator* it = setVar->createIterator();
	  for (int i = 0; i < nEleSet; i++) PrintVariables(((RooAbsPdf*)it->Next())->getVariables(),"vars");
	}
    }
  else
    {
      cout << "[ExtractYield::PrintVariables]\tWrong parameter: " << type << endl;
      exit (EXIT_FAILURE);
    }
}


void ClearVars (RooArgSet* vecVars)
{
  if (vecVars != NULL)
    {
      int nEle = vecVars->getSize();
      
      TIterator* it = vecVars->createIterator();
      for (int i = 0; i < nEle; i++) delete it->Next();
      
      vecVars->removeAll();
    }
}

void CloseAllAndQuit (TApplication* theApp, TFile* NtplFile)
{
  fileFitResults.close();

  if (NtplFile                   != NULL) NtplFile->Close("R");

  gROOT->CloseFiles();
  theApp->Terminate(0);
}


void SetStyle ()
{
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetTextFont(42);

  gStyle->SetOptFit(1112);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.11);
  gStyle->SetPadBottomMargin(0.12);

  gStyle->SetTitleFont(42,"x");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetTitleOffset(1.05,"x");
  gStyle->SetTitleOffset(0.95,"y");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"y");

  gStyle->SetLabelFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");

  TGaxis::SetMaxDigits(3);
  gStyle->SetStatY(0.9);
}


string MakeName (const RooDataSet* data, unsigned int ID)
{
  stringstream myString;


  myString.clear(); myString.str("");
  myString << data->GetName() << "_" << ID;


  return myString.str();
}


void DrawString (double Lumi, RooPlot* myFrame)
{
  stringstream myString;
  double scaleRespect2CMS = 0.75;


  myString.clear(); myString.str("");
  myString << "CMS";
  TLatex* LumiTex1 = new TLatex(0.18,0.9,myString.str().c_str());
  LumiTex1->SetTextFont(61);
  LumiTex1->SetTextSize(0.05);
  LumiTex1->SetTextColor(kBlack);
  LumiTex1->SetNDC(true);
  if (myFrame == NULL) LumiTex1->DrawLatex(0.18,0.9,myString.str().c_str());
  else
    {
      LumiTex1->Paint();
      myFrame->addObject(LumiTex1);
    }


  myString.clear(); myString.str("");
  myString << "#it{Preliminary}";
  TLatex* LumiTex2 = new TLatex(0.265,0.9,myString.str().c_str());
  LumiTex2->SetTextFont(42);
  LumiTex2->SetTextSize(0.05 * scaleRespect2CMS);
  LumiTex2->SetTextColor(kBlack);
  LumiTex2->SetNDC(true);
  if (myFrame == NULL) LumiTex2->DrawLatex(0.265,0.9,myString.str().c_str());
  else
    {
      LumiTex2->Paint();
      myFrame->addObject(LumiTex2);
    }


  myString.clear(); myString.str("");
  myString << Lumi <<  " fb#lower[0.4]{^{#font[122]{\55}1}} (8 TeV)";
  TLatex* LumiTex3 = new TLatex(0.8,0.9,myString.str().c_str());
  LumiTex3->SetTextFont(42);
  LumiTex3->SetTextSize(0.05 * scaleRespect2CMS);
  LumiTex3->SetTextColor(kBlack);
  LumiTex3->SetNDC(true);
  if (myFrame == NULL) LumiTex3->DrawLatex(0.8,0.9,myString.str().c_str());
  else
    {
      LumiTex3->Paint();
      myFrame->addObject(LumiTex3);
    }
}

void AddGaussConstraint (RooArgSet* vecConstr, RooAbsPdf* pdf, string varName)
{
  stringstream myString;


  myString.clear(); myString.str("");
  myString << varName;

  if (GetVar(pdf,myString.str().c_str())->isConstant() == false)
    {
      RooRealVar* varConstr = GetVar(pdf,myString.str().c_str());
      double mean  = GetVar(pdf,myString.str().c_str())->getVal();
      double sigma = (varConstr->getErrorHi() - varConstr->getErrorLo()) / 2.;
      
      myString << "_constr";
      
      RooGaussian* newConstr = new RooGaussian(myString.str().c_str(), myString.str().c_str(), *varConstr, RooConst(mean), RooConst(sigma));
      vecConstr->add(*newConstr);
    }
}

RooAbsPdf* MakeAngWithEffPDF (unsigned int q2BinIndx,RooRealVar* y, RooRealVar* z,RooRealVar* p, unsigned int FitType, bool useEffPDF, RooArgSet* VarsAng)
// ###################
// # y: cos(theta_l) #
// # z: cos(theta_K) #
// # p: phi
// ###################
{
  stringstream myString;
  double a,b,c;
  vector<RooRealVar*> vecParam;
  RooAbsPdf* AnglesPDF = NULL;

  if (
      (FitType == 206)
     ) 
    {
      // #####################################
      // # Make 3D angular*efficiency p.d.f. #
      // # For correctly tagged events       #
      // #####################################
      FlS = new RooRealVar("FlS","F_{L}", 0.0);
      P1S = new RooRealVar("P1S","P_{1}", -3.0, 3.0);
      P2S = new RooRealVar("P2S","P_{2}",  0.0);
      P3S = new RooRealVar("P3S","P_{3}", 0.0);

      P4pS = new RooRealVar("P4pS","P_{4p}", -3.0, 3.0);
      P5pS = new RooRealVar("P5pS","P_{5p}", -3.0, 3.0);
      P6pS = new RooRealVar("P6pS","P_{6p}", -3.0, 3.0);
      P8pS = new RooRealVar("P8pS","P_{8p}", -3.0, 3.0);

      VarsAng->add(*FlS);
      VarsAng->add(*P1S);
      VarsAng->add(*P2S);
      VarsAng->add(*P3S);
      VarsAng->add(*P4pS);
      VarsAng->add(*P5pS);
      VarsAng->add(*P6pS);
      VarsAng->add(*P8pS);
      VarsAng->add(*y);
      VarsAng->add(*z);
      VarsAng->add(*p);     
 
        myString.clear(); myString.str("");
        myString << "(9/(32 *" << Utility->PI << ") * ( 1/2 * (1-" << "FlS" << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * ( 1+ " << y->getPlotLabel() << "*" << y->getPlotLabel() << ") + 2 * " << "FlS" << "* " << z->getPlotLabel() << "*" << z->getPlotLabel() << " *( 1- " << y->getPlotLabel() << "*" << y->getPlotLabel() << ") +";
        myString << "0.5 * P1S" << " * (1-" << "FlS" << ") * (1- " << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * ( 1- " << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * cos (2 *"<< p->getPlotLabel()<< ") + ";
        myString << "P4pS" << "* 2 * cos(" << p->getPlotLabel() << ") * " << z->getPlotLabel() << "*" << y->getPlotLabel() << " * sqrt(" << "FlS" << "* (1-" << "FlS" << ") * ( 1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")) + " ;
        myString << "P5pS" << "* 2 * cos(" << p->getPlotLabel() << ") * " << z->getPlotLabel() << "* sqrt (" << "FlS" << "* (1-" << "FlS" << ") * ( 1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")) - " ;
        myString << "P6pS" << "* 2 * sin(" << p->getPlotLabel() << ") * " << z->getPlotLabel() << "* sqrt (" << "FlS" << "* (1-" << "FlS" << ") * ( 1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")) + " ;
        myString << "P8pS" << "* 2 * sin(" << p->getPlotLabel() << ") * " << z->getPlotLabel() << "*" << y->getPlotLabel() << " * sqrt(" << "FlS" << "* (1-" << "FlS" << ") * ( 1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")) + " ;
        myString << "P2S" << "* 2 * (1-" << "FlS" << ")*" << y->getPlotLabel() << "* (1-" <<  z->getPlotLabel() << "*" << z->getPlotLabel() << ") -" ;
        myString << "P3S" << " * sin( 2 *" << p->getPlotLabel() << ")  * ( 1-" << "FlS" << ") * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")))" ;

      cout << "\n[ExtractYield::MakeAngWithEffPDF]\t@@@ 3D angular p.d.f. @@@" << endl;
      cout << myString.str().c_str() << endl;
    
      AnglesPDF = new RooGenericPdf("AngleS",myString.str().c_str(),RooArgSet(*VarsAng));
    }

  vecParam.clear();
  return AnglesPDF;
}

unsigned int CopyFitResults (RooAbsPdf* pdf, unsigned int q2BinIndx, vector<vector<string>*>* fitParam)
// ####################################################################
// # Only for polynomial coeficients:                                 #
// # If errLo and errHi are 0.0 --> set poly. coefficient as constant #
// ####################################################################
{
  stringstream myString;
  stringstream myCoeff;
  double value, errLo, errHi;
 
  if (GetVar(pdf,"FlS") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(pdf,"FlS",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"FlS")->setConstant(false);
      cout << "Fl = " << FlS->getVal()  << endl;
      }
  if (GetVar(pdf,"P5pS") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("P5pS"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(pdf,"P5pS",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"P5pS")->setConstant(false);
      cout<<"P5p="<<myString.str()<<endl;
    } 
  if (GetVar(pdf,"P1S") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("P1S"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(pdf,"P1S",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"P1S")->setConstant(false);
      cout<<"P1="<<myString.str()<<endl;
     }
 if (GetVar(pdf,"P2S") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("P2S"))->operator[](q2BinIndx);
      SetValueAndErrors(pdf,"P2S",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"P2S")->setConstant(false);
      cout<<"P2="<<myString.str()<<endl;
    }
 if (GetVar(pdf,"P3S") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("P3S"))->operator[](q2BinIndx);
      SetValueAndErrors(pdf,"P3S",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"P3S")->setConstant(false);
      cout<<"P3="<<myString.str()<<endl;
    }
 if (GetVar(pdf,"P4pS") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("P4pS"))->operator[](q2BinIndx);
      SetValueAndErrors(pdf,"P4pS",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"P4pS")->setConstant(false);
      cout<<"P4p="<<myString.str()<<endl;
    }
 if (GetVar(pdf,"P6pS") != NULL)
    { 
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("P6pS"))->operator[](q2BinIndx);
      SetValueAndErrors(pdf,"P6pS",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"P6pS")->setConstant(false);
      cout<<"P6p="<<myString.str()<<endl;
    }
 if (GetVar(pdf,"P8pS") != NULL)
    { 
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("P8pS"))->operator[](q2BinIndx);
      SetValueAndErrors(pdf,"P8pS",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"P8pS")->setConstant(false);
      cout<<"P8p="<<myString.str()<<endl;
    }

  return value;
}


void MakeDatasets (B0KstMuMuTreeContent* NTuple, unsigned int FitType)
{
  stringstream myString;

  // ###########################
  // # Define useful variables #
  // ###########################
   B0mass  = new RooRealVar("B0mass","#font[12]{m}(#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{+}}}}#kern[-0.3]{#pi}#kern[-0.3]{#lower[0.6]{^{#font[122]{\55}}}}#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}})",Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str()),Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str()),"GeV");
  mumuMass           = new RooRealVar("mumuMass","#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}} inv. mass",0.0,6.0,"GeV");
  mumuMassE          = new RooRealVar("mumuMassE","#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}} inv. mass error",0.0,0.5,"GeV");
  ctK     = new RooRealVar("ctK","cos#theta_{K}",-1.0,1.0);
  ctL     = new RooRealVar("ctL","cos#theta_{L}",-1.00,1.0);
  phi     = new RooRealVar("phi","#phi", - Utility->PI ,Utility->PI);
  truthMatchSignal   = new RooRealVar("truthMatchSignal","Truth matching",0.0,1.0,"bool");
  rightFlavorTag     = new RooRealVar("rightFlavorTag","Right flavor tag",0.0,1.0,"bool");

  if ( (FitType == 206) )
    {
      RooArgSet Vars;
      Vars.add(*mumuMass);
      Vars.add(*ctK);
      Vars.add(*ctL);
      Vars.add(*phi);

      SingleCandNTuple           = new RooDataSet("SingleCandNTuple"          ,"SingleCandNTuple"          ,Vars);

      // #############################
      // # Load values from the tree #
      // #############################
      NTuple->ClearNTuple();
      NTuple->SetBranchAddresses(theTree);
      int nEntries = theTree->GetEntries();
      cout << "[ExtractYield::MakeDatasets]\t@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;
      for (int entry = 0; entry < nEntries; entry++)
	{
	  theTree->GetEntry(entry);
	 
	  if (
               (FitType == 206)
	     )
	    {
	      Vars.setRealValue("mumuMass",          NTuple->mumu_mass);
	      Vars.setRealValue("ctK",               NTuple->gen_cos_theta_k);
              Vars.setRealValue("ctL",               NTuple->gen_cos_theta_l);
              Vars.setRealValue("phi",               NTuple->gen_phi);

                         
	      // ########################
	      // # NTuple with all data #
	      // ########################
	      SingleCandNTuple->add(Vars);
            }
	}
      cout << "\n[ExtractYield::MakeDatasets]\t@@@ NTuple with all data @@@" << endl;
      SingleCandNTuple->Print("v");

      }
  // ####################################################
  // # Setting initial values for independent variables #
  // ####################################################
   B0mass->setVal(Utility->B0Mass);
   ctK ->setVal(0.0);
   ctL->setVal(0.0);
   phi->setVal(0.0);
 }



     //#############################
     //# GEN level fitting  PDF #
     //#############################
void InstantiateGen3AnglesFit    (RooAbsPdf** TotalPDF,
				bool useEffPDF,
				RooRealVar* y, RooRealVar* z,RooRealVar* p,
				unsigned int FitType,
				vector<vector<unsigned int>*>* configParam,
			        unsigned int q2BinIndx)
// #########################
// # y: angle cos(theta_l) #
// # z: angle cos(theta_K) #
// # p: angle phi          #
// #########################
{

  // ###################
  // # Local variables #
  // ###################
  stringstream myString;
  RooArgSet* VarsAng = new RooArgSet("VarsAng");

  if ((FitType == 206)) *TotalPDF = MakeAngWithEffPDF(q2BinIndx,y,z,p,FitType,useEffPDF,VarsAng);
  else
    {
      cout << "[ExtractYield::InstantiateGen3AnglesFit]\tIncorrect configuration sequence : useSignal = " << "useSignal" <<  endl;
      exit (EXIT_FAILURE);
    }

}
     //#############################
     //# Reco level fitting  PDF #
     //#############################


RooFitResult* Make3AnglesFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* y, RooRealVar* z,RooRealVar* p,  unsigned int FitType, RooArgSet* vecConstr, TCanvas* Canv, unsigned int ID)
{
  // ###################
  // # Local variables #
  // ###################
  RooFitResult* fitResult = NULL;
  stringstream myString;

  unsigned int nElements   = 0;
  unsigned int it          = 0;
  TString legNames[6]; 
  TLegend*   legY             = NULL;
  TLegend*   legZ             = NULL;
  TLegend*   legP             = NULL;

  if ((FitType == 206))
    {
      // ###################
      // # Make actual fit #
      // ###################
     
      fitResult = (*TotalPDF)->fitTo(*dataSet,Save(true),Minimizer(MINIMIZER),NumCPU(8));
      // ###################################################
      // # Set p.d.f. independent variables to known point #
      // ###################################################
      if (GetVar(*TotalPDF,y->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(y->getPlotLabel(),0.0);
      if (GetVar(*TotalPDF,z->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(z->getPlotLabel(),0.0);
      if (GetVar(*TotalPDF,p->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(p->getPlotLabel(),0.0);
      if (fitResult != NULL) fitResult->Print("v");
      
      // ###########################
      // # Costheta-l plot results #
      // ###########################
      if (PLOT==true ) // &&  CheckGoodFit(fitResult) == true)
      {
      Canv->cd(1);
      RooPlot* myFrameY = y->frame(NBINS);

      dataSet->plotOn(myFrameY, Name(MakeName(dataSet,ID).c_str()));
      if ((FitType == 206))  legNames[nElements++] = "GEN-MC";
      (*TotalPDF)->plotOn(myFrameY, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), Project(RooArgSet(*z,*p)));
      if ((FitType == 206))  legNames[nElements++] = "Gen p.d.f.";
      TPaveText* paveTextY = new TPaveText(0.12,0.82,0.4,0.86,"NDC");
      paveTextY->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
      CheckGoodFit(fitResult,paveTextY);
      paveTextY->SetTextAlign(11);
      paveTextY->SetBorderSize(0.0);
      paveTextY->SetFillStyle(0);
      paveTextY->SetTextSize(0.04);
      paveTextY->Paint();
      myFrameY->addObject(paveTextY);

      legY = new TLegend(0.78, 0.65, 0.97, 0.88, "");
      it = 0;
      for (unsigned int i = 0; i < 2*nElements; i++)
	{
	  TString objName = myFrameY->nameOf(i);
	  if ((objName == "") || (objName.Contains("paramBox") == true) || (objName.Contains("TPave") == true) || ((i > 0) && (objName == myFrameY->nameOf(i-1)))) continue;
	  TObject* obj = myFrameY->findObject(objName.Data());
	  legY->AddEntry(obj,legNames[it++],"PL");
	  legY->SetTextFont(42);
	}     
      legY->SetFillStyle(0);
      legY->SetFillColor(0);
      legY->SetTextSize(0.04);
      legY->SetBorderSize(0);
     
    myFrameY->Draw();
    legY->Draw("same");
     
      // ###########################
      // # Costheta-k plot results #
      // ###########################
      Canv->cd(2);
      RooPlot* myFrameZ = z->frame(NBINS);

      dataSet->plotOn(myFrameZ, Name(MakeName(dataSet,ID).c_str()));
      if ((FitType == 206))  legNames[nElements++] = "GEN-MC";
      (*TotalPDF)->plotOn(myFrameZ, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), Project(RooArgSet(*y,*p)));
      if ((FitType == 206))  legNames[nElements++] = "Gen p.d.f.";
 
      TPaveText* paveTextZ = new TPaveText(0.12,0.82,0.4,0.86,"NDC");
      paveTextZ->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
      paveTextZ->SetTextAlign(11);
      paveTextZ->SetBorderSize(0.0);
      paveTextZ->SetFillStyle(0);
      paveTextZ->SetTextSize(0.04);
      paveTextZ->Paint();
      myFrameZ->addObject(paveTextZ);

      legZ = new TLegend(0.78, 0.65, 0.97, 0.88, "");
      it = 0;
      for (unsigned int i = 0; i < 2*nElements; i++)
	{
	  TString objName = myFrameZ->nameOf(i);
	  if ((objName == "") || (objName.Contains("paramBox") == true) || (objName.Contains("TPave") == true) || ((i > 0) && (objName == myFrameZ->nameOf(i-1)))) continue;
	  TObject* obj = myFrameZ->findObject(objName.Data());
	  legZ->AddEntry(obj,legNames[it++],"PL");
	  legZ->SetTextFont(42);
	}     
      legZ->SetFillStyle(0);
      legZ->SetFillColor(0);
      legZ->SetTextSize(0.04);
      legZ->SetBorderSize(0);
     
    myFrameZ->Draw();
    legZ->Draw("same");


      // ###########################
      // # Costheta-l plot results #
      // ###########################
      Canv->cd(3);
      RooPlot* myFrameP = p->frame(NBINS);

      dataSet->plotOn(myFrameP, Name(MakeName(dataSet,ID).c_str()));
      if ((FitType == 206))  legNames[nElements++] = "GEN-MC";
      (*TotalPDF)->plotOn(myFrameP, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), Project(RooArgSet(*y,*z)));
      if ((FitType == 206))  legNames[nElements++] = "Gen p.d.f.";
      
      TPaveText* paveTextP = new TPaveText(0.12,0.82,0.4,0.86,"NDC");
      paveTextP->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameP->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
      paveTextP->SetTextAlign(11);
      paveTextP->SetBorderSize(0.0);
      paveTextP->SetFillStyle(0);
      paveTextP->SetTextSize(0.04);
      paveTextP->Paint();
      myFrameP->addObject(paveTextP);

      legP = new TLegend(0.78, 0.65, 0.97, 0.88, "");
      it = 0;
      for (unsigned int i = 0; i < 2*nElements; i++)
	{
	  TString objName = myFrameP->nameOf(i);
	  if ((objName == "") || (objName.Contains("paramBox") == true) || (objName.Contains("TPave") == true) || ((i > 0) && (objName == myFrameY->nameOf(i-1)))) continue;
	  TObject* obj = myFrameP->findObject(objName.Data());
	  legP->AddEntry(obj,legNames[it++],"PL");
	  legP->SetTextFont(42);
	}     
      legP->SetFillStyle(0);
      legP->SetFillColor(0);
      legP->SetTextSize(0.04);
      legP->SetBorderSize(0);
     
    myFrameP->Draw();
    legP->Draw("same");
    }

    Canv->Modified();
    Canv->Update();

    if (SAVEPLOT == true)
    {
     myString.clear(); myString.str("");
     myString << (*TotalPDF)->getPlotLabel() << "_Canv" << ID << ".pdf";
     Canv->Print(myString.str().c_str());

     myString.clear(); myString.str("");
     myString << (*TotalPDF)->getPlotLabel() << "_Canv" << ID << ".root";
     Canv->Print(myString.str().c_str());
    }

    } 
    return fitResult;
}


void Iterative3AnglesFitq2Bins (RooDataSet* dataSet,
				    bool useEffPDF,
				    RooRealVar* y, RooRealVar* z,RooRealVar* p,
				    int specBin,
				    unsigned int FitType,
				    vector<double>* q2Bins,
				    vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
				    RooArgSet* vecConstr,
                                    unsigned int ID)
{
  // ###################
  // # Local variables #
  // ###################
  stringstream myString;   
 
  RooAbsPdf*  TotalPDFq2Bins[q2Bins->size()-1];
  RooFitResult* fitResult;

  RooDataSet* dataSet_q2Bins[q2Bins->size()-1];
  for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
    {
      myString.clear(); myString.str("");
      myString << "(mumuMass*mumuMass) > " << q2Bins->operator[](i) << " && (mumuMass*mumuMass) <= " << q2Bins->operator[](i+1);
      cout << "\n[ExtractYield::IterativeAnglesFitq2Bins]\tCut string: " << myString.str() << endl;
      dataSet_q2Bins[i] = (RooDataSet*)dataSet->reduce(myString.str().c_str());
      cout << "[ExtractYield::IterativeAnglesFitq2Bins]\tNumber of events : " << dataSet_q2Bins[i]->sumEntries() << endl;

      TCanvas*    Gen[q2Bins->size()-1];
      Gen[i] = new TCanvas("Gen","Gen",10,10,700,500);
      Gen[i]->Divide(2,2);

      myString.clear(); myString.str("");
      myString << "GenTotalPDFq2Bin_" << i;
      InstantiateGen3AnglesFit(&TotalPDFq2Bins[i],useEffPDF,y,z,p,FitType,configParam,i);
       // #####################
      // # Initialize p.d.f. #
      // #####################
       CopyFitResults(TotalPDFq2Bins[i],i,fitParam);  

      // #####################
      // # Apply constraints #
      // #####################
      ClearVars(vecConstr);


      // ###################
      // # Perform the fit #
      // ##################
      fitResult = Make3AnglesFit(dataSet_q2Bins[i],&TotalPDFq2Bins[i],y,z,p,FitType,vecConstr,Gen[i],ID);
      if (CheckGoodFit(fitResult) == true) cout << "\n[ExtractYield::Iterative3AnglesFitq2Bins]\t@@@ Fit converged ! @@@" << endl;
      else                                 cout << "\n[ExtractYield::Iterative3AnglesFitq2Bins]\t@@@ Fit didn't converge ! @@@" << endl;

    }
}




// ==================
// ===> 4D MODEL <===
// ==================


int main(int argc, char** argv)
{
  if (argc >= 4)
    {
      // ##################
      // # Main variables #
      // ##################
      stringstream myString;
      
      string fileName           = "";
      string correct4Efficiency = "";
      string tmpFileName        = "";

      int specBin               = -1;
      unsigned int fileIndx     = 0;
      unsigned int FitType      = atoi(argv[1]);

      bool useEffPDF            = false;

      TFile* NtplFile           = NULL;

      if ((FitType == 206)
          && (argc >= 4))
	  
	{
	  ParameterFILE = Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str();


 	  // ###################
	  // # Read parameters #
 	  // ###################
	  Utility = new Utils(false);
	  Utility->ReadAllBins(ParameterFILE,&q2Bins);


	  // #################################
	  // # Check that FitType is correct #
	  // #################################
	  if (
	      (FitType == 206)
             )
	    {
	      fileName           = argv[2];
	      correct4Efficiency = argv[3];
	    }
          // ###################################
          // # Check that FitOption is correct #
          // ###################################    
          if ((correct4Efficiency != "noEffCorr") && (correct4Efficiency != "yesEffCorr"))
            {
              cout << "[ExtractYield::main]\tIncorrect option parameter " << correct4Efficiency << endl;
              exit (EXIT_FAILURE);
            }
  
          // #################################
          // # Read the q^2 bin and the rest #
          // #################################
          if (argc >= 5) specBin = atoi(argv[4]);
          if ((correct4Efficiency == "yesEffCorr")) useEffPDF = true;
          else if ((correct4Efficiency == "noEffCorr") || (correct4Efficiency == "yesEffCorr"))
            {
              if (argc >= 7)
                {
                 fileIndx = atoi(argv[6]);
                 if (argc == 8) tmpFileName = argv[7];
                }
            }

          cout << "\n[ExtractYield::main]\t@@@ Input variables from command line @@@" << endl;
          cout << "- input/outputFile.root = " << fileName.c_str() << endl;
          cout << "- correct4Efficiency = "    << correct4Efficiency << endl;
          cout << "- tmpFileName = "           << tmpFileName.c_str() << endl;
          cout << "- specBin = "               << specBin << endl;
          cout << "- fileIndx = "              << fileIndx << endl;
          cout << "- FitType = "               << FitType << endl;
          cout << "- useEffPDF = "             << useEffPDF << endl;
          cout << "- ParameterFILE = "         << ParameterFILE.c_str() << endl;

          cout << "\n[ExtractYield::main]\t@@@ Internal settings @@@" << endl;
          cout << "NBINS = "         << NBINS << endl;
          cout << "MULTYIELD = "     << MULTYIELD << endl;
          cout << "NCOEFFPOLYBKG = " << NCOEFFPOLYBKG << endl;

          cout << "\nPARAMETERFILEIN = " << PARAMETERFILEIN << endl;
          cout << "PARAMETERFILEOUT = "  << PARAMETERFILEOUT << endl;
 
          // ##########################
          // # Set histo layout style #
          // ##########################
          SetStyle();

          // ###################
          // # Read parameters #
          // ###################       
          Utility->ReadGenericParam(ParameterFILE);
          Utility->ReadFitStartingValues(ParameterFILE,&fitParam,&configParam,Utility->ParFileBlockN("fitValBins"));

          cout << "Use MINOS: "                          << Utility->GetGenericParam("UseMINOS").c_str() << " (0 = false; 1 = true)" << endl;

          // ###############################################################################################
          // # Read other parameters : this also allow to understand if the parameter file is well written #
          // ###############################################################################################
          LUMI = Utility->ReadLumi(ParameterFILE);
          if (Utility->WhatIsThis(ParameterFILE) == 0) cout << "\n[ExtractYield::main]\t@@@ I recognize that this is a DATA file @@@" << endl;
          else                                         cout << "\n[ExtractYield::main]\t@@@ I recognize that this is a Monte Carlo file @@@" << endl;


         // ###################
          // # Select fit type #
          // ###################
          if (
              (FitType == 206)
             )
            {
              NtplFile = new TFile(fileName.c_str(),"READ");
              theTree  = (TTree*) NtplFile->Get("tree");
              NTuple   = new B0KstMuMuTreeContent();
              NTuple->Init();

	      // #################
	      // # Make datasets #
	      // #################
	      cout << "\n[ExtractYield::main]\t@@@ Making datasets @@@" << endl;
	      MakeDatasets(NTuple,FitType);
              

	      // ##############################
	      // # Select the proper fit type #
	      // ##############################
	     
              // #############################
              // # 3D-fit P5P-Fl per q^2 bin #
	      // #############################
		  cout << "\n[ExtractYield::main]\t@@@ Now fit invariant mass, cos(theta_K) and cos(theta_l) per mumu q^2 bins @@@" << endl;

                  if ((FitType == 206))                      Iterative3AnglesFitq2Bins(SingleCandNTuple,
		  								       useEffPDF,
		  								       ctL,
		  								       ctK,
                                                                                       phi,
		  								       specBin,
		  								       FitType,
		  								       &q2Bins,
		  								       &configParam,&fitParam,
		  								       &vecConstr,
                                                                                       fileIndx);
		
	    }
            fileFitResults.close();           

	}
      else
	{
	  cout << "Wrong parameter: " << endl;
	  cout << "./ExtractYield [FitType] [input/output[if toy-MC]File.root] [noEffCorr yesEffCorr]" << endl;
	  cout << "               [q^2 bin to fit (0 - ...)]" << endl;

	  cout << "\n --> noEffCorr     = no eff. correction" << endl;
	  cout << " --> yesEffCorr    = use eff. correction" << endl;
	 
	  return EXIT_FAILURE;
	}
    }
  else
    {
          cout << "Parameter missing: " << endl;
	  cout << "./ExtractYield [FitType] [input/output[if toy-MC]File.root] [noEffCorr yesEffCorr]" << endl;
	  cout << "               [q^2 bin to fit (0 - ...)]" << endl;

	  cout << "\n --> noEffCorr     = no eff. correction" << endl;
	  cout << " --> yesEffCorr    = use eff. correction" << endl;
      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Signal @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << " FitType = 206: 3D  (B0Mass, cos(theta_K), cos(theta_l)) per q^2 bin on Gen-level" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      
      return EXIT_FAILURE;
    }

  return 0;
}
