// ################################################################################
// # Program to perform the full angular analysis of the decay B0 --> K*0 mu+ mu- #
// ################################################################################

#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <Riostream.h>
#include <TH3.h>
#include <TEnv.h>
#include <TSystem.h>
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
#include <RooEffProd.h>
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooHistFunc.h>
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
#include "RooBernsteinEffi.h"

#include <TBenchmark.h>
#include <time.h>
#include <ctime>
#include <fstream>
#include <cmath>
#include <cstdlib>
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
using std::fstream;
using std::ifstream;
using std::make_pair;
using namespace RooFit;


// ####################
// # Global constants #
// ####################
#define NBINS           20
#define MULTYIELD       1. // Multiplication factor to the number of entry in toy-MC
#define NCOEFFPOLYBKG   5  // Maximum number of coefficients (= degree) of the polynomial describing the combinatorial bkg

#define _USE_MATH_DEFINES

// ##########################################
// # Internal flags to control the workflow #
// ##########################################
#define EffV          false    //"Projection"  // choose efficiency function (Projection/Fitting) 
#define PLOT          true  //[true= plot the results]
#define SAVEPLOT      true   //2015-12-17
#define MINIMIZER     "Minuit" // Minimizer type for 3D MODEL actual fit ["Minuit"; "Minuit2"]

#define doTransf       false    // limits
#define ROOMINIMIZER   true     // if true: RooMinimizer, if false: Roofit
// ##################
// # External files #
// ##################
#define PARAMETERFILEIN  "/python/ParameterFile.txt"
#define PARAMETERFILEOUT "ParameterFileOut.txt"


// ############################################
// # Global variables from configuration file #
// ############################################
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
RooRealVar* phi;
RooRealVar* tagB0;
RooRealVar* genSignal;


// #################################
// # Variables and pdf for the fit #
// #################################

// ##################
// # Signal B0 mass #
// ##################

// #################
// # Signal angles #
// #################
RooRealVar* FlS;

RooRealVar* P4pS;
RooRealVar* P5pS;
RooRealVar* P6pS;
RooRealVar* P8pS;
RooRealVar* P1S;
RooRealVar* P2S;
RooRealVar* P3S;              

RooRealVar* S3S;
RooRealVar* S4S;
RooRealVar* S5S;
RooRealVar* AFBS;
RooRealVar* S7S;
RooRealVar* S8S;
RooRealVar* S9S;

RooAbsPdf*  AngleGood;
RooEffProd*  AngleS;

// ####################
// # Total signal pdf #
// ####################
RooAbsPdf* Signal;
RooAbsPdf* SignalT;

// ####################################
// # Combinatorial background B0 mass #
// ####################################

// #########################
// # Mistag signal B0 mass #
// #########################
RooRealVar* fracMisTag;
RooAbsPdf*  AngleMisTag;
// ##############################
// # Peaking background B0 mass #
// ##############################

// #####################
// # Background angles #
// #####################

// ########################
// # Total background pdf #
// ########################
RooRealVar* nMisTagFrac;

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

void SetStyle                   ();
string MakeName                 (const RooDataSet* data, int ID);
void DrawString                 (double Lumi, RooPlot* myFrame = NULL);
void AddGaussConstraint         (RooArgSet* vecConstr, RooAbsPdf* pdf, string varName);

double PrintFitResults          (RooAbsPdf* pdf, RooFitResult* fitResult);

RooAbsPdf* MakeAngWithEffPDF    (unsigned int q2BinIndx, RooRealVar* y, RooRealVar* z,RooRealVar* p, unsigned int FitType, bool useEffPDF, RooArgSet* VarsAng, unsigned int efftest);
unsigned int CopyFitResults     (RooAbsPdf* pdf, unsigned int q2BinIndx, int nscan, vector<vector<string>*>* fitParam,unsigned int countMisTag = 0,unsigned int countGoodTag =0);

void MakeDatasets               (B0KstMuMuTreeContent* NTuple, unsigned int FitType, int SampleType);

//###############
// 3-D models   #
//###############

void InstantiateGen3AnglesFit     (RooAbsPdf** TotalPDF,
                                   bool useEffPDF,
                                   RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                   string fitName, unsigned int FitType,
                                   vector<vector<string>*>* fitParam,
                                   unsigned int q2BinIndx, int efftest);

void InstantiateReco3AnglesFit (RooAbsPdf** TotalPDF,
                                bool useEffPDF,
                                RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                unsigned int FitType,
                                vector<vector<unsigned int>*>* configParam,
                                unsigned int q2BinIndx, int efftest);

RooFitResult* Make3AnglesFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* y, RooRealVar* z,RooRealVar* p,  unsigned int FitType, RooArgSet* vecConstr, TCanvas* Canv, int specBin);

void Iterative3AnglesFitq2Bins (RooDataSet* dataSet,
                                bool useEffPDF,
                                RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                int specBin,
                                unsigned int FitType,
                                vector<double>* q2Bins,
                                vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
                                RooArgSet* vecConstr, int efftest, int nscan);

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
            cout << "status from fitresult = " <<  fitResult->status() << endl;
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


string MakeName (const RooDataSet* data,  int ID)
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

double PrintFitResults (RooAbsPdf* pdf, RooFitResult* fitResult)
{
    double varVal, varValELo, varValEHi;

    cout << "[ExtractYield::PrintFitResults]\t Print Fitting Results" << endl;
    if (fitResult != NULL)
    {
        if (GetVar(pdf,"FlS") != NULL)
        {
            Utility->Transformer("FlS",doTransf,varVal,varValELo,varValEHi,fitResult,GetVar(pdf,"FlS"));
            cout << " # FlS +/- err" << endl;
            cout << varVal << " +/- " << (varValEHi - varValELo) / 2. << " (" << varValEHi << "/" << varValELo << ")" << endl;
        }
        if (GetVar(pdf,"P1S") != NULL)
        {
            cout << " # P1S +/- err" << endl;
            cout << GetVar(pdf,"P1S")->getVal() << " +/- " << GetVar(pdf,"P1S")->getError() <<  " (" << GetVar(pdf,"P1S")->getErrorHi() << "/" << GetVar(pdf,"P1S")->getErrorLo() << ")" << endl;
        }
        if (GetVar(pdf,"P2S") != NULL)
        {
            cout << " # P2S +/- err" << endl;
            cout << GetVar(pdf,"P2S")->getVal() << " +/- " << GetVar(pdf,"P2S")->getError() <<  " (" << GetVar(pdf,"P2S")->getErrorHi() << "/" << GetVar(pdf,"P2S")->getErrorLo() << ")" << endl;
        }
        if (GetVar(pdf,"P3S") != NULL)
        {
            cout << " # P3S +/- err" << endl;
            cout << GetVar(pdf,"P3S")->getVal() << " +/- " << GetVar(pdf,"P3S")->getError() <<  " (" << GetVar(pdf,"P3S")->getErrorHi() << "/" << GetVar(pdf,"P3S")->getErrorLo() << ")" << endl;
        }
        if (GetVar(pdf,"P4pS") != NULL)
        {
            cout << " # P4pS +/- err" << endl;
            cout << GetVar(pdf,"P4pS")->getVal() << " +/- " << GetVar(pdf,"P4pS")->getError() <<  " (" << GetVar(pdf,"P4pS")->getErrorHi() << "/" << GetVar(pdf,"P4pS")->getErrorLo() << ")" << endl;
        }
        if (GetVar(pdf,"P5pS") != NULL)
        {
            cout << " # P5pS +/- err" << endl;
            cout << GetVar(pdf,"P5pS")->getVal() << " +/- " << GetVar(pdf,"P5pS")->getError() <<  " (" << GetVar(pdf,"P5pS")->getErrorHi() << "/" << GetVar(pdf,"P5pS")->getErrorLo() << ")" << endl;
        }
        if (GetVar(pdf,"P6pS") != NULL)
        {
            cout << " # P6pS +/- err" << endl;
            cout << GetVar(pdf,"P6pS")->getVal() << " +/- " << GetVar(pdf,"P6pS")->getError() <<  " (" << GetVar(pdf,"P6pS")->getErrorHi() << "/" << GetVar(pdf,"P6pS")->getErrorLo() << ")" << endl;
        }
        if (GetVar(pdf,"P8pS") != NULL)
        {
            cout << " # P8pS +/- err" << endl;
            cout << GetVar(pdf,"P8pS")->getVal() << " +/- " << GetVar(pdf,"P8pS")->getError() <<  " (" << GetVar(pdf,"P8pS")->getErrorHi() << "/" << GetVar(pdf,"P8pS")->getErrorLo() << ")" << endl;
        }

        if (GetVar(pdf,"AFBS") != NULL)
        {
            Utility->Transformer("AFBS", doTransf, varVal,varValELo,varValEHi,fitResult,GetVar(pdf,"FlS"),GetVar(pdf,"AFBS"));
            cout << " # AFBS +/- err " <<endl;
            cout << varVal << " +/- " << (varValEHi - varValELo) / 2. << " (" << varValEHi << "/" << varValELo << ")" << endl;
        }
        if (GetVar(pdf,"S3S") != NULL)
        {
            cout << " # S3S +/- err" << endl;
            cout << GetVar(pdf,"S3S")->getVal() << " +/- " << GetVar(pdf,"S3S")->getError() <<  " (" << GetVar(pdf,"S3S")->getErrorHi() << "/" << GetVar(pdf,"S3S")->getErrorLo() << ")" << endl;
        }
        if (GetVar(pdf,"S4S") != NULL)
        {
            cout << " # S4S +/- err" << endl;
            cout << GetVar(pdf,"S4S")->getVal() << " +/- " << GetVar(pdf,"S4S")->getError() <<  " (" << GetVar(pdf,"S4S")->getErrorHi() << "/" << GetVar(pdf,"S4S")->getErrorLo() << ")" << endl;
        }
        if (GetVar(pdf,"S5S") != NULL)
        {
            cout << " # S5S +/- err" << endl;
            cout << GetVar(pdf,"S5S")->getVal() << " +/- " << GetVar(pdf,"S5S")->getError() <<  " (" << GetVar(pdf,"S5S")->getErrorHi() << "/" << GetVar(pdf,"S5S")->getErrorLo() << ")" << endl;
        }
        if (GetVar(pdf,"S7S") != NULL)
        {
            cout << " # S7S +/- err" << endl;
            cout << GetVar(pdf,"S7S")->getVal() << " +/- " << GetVar(pdf,"S7S")->getError() <<  " (" << GetVar(pdf,"S7S")->getErrorHi() << "/" << GetVar(pdf,"S7S")->getErrorLo() << ")" << endl;
        }
        if (GetVar(pdf,"S8S") != NULL)
        {
            cout << " # S8S +/- err" << endl;
            cout << GetVar(pdf,"S8S")->getVal() << " +/- " << GetVar(pdf,"S8S")->getError() <<  " (" << GetVar(pdf,"S8S")->getErrorHi() << "/" << GetVar(pdf,"S8S")->getErrorLo() << ")" << endl;
        }
        if (GetVar(pdf,"S9S") != NULL)
        {
            cout << " # S9S +/- err" << endl;
            cout << GetVar(pdf,"S9S")->getVal() << " +/- " << GetVar(pdf,"S9S")->getError() <<  " (" << GetVar(pdf,"S9S")->getErrorHi() << "/" << GetVar(pdf,"S9S")->getErrorLo() << ")" << endl;
        }
        if (GetVar(pdf,"S3S") != NULL)
        {
            cout << " # S3S +/- err" << endl;
            cout << GetVar(pdf,"S3S")->getVal() << " +/- " << GetVar(pdf,"S3S")->getError() <<  " (" << GetVar(pdf,"S3S")->getErrorHi() << "/" << GetVar(pdf,"S3S")->getErrorLo() << ")" << endl;
        }
        if (GetVar(pdf,"fracMisTag") != NULL)
        {
            cout << " # fracMisTag +/- err" << endl;
            cout << GetVar(pdf,"fracMisTag")->getVal() << " +/- " << GetVar(pdf,"fracMisTag")->getError() <<  " (" << GetVar(pdf,"fracMisTag")->getErrorHi() << "/" << GetVar(pdf,"fracMisTag")->getErrorLo() << ")" << endl;
        }
    }
    return 1;
}

RooAbsPdf* MakeAngWithEffPDF (unsigned int q2BinIndx, RooRealVar* y, RooRealVar* z,RooRealVar* p, unsigned int FitType, bool useEffPDF, RooArgSet* VarsAng, unsigned int efftest)
// ###################
// # y: cos(theta_l) #
// # z: cos(theta_K) #
// # p: phi
// ###################
{
    stringstream myString;
    RooAbsPdf* AnglesPDF = NULL;

    if ((FitType == 206) || (FitType == 106) || (FitType == 106 * 10))
    {
        // #####################################
        // # Make 3D angular*efficiency p.d.f. #
        // # For correctly tagged events       #
        // #####################################
        cout << "[ExtractYield::MakeAngWithEffPDF]\t@@@CP-averaged observables Pi(') @@@" << endl;
        FlS = new RooRealVar("FlS","F_{L}", 0.5, 0.0, 1.0);
        P1S = new RooRealVar("P1S","P_{1}", 0.0, -1.0, 1.0);
        P2S = new RooRealVar("P2S","P_{2}", 0.0);
        P3S = new RooRealVar("P3S","P_{3}", 0.0);

        P4pS = new RooRealVar("P4pS","P_{4p}", 0.0);
        P5pS = new RooRealVar("P5pS","P_{5p}", 0.0);
        P6pS = new RooRealVar("P6pS","P_{6p}", 0.0);
        P8pS = new RooRealVar("P8pS","P_{8p}", 0.0);

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

      if ((FitType == 106) || (FitType == 206) )
        {
          myString.clear(); myString.str("");
          if (atoi(Utility->GetGenericParam("UseSPwave").c_str()) == false) {
              // #####################
              // # P-wave decay rate #
              // #####################
              myString << "(9/(32 *" << Utility->PI << ") * ( 1/2 * (1-" << "FlS" << ") * (1-" << z->getPlotLabel()
                       << "*" << z->getPlotLabel() << ") * ( 1+ " << y->getPlotLabel() << "*" << y->getPlotLabel()
                       << ") + 2 * " << "FlS" << "* " << z->getPlotLabel() << "*" << z->getPlotLabel() << " *( 1- "
                       << y->getPlotLabel() << "*" << y->getPlotLabel() << ") +";
              myString << "0.5 * P1S" << " * (1-" << "FlS" << ") * (1- " << z->getPlotLabel() << "*"
                       << z->getPlotLabel() << ") * ( 1- " << y->getPlotLabel() << "*" << y->getPlotLabel()
                       << ") * cos (2 *" << p->getPlotLabel() << ") + ";
              myString << "P4pS" << "* 4 * cos(" << p->getPlotLabel() << ") * " << z->getPlotLabel() << "*"
                       << y->getPlotLabel() << " * sqrt(" << "FlS" << "* (1-" << "FlS" << ") * ( 1-"
                       << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * (1-" << z->getPlotLabel() << "*"
                       << z->getPlotLabel() << ")) + ";
              myString << "P5pS" << "* 2 * cos(" << p->getPlotLabel() << ") * " << z->getPlotLabel() << "* sqrt ("
                       << "FlS" << "* (1-" << "FlS" << ") * ( 1-" << y->getPlotLabel() << "*" << y->getPlotLabel()
                       << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")) + ";
              myString << "P6pS" << "* 2 * sin(" << p->getPlotLabel() << ") * " << z->getPlotLabel() << "* sqrt ("
                       << "FlS" << "* (1-" << "FlS" << ") * ( 1-" << y->getPlotLabel() << "*" << y->getPlotLabel()
                       << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")) + ";
              myString << "P8pS" << "* 4 * sin(" << p->getPlotLabel() << ") * " << z->getPlotLabel() << "*"
                       << y->getPlotLabel() << " * sqrt(" << "FlS" << "* (1-" << "FlS" << ") * ( 1-"
                       << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * (1-" << z->getPlotLabel() << "*"
                       << z->getPlotLabel() << ")) + ";
              myString << "P2S" << "* 2 * (1-" << "FlS" << ")*" << y->getPlotLabel() << "* (1-" << z->getPlotLabel()
                       << "*" << z->getPlotLabel() << ") -";
              myString << "P3S" << " * sin( 2 *" << p->getPlotLabel() << ")  * ( 1-" << "FlS" << ") * (1-"
                       << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * (1-" << z->getPlotLabel() << "*"
                       << z->getPlotLabel() << ")))";
          }
            cout << "\n[ExtractYield::MakeAngWithEffPDF]\t@@@ 3D angular p.d.f. Pi & Pi' @@@" << endl;
            cout << myString.str().c_str() << endl;
            if (useEffPDF == true){
                RooGenericPdf* _AnglesPDF = new RooGenericPdf("_AnglesPDF",myString.str().c_str(),RooArgSet(*VarsAng));
                // #############################
                // # Make 3D efficiency p.d.f. #
                // #############################
                RooAbsPdf* EffPdf_R = NULL;
                if (EffV == true){
                   RooAbsPdf* EffPdf_R = Utility->ReadRTEffPDF(q2BinIndx, efftest);
                   cout << "\n[ExtractYield::MakeAngWithEffPDF]\t@@@ The efficiency function from Alessio  @@@" << endl;
                  }
                else {
                   EffPdf_R =(RooAbsPdf*) Utility->GetFitReff(q2BinIndx);
              //   RooBernsteinEffi* EffPdf_R = Utility->GetFitReff(q2BinIndx);
                   cout << "\n[ExtractYield::MakeAngWithEffPDF]\t@@@ The efficiency function from Paolo @@@" << endl;
                  }
                AnglesPDF  = new RooEffProd("AngleS","Signal * Efficiency", *_AnglesPDF, *EffPdf_R); 
                cout << "\n[ExtractYield::MakeAngWithEffPDF]\t@@@ 3D angular*efficiency p.d.f. (P-wave) @@@" << endl;
                cout << myString.str().c_str() << endl;
             } 
           else 
           AnglesPDF = new RooGenericPdf("AngleS",myString.str().c_str(),RooArgSet(*VarsAng));
        }
      else if ( (FitType == 106 *10))
        { 
           myString.clear(); myString.str("");
           if (atoi(Utility->GetGenericParam("UseSPwave").c_str()) == false)
            {
              myString << "(9/(32 *" << Utility->PI << ") * ( 1/2 * (1-" << "FlS" << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * ( 1+ " << y->getPlotLabel() << "*" << y->getPlotLabel() << ") + 2 * " << "FlS" << "* " << z->getPlotLabel() << "*" << z->getPlotLabel() << " *( 1- " << y->getPlotLabel() << "*" << y->getPlotLabel() << ") +";
              myString << "0.5 * P1S" << " * (1-" << "FlS" << ") * (1- " << z->getPlotLabel() << "*" << z->getPlotLabel() << ") * ( 1- " << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * cos (2 *"<< p->getPlotLabel()<< ") + ";
              myString << "P4pS" << "* 4 * cos(" << p->getPlotLabel() << ") * " << z->getPlotLabel() << "*" << y->getPlotLabel() << " * sqrt(" << "FlS" << "* (1-" << "FlS" << ") * ( 1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")) + " ;
              myString << "P5pS" << "* 2 * cos(" << p->getPlotLabel() << ") * " << z->getPlotLabel() << "* sqrt (" << "FlS" << "* (1-" << "FlS" << ") * ( 1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")) - " ;
              myString << "P6pS" << "* 2 * sin(" << p->getPlotLabel() << ") * " << z->getPlotLabel() << "* sqrt (" << "FlS" << "* (1-" << "FlS" << ") * ( 1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")) - " ;
              myString << "P8pS" << "* 4 * sin(" << p->getPlotLabel() << ") * " << z->getPlotLabel() << "*" << y->getPlotLabel() << " * sqrt(" << "FlS" << "* (1-" << "FlS" << ") * ( 1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")) + " ;
              myString << "P2S" << "* 2 * (1-" << "FlS" << ")*" << y->getPlotLabel() << "* (1-" <<  z->getPlotLabel() << "*" << z->getPlotLabel() << ") + " ;
              myString << "P3S" << " * sin( 2 *" << p->getPlotLabel() << ")  * ( 1-" << "FlS" << ") * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")))" ;

             }
          cout << "\n[ExtractYield::MakeAngWithEffPDF]\t@@@ 3D angular p.d.f. Pi & Pi' @@@" << endl;
          cout << myString.str().c_str() << endl;
          if (useEffPDF == true)
            {
               RooGenericPdf* _AnglesPDF = new RooGenericPdf("_AnglesPDF",myString.str().c_str(),RooArgSet(*VarsAng));
               // #############################
               // # Make 3D efficiency p.d.f. #
               // #############################
               RooAbsPdf* EffPdf_W = Utility->ReadWTEffPDF(q2BinIndx, efftest);
               cout << "\n[ExtractYield::MakeAngWithEffPDF]\t@@@ The efficiency function from Alessio  @@@" << endl;
               AnglesPDF  = new RooEffProd("AngleS","Signal * Efficiency", *_AnglesPDF, *EffPdf_W);      
               cout << "\n[ExtractYield::MakeAngWithEffPDF]\t@@@ 3D angular*efficiency p.d.f. (P-wave) @@@" << endl;
               cout << myString.str().c_str() << endl;
            }
           else
           AnglesPDF = new RooGenericPdf("AngleS",myString.str().c_str(),RooArgSet(*VarsAng));
        }
  }
  return AnglesPDF;
}

unsigned int CopyFitResults (RooAbsPdf* pdf, unsigned int q2BinIndx, int nscan, vector<vector<string>*>* fitParam,unsigned int countMisTag, unsigned int countGoodTag)
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

      RooRealVar tmpVar("tmpVar","tmpVar",0.0);
      tmpVar.setVal(atof(fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[](q2BinIndx).c_str()));
     
      Utility->AntiTransformer("FlS",doTransf,value,errLo,errHi,&tmpVar);
      myString.clear(); myString.str("");
      myString << value << "   " << errLo << "   " << errHi;
      SetValueAndErrors(pdf,"FlS",1.0,&myString,&value,&errLo,&errHi);
      cout << "Fl = " << FlS->getVal()  << endl;
      }
  if (GetVar(pdf,"P5pS") != NULL)
    {
     if (nscan < 0)
       { 
         myString.clear(); myString.str("");
         myString << fitParam->operator[](Utility->GetFitParamIndx("P5pS"))->operator[](q2BinIndx).c_str();
         SetValueAndErrors(pdf,"P5pS",1.0,&myString,&value,&errLo,&errHi);
         GetVar(pdf,"P5pS")->setConstant(false);
       }
     else
       {
         myString.clear(); myString.str("");
         myString << RooRandom::uniform();
         SetValueAndErrors(pdf,"P5pS",1.0,&myString,&value,&errLo,&errHi);
         GetVar(pdf,"P5pS")->setConstant(false);
       }
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

 if (GetVar(pdf,"S3S") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("S3S"))->operator[](q2BinIndx);
      SetValueAndErrors(pdf,"S3S",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"S3S")->setConstant(false);
      cout<<"S3="<<myString.str()<<endl;
    }
 if (GetVar(pdf,"S4S") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("S4S"))->operator[](q2BinIndx);
      SetValueAndErrors(pdf,"S4S",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"S4S")->setConstant(false);
      cout<<"S4="<<myString.str()<<endl;
    }
 if (GetVar(pdf,"S5S") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("S5S"))->operator[](q2BinIndx);
      SetValueAndErrors(pdf,"S5S",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"S5S")->setConstant(false);
      cout<<"S5="<<myString.str()<<endl;
    }
 if (GetVar(pdf,"AFBS") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("AFBS"))->operator[](q2BinIndx);
      SetValueAndErrors(pdf,"AFBS",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"AFBS")->setConstant(false);

      RooRealVar tmpVar1("tmpVar1","tmpVar1",0.0);
      tmpVar1.setVal(atof(fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[](q2BinIndx).c_str()));

      RooRealVar tmpVar2("tmpVar2","tmpVar2",0.0);
      tmpVar2.setVal(atof(fitParam->operator[](Utility->GetFitParamIndx("AFBS"))->operator[](q2BinIndx).c_str()));

      Utility->AntiTransformer("AFBS",doTransf,value,errLo,errHi,&tmpVar1,&tmpVar2);
      myString.clear(); myString.str("");
      myString << value << "   " << errLo << "   " << errHi;
      SetValueAndErrors(pdf,"AFBS",1.0,&myString,&value,&errLo,&errHi);
      cout<<"AFB="<<myString.str()<<endl;
    }
 if (GetVar(pdf,"S7S") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("S7S"))->operator[](q2BinIndx);
      SetValueAndErrors(pdf,"S7S",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"S7S")->setConstant(false);
      cout<<"S7="<<myString.str()<<endl;
    }
 if (GetVar(pdf,"S8S") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("S8S"))->operator[](q2BinIndx);
      SetValueAndErrors(pdf,"S8S",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"S8S")->setConstant(false);
      cout<<"S8="<<myString.str()<<endl;
    }
 if (GetVar(pdf,"S9S") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("S9S"))->operator[](q2BinIndx);
      SetValueAndErrors(pdf,"S9S",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"S9S")->setConstant(false);
      cout<<"S9="<<myString.str()<<endl;
    }

  if (GetVar(pdf,"fracMisTag") != NULL)
    {
      if (atoi(Utility->GetGenericParam("CtrlMisTagWrkFlow").c_str()) != 3 )
        {
          myString.clear(); myString.str("");
          myString << fitParam->operator[](Utility->GetFitParamIndx("fracMisTag"))->operator[](q2BinIndx).c_str();
          SetValueAndErrors(pdf,"fracMisTag",1.0,&myString,&value,&errLo,&errHi);
         if ((atoi(Utility->GetGenericParam("CtrlMisTagWrkFlow").c_str()) == 0) || ((errLo == 0.0) && (errHi == 0.0))) GetVar(pdf,"fracMisTag")->setConstant(true);
         else                                                                                                          GetVar(pdf,"fracMisTag")->setConstant(false);
         }
      else
        {
          myString.clear(); myString.str("");
          myString << static_cast<double>(countMisTag) / static_cast<double>(countMisTag + countGoodTag);
          SetValueAndErrors(pdf,"fracMisTag",1.0,&myString,&value,&errLo,&errHi);
          GetVar(pdf,"fracMisTag")->setConstant(true);
        }
     cout << "fracMisTag =" << GetVar(pdf,"fracMisTag")->getVal() << endl;
    }

   return value;
}


void MakeDatasets (B0KstMuMuTreeContent* NTuple, unsigned int FitType, int SampleType)
{
  stringstream myString;

  // ###########################
  // # Define useful variables #
  // ###########################
  B0mass  = new RooRealVar("B0mass","#font[12]{m}(#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{+}}}}#kern[-0.3]{#pi}#kern[-0.3]{#lower[0.6]{^{#font[122]{\55}}}}#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}})",Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str()),Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str()),"GeV");
  mumuMass           = new RooRealVar("mumuMass","#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}} inv. mass",0.0,6.0,"GeV");
  mumuMassE          = new RooRealVar("mumuMassE","#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}} inv. mass error",0.0,0.5,"GeV");
  ctK     = new RooRealVar("ctK","cos#theta_{K}",-1.0,1.0);
  ctL     = new RooRealVar("ctL","cos#theta_{L}",-1.0,1.0);
  phi     = new RooRealVar("phi","#phi", - Utility->PI ,Utility->PI);
  tagB0   = new RooRealVar("tagB0","Tagged Events",-0.1,1.1);
  genSignal     = new RooRealVar("genSignal","gen Signal",0.0,3.0);

  if ( (FitType == 206) ||
       (FitType == 106)
     )
    {
      RooArgSet Vars;
      Vars.add(*B0mass);
      Vars.add(*mumuMass);
      Vars.add(*mumuMassE);
      Vars.add(*ctK);
      Vars.add(*ctL);
      Vars.add(*phi);
      Vars.add(*tagB0);
      Vars.add(*genSignal);

      SingleCandNTuple           = new RooDataSet("SingleCandNTuple"          ,"SingleCandNTuple"          ,Vars);
      SingleCandNTuple_RejectPsi = new RooDataSet("SingleCandNTuple_RejectPsi","SingleCandNTuple_RejectPsi",Vars);

      // #############################
      // # Load values from the tree #
      // #############################
      NTuple->ClearNTuple( SampleType );
      NTuple->SetBranchAddresses(theTree, SampleType);
      int nEntries = theTree->GetEntries();
      for (int entry = 0;  entry < 0.01* nEntries; entry++)
	{
	  theTree->GetEntry(entry);
	 
	  if (
              (FitType == 106) || (FitType == 206)
	     )
	    { 
              if ((FitType == 206) )
                {
                   Vars.setRealValue("mumuMass",          NTuple->mumu_mass);
                   Vars.setRealValue("ctK",               NTuple->gen_cos_theta_k);
                   Vars.setRealValue("ctL",               NTuple->gen_cos_theta_l);
                   Vars.setRealValue("phi",               NTuple->gen_phi);
                }
              else if (((FitType == 106))
                      )
                {
                   Vars.setRealValue("B0mass",            NTuple->tagged_mass);
        	   Vars.setRealValue("mumuMass",          NTuple->mumuMass);
                   Vars.setRealValue("mumuMassE",         NTuple->mumuMassE);
	           Vars.setRealValue("ctK",               NTuple->cos_theta_k);
                   Vars.setRealValue("ctL",               NTuple->cos_theta_l);
                   Vars.setRealValue("phi",               NTuple->phi_kst_mumu);
                   Vars.setRealValue("tagB0",             NTuple->tagB0);
                   Vars.setRealValue("genSignal",         NTuple->genSignal);                         
                }

              // ########################
              // # NTuple with all data #
              // ########################
	      SingleCandNTuple->add(Vars);
          if (
               (((FitType == 106)) &&
                (((strcmp(CTRLfitWRKflow.c_str(),"trueGoodtag") == 0) && ((NTuple->tagB0 ==1 && NTuple->genSignal ==1) || (NTuple->tagB0 == 0 && NTuple->genSignal ==2))) ||
                 ((strcmp(CTRLfitWRKflow.c_str(),"trueMistag")  == 0) && (!((NTuple->tagB0 ==1 && NTuple->genSignal ==1) || (NTuple->tagB0 ==0 && NTuple->genSignal ==2))))||
                  (strcmp(CTRLfitWRKflow.c_str(),"trueAll")     == 0) ||
                  (strcmp(CTRLfitWRKflow.c_str(),"allEvts")     == 0)) &&
                (Utility->PsiRejection(NTuple->tagged_mass,NTuple->mumuMass,NTuple->mumuMassE,"rejectPsi",true) == true))
             )
           SingleCandNTuple_RejectPsi->add(Vars);             
         }	
       }
      cout << "\n[ExtractYield::MakeDatasets]\t@@@ NTuple with all data @@@" << endl;
      SingleCandNTuple->Print("v");

      cout << "[ExtractYield::MakeDatasets]\t@@@ NTuple without J/psi and psi(2S) regions @@@" << endl;
      SingleCandNTuple_RejectPsi->Print("v");
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
			                      unsigned int q2BinIndx,  int efftest)
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

  if ((FitType == 206)) *TotalPDF = MakeAngWithEffPDF(q2BinIndx,y,z,p,FitType,useEffPDF,VarsAng, efftest);
  else
    {
      cout << "[ExtractYield::InstantiateGen3AnglesFit]\tIncorrect configuration sequence : useSignal = " << "useSignal" <<  endl;
      exit (EXIT_FAILURE);
    }

}
     //#############################
     //# Reco level fitting  PDF #
     //#############################
void InstantiateReco3AnglesFit (RooAbsPdf** TotalPDF,
                                bool useEffPDF,
                                RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                unsigned int FitType,
                                vector<vector<unsigned int>*>* configParam,
                                unsigned int q2BinIndx,  int efftest )
// #########################
// # y: angle cos(theta_l) #
// # z: angle cos(theta_K) #
// # p: angle phi          #
// #########################
{
  // ################################
  //  # Read configuration variables #
  //  ################################
  unsigned int useSignal = configParam->operator[](Utility->GetConfigParamIndx("SigType"))->operator[](q2BinIndx);

  // ###################
  // # Local variables #
  // ###################
  stringstream myString;
  RooArgSet* VarsAng = new RooArgSet("VarsAng");

  if ((FitType == 106)){
      if ((strcmp(CTRLfitWRKflow.c_str(),"trueAll") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"allEvts") == 0) ){
          AngleGood = MakeAngWithEffPDF(q2BinIndx,y,z,p,FitType,useEffPDF,VarsAng,efftest);
          AngleMisTag = MakeAngWithEffPDF(q2BinIndx,y,z,p,FitType*10,useEffPDF,VarsAng,efftest);
          fracMisTag = new RooRealVar("fracMisTag","Fraction of mistag",0.0,0.0,1.0);
          *TotalPDF = new RooAddPdf("SignalRECO","AngleGood + nMisTag * (AngleMisTag)",RooArgSet(*AngleMisTag,*AngleGood),*fracMisTag);
      }
      else if (strcmp(CTRLfitWRKflow.c_str(),"trueGoodtag") == 0)
          *TotalPDF = MakeAngWithEffPDF(q2BinIndx,y,z,p,FitType,useEffPDF,VarsAng,efftest);
      else if (strcmp(CTRLfitWRKflow.c_str(),"trueMistag")  == 0)
          *TotalPDF = MakeAngWithEffPDF(q2BinIndx,y,z,p,FitType*10,useEffPDF,VarsAng,efftest);
  }
  else {
      cout << "[ExtractYield::InstantiateGen3AnglesFit]\tIncorrect configuration sequence : useSignal = " << "useSignal" <<  endl;
      exit (EXIT_FAILURE);
       }
}


RooFitResult* Make3AnglesFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* y, RooRealVar* z,RooRealVar* p,  unsigned int FitType, RooArgSet* vecConstr, TCanvas* Canv, int specBin)
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

  if ((FitType == 206) || (FitType == 106))
    {
      clock_t fitstart = clock();
      TBenchmark* Fittime = new TBenchmark();
      Fittime->Start("test");
      // ###################
      // # Make actual fit #
      // ###################
      if (ROOMINIMIZER == true )
       {
         clock_t nlls = clock();
         TBenchmark* nll = new TBenchmark();
         nll->Start("nll");
         RooAbsReal* MC_nll = (*TotalPDF)->createNLL(*dataSet,RooFit::NumCPU(10));    
         RooMinuit Minuit(*MC_nll) ;
         Minuit.setPrintLevel(-1);
         clock_t nlle = clock();
         double nlltime = (double)(nlle - nlls)/ CLOCKS_PER_SEC;
         cout <<  "[ExtractYield::Make3AnglesFit]\t @@@ Time for nll clock_t \t" << nlltime << "s\t @@@" << endl;
         nll->Stop("nll");
         nll->Show("nll");

         clock_t migrads = clock();
         TBenchmark* migrad = new TBenchmark();
         migrad->Start("migrad");
         Minuit.migrad() ;
         clock_t migrade = clock();
         double migradtime = (double)(migrade - migrads)/ CLOCKS_PER_SEC;
         cout <<  "[ExtractYield::Make3AnglesFit]\t @@@ Time for Migrad clock_t \t" << migradtime << "s\t @@@" << endl;
         migrad->Stop("migrad");
         migrad->Show("migrad");
        
         if (Minuit.migrad() == 0) {
         clock_t hesses = clock();
         TBenchmark* hesse = new TBenchmark();
         hesse->Start("hesse");
         Minuit.hesse();
         clock_t hessee = clock();
         double hessetime = (double)(hessee - hesses)/ CLOCKS_PER_SEC;
         cout <<  "[ExtractYield::Make3AnglesFit]\t @@@ Time for Hesse clock_t \t" << hessetime << "s\t @@@" << endl;
         hesse->Stop("hesse");
         hesse->Show("hesse");
         fitResult = Minuit.save(); 
         }
        else exit(EXIT_FAILURE);
       }
      else 
	fitResult = (*TotalPDF)->fitTo(*dataSet,Save(true),Minos(0),Minimizer(MINIMIZER),NumCPU(8));
      cout <<  "[ExtractYield::Make3AnglesFit]\t @@@ status = " << fitResult->status() << endl;
      clock_t fitend = clock();
      double fittime = (double) (fitend - fitstart)/ CLOCKS_PER_SEC;
       cout <<  "[ExtractYield::Make3AnglesFit]\t @@@ Time for fitting \t" << fittime << "s\t @@@" << endl; 
       
       Fittime->Stop("test");
       Fittime->Show("test");
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
       clock_t plotstart = clock();
      if ((PLOT==true ) &&  CheckGoodFit(fitResult) == true)
      {
      Canv->cd(1);
      RooPlot* myFrameY = y->frame(NBINS);

      dataSet->plotOn(myFrameY, Name(MakeName(dataSet,specBin).c_str()));
      if ((FitType == 206) )  legNames[nElements++] = "GEN-MC";
      else if ((FitType == 106) )  legNames[nElements++] = "Reco-MC";
      (*TotalPDF)->plotOn(myFrameY, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), Project(RooArgSet(*z,*p)));
      if ((FitType == 206) )  legNames[nElements++] = "Gen p.d.f.";
      else if ((FitType == 106))  legNames[nElements++] = "Reco p.d.f.";
      TPaveText* paveTextY = new TPaveText(0.12,0.78,0.4,0.86,"NDC");
      paveTextY->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,specBin).c_str())));
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
     
      myFrameY->SetMaximum(1.4*myFrameY->GetMaximum());
      myFrameY->Draw();
      legY->Draw("same");
     
      // ###########################
      // # Costheta-k plot results #
      // ###########################
      Canv->cd(2);
      RooPlot* myFrameZ = z->frame(NBINS);

      dataSet->plotOn(myFrameZ, Name(MakeName(dataSet,specBin).c_str()));
      if ((FitType == 206))  legNames[nElements++] = "GEN-MC";
      else if ((FitType == 106) )  legNames[nElements++] = "Reco-MC";
      (*TotalPDF)->plotOn(myFrameZ, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), Project(RooArgSet(*y,*p)));
      if ((FitType == 206))  legNames[nElements++] = "Gen p.d.f.";
      else if ((FitType == 106))  legNames[nElements++] = "Reco p.d.f";
 
      TPaveText* paveTextZ = new TPaveText(0.12,0.82,0.4,0.86,"NDC");
      paveTextZ->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,specBin).c_str())));
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
     
      myFrameZ->SetMaximum(1.4*myFrameZ->GetMaximum());
      myFrameZ->Draw();
      legZ->Draw("same");


      // ###########################
      // # Costheta-l plot results #
      // ###########################
      Canv->cd(3);
      RooPlot* myFrameP = p->frame(NBINS);

      dataSet->plotOn(myFrameP, Name(MakeName(dataSet,specBin).c_str()));
      if ((FitType == 206))  legNames[nElements++] = "GEN-MC";
      else if ((FitType == 106) )  legNames[nElements++] = "Reco-MC";
      (*TotalPDF)->plotOn(myFrameP, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), Project(RooArgSet(*y,*z)));
      if ((FitType == 206) )  legNames[nElements++] = "Gen p.d.f.";
      else if ((FitType == 106))  legNames[nElements++] = "Reco p.d.f";

      TPaveText* paveTextP = new TPaveText(0.12,0.82,0.4,0.86,"NDC");
      paveTextP->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameP->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,specBin).c_str())));
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
     
      myFrameP->SetMaximum(1.4*myFrameP->GetMaximum());
      myFrameP->Draw();
      legP->Draw("same");
    }

    Canv->Modified();
    Canv->Update();

    if (SAVEPLOT == true)
    {
     myString.clear(); myString.str("");
     myString << (*TotalPDF)->getPlotLabel() << "_Canv" << specBin << ".pdf";
     Canv->Print(myString.str().c_str());

     myString.clear(); myString.str("");
     myString << (*TotalPDF)->getPlotLabel() << "_Canv" << specBin  << ".root";
     Canv->Print(myString.str().c_str());
    }
   clock_t plotend = clock(); 
   double plottime = (double) (plotend - plotstart)/ CLOCKS_PER_SEC;
   cout <<  "[ExtractYield::Make3AnglesFit]\t @@@ Time for plotting \t" << plottime << "s\t @@@" << endl;
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
				    RooArgSet* vecConstr, int efftest,
                                    int nscan)
{
  // ###################
  // # Local variables #
  // ###################
  stringstream myString;   

  if ( (FitType == 106) || (FitType == 206) )
    {
      RooAbsPdf*  TotalPDFq2Bins[q2Bins->size()-1];
      RooFitResult* fitResult = NULL;

      RooDataSet* dataSet_q2Bins[q2Bins->size()-1];
      for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
       {
         myString.clear(); myString.str("");
         myString << "(mumuMass*mumuMass) > " << q2Bins->operator[](i) << " && (mumuMass*mumuMass) <= " << q2Bins->operator[](i+1);
         cout << "\n[ExtractYield::IterativeAnglesFitq2Bins]\tCut string: " << myString.str() << endl;
         dataSet_q2Bins[i] = (RooDataSet*)dataSet->reduce(myString.str().c_str());
         cout << "[ExtractYield::IterativeAnglesFitq2Bins]\tNumber of events : " << dataSet_q2Bins[i]->sumEntries() << endl;

      if ( FitType == 206)
        {
          TCanvas*    Gen[q2Bins->size()-1];
          Gen[i] = new TCanvas("Gen","Gen",10,10,700,500);
          Gen[i]->Divide(2,2);

          myString.clear(); myString.str("");
          myString << "GenTotalPDFq2Bin_" << i;
          InstantiateGen3AnglesFit(&TotalPDFq2Bins[i],useEffPDF,y,z,p,FitType,i,efftest);

          // ########################################
          // # Initialize p.d.f & Apply constraints #
          // ########################################
          clock_t copystart = clock();
          CopyFitResults(TotalPDFq2Bins[i],i,nscan,fitParam,0,0);
          ClearVars(vecConstr);
          clock_t copyend = clock();
          double  copytime = (double) (copyend - copystart)/ CLOCKS_PER_SEC;
          cout <<  "[ExtractYield::Iterative3AnglesFitq2Bins]\t @@@ Time for Copying \t" << copytime << "s\t @@@" << endl;
 
          fitResult = Make3AnglesFit(dataSet_q2Bins[i],&TotalPDFq2Bins[i],y,z,p,FitType,vecConstr,Gen[i],specBin);
          if (CheckGoodFit(fitResult) == true) cout << "\n[ExtractYield::Iterative3AnglesFitq2Bins]\t@@@ Fit converged ! @@@" << endl;
          else                                 cout << "\n[ExtractYield::Iterative3AnglesFitq2Bins]\t@@@ Fit didn't converge ! @@@" << endl;
                
        }
      else if (  FitType == 106)
        {
          TCanvas*    Reco[q2Bins->size()-1]; 
          Reco[i] = new TCanvas("Reco","Reco",10,10,700,500);
          Reco[i]->Divide(2,2);

          unsigned int countMisTag  = 0;
          unsigned int countGoodTag = 0;
          for (int j = 0; j < static_cast<int>(dataSet_q2Bins[i]->sumEntries()); j++)
           {
             if ( ((dataSet_q2Bins[i]->get(j)->getRealValue("tagB0") == 1) && (dataSet_q2Bins[i]->get(j)->getRealValue("genSignal") == 1)) || ((dataSet_q2Bins[i]->get(j)->getRealValue("tagB0") == 0) && (dataSet_q2Bins[i]->get(j)->getRealValue("genSignal") == 2)))   
             countGoodTag++;
             else
              countMisTag++;
           } 
          cout << "[ExtractYield::IterativeMassFitq2Bins]\tDynamic mis-tag fraction : " << static_cast<double>(countMisTag) / static_cast<double>(countMisTag + countGoodTag) << " = (" << countMisTag << "/(" << countMisTag << "+" << countGoodTag << "))" << endl;
          cout << countGoodTag << endl;
          cout << countMisTag << endl; 
          myString.clear(); myString.str("");
          myString << "RecoTotalPDFq2Bin_" << i;
          InstantiateReco3AnglesFit(&TotalPDFq2Bins[i],useEffPDF,y,z,p,FitType,configParam,i,efftest);
      // ########################################
      // # Initialize p.d.f. & Apply constraints#
      // ########################################
       if (nscan < 0) CopyFitResults(TotalPDFq2Bins[i],i,nscan,fitParam,countMisTag,countGoodTag);
       else 
         {
           cout << "\n[ExtractYield::IterativeMass3AnglesFitq2Bins]\t@@@ I am scanning the fitting values @@@" << endl;
           for (int scan = 1; scan < nscan + 1; scan++ )
             {
               CopyFitResults(TotalPDFq2Bins[i],i,nscan,fitParam,countMisTag,countGoodTag);
             }
         }
       ClearVars(vecConstr);

      // ###################
      // # Perform the fit #
      // ##################
      fitResult = Make3AnglesFit(dataSet_q2Bins[i],&TotalPDFq2Bins[i],y,z,p,FitType,vecConstr,Reco[i],specBin);
      if (CheckGoodFit(fitResult) == true) cout << "\n[ExtractYield::Iterative3AnglesFitq2Bins]\t@@@ Fit converged ! @@@" << endl;
      else                                 cout << "\n[ExtractYield::Iterative3AnglesFitq2Bins]\t@@@ Fit didn't converge ! @@@" << endl;
        }
      PrintFitResults(*&TotalPDFq2Bins[i], fitResult);
    }
  }
}



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
      unsigned int FitType      = atoi(argv[1]);
      int SampleType            = -1;

      int efftest               = -1;
      int nscan                 = -1;

      bool useEffPDF            = false;

      TFile* NtplFile           = NULL;

      if (((FitType == 106) || (FitType == 206))
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
	  if ((FitType == 106) || (FitType == 206)
             )
	    {
	      fileName           = argv[2];
	      correct4Efficiency = argv[3];
	    }
          if ((FitType == 106))
          SampleType         = 1;
          else if ((FitType == 206) )
          SampleType         = 0;
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

          // ################################
          // # test different eff Functions #
          // ################################
          if (argc >= 6) efftest = atoi(argv[5]);
          if (efftest < 3 || efftest >8) 
             {
              cout << "[ExtractYield::main]\tIncorrect option parameter " << efftest << endl;
              exit (EXIT_FAILURE);
             }

          // ######################
          // # scan initial value #   
          // ######################
          if (argc >= 7) nscan = atoi(argv[6]);
          if (nscan <0 || nscan >100)
             {
               cout << "[ExtractYield::main]\tIncorrect option parameter " << nscan << endl;
               exit (EXIT_FAILURE);
             }
          if ((correct4Efficiency == "yesEffCorr")) useEffPDF = true;

          cout << "\n[ExtractYield::main]\t@@@ Input variables from command line @@@" << endl;
          cout << "- input/outputFile.root = " << fileName.c_str() << endl;
          cout << "- correct4Efficiency = "    << correct4Efficiency << endl;
          cout << "- specBin = "               << specBin << endl;
          cout << "- FitType = "               << FitType << endl;
          cout << "- useEffPDF = "             << useEffPDF << endl;
          cout << "- SampleType = "            << SampleType << endl;
          cout << "- nscan = "                 << nscan << endl;       
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

          CTRLfitWRKflow = Utility->GetGenericParam("CtrlFitWrkFlow");

          cout << "Use MINOS: "                          << Utility->GetGenericParam("UseMINOS").c_str() << " (0 = false; 1 = true)" << endl;
          cout << "Control fit workflow: "               << CTRLfitWRKflow.c_str() << endl;
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
              (FitType == 106) || (FitType == 206) 
             )
            {
              NtplFile = new TFile(fileName.c_str(),"READ");
              if ( (FitType == 106))      theTree  = (TTree*) NtplFile->Get("ntuple");
              else if ( (FitType == 206))      theTree  = (TTree*) NtplFile->Get("tree");
              NTuple   = new B0KstMuMuTreeContent( SampleType);
              NTuple->Init(SampleType);

	      // #################
	      // # Make datasets #
	      // #################
	      cout << "\n[ExtractYield::main]\t@@@ Making datasets @@@" << endl;
              clock_t start = clock();
	      MakeDatasets(NTuple,FitType,SampleType);
              clock_t end = clock();
              double time = (double) (end - start)/ CLOCKS_PER_SEC;
              cout << "\n[ExtractYield::main]\t@@@ Time Reading tree is\n" << time <<"s @@@" << endl;
         
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
                                                                                    &vecConstr, efftest,
                                                                                    nscan);
               
	        else if ((FitType == 106))                Iterative3AnglesFitq2Bins(SingleCandNTuple_RejectPsi,
	                                                                            useEffPDF,
	                                                                            ctL,
	                                                                            ctK,
	                                                                            phi,
	                                                                            specBin,
	                                                                            FitType,
	                                                                            &q2Bins,
	                                                                            &configParam,&fitParam,
	                                                                            &vecConstr, efftest,
	                                                                            nscan);
	    }

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
      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  Signa  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "FitType = 206: 3D P5p-Fl (cos(theta_K), cos(theta_l)) per q^2 bin on Gen-level" << endl;
      cout << "FitType = 205: 3D Si-AFB (cos(theta_K), cos(theta_l)) per q^2 bin on Gen-level" << endl;
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      
      return EXIT_FAILURE;
    }

  return 0;
}
