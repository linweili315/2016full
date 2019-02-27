#ifndef UTILS_H
#define UTILS_H

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TH3.h>
#include <TMatrixTSym.h>
#include <TGraphAsymmErrors.h>

#if ROOFIT

#include <RooDataHist.h>
#include <RooAbsReal.h>
#include <TH3.h>
#include <TTree.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include  <TStopwatch.h>
#include <RooRealProxy.h>
#include <RooCategoryProxy.h>
#include <RooAbsCategory.h>

#include  <TStopwatch.h>
#include <RooAbsPdf.h>
#include <RooHistPdf.h>
#include <RooFitResult.h>
#endif

#include <string>
#include <vector>
#include "B0KstMuMuTreeContent.h"


class Utils
{
 public:
  
  Utils(bool rightFlavorTag = true);
  ~Utils() {};

  // #################################
  // # Data structure for efficiency #
  // #################################
  struct _effStruct
  {
    // ###########################################################
    // # Numerators and Denominators = number of events * weight #
    // ###########################################################
    double *Num1, *Num2;
    double *Den1, *Den2;
    
    // ###################################################
    // # Poissonian errors = number of events * weight^2 #
    // ###################################################
    double *Err2PoisNum1, *Err2PoisNum2;
    double *Err2PoisDen1, *Err2PoisDen2;
    
    // #################
    // # Weight errors #
    // #################
    double *Err2WeigNum1, *Err2WeigNum2;
    double *Err2WeigDen1, *Err2WeigDen2;

    // #########################
    // # Correspondence :      #
    // # GenFilter     <--> N1 #
    // # SingleCand    <--> N2 #
    // # GenNoFilter   <--> D1 #
    // # AllCandFilter <--> D2 #
    // #########################
  };
  typedef struct _effStruct effStruct;

  struct _effValue
  {
    double Num1, Num2;
    double Den1, Den2;
    
    // Poissonian errors
    double Err2PoisNum1, Err2PoisNum2;
    double Err2PoisDen1, Err2PoisDen2;
    
    // Weight errors
    double Err2WeigNum1, Err2WeigNum2;
    double Err2WeigDen1, Err2WeigDen2;
  };
  typedef struct _effValue effValue;

  void ReadAllBins  (std::string fileName, std::vector<double>* q2Bins);
  void Readq2Bins   (std::string fileName, std::vector<double>* q2Bins);
  double ReadLumi                   (std::string fileName);

 
  void ReadFitStartingValues             (std::string fileName, std::vector<std::vector<std::string>*>* vecParam, std::vector<std::vector<unsigned int>*>* configParam, const unsigned int dataBlockN);
  void ReadParVsq2Bins                   (std::string fileName, std::string praName, std::vector<std::string>** vecParam); 

  int WhatIsThis (std::string fileName);
  void SaveFitValues (std::string fileName, std::vector<std::string>* vecParStr, int indx, std::string howOpen, std::string str = "");
  
  std::string MakeAnalysisPATH (std::string relativePath);
  unsigned int ParFileBlockN   (std::string blockName);

  unsigned int GetFitParamIndx    (std::string varName);
  unsigned int GetConfigParamIndx (std::string varName);

  bool PsiRejection   (double myB0Mass, double myMuMuMass, double myMuMuMassE, std::string seleType, bool B0andPsiCut = false);

  void   ReadGenericParam     (std::string fileName);
  bool   SetGenericParam      (std::string parName, std::string val);
  std::string GetGenericParam (std::string parName);

  double* MakeBinning (std::vector<double>* STLvec);
  #if ROOFIT
  std::string Transformer (std::string varName, bool doIt, double& varValOut, double& varValOutELo, double& varValOutEHi, RooFitResult* fitResult = NULL, RooRealVar* varValIn1 = NULL, RooRealVar* varValIn2 = NULL); 
  void AntiTransformer (std::string varName, bool doIt, double& varValOut, double& varValOutELo, double& varValOutEHi, RooRealVar* varValIn1 = NULL, RooRealVar* varValIn2 = NULL);
  RooAbsPdf* ReadRTEffPDF ( unsigned int q2Indx, unsigned int test);
  RooAbsPdf* ReadWTEffPDF ( unsigned int q2Indx, unsigned int test);
  #endif
  double muonMass;
  double pionMass;
  double kaonMass;
  double kstMass;
  double B0Mass;
  double JPsiMass;
  double PsiPMass;

  double JPsiBF;
  double JPsiKpiBF;
  double KstMuMuBF;
  double KstKpiMuMuBF;
  double PsiPBF;
  double PsiPKpiBF;

  double muonMassErr;
  double pionMassErr;
  double kaonMassErr;
  double B0MassErr;
  double kstSigma;

  double PI;

  unsigned int NcoeffThetaL;
  unsigned int NcoeffThetaK;
  unsigned int NcoeffPhi;

  bool RIGHTflavorTAG;

  int B0ToKstMuMu;
  int B0ToJPsiKst;
  int B0ToPsi2SKst;

  unsigned int nFitParam;
  unsigned int nConfigParam;


 private:

  TF1* KstMassShape;

  std::vector<std::string> HLTpath;
  std::vector<double> VecHLTCutVar1;
  std::vector<double> VecHLTCutVar2;
  std::vector<double> VecHLTentries;
  std::vector<std::string> TrigTable;

  std::vector<double> PreCuts;
  std::vector<double> SeleCuts;
  std::vector<std::string> GenericPars;

  double ProbThreshold;
  double scrambleFraction;

  unsigned int nFitObserv;

};
#endif
