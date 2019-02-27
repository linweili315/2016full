#include <TTree.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TH1.h>
#include <Math/Functor.h>
#include <TMath.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>

#include <sstream>
#include <time.h>
#include <ctime>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>

using std::cout;
using std::endl;
using std::stringstream;
using std::string;
using namespace RooFit;

// #############
// # Variables #
// #############
double q2Min = 0.;
double q2Max = 0.;
int q2Bin = -1;
double mumuMass2 = -1.;

double wight ;
double M_6s ;
double M_6c ;
double f_1s;
double f_3 ;
double f_4 ;
double f_5;
double f_6s ;
double f_6c;
double f_7 ;
double f_8 ;
double f_9 ;

double Mi ;
// ##############
// # Parameters #
// ##############
double value;
double error;
double Fl = 0.;
double AFB = 0.;
double S3 = 0.;
double S4 = 0.;
double S5 = 0.;
double S6 = 0.;
double S7 = 0.;
double S8 = 0.;
double S9 = 0.;

double P1 = 0.;
double P2 = 0.;
double P3 = 0.;
double P4p = 0.;
double P5p = 0.;
double P6p = 0.;
double P8p = 0.;

// ####################
// # angular varables #
// ####################
// #################
// # Reco-Level MC #
// #################
double cos_theta_k = 0;
TBranch *b_cos_theta_k = 0;

double cos_theta_l = 0;
TBranch *b_cos_theta_l = 0;

double phi_kst_mumu = 0;
TBranch *b_phi_kst_mumu = 0;

float mumuMass = 0;
TBranch *b_mumuMass = 0;

float tagB0 = 0;
TBranch *b_tagB0 = 0;

float genSignal = 0;
TBranch *b_genSignal = 0;

// ################
// # gen-Level MC #
// ################
float gen_cos_theta_k = 0;
TBranch *b_gen_cos_theta_k = 0;

float gen_cos_theta_l = 0;
TBranch *b_gen_cos_theta_l = 0;

float gen_phi = 0;
TBranch *b_gen_phi = 0;

float mumu_mass = 0;
TBranch *b_mumu_mass = 0;

void quzhi() {
    char a = '0' + q2Bin;
    switch (a) {
        case '0' :
            q2Min = 1.0;
            q2Max = 2.0;
            q2Bin = 0;
            break;

        case '1' :
            q2Min = 2.0;
            q2Max = 4.3;
            q2Bin = 1;
            break;

        case '2' :
            q2Min = 4.3;
            q2Max = 6.0;
            q2Bin = 2;
            break;

        case '3':
            q2Min = 6.0;
            q2Max = 8.68;
            q2Bin = 3;
            break;

        case '4':
            q2Min = 8.68;
            q2Max = 10.09;
            q2Bin = 4;
            break;

        case '5':
            q2Min = 10.09;
            q2Max = 12.86;
            q2Bin = 5;
            break;

        case '6':
            q2Min = 12.86;
            q2Max = 14.18;
            q2Bin = 6;
            break;

        case '7':
            q2Min = 14.18;
            q2Max = 16.0;
            q2Bin = 7;
            break;

        case '8':
            q2Min = 16.0;
            q2Max = 19.0;
            q2Bin = 8;
            break;

        default:
            break;
    }
}

void GenCalValue(int q2Bin, double q2Min, double q2Max, string Paratype) 
{
    value = 0.0;
    error= 0.0;
    M_6s = 0.;
    M_6c = 0.;
    f_1s = 0.;
    f_3 = 0.;
    f_4 = 0.;
    f_5 = 0.;
    f_6s = 0.;
    f_6c = 0.;
    f_7 = 0.;
    f_8 = 0.;
    f_9 = 0.;
   
    TH1D *f1s = new TH1D("f1s", "f1s", 100, 0.0, 1.0);
    TH1D *m6s = new TH1D("m6s", "m6s", 200, -5.0, 5.0);
    TH1D *m6c = new TH1D("m6c", "m6c", 200, -5.0, 5.0);
    TH1D *f3 = new TH1D("f3", "f3", 100, -5.0, 5.0);
    TH1D *f4 = new TH1D("f4", "f4", 200, -5.0, 5.0);
    TH1D *f5 = new TH1D("f5", "f5", 200, -5.0, 5.0);
    TH1D *f7 = new TH1D("f7", "f7", 200, -5.0, 5.0);
    TH1D *f8 = new TH1D("f8", "f8", 200, -5.0, 5.0);
    TH1D *f9 = new TH1D("f9", "f9", 200, -5.0, 5.0);
    // ############################
    // # Read tree and set Branch #
    // ############################
    TFile *f = new TFile("gen_B0_miniaodWithoutGenCuts.root");
    cout << "\n[Moment::GenCalValue]\tTry to open\t" << "gen_B0_miniaodWithoutGenCuts.root" << endl;

    TTree *t = (TTree *) f->Get("tree");
    t->SetBranchAddress("gen_cos_theta_k", &gen_cos_theta_k, &b_gen_cos_theta_k);
    t->SetBranchAddress("gen_cos_theta_l", &gen_cos_theta_l, &b_gen_cos_theta_l);
    t->SetBranchAddress("gen_phi", &gen_phi, &b_gen_phi);
    t->SetBranchAddress("mumu_mass", &mumu_mass, &b_mumu_mass);

    Int_t entries = (Int_t) t->GetEntries();
    cout << "\n[Moment::GenCalValue]\tTotal number of events in the tree: " << entries << " @@@" << endl;
    for (Int_t i = 0; i < entries; i++)
      {
        t->GetEntry(i);
        mumuMass2 = mumu_mass * mumu_mass;

       // ###############################
       // # define orthogonal functions #
       // ###############################
        f_1s = 1 - gen_cos_theta_k * gen_cos_theta_k;
        M_6s = (1 - gen_cos_theta_k * gen_cos_theta_k) * gen_cos_theta_l;
        M_6c = gen_cos_theta_k * gen_cos_theta_k * gen_cos_theta_l;
        f_3 = (1 - gen_cos_theta_k * gen_cos_theta_k) * (1 - gen_cos_theta_l * gen_cos_theta_l) * TMath::Cos(2 * gen_phi);
        f_4 = 4 * gen_cos_theta_k * gen_cos_theta_l * TMath::Cos(gen_phi) * TMath::Sqrt((1 - gen_cos_theta_k * gen_cos_theta_k) * (1 - gen_cos_theta_l * gen_cos_theta_l));
        f_5 = 2 * gen_cos_theta_k * TMath::Cos(gen_phi)* TMath::Sqrt((1 - gen_cos_theta_k * gen_cos_theta_k) * (1 - gen_cos_theta_l * gen_cos_theta_l));
        f_7 = 2 * gen_cos_theta_k * TMath::Sin(gen_phi) * TMath::Sqrt((1 - gen_cos_theta_k * gen_cos_theta_k) * (1 - gen_cos_theta_l * gen_cos_theta_l));
        f_8 = 4 * gen_cos_theta_k * gen_cos_theta_l * TMath::Sin(gen_phi) * TMath::Sqrt((1 - gen_cos_theta_k * gen_cos_theta_k) * (1 - gen_cos_theta_l * gen_cos_theta_l));
        f_9 = (1 - gen_cos_theta_k * gen_cos_theta_k) * (1 - gen_cos_theta_l * gen_cos_theta_l) * TMath::Sin(2 * gen_phi);
      // ##################################
      // # begin to compute the variables #
      // ##################################
      if (mumuMass2 > q2Min && mumuMass2 < q2Max) 
         {
            f1s->Fill(f_1s);
            m6s->Fill(M_6s);
            m6c->Fill(M_6c);
            f3->Fill(f_3);
            f4->Fill(f_4);
            f5->Fill(f_5);
            f7->Fill(f_7);
            f8->Fill(f_8);
            f9->Fill(f_9);       
         }
      } // end for 
    if  (Paratype == "FlS") 
      {
        value = 2.0 - 2.5 * f1s->GetMean();
        error = 2.5 * f1s->GetRMS() /TMath::Sqrt(f1s->GetEntries() - 1);
      }
    else if  (Paratype == "AFBS")
      {
        value = 3.0/ 4.0 * (3.0 * m6s->GetMean() - 2.0 * m6c->GetMean());
        error = 3.0 / 4.0 * TMath::Sqrt(4.0 * TMath::Power(m6c->GetRMS(), 2) / (m6c->GetEntries() - 1) +
                                        9.0 * TMath::Power(m6s->GetRMS(), 2) / (m6s->GetEntries() - 1));
      }
    else if  (Paratype == "S3S")
      {
        value = 25.0 / 8.0 * f3->GetMean();
        error = 25.0 / 8.0 * f3->GetRMS() /TMath::Sqrt(f3->GetEntries() - 1);
      }
    else if  (Paratype == "S4S")
      {
        value = 25.0 / 8.0 * f4->GetMean();
        error = 25.0 / 8.0 * f4->GetRMS() /TMath::Sqrt(f4->GetEntries() - 1);
      }
    else if (Paratype == "S5S")
      {
        value = 2.5 * f5->GetMean();
        error = 2.5 * f5->GetRMS() /TMath::Sqrt(f5->GetEntries() - 1);
      }
    else if (Paratype == "S7S")
      {
        value = 2.5 * f7->GetMean();
        error = 2.5 * f7->GetRMS() /TMath::Sqrt(f7->GetEntries() - 1);
      }
    else if  (Paratype == "S8S")
      {
        value = 25.0 / 8.0 * f8->GetMean();
        error = 25.0 / 8.0 * f8->GetRMS() /TMath::Sqrt(f8->GetEntries() - 1);
      }
    else if  (Paratype == "S9S")
      {
        value = 25.0 / 8.0 * f9->GetMean();
        error = 25.0 / 8.0 * f9->GetRMS() /TMath::Sqrt(f9->GetEntries() - 1);
      }
    else if  (Paratype == "P1S")
      {
       value = 25./4. * f3->GetMean() / ( 1. - (2.0 - 2.5 * f1s->GetMean()) );
       error = TMath::Sqrt(TMath::Power(2/(1.0 - (2.0 - 2.5 * f1s->GetMean())) * (25/8 * f3->GetRMS() /TMath::Sqrt(f3->GetEntries() - 1)), 2 ) + TMath::Power( (25/4 * f3->GetMean() / TMath::Power( (1.0 - (2.0 - 2.5 * f1s->GetMean())), 2) ) * (2.5 * f1s->GetRMS() /TMath::Sqrt(f1s->GetEntries() - 1)), 2 ) );
      }
    else if  (Paratype == "P2S")
      {
       value = 0.5 * m6s->GetMean() / ( 1. - (2.0 - 2.5 * f1s->GetMean()) );
       error = TMath::Sqrt( TMath::Power( 1/(2 * (1- (2.0 - 2.5 * f1s->GetMean()))) ,2) * (4 * TMath::Power(m6c->GetRMS(), 2) / (m6c->GetEntries() - 1) + 9 * TMath::Power(m6s->GetRMS(), 2) / (m6s->GetEntries() - 1)) + TMath::Power( (3.0 * m6s->GetMean() - 2.0 * m6c->GetMean())/( 2 * TMath::Power( (1-(2.0 - 2.5 * f1s->GetMean())) ,2 )) ,2) * (TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1)) );
      }
    else if (Paratype == "P3S")
      {
       value = - 25.0/8.0 * f9->GetMean() / (1.0 - (2.0 - 2.5 * f1s->GetMean()));
       double a = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
       double b = TMath::Power( 25/8 * f9->GetRMS() ,2 ) / (f9->GetEntries() - 1);
      error = TMath::Sqrt(TMath::Power( 1.0 / (1.0 - (2.0 - 2.5 * f1s->GetMean())) ,2 ) * b + TMath::Power( 25/8 * f9->GetMean() / TMath::Power((1.0 - (2.0 - 2.5 * f1s->GetMean())) ,2) ,2) * a );
      }
    else if (Paratype == "P4pS")
      {
       value = 4 * 25/8 * f4->GetMean() /TMath::Sqrt((2.0 - 2.5 * f1s->GetMean()) * ( 1.0 - (2.0 - 2.5 * f1s->GetMean())) );
       double a = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
       double b = TMath::Power( 25/8 * f4->GetRMS() ,2 ) / (f4->GetEntries() - 1);
       error = 4 * TMath::Sqrt( 1/ ((2.0 - 2.5 * f1s->GetMean()) * ( 1- (2.0 - 2.5 * f1s->GetMean()))) * b + TMath::Power((1 - 2 * (2.0 - 2.5 * f1s->GetMean())) * 25/8 * f4->GetMean() ,2) * a / ( 4 * TMath::Power((1 - (2.0 - 2.5 * f1s->GetMean())) ,3)) );
      }
    else if (Paratype == "P5pS")
      {
        double a = 2.0 - 2.5 * f1s->GetMean();
        value = 2.5 * f5->GetMean() /TMath::Sqrt(a * ( 1 - a) );
        double b = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
        double c = TMath::Power( 2.5 * f5->GetRMS() ,2 ) / (f5->GetEntries() - 1);
        error = TMath::Sqrt( 1/ (a * ( 1- a)) * c + TMath::Power((1 - 2 * a) * 2.5 * f5->GetMean() ,2) * b /( 4 * TMath::Power((1 - a) ,3)) );
      }
    else if (Paratype == "P6pS")
      {
        double a = 2.0 - 2.5 * f1s->GetMean();
        value = - 2.5 * f7->GetMean() /TMath::Sqrt(a * ( 1 - a) );
        double b = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
        double c = TMath::Power( 2.5 * f7->GetRMS() ,2 ) / (f7->GetEntries() - 1);
        error = TMath::Sqrt( 1/ (a * ( 1- a)) * c + TMath::Power((1 - 2 * a) * 2.5 * f7->GetMean() ,2) * b /(4 * TMath::Power((1 - a) ,3)) ); 
       }
    else if (Paratype == "P8pS")
      {
        double a = 2.0 - 2.5 * f1s->GetMean();
        value = 25/8 * f8->GetMean() /TMath::Sqrt(a * ( 1 - a) );
        double b = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
        double c = TMath::Power( 25/8 * f8->GetRMS() ,2 ) / (f8->GetEntries() - 1);
        error = TMath::Sqrt( 1/ (a * ( 1- a)) * c + TMath::Power((1 - 2 * a) * 25/8 * f8->GetMean() ,2) * b / ( 4 * TMath::Power((1 - a) ,3)) );
      }
    cout << "\n[Moment::GenCalValue]\t q2bin = " << q2Bin << ":" << Paratype << " = " << value << "+/- " << error  << endl;
}

void ReCalValue(int q2Bin, double q2Min, double q2Max, string TagType, string Paratype) {

    wight = 0.0;
    M_6s = 0.;
    M_6c = 0.;
    f_1s = 0.;
    f_3 = 0.;
    f_4 = 0.;
    f_5 = 0.;
    f_6s = 0.;
    f_6c = 0.;
    f_7 = 0.;
    f_8 = 0.;
    f_9 = 0.;
    // ############################
    // # Read Effcicency Function #
    // ############################
    TH1D *f1s = new TH1D("f1s", "f1s", 100, -100., 100.0);
    TH1D *m6s = new TH1D("m6s", "m6s", 200, -100.0, 100.0);
    TH1D *m6c = new TH1D("m6c", "m6c", 200, -100.0, 100.0);
    TH1D *f3 = new TH1D("f3", "f3", 200, -100.0, 100.0);
    TH1D *f4 = new TH1D("f4", "f4", 200, -100.0, 100.0);
    TH1D *f5 = new TH1D("f5", "f5", 200, -100.0, 100.0);
    TH1D *f7 = new TH1D("f7", "f7", 200, -100.0, 100.0);
    TH1D *f8 = new TH1D("f8", "f8", 200, -100.0, 100.0);
    TH1D *f9 = new TH1D("f9", "f9", 200, -100.0, 100.0);

    int test = 4;  // you can choose the Maximum order
    stringstream myString;
    myString.clear();
    myString.str("");
    myString << "effProjection_sh" << test << "o_b" << q2Bin << "ct_25_25_25.root";
    cout << "\n[Moment::ReCalValue]\tTry to open" << myString.str().c_str() << endl;

    TFile *Eff = new TFile(myString.str().c_str(), "READ");
    RooWorkspace *w = (RooWorkspace *) Eff->Get("ws");
    RooAbsPdf *EffPdf = (RooAbsPdf *) w->function("projectedFunc");
    RooArgSet *param = EffPdf->getVariables();
   

    // ##################
    // # Read Data & MC #
    // ##################
    cout << "\n[Moment::CalValue]\t @@@ Making datasets @@@ " << endl;
    TFile *f = new TFile("2016MC_RECO_p1p2p3_newtag_LMNR_addW_add4BDT_addvars_bestBDTv4.root");
    TTree *t = (TTree *) f->Get("ntuple");
    t->SetBranchAddress("cos_theta_k", &cos_theta_k, &b_cos_theta_k);
    t->SetBranchAddress("cos_theta_l", &cos_theta_l, &b_cos_theta_l);
    t->SetBranchAddress("phi_kst_mumu", &phi_kst_mumu, &b_phi_kst_mumu);
    t->SetBranchAddress("mumuMass", &mumuMass, &b_mumuMass);
    t->SetBranchAddress("tagB0", &tagB0, &b_tagB0);
    t->SetBranchAddress("genSignal", &genSignal, &b_genSignal);

    Int_t entries = (Int_t) t->GetEntries();
    cout << "\n[Moment::CalValue]\tTotal number of events in the tree: " << entries << " @@@" << endl;
    for (Int_t i = 0; i < entries; i++) {
        t->GetEntry(i);
        if (((TagType == "good") && ((tagB0 == 1 && genSignal == 1) || (tagB0 == 0 && genSignal == 2))) ||
            ((TagType == "mis") && (!((tagB0 == 1 && genSignal == 1) || (tagB0 == 0 && genSignal == 2))))) {
            mumuMass2 = mumuMass * mumuMass;

            // ##################
            // # complete wight #
            // ##################
            param->setRealValue("ctK", cos_theta_k);
            param->setRealValue("ctL", cos_theta_l);
            param->setRealValue("phi", phi_kst_mumu);

            // ###############################
            // # define orthogonal functions #
            // ###############################
            f_1s = 1 - cos_theta_k * cos_theta_k;
            M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
            M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
            f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
            f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);
            if (EffPdf->getValV() < 0.0)  
            cout  << "Eff = " << EffPdf->getValV() << " " << cos_theta_k << "" << cos_theta_l << "" << phi_kst_mumu<< endl;
            if (mumuMass2 > q2Min && mumuMass2 < q2Max) {
                wight = std::max(1 / (EffPdf->getVal()),0.);
                f1s->Fill(f_1s, wight);
                m6s->Fill(M_6s, wight);
                m6c->Fill(M_6c, wight);
                f3->Fill(f_3, wight);
                f4->Fill(f_4, wight);
                f5->Fill(f_5, wight);
                f7->Fill(f_7, wight);
                f8->Fill(f_8, wight);
                f9->Fill(f_9, wight);
       
             }
        }
    }

    if  (Paratype == "FlS")
      {
        cout << f1s->GetMean() << endl;
        value = 2.0 - 2.5 * f1s->GetMean();
        error = 2.5 * f1s->GetRMS() /TMath::Sqrt(f1s->GetEntries() - 1);
      }
    else if  (Paratype == "AFBS")
      {
        value = 3.0/ 4.0 * (3.0 * m6s->GetMean() - 2.0 * m6c->GetMean());
        error = 3.0 / 4.0 * TMath::Sqrt(4.0 * TMath::Power(m6c->GetRMS(), 2) / (m6c->GetEntries() - 1) +
                                        9.0 * TMath::Power(m6s->GetRMS(), 2) / (m6s->GetEntries() - 1));
      }
    else if  (Paratype == "S3S")
      {
        value = 25.0 / 8.0 * f3->GetMean();
        error = 25.0 / 8.0 * f3->GetRMS() /TMath::Sqrt(f3->GetEntries() - 1);
      }
    else if  (Paratype == "S4S")
      {
        value = 25.0 / 8.0 * f4->GetMean();
        error = 25.0 / 8.0 * f4->GetRMS() /TMath::Sqrt(f4->GetEntries() - 1);
      }
    else if (Paratype == "S5S")
      {
        value = 2.5 * f5->GetMean();
        error = 2.5 * f5->GetRMS() /TMath::Sqrt(f5->GetEntries() - 1);
      }
    else if (Paratype == "S7S")
      {
        value = 2.5 * f7->GetMean();
        error = 2.5 * f7->GetRMS() /TMath::Sqrt(f7->GetEntries() - 1);
      }
    else if  (Paratype == "S8S")
      {
        value = 25.0 / 8.0 * f8->GetMean();
        error = 25.0 / 8.0 * f8->GetRMS() /TMath::Sqrt(f8->GetEntries() - 1);
      }
    else if  (Paratype == "S9S")
      {
        value = 25.0 / 8.0 * f9->GetMean();
        error = 25.0 / 8.0 * f9->GetRMS() /TMath::Sqrt(f9->GetEntries() - 1);
      }
    else if  (Paratype == "P1S")
      {
       value = 25./4. * f3->GetMean() / ( 1. - (2.0 - 2.5 * f1s->GetMean()) );
      }
    else if  (Paratype == "P2S")
      {
       value = 0.5 * m6s->GetMean() / ( 1. - (2.0 - 2.5 * f1s->GetMean()) );
      }
    else if (Paratype == "P3S")
      {
       value = - 25.0/8.0 * f9->GetMean() / (1.0 - (2.0 - 2.5 * f1s->GetMean()));
      }
    else if (Paratype == "P4pS")
      {
       value = 4 * 25/8 * f4->GetMean() /TMath::Sqrt((2.0 - 2.5 * f1s->GetMean()) * ( 1.0 - (2.0 - 2.5 * f1s->GetMean())) );
      }
    else if (Paratype == "P5pS")
      {
        double a = 2.0 - 2.5 * f1s->GetMean();
        value = 2.5 * f5->GetMean() /TMath::Sqrt(a * ( 1 - a) );
      }
    else if (Paratype == "P6pS")
      {
        double a = 2.0 - 2.5 * f1s->GetMean();
        value = - 2.5 * f7->GetMean() /TMath::Sqrt(a * ( 1 - a) );
       }
    else if (Paratype == "P8pS")
      {
        double a = 2.0 - 2.5 * f1s->GetMean();
        value = 25/8 * f8->GetMean() /TMath::Sqrt(a * ( 1 - a) );
      }
 cout << "\n[Moment::RecoCalValue]\t q2bin = " << q2Bin << ":" << Paratype << " = " << value << "+/- " << error  << endl;


}

int main(int argc, char **argv) {

    string Paratype = "";
    string TagType = "";
    string SampleType = "";

    // ###################################
    // # Check that Parameter is correct #
    // ###################################
    Paratype = argv[1];
    if ((Paratype != "FlS") && (Paratype != "AFBS") && (Paratype != "P1S") && (Paratype != "P2S") &&
        (Paratype != "P3S") && (Paratype != "P4pS") && (Paratype != "P5pS") && (Paratype != "P6pS") &&
        (Paratype != "P8pS") && (Paratype != "S3S") && (Paratype != "S4S") && (Paratype != "S5S") &&
        (Paratype != "S6S") && (Paratype != "S7S") && (Paratype != "S8S") && (Paratype != "S9S")) {
        cout << "\n[Moment::main]\tIncorrect Parameters\t" << Paratype << endl;
    }
    cout << "\n[Moment::main]\tParameter Type = \t" << Paratype << endl;

    // ####################
    // # select Data & MC #
    // ####################
    SampleType = argv[2];
    if ((SampleType != "gen") && (SampleType != "reco")) {
        cout << "\n[Moment::main]\tIncorrect Sample Type\t" << SampleType << endl;
    }
    cout << "\n[Moment::main]\tSample Type = \t" << SampleType << endl;
    if (SampleType == "reco")
      {
        // ##################
        // # Check Tag Type #
        // ##################
        TagType = argv[3];
        if ((TagType != "good") && (TagType != "mis")) {
        cout << "\n[Moment::main]\tIncorrect Events Tag Type\t" << TagType << endl;
        }
        cout << "\n[Moment::main]\tEvents Tag Type = \t" << TagType << endl;
      }

    // ###########################
    // # calculate final results #
    // ###########################
    for (q2Bin = 0; q2Bin < 9; q2Bin++) {
        quzhi();
        if (SampleType == "gen") 
            GenCalValue(q2Bin, q2Min, q2Max, Paratype);
        else if (SampleType == "reco") 
            ReCalValue(q2Bin, q2Min, q2Max, TagType, Paratype);
    }
}
