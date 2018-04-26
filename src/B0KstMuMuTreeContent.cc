#include "../interface/B0KstMuMuTreeContent.h"
#include <iostream>

B0KstMuMuTreeContent::B0KstMuMuTreeContent ()
{
  // ### B0 Parameter ### 
  b0_id    = 0;
  b0_pt    = 0;
  b0_eta   = 0;
  b0_phi   = 0;

  // ### Kst Parameter ###
  kst_pt   = 0;
  kst_eta  = 0;
  kst_phi  = 0;

  // ### kion parameter ###
  k_pt     = 0;
  k_eta    = 0;
  k_phi    = 0;

  // ### pi parameter ###
  pi_pt    = 0;
  pi_eta   = 0;
  pi_phi   = 0;

  // ### mu- parameter ###
  mum_pt   = 0;
  mum_eta  = 0;
  mum_phi  = 0;

  // ### mu+ parameter ###
  mup_pt   = 0;
  mup_eta  = 0;
  mup_phi  = 0;

  // ### dimuon mass ###
  mumu_mass = 0;

  // angular parameter ###
  gen_cos_theta_l  = 0;
  gen_cos_theta_k  = 0;
  gen_phi          = 0;
}

void B0KstMuMuTreeContent::Init ()
{
  mumu_mass = 0;

  gen_cos_theta_l  = 0;
  gen_cos_theta_k  = 0;
  gen_phi          = 0;
}
void B0KstMuMuTreeContent::ClearNTuple  ()
{

  // ### B0 Parameter ###
  b0_id    = 0;
  b0_pt    = 0;
  b0_eta   = 0;
  b0_phi   = 0;

  // ### Kst Parameter ###
  kst_pt   = 0;
  kst_eta  = 0;
  kst_phi  = 0;

  // ### kion parameter ###
  k_pt     = 0;
  k_eta    = 0;
  k_phi    = 0;

  // ### pi parameter ###
  pi_pt    = 0;
  pi_eta   = 0;
  pi_phi   = 0;

  // ### mu- parameter ###
  mum_pt   = 0;
  mum_eta  = 0;
  mum_phi  = 0;

  // ### mu+ parameter ###
  mup_pt   = 0;
  mup_eta  = 0;
  mup_phi  = 0;

  // ### dimuon mass ###
  mumu_mass = 0;

  // angular parameter ###
  gen_cos_theta_l  = 0;
  gen_cos_theta_k  = 0;
  gen_phi          = 0;
}

void B0KstMuMuTreeContent::MakeTreeBranches (TTree* theTree)
{
  // ### B0 parameter ###
  theTree->Branch("b0_id",    &b0_id,    "b0_id/F");
  theTree->Branch("b0_pt",    &b0_pt,    "b0_pt/F");
  theTree->Branch("b0_eta",   &b0_eta,   "b0_eta/F");
  theTree->Branch("b0_phi",   &b0_phi,   "b0_phi/F"); 

  // ### Kst parameter ###
  theTree->Branch("kst_pt",    &kst_pt,    "kst_pt/F");
  theTree->Branch("kst_eta",   &kst_eta,   "kst_eta/F");
  theTree->Branch("kst_phi",   &kst_phi,   "kst_phi/F");

  // ### ki parameter ###
  theTree->Branch("k_pt",    &k_pt,    "k_pt/F");
  theTree->Branch("k_eta",   &k_eta,   "k_eta/F");
  theTree->Branch("k_phi",   &k_phi,   "k_phi/F");

  // ### pi parameter ###
  theTree->Branch("pi_pt",    &pi_pt,    "pi_pt/F");
  theTree->Branch("pi_eta",   &pi_eta,   "pi_eta/F");
  theTree->Branch("pi_phi",   &pi_phi,   "pi_phi/F");

  // ### mu- papameter ###
  theTree->Branch("mum_pt",    &mum_pt,    "mum_pt/F");
  theTree->Branch("mum_eta",   &mum_eta,   "mum_eta/F");
  theTree->Branch("mum_phi",   &mum_phi,   "mum_phi/F");

  // ### mu+ parameter ###
  theTree->Branch("mup_pt",    &mup_pt,    "mup_pt/F");
  theTree->Branch("mup_eta",   &mup_eta,   "mup_eta/F");
  theTree->Branch("mup_phi",   &mup_phi,   "mup_phi/F");

  // ### dumion mass ###
  theTree->Branch("mumu_mass",    &mumu_mass,    "mumu_mass/F");

  // ### angular parameter ###
  theTree->Branch("gen_cos_theta_l",    &gen_cos_theta_l,   "gen_cos_theta_l/F");
  theTree->Branch("gen_cos_theta_k",    &gen_cos_theta_k,   "gen_cos_theta_k/F");
  theTree->Branch("gen_phi",            &gen_phi,           "gen_phi/F");
}

void B0KstMuMuTreeContent::SetBranchAddresses (TTree* theTree)
{
  theTree->SetBranchAddress("b0_id", &b0_id);
  theTree->SetBranchAddress("b0_pt", &b0_pt);
  theTree->SetBranchAddress("b0_eta", &b0_eta);
  theTree->SetBranchAddress("b0_phi", &b0_phi);

  theTree->SetBranchAddress("kst_pt", &kst_pt);
  theTree->SetBranchAddress("kst_eta", &kst_eta);
  theTree->SetBranchAddress("kst_phi", &kst_phi);

  theTree->SetBranchAddress("k_pt", &k_pt);
  theTree->SetBranchAddress("k_eta", &k_eta);
  theTree->SetBranchAddress("k_phi", &k_phi);

  theTree->SetBranchAddress("pi_pt", &pi_pt);
  theTree->SetBranchAddress("pi_eta", &pi_eta);
  theTree->SetBranchAddress("pi_phi", &pi_phi);

  theTree->SetBranchAddress("mum_pt", &mum_pt);
  theTree->SetBranchAddress("mum_eta", &mum_eta);
  theTree->SetBranchAddress("mum_phi", &mum_phi);
 
  theTree->SetBranchAddress("mup_pt", &mup_pt);
  theTree->SetBranchAddress("mup_eta", &mup_eta);
  theTree->SetBranchAddress("mup_phi", &mup_phi);

  theTree->SetBranchAddress("mumu_mass", &mumu_mass);

  theTree->SetBranchAddress("gen_cos_theta_l", &gen_cos_theta_l);
  theTree->SetBranchAddress("gen_cos_theta_k", &gen_cos_theta_k);
  theTree->SetBranchAddress("gen_phi", &gen_phi);
}

void B0KstMuMuTreeContent::CopyCandidate (B0KstMuMuTreeContent* NTupleIn)
{ 
  CopyScalars(NTupleIn);
}

void B0KstMuMuTreeContent::CopyScalars (B0KstMuMuTreeContent* NTupleIn)
{
  b0_id               = NTupleIn->b0_id;
  b0_pt               = NTupleIn->b0_pt;
  b0_eta              = NTupleIn->b0_eta;
  b0_phi              = NTupleIn->b0_phi;
 
  kst_pt              = NTupleIn->kst_pt;
  kst_eta             = NTupleIn->kst_eta;
  kst_phi             = NTupleIn->kst_phi;

  k_pt                = NTupleIn->k_pt;
  k_eta               = NTupleIn->k_eta;
  k_phi               = NTupleIn->k_phi;

  pi_pt               = NTupleIn->pi_pt;
  pi_eta              = NTupleIn->pi_eta;
  pi_phi              = NTupleIn->pi_phi;

  mum_pt              = NTupleIn->mum_pt;
  mum_eta             = NTupleIn->mum_eta;
  mum_phi             = NTupleIn->mum_phi;

  mup_pt              = NTupleIn->mup_pt;
  mup_eta             = NTupleIn->mup_eta;
  mup_phi             = NTupleIn->mup_phi;

  mumu_mass           = NTupleIn->mumu_mass;

  gen_cos_theta_l     = NTupleIn->gen_cos_theta_l;
  gen_cos_theta_k     = NTupleIn->gen_cos_theta_k;
  gen_phi             = NTupleIn->gen_phi;

}






