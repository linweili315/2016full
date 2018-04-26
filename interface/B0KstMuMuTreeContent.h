#ifndef B0KSTMUMUTREECONTENT_H
#define B0KSTMUMUTREECONTENT_H

#include "TTree.h"


class B0KstMuMuTreeContent
{
 public:
  
  B0KstMuMuTreeContent ();

  void Init ();
  void ClearNTuple ();
  void MakeTreeBranches (TTree* theTree);
  void SetBranchAddresses (TTree* theTree);
  void CopyCandidate (B0KstMuMuTreeContent* NTupleIn);

  // ################
  // # B0 parameter # 
  // ################
  float b0_id, b0_pt, b0_eta, b0_phi; 
 
  // ### Kst Parameters ###
  float kst_pt, kst_eta, kst_phi;

  // ### ki Parameters ###
  float k_pt, k_eta, k_phi;
 
  // ### pi Parameters ###
  float pi_pt, pi_eta, pi_phi;

  // ### mu- Parameters ###
  float mum_pt, mum_eta, mum_phi;

  // ### mu+ Parameters ###
  float mup_pt, mup_eta, mup_phi;

  // ### dimuon mass ###
  float mumu_mass;

  // ### angular parameters ###
  float gen_cos_theta_l, gen_cos_theta_k, gen_phi;

 private:

  void CopyScalars (B0KstMuMuTreeContent* NTupleIn);
};

#endif
