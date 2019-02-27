#ifndef B0KSTMUMUTREECONTENT_H
#define B0KSTMUMUTREECONTENT_H

#include "TTree.h"


class B0KstMuMuTreeContent
{
 public:
  
  B0KstMuMuTreeContent (int SampleType);

  void Init ( int SampleType );
  void ClearNTuple ( int SampleType );
  void MakeTreeBranches (TTree* theTree, int SampleType );
  void SetBranchAddresses (TTree* theTree, int SampleType );
  void CopyCandidate (B0KstMuMuTreeContent* NTupleIn, int SampleType);

  // #################
  // # gen parameter # 
  // #################
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

  // #############################
  // # reco MC & DATA parameters #
  // #############################
  
  float runN, eventN, recoVtxN, evWeight, evWeightE2, trueNumInteractionsMC, trig;

  // ### B0 parameters ###
  float bsX, bsY, bMass, bMassE, bBarMass, bBarMassE, bVtxCL, bCosAlphaBS, bCosAlphaBSE;
  float bLBS, bLBSE, bDCABS, bDCABSE, bctauPVBS, bctauPVBSE;

  // ### Kst parameters ###
  float kstMass, kstMassE, kstBarMass, kstBarMassE, kstVtxCL, kkMass;
  float kstTrkmCL, kstTrkmDCABS, kstTrkmDCABSE, kstTrkmdxyVtx, kstTrkmdzVtx, kstTrkmMinIP2D, kstTrkmMinIP2DE;
  float kstTrkpCL, kstTrkpDCABS, kstTrkpDCABSE, kstTrkpdxyVtx, kstTrkpdzVtx, kstTrkpMinIP2D, kstTrkpMinIP2DE;

  // ### dimuon parameters ###
  float mumuMass, mumuMassE, mumuVtxCL, mumuCosAlphaBS, mumuCosAlphaBSE, mumuLBS, mumuLBSE, mumuDCA;

  // ### mu- parameters ###
  float mumNormChi2, mumdxyVtx, mumdzVtx, mumMinIP2D, mumMinIP2DE, mumGlobalMuon, mumTrackerMuon, mumStandAloneMuon, mumTMOneStationTight;

  // ### mu+ parameters ###
  float mupNormChi2, mupdxyVtx, mupdzVtx, mupMinIP2D, mupMinIP2DE, mupGlobalMuon, mupTrackerMuon, mupStandAloneMuon, mupTMOneStationTight;

  float tagB0, truthMatchMum, truthMatchMup, truthMatchTrkm, truthMatchTrkp, mumDeltaRwithMC, mupDeltaRwithMC, kstTrkpDeltaRwithMC, kstTrkmDeltaRwithMC;
  float genSignal;
  float kstTrkmGlobalMuon, kstTrkmTrackerMuon, kstTrkmTMOneStationTight, kstTrkmTMOneStationLoose, kstTrkmTrackerMuonArbitrated;
  float kstTrkpGlobalMuon, kstTrkpTrackerMuon, kstTrkpTMOneStationTight, kstTrkpTMOneStationLoose, kstTrkpTrackerMuonArbitrated;

  float bPt, kstPt, mumuPt, mumPt, mupPt, kstTrkmPt, kstTrkpPt;
  float bPhi, kstPhi, mumuPhi, mumPhi, mupPhi, kstTrkmPhi, kstTrkpPhi;
  float bEta, kstEta, mumuEta, mumEta, mupEta, kstTrkmEta, kstTrkpEta;
  double iso_mum, iso_mup, iso_trkm, iso_trkp, sum_iso;
  float kstarmass;
  bool pass_preselection;
  float bdt_prob_v1, bdt_prob_v2, bdt_prob_v3, bdt_prob_v4;
  double cos_theta_l, cos_theta_k, phi_kst_mumu;
  float tagged_mass;

 private:

  void CopyScalars (B0KstMuMuTreeContent* NTupleIn, int SampleType);
};

#endif
