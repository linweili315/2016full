#include "../interface/B0KstMuMuTreeContent.h"
#include <iostream>

// ###################################
// # Gen-level MC, SampleType == 0   #
// # Reco MC & DATA, SampleType == 1 #
// ###################################


B0KstMuMuTreeContent::B0KstMuMuTreeContent ( int SampleType )
{
  if (SampleType == 0)
   {
    // ########################
    // # Gen-Level parameters #
    // ########################

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
  else if (SampleType == 1)
    {
     // #############################
     // # reco MC & DATA parameters #
     // #############################
     runN = 0;
     eventN = 0;
     recoVtxN = 0;
     evWeight = 0;
     evWeightE2 = 0;
     trig = 0;
     bsX = 0;
     bsY = 0;
     bMass = 0;
     bMassE = 0;
     bBarMass = 0;
     bBarMassE = 0;
     bVtxCL = 0;
     bCosAlphaBS = 0;
     bCosAlphaBSE = 0;
     bLBS = 0;
     bLBSE = 0;
     bDCABS = 0;
     bDCABSE = 0;
     bctauPVBS = 0;
     bctauPVBSE = 0;
     kstMass = 0;
     kstMassE = 0;
     kstBarMass = 0;
     kstBarMassE = 0;
     kstVtxCL = 0;
     kkMass = 0;
     mumuMass = 0;
     mumuMassE = 0;
     mumuVtxCL = 0;
     mumuCosAlphaBS = 0;
     mumuCosAlphaBSE = 0;
     mumuLBS = 0;
     mumuLBSE = 0;

     mumuDCA = 0;
     mumNormChi2 = 0;
     mumdxyVtx = 0;
     mumdzVtx = 0;
     mumMinIP2D = 0;
     mumMinIP2DE = 0;
     mupNormChi2 = 0;
     mupdxyVtx = 0;
     mupdzVtx = 0;
     mupMinIP2D = 0;
     mupMinIP2DE = 0;

     kstTrkmCL = 0;
     kstTrkmDCABS = 0;
     kstTrkmDCABSE = 0;
     kstTrkmdxyVtx = 0;
     kstTrkmdzVtx = 0;
     kstTrkmMinIP2D = 0;
     kstTrkmMinIP2DE = 0;
     kstTrkpCL = 0;
     kstTrkpDCABS = 0;
     kstTrkpDCABSE = 0;
     kstTrkpdxyVtx = 0;
     kstTrkpdzVtx = 0;
     kstTrkpMinIP2D = 0;
     kstTrkpMinIP2DE = 0;

     tagB0 = 0;
     genSignal = 0;
     mumGlobalMuon = 0;
     mumTrackerMuon = 0;
     mumStandAloneMuon = 0;
     mumTMOneStationTight = 0;
     mupGlobalMuon = 0;
     mupTrackerMuon = 0;
     mupStandAloneMuon = 0;
     mupTMOneStationTight = 0;

     kstTrkmGlobalMuon = 0;
     kstTrkmTrackerMuon = 0;
     kstTrkmTMOneStationTight = 0;
     kstTrkmTMOneStationLoose = 0;
     kstTrkmTrackerMuonArbitrated = 0;
     kstTrkpGlobalMuon = 0;
     kstTrkpTrackerMuon = 0;
     kstTrkpTMOneStationTight = 0;
     kstTrkpTMOneStationLoose = 0;
     kstTrkpTrackerMuonArbitrated = 0;
  
     bPt = 0;
     kstPt = 0;
     mumuPt = 0;
     mumPt = 0;
     mupPt = 0;
     kstTrkmPt = 0;
     kstTrkpPt = 0;
     bPhi = 0;
     kstPhi = 0;
     mumuPhi = 0;
     mumPhi = 0;
     mupPhi = 0;
     kstTrkmPhi = 0;
     kstTrkpPhi = 0;
     bEta = 0;
     kstEta = 0;
     mumuEta = 0;
     mumEta = 0;
     mupEta = 0;
     kstTrkmEta = 0;
     kstTrkpEta = 0;
     iso_mum = 0;
     iso_mup = 0;
     iso_trkm = 0;
     iso_trkp = 0;
     sum_iso = 0;
     kstarmass = 0;
     pass_preselection = NULL;
     bdt_prob_v1 = 0;
     bdt_prob_v2 = 0;
     bdt_prob_v3 = 0;
     bdt_prob_v4 = 0;
     cos_theta_l = 0;
     cos_theta_k = 0;
     phi_kst_mumu = 0;
     tagged_mass = 0;
   }
}

void B0KstMuMuTreeContent::Init ( int SampleType )
{
 if (SampleType == 0)
   {
     // ########################
     // # Gen-level parameters #
     // ########################
     mumu_mass = 0;

     gen_cos_theta_l  = 0;
     gen_cos_theta_k  = 0;
     gen_phi          = 0;
    }
  else if (SampleType == 1 )
    {
      // #############################
      // # Reco MC & DATA parameters #
      // #############################
      // ### B0 ###
      bsX = 0;
      bsY = 0;
      bMass = 0;
      bMassE = 0;
      bBarMass = 0;
      bBarMassE = 0;

      bPt = 0;
      bPhi = 0;
      bEta = 0;
 
      // ### Kst ###
      kstMass = 0;
      kstMassE = 0;
      kstBarMass = 0;
      kstBarMassE = 0;
   
      kstPt = 0;
      kstPhi = 0;
      kstEta = 0;

      kkMass = 0;
 
      // ### dimuon ###
      mumuMass = 0;
      mumuMassE = 0;
  
      mumuPt = 0;
      mumuPhi = 0;
      mumuEta = 0;
 
      // ### mum ###
      mumPt = 0;
      mumPhi = 0;
      mumEta = 0;

      // ### mup ###
      mupPt = 0;
      mupPhi = 0;
      mupEta = 0; 

      // ### tag information ###
      tagB0 = 0;
      genSignal = 0;

       // ### Three angulars ###
      cos_theta_l = 0;
      cos_theta_k = 0;
      phi_kst_mumu = 0;

      kstarmass = 0;
      tagged_mass = 0; 
  }

}
void B0KstMuMuTreeContent::ClearNTuple  (int SampleType )
{
  if (SampleType == 0)
    {
       // ########################
       // # Gen-Level Parameters #
       // ########################
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
  
  // #############################
  // # reco MC & DATA Parameters #
  // #############################
  else if (SampleType ==1)
    {
       // ### events information ###
       runN = 0;
       eventN = 0;
       recoVtxN = 0;
       evWeight = 0;
       evWeightE2 = 0;
       trig = 0;
 
       // ### B0 information ###
       bsX = 0;
       bsY = 0;
       bMass = 0;
       bMassE = 0;
       bBarMass = 0;
       bBarMassE = 0;
       bVtxCL = 0;
       bCosAlphaBS = 0;
       bCosAlphaBSE = 0;
       bLBS = 0;
       bLBSE = 0;
       bDCABS = 0;
       bDCABSE = 0;
       bctauPVBS = 0;
       bctauPVBSE = 0;
  
       tagB0 = 0;
       genSignal = 0; 
      // ### Kst information ###
       kstMass = 0;
       kstMassE = 0;
       kstBarMass = 0;
       kstBarMassE = 0;
       kstVtxCL = 0;
       kkMass = 0;

       kstTrkmCL = 0;
       kstTrkmDCABS = 0;
       kstTrkmDCABSE = 0;
       kstTrkmdxyVtx = 0;
       kstTrkmdzVtx = 0;
       kstTrkmMinIP2D = 0;
       kstTrkmMinIP2DE = 0;
       kstTrkpCL = 0;
       kstTrkpDCABS = 0;
       kstTrkpDCABSE = 0;
       kstTrkpdxyVtx = 0;
       kstTrkpdzVtx = 0;
       kstTrkpMinIP2D = 0;
       kstTrkpMinIP2DE = 0;

       // ### dimuon information ###
       mumuMass = 0;
       mumuMassE = 0;
       mumuVtxCL = 0;
       mumuCosAlphaBS = 0;
       mumuCosAlphaBSE = 0;
       mumuLBS = 0;
       mumuLBSE = 0;
       mumuDCA = 0;

       // ### mum information ###
       mumNormChi2 = 0;
       mumdxyVtx = 0;
       mumdzVtx = 0;
       mumMinIP2D = 0;
       mumMinIP2DE = 0;
 
       mumGlobalMuon = 0;
       mumTrackerMuon = 0;
       mumStandAloneMuon = 0;
       mumTMOneStationTight = 0;

       // ### mup information ###
       mupNormChi2 = 0;
       mupdxyVtx = 0;
       mupdzVtx = 0;
       mupMinIP2D = 0;
       mupMinIP2DE = 0;

       mupGlobalMuon = 0;
       mupTrackerMuon = 0;
       mupStandAloneMuon = 0;
       mupTMOneStationTight = 0;

       kstTrkmGlobalMuon = 0;
       kstTrkmTrackerMuon = 0;
       kstTrkmTMOneStationTight = 0;
       kstTrkmTMOneStationLoose = 0;
       kstTrkmTrackerMuonArbitrated = 0;
       kstTrkpGlobalMuon = 0;
       kstTrkpTrackerMuon = 0;
       kstTrkpTMOneStationTight = 0;
       kstTrkpTMOneStationLoose = 0;
       kstTrkpTrackerMuonArbitrated = 0;

       bPt = 0;
       kstPt = 0;
       mumuPt = 0;
       mumPt = 0;
       mupPt = 0;
       kstTrkmPt = 0;
       kstTrkpPt = 0;
       bPhi = 0;
       kstPhi = 0;
       mumuPhi = 0;
       mumPhi = 0;
       mupPhi = 0;
       kstTrkmPhi = 0;
       kstTrkpPhi = 0;
       bEta = 0;
       kstEta = 0;
       mumuEta = 0;
       mumEta = 0;
       mupEta = 0;
       kstTrkmEta = 0;
       kstTrkpEta = 0;
       iso_mum = 0;
       iso_mup = 0;
       iso_trkm = 0;
       iso_trkp = 0;
       sum_iso = 0;
       kstarmass = 0;
       pass_preselection = NULL;
       bdt_prob_v1 = 0;
       bdt_prob_v2 = 0;
       bdt_prob_v3 = 0;
       bdt_prob_v4 = 0;
       cos_theta_l = 0;
       cos_theta_k = 0;
       phi_kst_mumu = 0;
       tagged_mass = 0;
   }
}

void B0KstMuMuTreeContent::MakeTreeBranches (TTree* theTree, int SampleType)
{
 // ##########
 // # Gen MC #
 // ########## 
 if (SampleType == 0)
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
 // ##################
 // # Reco MC & DATA #
 // ##################
 else if (SampleType == 1)
  {
   theTree->Branch("runN", &runN, "runN/F");
   theTree->Branch("eventN", &eventN, "eventN/F");
   theTree->Branch("recoVtxN", &recoVtxN, "recoVtxN/F");
   theTree->Branch("evWeight", &evWeight, "evWeight/F");
   theTree->Branch("evWeightE2", &evWeightE2, "evWeightE2/F");
   theTree->Branch("trig", &trig, "trig/F");
   theTree->Branch("bsX", &bsX, "bsX/F");
   theTree->Branch("bsY", &bsY, "bsY/F");
   theTree->Branch("bMass", &bMass, "bMass/F");
   theTree->Branch("bMassE", &bMassE,"bMassE/F");
   theTree->Branch("bBarMass", &bBarMass, "bBarMass/F");
   theTree->Branch("bBarMassE", &bBarMassE, "bBarMassE/F");
   theTree->Branch("bVtxCL", &bVtxCL, "bVtxCL/F");
   theTree->Branch("bCosAlphaBS", &bCosAlphaBS, "bCosAlphaBS/F");
   theTree->Branch("bCosAlphaBSE", &bCosAlphaBSE, "bCosAlphaBSE/F");
   theTree->Branch("bLBS", &bLBS, "bLBS/F");
   theTree->Branch("bLBSE", &bLBSE, "bLBSE/F");
   theTree->Branch("bDCABS", &bDCABS, "bDCABS/F");
   theTree->Branch("bDCABSE", &bDCABSE, "bDCABSE/F");
   theTree->Branch("bctauPVBS", &bctauPVBS, "bctauPVBS/F");
   theTree->Branch("bctauPVBSE", &bctauPVBSE, "bctauPVBSE/F");
   theTree->Branch("kstMass", &kstMass, "kstMass/F");
   theTree->Branch("kstMassE", &kstMassE, "kstMassE/F");
   theTree->Branch("kstBarMass", &kstBarMass, "kstBarMass/F");
   theTree->Branch("kstBarMassE", &kstBarMassE, "kstBarMassE/F");
   theTree->Branch("kstVtxCL", &kstVtxCL, "kstVtxCL/F");
   theTree->Branch("kkMass", &kkMass, "kkMass/F");
   theTree->Branch("mumuMass", &mumuMass, "mumuMass/F");
   theTree->Branch("mumuMassE", &mumuMassE, "mumuMassE/F");

   theTree->Branch("mumuVtxCL", &mumuVtxCL, "mumuVtxCL/F");
   theTree->Branch("mumuCosAlphaBS", &mumuCosAlphaBS, "mumuCosAlphaBS/F");
   theTree->Branch("mumuCosAlphaBSE", &mumuCosAlphaBSE, "mumuCosAlphaBSE/F");
   theTree->Branch("mumuLBS", &mumuLBS, "mumuLBS/F");
   theTree->Branch("mumuLBSE", &mumuLBSE, "mumuLBSE/F");
   theTree->Branch("mumuDCA", &mumuDCA, "mumuDCA/F");
   theTree->Branch("mumNormChi2", &mumNormChi2, "mumNormChi2/F");
   theTree->Branch("mumdxyVtx", &mumdxyVtx, "mumdxyVtx/F");
   theTree->Branch("mumdzVtx", &mumdzVtx, "mumdzVtx/F");
   theTree->Branch("mumMinIP2D", &mumMinIP2D, "mumMinIP2D/F");
   theTree->Branch("mumMinIP2DE", &mumMinIP2DE, "mumMinIP2DE/F");
   theTree->Branch("mupNormChi2", &mupNormChi2, "mupNormChi2/F");
   theTree->Branch("mupdxyVtx", &mupdxyVtx, "mupdxyVtx/F");
   theTree->Branch("mupdzVtx", &mupdzVtx, "mupdzVtx/F");
   theTree->Branch("mupMinIP2D", &mupMinIP2D, "mupMinIP2D/F");
   theTree->Branch("mumMinIP2DE", &mumMinIP2DE, "mupMinIP2DE/F");

   theTree->Branch("mupNormChi2", &mupNormChi2, "mupNormChi2/F");
   theTree->Branch("mupdxyVtx", &mupdxyVtx, "mupdxyVtx/F");
   theTree->Branch("mupdzVtx", &mupdzVtx, "mupdzVtx/F");
   theTree->Branch("mupMinIP2D", &mupMinIP2D, "mupMinIP2D/F");
   theTree->Branch("mupMinIP2DE", &mupMinIP2DE, "mupMinIP2DE/F");
   theTree->Branch("kstTrkmCL", &kstTrkmCL, "kstTrkmCL/F");
   theTree->Branch("kstTrkmDCABS", &kstTrkmDCABS, "kstTrkmDCABS/F");
   theTree->Branch("kstTrkmDCABSE", &kstTrkmDCABSE, "kstTrkmDCABSE/F");
   theTree->Branch("kstTrkmdxyVtx", &kstTrkmdxyVtx, "kstTrkmdxyVtx/F");
   theTree->Branch("kstTrkmdzVtx", &kstTrkmdzVtx, "kstTrkmdzVtx/F");
   theTree->Branch("kstTrkmMinIP2D", &kstTrkmMinIP2D, "kstTrkmMinIP2D/F");
   theTree->Branch("kstTrkmMinIP2DE", &kstTrkmMinIP2DE, "kstTrkmMinIP2DE/F");
   theTree->Branch("kstTrkpCL", &kstTrkpCL, "kstTrkpCL/F");
   theTree->Branch("kstTrkpDCABS", &kstTrkpDCABS, "kstTrkpDCABS/F");
   theTree->Branch("kstTrkpDCABSE", &kstTrkpDCABSE, "kstTrkpDCABSE/F");
   theTree->Branch("kstTrkpdxyVtx", &kstTrkpdxyVtx, "kstTrkpdxyVtx/F");
   theTree->Branch("kstTrkpdzVtx", &kstTrkpdzVtx, "kstTrkpdzVtx/F");
   theTree->Branch("kstTrkpMinIP2D", &kstTrkpMinIP2D, "kstTrkpMinIP2D/F");
   theTree->Branch("kstTrkpMinIP2DE", &kstTrkpMinIP2DE, "kstTrkpMinIP2DE/F");
   theTree->Branch("tagB0", &tagB0, "tagB0/F");
   theTree->Branch("genSiganl", &genSignal, "genSiganl/F");
   theTree->Branch("mumGlobalMuon", &mumGlobalMuon, "mumGlobalMuon/F");
   theTree->Branch("mumTrackerMuon", &mumTrackerMuon, "mumTrackerMuon/F");
   theTree->Branch("mumStandAloneMuon", &mumStandAloneMuon, "mumStandAloneMuon/F");
   theTree->Branch("mumTMOneStationTight", &mumTMOneStationTight, "mumTMOneStationTight/F");
   theTree->Branch("mupGlobalMuon", &mupGlobalMuon, "mupGlobalMuon/F");
   theTree->Branch("mupTrackerMuon", &mupTrackerMuon, "mupTrackerMuon/F");
   theTree->Branch("mupStandAloneMuon", &mupStandAloneMuon, "mupStandAloneMuon/F");
   theTree->Branch("mupTMOneStationTight", &mupTMOneStationTight, "mupTMOneStationTight/F");
   theTree->Branch("kstTrkmGlobalMuon", &kstTrkmGlobalMuon, "kstTrkmGlobalMuon/F");
   theTree->Branch("kstTrkmTrackerMuon", &kstTrkmTrackerMuon, "kstTrkmTrackerMuon/F");
   theTree->Branch("kstTrkmTMOneStationTight", &kstTrkmTMOneStationTight, "kstTrkmTMOneStationTight/F");
   theTree->Branch("kstTrkmTMOneStationLoose", &kstTrkmTMOneStationLoose, "kstTrkmTMOneStationLoose/F");
   theTree->Branch("kstTrkmTrackerMuonArbitrated", &kstTrkmTrackerMuonArbitrated, "kstTrkmTrackerMuonArbitrated/F");

   theTree->Branch("kstTrkpGlobalMuon", &kstTrkpGlobalMuon, "kstTrkpGlobalMuon/F");
   theTree->Branch("kstTrkpTrackerMuon", &kstTrkpTrackerMuon, "kstTrkpTrackerMuon/F");
   theTree->Branch("kstTrkpTMOneStationTight", &kstTrkpTMOneStationTight, "kstTrkpTMOneStationTight/F");
   theTree->Branch("kstTrkpTMOneStationLoose", &kstTrkpTMOneStationLoose, "kstTrkpTMOneStationLoose/F");
   theTree->Branch("kstTrkpTrackerMuonArbitrated", &kstTrkpTrackerMuonArbitrated, "kstTrkpTrackerMuonArbitrated/F");
   theTree->Branch("bPt", &bPt, "bPt/F");
   theTree->Branch("kstPt", &kstPt, "kstPt/F");
   theTree->Branch("mumuPt", &mumuPt, "mumuPt/F");
   theTree->Branch("mumPt", &mumPt, "mumPt/F");
   theTree->Branch("mupPt", &mupPt, "mupPt/F");
   theTree->Branch("kstTrkmPt", &kstTrkmPt, "kstTrkmPt/F");
   theTree->Branch("kstTrkpPt", &kstTrkpPt, "kstTrkpPt/F");

   theTree->Branch("bPhi", &bPhi, "bPhi/F");
   theTree->Branch("kstPhi", &kstPhi, "kstPhi/F");
   theTree->Branch("mumuPhi", &mumuPhi, "mumuPhi/F");
   theTree->Branch("mumPhi", &mumPhi, "mumPhi/F");
   theTree->Branch("mupPhi", &mupPhi, "mupPhi/F");
   theTree->Branch("kstTrkmPhi", &kstTrkmPhi, "kstTrkmPhi/F");
   theTree->Branch("kstTrkpPhi", &kstTrkpPhi, "kstTrkpPhi/F");
   theTree->Branch("bEta", &bEta, "bEta/F");
   theTree->Branch("kstEta", &kstEta, "kstEta/F");
   theTree->Branch("mumuEta", &mumuEta, "mumuEta/F");
   theTree->Branch("mumEta", &mumEta, "mumEta/F");
   theTree->Branch("mupEta", &mupEta, "mupEta/F");
   theTree->Branch("kstTrkmEta", &kstTrkmEta, "kstTrkmEta/F");
   theTree->Branch("kstTrkpEta", &kstTrkpEta, "kstTrkpEta/F");

   theTree->Branch("iso_mum", &iso_mum, "iso_mum/F");
   theTree->Branch("iso_mup", &iso_mup, "iso_mup/F");
   theTree->Branch("iso_trkm", &iso_trkm, "iso_trkm/F");
   theTree->Branch("iso_trkp", &iso_trkp, "iso_trkp/F");
   theTree->Branch("sum_iso", &sum_iso, "sum_iso/F");
   theTree->Branch("kstarmass", &kstarmass, "kstarmass/F");
   theTree->Branch("pass_preselection", &pass_preselection);
   theTree->Branch("bdt_prob_v1", &bdt_prob_v1, "bdt_prob_v1/F");
   theTree->Branch("bdt_prob_v2", &bdt_prob_v2, "bdt_prob_v2/F");
   theTree->Branch("bdt_prob_v3", &bdt_prob_v3, "bdt_prob_v3/F");
   theTree->Branch("bdt_prob_v4", &bdt_prob_v4, "bdt_prob_v4/F");
   theTree->Branch("cos_theta_l", &cos_theta_l, "cos_theta_l/D");
   theTree->Branch("cos_theta_k", &cos_theta_k, "cos_theta_k/D");
   theTree->Branch("phi_kst_mumu", &phi_kst_mumu, "phi_kst_mumu/D");
   theTree->Branch("tagged_mass", &tagged_mass, "tagged_mass/F");

  }
}

void B0KstMuMuTreeContent::SetBranchAddresses (TTree* theTree, int SampleType)
{
 // #############
 // # Gen-Level #
 // ############
 if (SampleType == 0)
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

  // ##################
  // # reco MC & DATA #
  // ##################
  else if (SampleType == 1)
   {
  theTree->SetBranchAddress("runN", &runN);
  theTree->SetBranchAddress("eventN", &eventN);
  theTree->SetBranchAddress("recoVtxN", &recoVtxN);
  theTree->SetBranchAddress("evWeight", &evWeight);
  theTree->SetBranchAddress("evWeightE2", &evWeightE2);
  theTree->SetBranchAddress("trig", &trig);
  theTree->SetBranchAddress("bsX", &bsX);
  theTree->SetBranchAddress("bsY", &bsY);
  theTree->SetBranchAddress("bMass", &bMass);
  theTree->SetBranchAddress("bMassE", &bMassE);
  theTree->SetBranchAddress("bBarMass", &bBarMass);
  theTree->SetBranchAddress("bBarMassE", &bBarMassE);
  theTree->SetBranchAddress("bVtxCL", &bVtxCL);
  theTree->SetBranchAddress("bCosAlphaBS", &bCosAlphaBS);
  theTree->SetBranchAddress("bCosAlphaBSE", &bCosAlphaBSE);
  theTree->SetBranchAddress("bLBS", &bLBS);
  theTree->SetBranchAddress("bLBSE", &bLBSE);
  theTree->SetBranchAddress("bDCABS", &bDCABS);
  theTree->SetBranchAddress("bDCABSE", &bDCABSE);
  theTree->SetBranchAddress("bctauPVBS", &bctauPVBS);
  theTree->SetBranchAddress("bctauPVBSE", &bctauPVBSE);
  theTree->SetBranchAddress("kstMass", &kstMass);
  theTree->SetBranchAddress("kstMassE", &kstMassE);
  theTree->SetBranchAddress("kstBarMass", &kstBarMass);
  theTree->SetBranchAddress("kstBarMassE", &kstBarMassE);
  theTree->SetBranchAddress("kstVtxCL", &kstVtxCL);
  theTree->SetBranchAddress("kkMass", &kkMass);
  theTree->SetBranchAddress("mumuMass", &mumuMass);
  theTree->SetBranchAddress("mumuMassE", &mumuMassE);
  theTree->SetBranchAddress("mumuVtxCL", &mumuVtxCL);
  theTree->SetBranchAddress("mumuCosAlphaBS", &mumuCosAlphaBS);
  theTree->SetBranchAddress("mumuCosAlphaBSE", &mumuCosAlphaBSE);
  theTree->SetBranchAddress("mumuLBS", &mumuLBS);
  theTree->SetBranchAddress("mumuLBSE", &mumuLBSE);
  theTree->SetBranchAddress("mumuDCA", &mumuDCA);
  theTree->SetBranchAddress("mumNormChi2", &mumNormChi2);
  theTree->SetBranchAddress("mumdxyVtx", &mumdxyVtx);
  theTree->SetBranchAddress("mumdzVtx", &mumdzVtx);
  theTree->SetBranchAddress("mumMinIP2D", &mumMinIP2D);
  theTree->SetBranchAddress("mumMinIP2DE", &mumMinIP2DE);
  theTree->SetBranchAddress("mupNormChi2", &mupNormChi2);
  theTree->SetBranchAddress("mupdxyVtx", &mupdxyVtx);
  theTree->SetBranchAddress("mupdzVtx", &mupdzVtx);
  theTree->SetBranchAddress("mupMinIP2D", &mupMinIP2D);
  theTree->SetBranchAddress("mupMinIP2DE", &mupMinIP2DE);
  theTree->SetBranchAddress("kstTrkmCL", &kstTrkmCL);
  theTree->SetBranchAddress("kstTrkmDCABS", &kstTrkmDCABS);
  theTree->SetBranchAddress("kstTrkmDCABSE", &kstTrkmDCABSE);
  theTree->SetBranchAddress("kstTrkmdxyVtx", &kstTrkmdxyVtx);
  theTree->SetBranchAddress("kstTrkmdzVtx", &kstTrkmdzVtx);
  theTree->SetBranchAddress("kstTrkmMinIP2D", &kstTrkmMinIP2D);
  theTree->SetBranchAddress("kstTrkmMinIP2DE", &kstTrkmMinIP2DE);
  theTree->SetBranchAddress("kstTrkpCL", &kstTrkpCL);
  theTree->SetBranchAddress("kstTrkpDCABS", &kstTrkpDCABS);
  theTree->SetBranchAddress("kstTrkpDCABSE", &kstTrkpDCABSE);
  theTree->SetBranchAddress("kstTrkpdxyVtx", &kstTrkpdxyVtx);
  theTree->SetBranchAddress("kstTrkpdzVtx", &kstTrkpdzVtx);
  theTree->SetBranchAddress("kstTrkpMinIP2D", &kstTrkpMinIP2D);
  theTree->SetBranchAddress("kstTrkpMinIP2DE", &kstTrkpMinIP2DE);
  theTree->SetBranchAddress("tagB0", &tagB0);
  theTree->SetBranchAddress("genSignal", &genSignal);
  theTree->SetBranchAddress("mumGlobalMuon", &mumGlobalMuon);
  theTree->SetBranchAddress("mumTrackerMuon", &mumTrackerMuon);
  theTree->SetBranchAddress("mumStandAloneMuon", &mumStandAloneMuon);
  theTree->SetBranchAddress("mumTMOneStationTight", &mumTMOneStationTight);
  theTree->SetBranchAddress("mupGlobalMuon", &mupGlobalMuon);
  theTree->SetBranchAddress("mupTrackerMuon", &mupTrackerMuon);
  theTree->SetBranchAddress("mupStandAloneMuon", &mupStandAloneMuon);
  theTree->SetBranchAddress("mupTMOneStationTight", &mupTMOneStationTight);
  theTree->SetBranchAddress("kstTrkmGlobalMuon", &kstTrkmGlobalMuon);
  theTree->SetBranchAddress("kstTrkmTrackerMuon", &kstTrkmTrackerMuon);
  theTree->SetBranchAddress("kstTrkmTMOneStationTight", &kstTrkmTMOneStationTight);
  theTree->SetBranchAddress("kstTrkmTMOneStationLoose", &kstTrkmTMOneStationLoose);
  theTree->SetBranchAddress("kstTrkmTrackerMuonArbitrated", &kstTrkmTrackerMuonArbitrated);
  theTree->SetBranchAddress("kstTrkpGlobalMuon", &kstTrkpGlobalMuon);
  theTree->SetBranchAddress("kstTrkpTrackerMuon", &kstTrkpTrackerMuon);
  theTree->SetBranchAddress("kstTrkpTMOneStationTight", &kstTrkpTMOneStationTight);
  theTree->SetBranchAddress("kstTrkpTMOneStationLoose", &kstTrkpTMOneStationLoose);
  theTree->SetBranchAddress("kstTrkpTrackerMuonArbitrated", &kstTrkpTrackerMuonArbitrated);
  theTree->SetBranchAddress("bPt", &bPt);
  theTree->SetBranchAddress("kstPt", &kstPt);
  theTree->SetBranchAddress("mumuPt", &mumuPt);
  theTree->SetBranchAddress("mumPt", &mumPt);
  theTree->SetBranchAddress("mupPt", &mupPt);
  theTree->SetBranchAddress("kstTrkmPt", &kstTrkmPt);
  theTree->SetBranchAddress("kstTrkpPt", &kstTrkpPt);
  theTree->SetBranchAddress("bPhi", &bPhi);
  theTree->SetBranchAddress("kstPhi", &kstPhi);
  theTree->SetBranchAddress("mumuPhi", &mumuPhi);
  theTree->SetBranchAddress("mumPhi", &mumPhi);
  theTree->SetBranchAddress("mupPhi", &mupPhi);
  theTree->SetBranchAddress("kstTrkmPhi", &kstTrkmPhi);
  theTree->SetBranchAddress("kstTrkpPhi", &kstTrkpPhi);
  theTree->SetBranchAddress("bEta", &bEta);
  theTree->SetBranchAddress("kstEta", &kstEta);
  theTree->SetBranchAddress("mumuEta", &mumuEta);
  theTree->SetBranchAddress("mumEta", &mumEta);
  theTree->SetBranchAddress("mupEta", &mupEta);
  theTree->SetBranchAddress("kstTrkmEta", &kstTrkmEta);
  theTree->SetBranchAddress("kstTrkpEta", &kstTrkpEta);
  theTree->SetBranchAddress("iso_mum", &iso_mum);
  theTree->SetBranchAddress("iso_mup", &iso_mup);
  theTree->SetBranchAddress("iso_trkm", &iso_trkm);
  theTree->SetBranchAddress("iso_trkp", &iso_trkp);
  theTree->SetBranchAddress("sum_iso", &sum_iso);
  theTree->SetBranchAddress("kstarmass", &kstarmass);
  theTree->SetBranchAddress("pass_preselection", &pass_preselection);
  theTree->SetBranchAddress("bdt_prob_v1", &bdt_prob_v1);
  theTree->SetBranchAddress("bdt_prob_v2", &bdt_prob_v2);
  theTree->SetBranchAddress("bdt_prob_v3", &bdt_prob_v3);
  theTree->SetBranchAddress("bdt_prob_v4", &bdt_prob_v4);
  theTree->SetBranchAddress("cos_theta_l", &cos_theta_l);
  theTree->SetBranchAddress("cos_theta_k", &cos_theta_k);
  theTree->SetBranchAddress("phi_kst_mumu", &phi_kst_mumu);
  theTree->SetBranchAddress("tagged_mass", &tagged_mass);
 } 
}
void B0KstMuMuTreeContent::CopyCandidate (B0KstMuMuTreeContent* NTupleIn, int SampleType)
{ 
  CopyScalars(NTupleIn, SampleType);
}

void B0KstMuMuTreeContent::CopyScalars (B0KstMuMuTreeContent* NTupleIn, int SampleType )
{
 // ################
 // # Gen-Level MC #
 // ################
 if (SampleType == 0)
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

  // ##################
  // # reco MC & DATA #
  // ##################
  else if (SampleType == 1)
   {
  
  runN                = NTupleIn->runN;
  eventN              = NTupleIn->eventN;
  recoVtxN            = NTupleIn->recoVtxN;
  evWeight            = NTupleIn->evWeight;
  evWeightE2          = NTupleIn->evWeightE2;
  trig                = NTupleIn->trig;
  bsX                 = NTupleIn->bsX;
  bsY                 = NTupleIn->bsY;
  bMass               = NTupleIn->bMass;
  bMassE              = NTupleIn->bMassE;
  bBarMass            = NTupleIn->bBarMass;
  bBarMassE           = NTupleIn->bBarMassE;
  bVtxCL              = NTupleIn->bVtxCL;
  bCosAlphaBS         = NTupleIn->bCosAlphaBS;
  bCosAlphaBSE        = NTupleIn->bCosAlphaBSE;
  bLBS                = NTupleIn->bLBS;
  bLBSE               = NTupleIn->bLBSE;
  bDCABS              = NTupleIn->bDCABS;
  bDCABSE             = NTupleIn->bDCABSE;
  bctauPVBS           = NTupleIn->bctauPVBS;
  bctauPVBSE          = NTupleIn->bctauPVBSE;
  kstMass             = NTupleIn->kstMass;
  kstMassE            = NTupleIn->kstMassE;
  kstBarMass          = NTupleIn->kstBarMass;
  kstBarMassE         = NTupleIn->kstBarMassE;
  kstVtxCL            = NTupleIn->kstVtxCL;
  kkMass              = NTupleIn->kkMass;
  mumuMass            = NTupleIn->mumuMass;
  mumuMassE           = NTupleIn->mumuMassE;
  mumuVtxCL           = NTupleIn->mumuVtxCL;
  mumuCosAlphaBS      = NTupleIn->mumuCosAlphaBS;
  mumuCosAlphaBSE     = NTupleIn->mumuCosAlphaBSE;
  mumuLBS             = NTupleIn->mumuLBS;
  mumuLBSE            = NTupleIn->mumuLBSE;

  mumuDCA             = NTupleIn->mumuDCA;
  mumNormChi2         = NTupleIn->mumNormChi2;
  mumdxyVtx           = NTupleIn->mumdxyVtx;
  mumdzVtx            = NTupleIn->mumdzVtx;
  mumMinIP2D          = NTupleIn->mumMinIP2D;
  mumMinIP2DE         = NTupleIn->mumMinIP2DE;
  mupNormChi2         = NTupleIn->mupNormChi2;
  mupdxyVtx           = NTupleIn->mupdxyVtx;
  mupdzVtx            = NTupleIn->mupdzVtx;
  mupMinIP2D          = NTupleIn->mupMinIP2D;
  mupMinIP2DE         = NTupleIn->mupMinIP2DE;

  kstTrkmCL           = NTupleIn->kstTrkmCL;
  kstTrkmDCABS        = NTupleIn->kstTrkmDCABS;
  kstTrkmDCABSE       = NTupleIn->kstTrkmDCABSE;
  kstTrkmdxyVtx       = NTupleIn->kstTrkmdxyVtx;
  kstTrkmdzVtx        = NTupleIn->kstTrkmdzVtx;
  kstTrkmMinIP2D      = NTupleIn->kstTrkmMinIP2D;
  kstTrkmMinIP2DE     = NTupleIn->kstTrkmMinIP2DE;
  kstTrkpCL           = NTupleIn->kstTrkpCL;
  kstTrkpDCABS        = NTupleIn->kstTrkpDCABS;
  kstTrkpDCABSE       = NTupleIn->kstTrkpDCABSE;
  kstTrkpdxyVtx       = NTupleIn->kstTrkpdxyVtx;
  kstTrkpdzVtx        = NTupleIn->kstTrkpdzVtx;
  kstTrkpMinIP2D      = NTupleIn->kstTrkpMinIP2D;
  kstTrkpMinIP2DE     = NTupleIn->kstTrkpMinIP2DE;

  tagB0               = NTupleIn->tagB0;
  genSignal           = NTupleIn->genSignal;
  mumGlobalMuon       = NTupleIn->mumGlobalMuon;
  mumTrackerMuon      = NTupleIn->mumTrackerMuon;
  mumStandAloneMuon   = NTupleIn->mumStandAloneMuon;
  mumTMOneStationTight = NTupleIn->mumTMOneStationTight;
  mupGlobalMuon       = NTupleIn->mupGlobalMuon;
  mupTrackerMuon      = NTupleIn->mupTrackerMuon;
  mupStandAloneMuon   = NTupleIn->mupStandAloneMuon;
  mupTMOneStationTight = NTupleIn->mupTMOneStationTight;

  kstTrkmGlobalMuon   = NTupleIn->kstTrkmGlobalMuon;
  kstTrkmTrackerMuon  = NTupleIn->kstTrkmTrackerMuon;
  kstTrkmTMOneStationTight = NTupleIn->kstTrkmTMOneStationTight;
  kstTrkmTMOneStationLoose = NTupleIn->kstTrkmTMOneStationLoose;
  kstTrkmTrackerMuonArbitrated = NTupleIn->kstTrkmTrackerMuonArbitrated;
  kstTrkpGlobalMuon    = NTupleIn->kstTrkpGlobalMuon;
  kstTrkpTrackerMuon   = NTupleIn->kstTrkpTrackerMuon;
  kstTrkpTMOneStationTight = NTupleIn->kstTrkpTMOneStationTight;
  kstTrkpTMOneStationLoose = NTupleIn->kstTrkpTMOneStationLoose;
  kstTrkpTrackerMuonArbitrated = NTupleIn->kstTrkpTrackerMuonArbitrated;

  bPt                  = NTupleIn->bPt;
  kstPt                = NTupleIn->kstPt;
  mumuPt               = NTupleIn->mumuPt;
  mumPt                = NTupleIn->mumPt;
  mupPt                = NTupleIn->mupPt;
  kstTrkmPt            = NTupleIn->kstTrkmPt;
  kstTrkpPt            = NTupleIn->kstTrkpPt;
  bPhi                 = NTupleIn->bPhi;
  kstPhi               = NTupleIn->kstPhi;
  mumuPhi              = NTupleIn->mumuPhi;
  mumPhi               = NTupleIn->mumPhi;
  mupPhi               = NTupleIn->mupPhi;
  kstTrkmPhi           = NTupleIn->kstTrkmPhi;
  kstTrkpPhi           = NTupleIn->kstTrkpPhi;
  bEta                 = NTupleIn->bEta;
  kstEta               = NTupleIn->kstEta;
  mumuEta              = NTupleIn->mumuEta;
  mumEta               = NTupleIn->mumEta;
  mupEta               = NTupleIn->mupEta;
  kstTrkmEta           = NTupleIn->kstTrkmEta;
  kstTrkpEta           = NTupleIn->kstTrkpEta;
  iso_mum              = NTupleIn->iso_mum;
  iso_mup              = NTupleIn->iso_mup;
  iso_trkm             = NTupleIn->iso_trkm;
  iso_trkp             = NTupleIn->iso_trkp;
  sum_iso              = NTupleIn->sum_iso;
  kstarmass            = NTupleIn->kstarmass;
  pass_preselection    = NTupleIn->pass_preselection;
  bdt_prob_v1          = NTupleIn->bdt_prob_v1;
  bdt_prob_v2          = NTupleIn->bdt_prob_v2;
  bdt_prob_v3          = NTupleIn->bdt_prob_v3;
  bdt_prob_v4          = NTupleIn->bdt_prob_v4;
  cos_theta_l          = NTupleIn->cos_theta_l;
  cos_theta_k          = NTupleIn->cos_theta_k;
  phi_kst_mumu         = NTupleIn->phi_kst_mumu;
  tagged_mass          = NTupleIn->tagged_mass;
  }
}






