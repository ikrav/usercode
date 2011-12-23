#ifndef EWKANA_NTUPLER_EWKANADEFS_HH 
#define EWKANA_NTUPLER_EWKANADEFS_HH

enum EMuType 
{ 
  kGlobal     = 1, 
  kTracker    = 2, 
  kStandalone = 4
};

enum EQualityBit
{ 
  // descriptions from DataFormats/MuonReco/interface/MuonSelectors.h
  kAll  			    = 0x000001,  // dummy options - always true
  kAllGlobalMuons		    = 0x000002,  // checks isGlobalMuon flag
  kAllStandAloneMuons		    = 0x000004,  // checks isStandAloneMuon flag
  kAllTrackerMuons		    = 0x000008,  // checks isTrackerMuon flag
  kTrackerMuonArbitrated	    = 0x000010,  // resolve ambiguity of sharing segments
  kAllArbitrated		    = 0x000020,  // all muons with the tracker muon arbitrated
  kGlobalMuonPromptTight	    = 0x000040,  // global muons with tighter fit requirements
  kTMLastStationLoose		    = 0x000080,  // penetration depth loose selector
  kTMLastStationTight		    = 0x000100,  // penetration depth tight selector
  kTM2DCompatibilityLoose	    = 0x000200,  // likelihood based loose selector
  kTM2DCompatibilityTight	    = 0x000400,  // likelihood based tight selector
  kTMOneStationLoose		    = 0x000800,  // require one well matched segment
  kTMOneStationTight		    = 0x001000,  // require one well matched segment
  kTMLastStationOptimizedLowPtLoose = 0x002000,  // combination of TMLastStation and TMOneStation
  kTMLastStationOptimizedLowPtTight = 0x004000,  // combination of TMLastStation and TMOneStation
  kGMTkChiCompatibility 	    = 0x008000,  // require tk stub have good chi2 relative to glb track
  kGMStaChiCompatibility	    = 0x010000,  // require sta stub have good chi2 compatibility relative to glb track
  kGMTkKinkTight		    = 0x020000,  // require a small kink value in the tracker stub
  kTMLastStationAngLoose	    = 0x040000,  // TMLastStationLoose with additional angular cuts
  kTMLastStationAngTight	    = 0x080000,  // TMLastStationTight with additional angular cuts
  kTMOneStationAngLoose 	    = 0x100000,  // TMOneStationLoose with additional angular cuts
  kTMOneStationAngTight 	    = 0x200000,  // TMOneStationTight with additional angular cuts
  //The two algorithms that follow are identical to what were known as
  //TMLastStationOptimizedLowPt* (sans the Barrel) as late as revision
  //1.7 of this file. The names were changed because indeed the low pt
  //optimization applies only to the barrel region, whereas the sel-
  //ectors above are more efficient at low pt in the endcaps, which is
  //what we feel is more suggestive of the algorithm name. This will be
  //less confusing for future generations of CMS members, I hope...
  kTMLastStationOptimizedBarrelLowPtLoose = 0x400000,  // combination of TMLastStation and TMOneStation but with low pT optimization in barrel only
  kTMLastStationOptimizedBarrelLowPtTight = 0x800000   // combination of TMLastStation and TMOneStation but with low pT optimization in barrel only
}; 

enum ETriggerBit
{
  // MuEG
  kHLT_Mu17_Ele8_CaloIdL          = 1UL<<1, // MC, data
  kHLT_Mu8_Ele17_CaloIdL          = 1UL<<2, // MC, data
  kHLT_Mu8_Photon20_CaloIdVT_IsoT = 1UL<<3, // MC, data
  kHLT_Mu15_Photon20_CaloIdL      = 1UL<<4, // MC, data

  // DoubleElectron
  kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL                                   = 1UL<<5,  // MC, data
  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL = 1UL<<6,  // MC, data
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30                              = 1UL<<7,  // MC, data
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30                             = 1UL<<8,  // data
  kHLT_Ele32_CaloIdL_CaloIsoVL_SC17                                                     = 1UL<<9,  // MC, data
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17                                       = 1UL<<10, // data
  kHLT_Ele8                                                                             = 1UL<<11, // MC, data
  kHLT_Ele8_CaloIdL_TrkIdVL                                                             = 1UL<<12, // MC, data
  kHLT_Ele8_CaloIdL_CaloIsoVL                                                           = 1UL<<13, // MC, data
  kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL                                          = 1UL<<14, // data
  kHLT_Ele17_CaloIdL_CaloIsoVL                                                          = 1UL<<15, // MC, data
  kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40                                                     = 1UL<<16, // MC, data
  kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL                                    = 1UL<<17, // MC, data

  // SingleElectron
  kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT = 1UL<<18, // MC, data
  kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT = 1UL<<19  // MC, data  
};

enum ETriggerObjBit
{
  // MuEG 
  kHLT_Mu17_Ele8_CaloIdL_MuObj          = 1UL<<0,
  kHLT_Mu17_Ele8_CaloIdL_EGObj          = 1UL<<1,  
  kHLT_Mu8_Ele17_CaloIdL_MuObj          = 1UL<<2,
  kHLT_Mu8_Ele17_CaloIdL_EGObj          = 1UL<<3,
  kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj = 1UL<<4,
  kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj = 1UL<<5,
  kHLT_Mu15_Photon20_CaloIdL_MuObj      = 1UL<<6,
  kHLT_Mu15_Photon20_CaloIdL_EGObj      = 1UL<<7,

  // DoubleElectron
  kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj                                   = 1UL<<8,
  kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj                                   = 1UL<<9,
  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj = 1UL<<10,
  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj = 1UL<<11,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj                               = 1UL<<12,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj                                = 1UL<<13,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj                             = 1UL<<14,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele2Obj                             = 1UL<<15,
  kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj                                                      = 1UL<<16,
  kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_SCObj                                                       = 1UL<<17,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj                                        = 1UL<<18,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_SCObj                                         = 1UL<<19,
  kHLT_Ele8_EleObj                                                                              = 1UL<<20,
  kHLT_Ele8_CaloIdL_TrkIdVL_EleObj                                                              = 1UL<<21,
  kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj                                                            = 1UL<<22,
  kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_EleObj                                           = 1UL<<23,
  kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj                                                           = 1UL<<24,
  kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj                                                      = 1UL<<25,
  kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_JetObj                                                      = 1UL<<26,
  kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj                                     = 1UL<<27,
  kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj                                     = 1UL<<28,
  
  // SingleElectron
  kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj = 1UL<<29,
  kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj = 1UL<<30
};
#endif
