    // Auxiliary event info data object
    struct TEventInfo
    {
      UInt_t eventId;      // event ID
      UInt_t nGenZ;        // number of generated Z in event
      UInt_t nDimuons;     // number of reconstructed dimuons
      UInt_t nMuons;       // number of reconstructed muons selected
      UInt_t nJets;        // number of reconstructed jets      
      UInt_t triggerBits;  // HLT trigger bits
    };
    
    // Generator level info data object
    struct TGenInfo
    {
      UInt_t eventId;                // event ID
      Float_t zmass, zpt, zy, zphi;  // kinematics
      Float_t mass, pt, y, phi;      // dimuon kinematics
      Float_t pt_1, eta_1, phi_1;    // muon kinematics
      Float_t pt_2, eta_2, phi_2;    // anti-muon kinematics
    };
    
    // Reco muon/track data object
    struct TMuon
    {
      UInt_t eventId;              // event ID
      Float_t pt, eta, phi;        // kinematics
      Float_t trkIso03, trkIso05;  // track isolation
      Float_t emIso03, emIso05;    // ECAL-based isolation
      Float_t hadIso03, hadIso05;  // HCAL-based isolation
      Float_t hoIso03, hoIso05;    // HO-based isolation
      Float_t d0, d0Err;           // impact parameter
      Float_t caloComp;            // measure of compatibility with a MIP through calorimeters
      Float_t segComp;             // measure of compatibility with
      Float_t nchi2;               // track chi^2/ndf 
      Int_t q;                     // charge
      UInt_t isGoodBits;           // bits for various muon quality criteria
      UInt_t typeBits;             // global muon, tracker muon, or standalone muon
      UInt_t nTkHits;              // number of inner tracker hits
      UInt_t hltMatchBits;         // bits for matching with HLT primitives 
      UInt_t isUsed;               // used in a dimuon combination?     
    };
    
    // Reco dimuon data object
    struct TDimuon
    {
      UInt_t eventId;            // event ID
      Float_t mass, pt, y, phi;  // dimuon kinematics
      
      // leading pT muon
      Float_t pt_1, eta_1, phi_1;
      Float_t trkIso03_1, trkIso05_1;
      Float_t emIso03_1, emIso05_1;
      Float_t hadIso03_1, hadIso05_1;
      Float_t hoIso03_1, hoIso05_1;
      Float_t d0_1, d0Err_1;
      Float_t caloComp_1;
      Float_t segComp_1;
      Float_t nchi2_1;
      Int_t q_1;
      UInt_t isGoodBits_1;
      UInt_t typeBits_1;
      UInt_t nTkHits_1;
      UInt_t hltMatchBits_1;
             
      // lagging pT muon
      Float_t pt_2, eta_2, phi_2;
      Float_t trkIso03_2, trkIso05_2;
      Float_t emIso03_2, emIso05_2;
      Float_t hadIso03_2, hadIso05_2;
      Float_t hoIso03_2, hoIso05_2;
      Float_t d0_2, d0Err_2;
      Float_t caloComp_2;
      Float_t segComp_2;
      Float_t nchi2_2;
      Int_t q_2;
      UInt_t isGoodBits_2;
      UInt_t typeBits_2;
      UInt_t nTkHits_2;
      UInt_t hltMatchBits_2;
    };    
    
    enum ETriggerBit
    {
      kHLT_Mu9   = 1,
      kHLT_Mu11  = 2,
      kHLT_Mu15  = 4,
      kHLT_Jet30 = 8
    };
    
    enum EMuType 
    { 
      kGlobal     = 1, 
      kTracker    = 2, 
      kStandalone = 4
    };
    
    enum ESelBit
    {
      kAllArbitrated          =   1,  // All arbitration (DT/CSC/RPC Hits) put on at least one segments given a global Muon
      kPromptTight            =   2,  // Standard global muon identification
      kTMLastStationLoose     =   4,  // Loose matching requirements on lastmost muon station of reco 
      kTMLastStationTight     =   8,  // Tight matching requirements on lastmost muon station of reco 
      kTMOneStationLoose      =  16,  // Loose matching requirements on at least one muon station of reco 
      kTMOneStationTight      =  32,  // Tight matching requirements on at least one muon station of reco 
      kTM2DCompatibilityLoose =  64,  // Loose requirement on sum of compatabiliity variables ===> 1.2 Segment compatibility + 0.8 calo compatibility > 0.8
      kTM2DCompatibilityTight = 128   // Tight requirement on sum of compatabiliity variables ===> 1.2 Segment compatibility + 0.8 calo compatibility > 1.2
    };
