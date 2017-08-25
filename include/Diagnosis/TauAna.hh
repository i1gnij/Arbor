#ifndef _TauAna_hh_
#define _TauAna_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>

class TTree;

class TauAna  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new TauAna ; }

		TauAna();

		~TauAna() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::string _colAdcVals;
//		TFile *tree_file;

		int _overwrite;
		TTree *_outputTree;
		TTree *_outputTau;
		TTree *_outputMCTau;
		TTree *_outputEvt;
		unsigned int _eventNr;
		int _Num, _evtN; 
		float _Charge, _En, _VisEn, _VisEnCh, _Reco_Charge, _Reco_VisEn, _Reco_VisEnCh, _DeltaTheta, _DeltaPhi, _DeltaR, _DeltaPhiVis, _DeltaThetaVis, _cone, _En_tr_all, _TauEnergy, _TauCharge, _TauMass, _MCTauEnergy, _MCTauCharge, _MCTauMass,_TauRecoE, _cosTh, _cosThR;
		float _MCTauPMass,_TauPMass,_MCTauPEn,_TauPEn,_MCTauPP[3],_TauPP[3],_MCTauMMass,_TauMMass,_MCTauMEn,_TauMEn,_MCTauMP[3],_TauMP[3],_tauAngle;		
		int _NMuon, _NEle, _NPion, _NKaon, _NPhoton, _NUndef, _Reco_NMuon, _Reco_NEle, _Reco_NPion, _Reco_NPhoton, _Reco_NUndef, _TauNumber, _MCTauM,_MCTauP, _MCTauNumber, _TauReco, _Tag_e,_Tag_mu, _Tag_pi,_Tag_photon, _Ancharged,_Bncharged,_Anphoton,_Bnphoton,_Ane,_Anmu,_Anpi,_Bne,_Bnmu,_Bnpi,_dimu; 
		float _Aimp, _Bimp, _LD0, _LZ0, _NLD0, _NLZ0;
		float _TauP[3];
		float _MC_TauP[3];
		float _VisP[3];
		float _P[3];
		int _HDecay[999];
		int _HDDecay[999];
		float _Reco_VisP[3];
		float _Mass, _MCMass; 
		float _AvisEn,_BvisEn,_AMass,_BMass,_Acone1,_Acone2,_Acone3,_Bcone1,_Bcone2,_Bcone3,_recoilM,_Acharged,_Bcharged,_ATauM,_BTauM;
		float _AvisP[3];
		float _BvisP[3];
		TH1F *_photonmass, *_photonmass2, *_photonmass3; 
		float _Reco_recoilM,_ArecoEn,_BrecoEn,_ArecoM,_BrecoM;
		int _ArecoNch,_BrecoNch,_ArecoNph,_BrecoNph,_Arecoe,_Arecomu,_Arecopi,_Brecoe,_Brecomu,_Brecopi;
		float _ArecoP[3];
                float _BrecoP[3];
		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


