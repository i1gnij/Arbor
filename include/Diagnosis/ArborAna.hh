#ifndef _ArborAna_hh_
#define _ArborAna_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>

class TTree;

class ArborAna  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new ArborAna ; }

		ArborAna();

		~ArborAna() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::ostream *_output;
		TTree *_outputTree;
		TTree *_outputPFO; 
		std::vector<std::string> _inputArborParticle;

		int _overwrite;
		int _Num;
		unsigned int _eventNr;

		float _ThetaTau, _PhiTau, _ThetaVis, _PhiVis, _EnTau, _MassTau, _InvMassTau; 
		float _Theta_Calo, _Theta_Clu, _Phi_Calo, _Phi_Clu; 
		int _NChargeMCP, _NNeutralMCP, _NMuon, _NPion, _NPion0, _NEM, _NPhoton, _NNeutrino;
		float _EChargeMCP, _ENeutralMCP, _ENeutrino; 
		int _NFS_Charged, _NFS_Neutral; 
		float _CE_sR, _CE_lR, _NE_sR, _NE_lR, _InvE_sR, _InvE_lR; 
		float _EnSeg[10], _EnIsoFrag[2], _EnIsoHit[2];

		int _NChargeArbor, _NNeutralArbor, _NEffChargeArbor, _NEffNeutralArbor;
		float _E_noCluTrk, _HQ_Charged;
		float _EChargeArbor, _ENeutralArbor;
		float _TotalEVis, _InvMassVis, _ThetaReco, _PhiReco;
		float _TotalPVis[3], _CoM, _QuarkP[3], _CosTheta; 
		float _Charge, _TrkEn, _CluEn; 
		int _NClu, _Index; 
		float _PandoraTotalEn, _PandoraChEn, _PandoraNeEn; 
		float _L_MCPEn, _L_ChargeArborEn; 
		int _L_MCPID; 
		float _LeadingCMCP[3], _L_ChargeArborP[3];
};

#endif


