#ifndef _TreeAna_hh_
#define _TreeAna_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Cluster.h>
#include <EVENT/Track.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>

class TTree;

class TreeAna  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new TreeAna ; }

		TreeAna();

		~TreeAna() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::ostream *_output;
		std::vector<std::string> _inputArborCollections;
		std::vector<std::string> _inputArborParticle;

		TTree *_outputTree;
		TTree *_outputArbor; 
		TTree *_outputEvt; 
		int _overwrite;
		int _Num;
		unsigned int _eventNr;

		float _MCPEn;
		float _MCPP[3];
		float _PosClu[3];
		int _nLECAL, _nLHCAL; 
		int _NNeC, _NChC; 

		float _En;
		float _P[3];
		int _ArborPID, _MCPID, _MCPIDOri; 
		int _NTrack, _NClu, _NH_ECAL, _NH_HCAL, _NLEcal, _NLHcal;
		float _E_Trk, _E_Clu, _EE_Clu, _EH_Clu, _FD_ECAL, _FD_HCAL, _FD_all, _MinDisHel, _MaxDisHel, _minDepth, _maxDepth; 
		float _FD[8];
		int _NH[8];
		int _NL[8];


		float _ECALEn, _HCALEn, _ECALFD, _HCALFD, _CluEn, _CluFD;
		int _CluSize, _ECALSize, _HCALSize; 
		int _NTrk; 		
	
		//int _nTrk, _nClu, _nTrkHit, _nCluHit; 
		int _Index, _type, _nHClu; 
		float _Depth, _LeadDepth; // _FD;
		float _Pos[3], _DirB[3], _DirF[3];
		float _DirAngle1, _DirAngle2, _DirAngle3; 	 	
		int _Size;
		float _Energy, _MinDepth; 

		float _TotalEn, _TotalP[3], _ChEn, _NeEn; 
		float _LCEn, _LCFD, _LCFD_E, _LCFD_H, _LCHCALEn, _LCECALEn;
		int _LCSize, _LCHCALSize, _LCECALSize, _LCnLECAL, _LCnLHCAL; 
		float _TotalClEn;
		int _LCNH[4];
		float _LCFD_E3[3], _LCFD_H3[3];

		int _nLFD01, _nLNH20; 

		//float _MCPEn, _TrkEn, _CluEn, _Energy, _Charge;
		//float _MCP[3], _TrkVtxP[3], _TrkEndP[3], _P[3], _EnSeg[10];
		//float _TotalEn, _TotalChEn, _TotalNeuEn; 
		//int _NArborC, _NArborN; 
		//float _Depth, _LeadCDepth, _LeadCEn, _FD, _SumDE, _LeadDE; 
		//int _LeadCHits; 
};

#endif


