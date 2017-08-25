#ifndef _ArborPID_hh_
#define _ArborPID_hh_

//#include <RConfigure.h>
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
#include <TH3.h>
class TTree;

class ArborPID  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new ArborPID ; }

		ArborPID();

		~ArborPID() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::vector<std::string> _hcalCollections; 
		std::vector<std::string> _branchCollections; 
		std::ostream _output; 

		int _overwrite;
		TTree *_outputTree;

		TTree *_outputHits;
		float _HitX, _HitY, _HitZ; 
		int _nLayer, _DepF; 	
		int _CPid; 

		int _Num;
		unsigned int _eventNr;
		int _MCPID; 
		float _MCEnergy;
		float _MCMomentum[3];
		int _RecoPID;
		float _RecoEnergy;
		float _RecoMomentum[3];

		float _FDTotal;
		float _FDEcal, _FDHcal;
		float _FD[6];
		int _NH[6]; 
		int _NHPL[80];
		int _CluSize;
		float _CluEn[6]; 
		float _CluTotalEn; 
		int _NBush; 
		float _LCHitWeight, _LCEnWeight; 
		float _BushDepth;

		std::vector<CalorimeterHit*> EcalStart;
		std::vector<CalorimeterHit*> EcalMiddle;
		std::vector<CalorimeterHit*> EcalEnd;
		std::vector<CalorimeterHit*> HcalStart;
		std::vector<CalorimeterHit*> HcalMiddle;
		std::vector<CalorimeterHit*> HcalEnd;
		std::vector<CalorimeterHit*> EcalAll;
		std::vector<CalorimeterHit*> HcalAll;

		/*
		   int _NHitTotal; 
		   int _NHitArbor;
		   int _NumBranch; 
		   float _SumLength; 
		   int _MCPID[2];
		   float _MCPEn[2];
		   float _MCPDis; 

		   int _BranchSize; 
		   int _BrNHMC1; 
		   int _BrNHMC2; 
		   int _BrNHUndef; 
		   float _BranchPurity; 
		   int _LeadMCPID;
		   int _ClusterCollID; 
		   int _ClusterID;

		   int _NBushC, _NBushN, _NBushD;
		   int _NHitC, _NHitN, _NHitDrop;

		   std::vector<int> _NHitBranch; 
		   std::vector<float> _MCPEnVec;
		   std::vector<int> _MCPIDVec;
		   std::vector<float> _ClusterLength;
		   std::vector<int> SortMCPIndex;
		   std::vector<int> SortClusterIndex;
		   std::vector<int> _NHitBranchMCP1;
		   std::vector<int> _NHitBranchMCP2;
		   std::vector<int> _NHitBranchRem; 
		   std::vector<int> _ArborBranchPurity; 

		   float _BranchLength[10];
		   int _NHitBranMC[2][10];
		   float _DHCALFirstThreshold; 
		   int _WrongNHits; 
		   int _NHMCPSize; 

		   std::string _fileName;
		   std::ostream *_output;
		   std::string _histFileName;
		 */
};

#endif


