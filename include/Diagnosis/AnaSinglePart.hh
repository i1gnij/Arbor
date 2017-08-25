#ifndef _AnaSinglePart_hh_
#define _AnaSinglePart_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/ReconstructedParticle.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>
#include <TVector3.h>

class TTree;

class AnaSinglePart  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new AnaSinglePart ; }

		AnaSinglePart();

		~AnaSinglePart() {};

		void init();

		void processEvent( LCEvent * evtP );
            
                
		void end();

	protected:
		std::string _treeFileName;
		std::ostream *_output;
		int _overwrite;
		TTree *_outputEvt;
		TTree *_outputPFO;
                int _eventNr, _Num, _filenum;

/////////////////////Evt///////////////////////////
		int _nMCP, _nSelMCP12B2;
		int _MCPOID;
		float _MCPOTheta, _MCPOPhi, _MCPOEn;

		float _THEn;

		int _nClu;
		float _LCEn, _TCEn;

		int _nPFO_a, _nCHPFOs_a, _nNEPFOs_a, _nCMIPPFOs_a, _nCEMPFOs_a, _nCHADPFOs_a, _nNEMPFOs_a, _nNHADPFOs_a;
		float _leadingPFOEn, _leadingNePFOEn, _leadingChPFOEn;
		float _LChPFOTheta, _LNePFOTheta, _LChPFOPhi, _LNePFOPhi;
		float _TotalRecoP4_a[4];
		

/////////////////////PFO///////////////////////////

		int _Charge_a, _Type_a, _MCTType_a;
		float _PPFO_a[4], _PFOCluEn_a;
		float _PFOTheta, _PFOPhi;
		float _trkD0, _trkZ0;
		float _trkSP[3], _trkEP[3];
		int  _NHTrk;


};

#endif




