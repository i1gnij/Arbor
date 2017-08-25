#ifndef _HitEff_hh_
#define _HitEff_hh_

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

class HitEff  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new HitEff ; }

		HitEff();

		~HitEff() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::vector<std::string> _EcalInputHits; 
		std::vector<std::string> _HcalInputHits; 
		std::ostream _output; 

		int _overwrite;
		TTree *_outputTree;	// Used only for single particle analysis
		//TTree *_outputHits;
		
		int _Num, _eventNr; 
		int _MCPID;
		float _MCPEn;
		float _MCP[3];
		
		int _NBush_E, _NBush_H; 
		int _NHit_E, _NHit_H;
		int _NLimHit_E, _NLimHit_H;
		int _NCluHit_E, _NCluHit_H; 
		int _NLCHit_E, _NLCHit_H; 
		
		float _En_E, _En_H; 
		float _LimEn_E, _LimEn_H;
		float _CluEn_E, _CluEn_H;
		float _LCEn_E, _LCEn_H; 

};

#endif


