#ifndef _AnaTreeNew_hh_
#define _AnaTreeNew_hh_

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

class AnaTreeNew  : public marlin::Processor
{

	public:

		Processor*  newProcessor() { return new AnaTreeNew ; }

		AnaTreeNew();

		~AnaTreeNew() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();


	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::ostream *_output;
		std::vector<std::string> _inputArborCollections;
		std::vector<std::string> _inputPandoraCollections;
		std::vector<std::string> _inputMCParticle;

		TTree *_outputTree;
		TTree *_outputAPFO; 
		TTree *_outputPPFO; 
		int _overwrite;
		int _Num;
		unsigned int _eventNr;
		int MCCharge;


		int nCHMCP, nCHPFO_a, nCHPFO_p, nNEMCP, nNEPFO_a, nNEPFO_p;
		float thetaMCP, phiMCP, energyMCP, PMCP, thetaPdP, phiPdP, energyPdP, thetaAbP, phiAbP, energyAbP;
		int MCPDGMother, PdPID, AbPID;

		float AbD0, AbZ0, AbLCEn, Ab_dpp, Ab_dEENeu, PdD0, PdZ0, PdLCEn, Pd_dpp, Pd_dEENeu;
		int AbPIDMatrix[5][5], PdPIDMatrix[5][5];
		int TypeA, ChargeA, TrackHitA;
		float EnergyA, CluEnA, PA[3], StartPosA[3], EndPosA[3], CluEnComA[2];
		
		float _LCEn, _TCEn;

		int TypeP, ChargeP, TrackHitP;
		float EnergyP, CluEnP, PP[3], StartPosP[3], EndPosP[3], CluEnComP[2];
		int nRecMu_p, nRecEle_p, nRecHad_p, nRecGamma_p, nRecOtherNeu_p, nRecMu_a, nRecEle_a, nRecHad_a, nRecGamma_a, nRecOtherNeu_a;

};


#endif
