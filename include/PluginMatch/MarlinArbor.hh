#ifndef _MarlinArbor_hh_
#define _MarlinArbor_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/LCCollection.h>
#include <TNtuple.h>
#include <TObject.h>
#include <Arbor.hh>

#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

class TTree;

class MarlinArbor  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new MarlinArbor ; }

		MarlinArbor();

		~MarlinArbor() {};

		void init();

		void HitsPreparation(); // LCEvent * evtP);

		void LinkVisulization( LCEvent * evtPP, std::string Name, std::vector<CalorimeterHit*> Hits, linkcoll inputLinks );

		// void ClusterBuilding( LCEvent * evtPP, std::string Name, std::vector<CalorimeterHit*> Hits, std::vector< std::vector<int> > BranchOrder );

		void MakeIsoHits( LCEvent * evtPP, std::vector<CalorimeterHit*> inputCaloHits, std::string outputBushCollection );

		void MakeBush( LCEvent * evtPP, std::vector<std::string> inputTreeCollections, std::string outputBushCollection );

		void MakeMergedBush( LCEvent * evtPP, std::vector<std::string> inputTreeCollections, std::string outputBushCollection, float TreeSeedDisThreshold );

		void processEvent( LCEvent * evtP );

		void end();

		static marlin::StringParameters* getParameters() { return arborParas; }

	protected:
		std::string _colName;
		std::vector<std::string> _CalCollections;
		std::vector<std::string> _SimCalCollections;
		std::vector<std::string> _garlicCollections;
		std::vector<std::string> _endHitTrackCollections;
		std::vector<std::string> _EcalPreShowerCollections;
		std::vector<std::string> _EcalCalCollections;
		std::vector<std::string> _HcalCalCollections;
        std::vector<float> _cepc_thresholds; 
        std::vector<float> _arbor_sel_params;

        std::string _DisABFileName;

        //		TTree *_outputTree;
        //       TFile *disABFile;
        //      TTree *disABTree;
        int _EH; 
        float _HitPos[3];
        float _BushP[3];
        float _CloseDis; 
        float _HitEnergy;

        int _CellSize; 
        int _CaloTrackLengthCut;

        int _Num, _Seg, _eventNr;
        int numElements;

        float _DHCALFirstThreshold; 
        float _InitLinkDisThreshold;

        bool _DHCALSimuDigiMode; 
        bool _FlagInputSimHit;
        bool _FlagMutePhoton;
        bool _FlagMuteChargeParticle;
        bool _FlagMuteGarlicHits;
        bool _FlagUseTrackerEndHit; 

        TH2F *_h1, *_h2, *_h7; 
        TH1F *_h3, *_h4, *_h5, *_h6; 
        std::ostream _output;

    private:

        static marlin::StringParameters* arborParas;
};


#endif


