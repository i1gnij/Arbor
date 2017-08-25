#include <MarlinArbor.hh>
#include <ArborTool.hh>
#include <ArborToolLCIO.hh>
#include <ArborHit.hh>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <EVENT/LCGenericObject.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/SimCalorimeterHitImpl.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/ClusterImpl.h>
#include "UTIL/CellIDDecoder.h"
#include "UTIL/Operators.h"

#include <values.h>
#include <string>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <Rtypes.h>
#include <sstream>
#include <set>
#include <TVector3.h>
#include <vector>
#include <algorithm>

#include "DetectorPos.hh"

using namespace std;
using namespace lcio ;
using namespace marlin ;

marlin::StringParameters* MarlinArbor::arborParas = 0;

const string AHCALCellIDDecoder = "M:3,S-1:4,I:9,J:9,K-1:7";    //Same for Ecal and AHCAL
//const string ECALCellIDDecoder  = "M:3,S-1:4,I:9,J:9,K-1:6";	//Need to verify that for DHCAL
const string DHCALCellIDDecoder = "M:3,S-1:4,I:9,J:9,K-1:7";
const string ECALCellIDDecoder = "S-1:3,M:3,K-1:6,I:16,GRZone:3,J:32:16";
CellIDDecoder<CalorimeterHit> idDecoder_ECAL(ECALCellIDDecoder);
CellIDDecoder<CalorimeterHit> idDecoder_DHCAL(DHCALCellIDDecoder);

extern linkcoll InitLinks; 
extern linkcoll IterLinks_1; 
extern linkcoll IterLinks; 
extern linkcoll links_debug; 
extern branchcoll Trees; 
extern std::vector<int> IsoHitsIndex;

std::vector<std::string> CaloHitCollections;

MarlinArbor aMarlinArbor ;
MarlinArbor::MarlinArbor()
        : Processor("MarlinArbor")
          ,
          _output(0)
{
        _description = " Tree Algorithm in Marlin" ;

	std::vector<std::string> EcalCaloHitCollections;
	EcalCaloHitCollections.push_back(std::string("ECALBarrel"));
	EcalCaloHitCollections.push_back(std::string("ECALEndcap"));
	registerInputCollections( LCIO::CALORIMETERHIT,      //adding a flag to Calo /SimCalo
			"EcalHitCollections" ,
			"Hit Collection Names" ,
			_EcalCalCollections ,
			EcalCaloHitCollections);
	
	std::vector<std::string> HcalCaloHitCollections;
        HcalCaloHitCollections.push_back(std::string("HCALBarrel"));
        registerInputCollections( LCIO::CALORIMETERHIT,      //adding a flag to Calo /SimCalo
                        "HcalHitCollections" ,
                        "Hit Collection Names" ,
                        _HcalCalCollections ,
                        HcalCaloHitCollections);

	std::vector<float> cepc_thresholds;
	cepc_thresholds.push_back(10);
	cepc_thresholds.push_back(90);
	cepc_thresholds.push_back(50);
	cepc_thresholds.push_back(7.5);
	registerProcessorParameter( "ThresholdsforArborBuilding" ,
                        "Thresholds: EE, EH, HH, EE_Seed" ,
                        _cepc_thresholds ,
                        cepc_thresholds);
    //jing Sep 28
    std::string DisABFileName = "hits_th.root";
    registerProcessorParameter("DisABFileName",
            "root file name for dis ab",
            _DisABFileName,
            DisABFileName);

	std::vector<float> arbor_sel_params;
	cepc_thresholds.push_back(6.744);
	cepc_thresholds.push_back(55.72);
	cepc_thresholds.push_back(2.554);
	cepc_thresholds.push_back(32.06);
	cepc_thresholds.push_back(174.);
	cepc_thresholds.push_back(39.7);
	cepc_thresholds.push_back(1.198);
	cepc_thresholds.push_back(25.66);
	registerProcessorParameter( "LinkParamforArborBuilding" ,
                        "Thresholds: MeanX, RMSX, MeanY, RMSY, MeanX2, RMSX2, MeanY2, RMSY2" ,
                        _arbor_sel_params,
                        arbor_sel_params);
}

void MarlinArbor::init() {

    CreateROOT(_DisABFileName);
	printParameters();

    //jing Sep 28
    //disABFile = new TFile(_DisABFileName.c_str(),"RECREATE");
    //disABTree = new TTree("DisAB","DisAB");


	/*
	   for(unsigned int i0 = 0; i0 < _EcalPreShowerCollections.size();i0++)
	   {
	   CaloHitCollections.push_back(_EcalCalCollections[i0].c_str());
	   }
	   */
	for(unsigned int i1 = 0; i1 < _EcalCalCollections.size(); i1++)
	{
		CaloHitCollections.push_back(_EcalCalCollections[i1].c_str());
	}
	for(unsigned int i2 = 0; i2 < _HcalCalCollections.size(); i2++)
	{
		CaloHitCollections.push_back(_HcalCalCollections[i2].c_str());
	}
	CaloHitCollections.push_back("ECALPSHitCollection");	// Ok for convieient

}

void MarlinArbor::HitsPreparation()
{
	cout<<"Start to prepare Hits"<<endl;
}

void MarlinArbor::LinkVisulization( LCEvent * evtPP, std::string Name, std::vector<CalorimeterHit*> Hits, linkcoll inputLinks ) 
{

	LCCollectionVec *colllink = new LCCollectionVec(LCIO::LCRELATION);
	LCFlagImpl linkflag;
	linkflag.setBit(LCIO::CHBIT_LONG);
	colllink->setFlag(linkflag.getFlag());

	int NLink = inputLinks.size();
	std::pair<int, int> currlink; 
	for(int i0 = 0; i0 < NLink; i0++)
	{
		currlink = inputLinks[i0];

		CalorimeterHit* a_hit = Hits[currlink.first];
		CalorimeterHit* b_hit = Hits[currlink.second];

		LCRelationImpl *a_link = new LCRelationImpl(a_hit, b_hit);
		colllink->addElement( a_link );
	}

	evtPP->addCollection(colllink, Name);
}

void MarlinArbor::MakeBush( LCEvent * evtPP, std::vector<std::string> inputTreeCollections, std::string outputBushCollection )
{
	LCCollection *currbushcoll = new LCCollectionVec(LCIO::CLUSTER);
	LCFlagImpl bushflag;
	bushflag.setBit(LCIO::CHBIT_LONG);
	currbushcoll->setFlag(bushflag.getFlag());	

	int Ncoll = inputTreeCollections.size();
	int NTree = 0;
	for(int i = 0; i < Ncoll; i++)
	{
		LCCollection * a_TreeColl = evtPP ->getCollection(inputTreeCollections[i].c_str());
		NTree = a_TreeColl->getNumberOfElements();
		for(int j = 0; j < NTree; j++)
		{
			Cluster * a_tree = dynamic_cast<Cluster*>(a_TreeColl->getElementAt(j));
			ClusterImpl *a_bush = NaiveCluImpl(a_tree);
			currbushcoll->addElement(a_bush);
		}
	}

	evtPP->addCollection( currbushcoll, outputBushCollection );
}

void MarlinArbor::MakeIsoHits( LCEvent * evtPP, std::vector<CalorimeterHit*> inputCaloHits, std::string outputBushCollection )
{
	LCCollection * isohitcoll = new LCCollectionVec(LCIO::CALORIMETERHIT);

	string initString = "M:3,S-1:3,I:9,J:9,K-1:6";          //Need to verify
	isohitcoll->parameters().setValue(LCIO::CellIDEncoding, initString);

	LCFlagImpl flag;
	flag.setBit(LCIO::CHBIT_LONG);                  
	flag.setBit(LCIO::CHBIT_ID1);
	flag.setBit(LCIO::RCHBIT_ENERGY_ERROR);
	isohitcoll->setFlag(flag.getFlag());

	int nhit = inputCaloHits.size();

	for(int i = 0; i < nhit; i++)
	{
		CalorimeterHit* a_hit = inputCaloHits[i];
		CalorimeterHitImpl * collhit = new CalorimeterHitImpl();
		collhit->setPosition(a_hit->getPosition());
		collhit->setCellID0(a_hit->getCellID0());
		collhit->setCellID1(a_hit->getCellID1());
		collhit->setEnergy(a_hit->getEnergy());
		isohitcoll->addElement(collhit);
	}

	evtPP->addCollection(isohitcoll, outputBushCollection);
}

void MarlinArbor::processEvent( LCEvent * evtP )
{
	if(evtP)
	{

		int EvtNr = evtP->getEventNumber();
		_eventNr = EvtNr; 

		std::cout<<EvtNr<<" events processed"<<std::endl;

		MarlinArbor::HitsPreparation();	//Absorb isolated hits; 

		TVector3 currHitPos;

		std::vector< TVector3 > inputHitsPos;
		std::vector< ArborHit > inputABHit; 
		std::vector< CalorimeterHit* > inputHits;  
		std::vector< std::vector<int> > Sequence; 
		int LayerNum = 0; 
		int StaveNum = 0; 
		int SubDId = -10; 
		float Depth = 0; 
		int KShift = 0; 
		//int NTrkHit = 0;
		TVector3 TrkEndPointPos; 
		//float TrkEP[3];
		/*
		   try{
		   LCCollection * TrackColl = evtP ->getCollection("MarlinTrkTracks");
		   LCCollection * trkconvhitcoll = new LCCollectionVec(LCIO::CALORIMETERHIT);
		   LCFlagImpl bushflag;
		   bushflag.setBit(LCIO::CHBIT_LONG);
		   trkconvhitcoll->setFlag(bushflag.getFlag());

		   int NTrk = TrackColl->getNumberOfElements();
		   for(int i0 = 0; i0 < NTrk; i0++)	//Here should put the limit on track energy to ensure safety; or, calculate the reference direction
		   {
		   Track * a_Trk = dynamic_cast<Track*>( TrackColl->getElementAt( i0 ) );
		   NTrkHit = a_Trk->getTrackerHits().size();
		   TrkEndPointPos = (a_Trk->getTrackerHits()[NTrkHit - 1])->getPosition();
		   if( (TrkEndPointPos.Perp() > 1600 || abs(TrkEndPointPos[2]) > 2100 ) && 1 == 0 )	//tmply mute
		   {
		   ArborHit a_abhit(TrkEndPointPos, -2, 0, 0, 0, 0);
		   inputABHit.push_back(a_abhit);
		   CalorimeterHitImpl * a_hit = new CalorimeterHitImpl; 
		   a_hit->setEnergy(0);
		   TrkEP[0] = TrkEndPointPos.X();
		   TrkEP[1] = TrkEndPointPos.Y();
		   TrkEP[2] = TrkEndPointPos.Z();
		   a_hit->setPosition(TrkEP);
		   trkconvhitcoll->addElement(a_hit);
		   inputHits.push_back(a_hit);
		   }
		   }

		   evtP->addCollection(trkconvhitcoll, "TrkConvCaloHit");

		   }catch(lcio::DataNotAvailableException zero){}	//PreShower should also include it...
		   */
		std::vector<CalorimeterHit*> IsoHits;
		// std::vector<CalorimeterHit*> Ecal/HcalIsoHits; 

		for(unsigned int i1 = 0; i1 < CaloHitCollections.size(); i1++)
		{
			try{

				KShift = 0;
				SubDId = -1; 

				if( i1 < _EcalCalCollections.size() )
					SubDId = 1; 
				else if( i1 < _EcalCalCollections.size() + _HcalCalCollections.size() )
					SubDId = 2;
				else
					SubDId = 3; 

				if(i1 >  _EcalCalCollections.size() - 1)
					KShift = 100; 
				else if( i1 == CaloHitCollections.size() - 2)	//HCAL Ring
					KShift = 50;

				LCCollection * CaloHitColl = evtP ->getCollection(CaloHitCollections[i1].c_str());
				{
					int NHitsCurrCol = CaloHitColl->getNumberOfElements();
					CellIDDecoder<CalorimeterHit> idDecoder(CaloHitColl);
					cout<<"Sub Collection "<<CaloHitCollections[i1].c_str()<<endl; 
					for(int i2 = 0; i2 < NHitsCurrCol; i2++)
					{
						CalorimeterHit * a_hit = dynamic_cast<CalorimeterHit*>(CaloHitColl->getElementAt(i2));
						currHitPos = a_hit->getPosition();
						Depth = DisSeedSurface(currHitPos);
						if( CaloHitCollections[i1] == "LCAL")
						{
							LayerNum = idDecoder(a_hit)["K"];
							StaveNum = 0;
						}
						else
						{
							LayerNum = idDecoder(a_hit)["K-1"] + 1 + KShift;
							StaveNum = idDecoder(a_hit)["S-1"] + 1;
						}
						ArborHit a_abhit(currHitPos, LayerNum, 0, Depth, StaveNum, SubDId, a_hit->getEnergy());
						// cout<<"Time "<<a_hit->getTime()<<endl; 
						inputABHit.push_back(a_abhit);
						inputHits.push_back(a_hit);
					}
				}

				// cout<<i1<<"  Stat  "<<SubDId<<" ~~~ "<<inputABHit.size()<<endl; 

			}catch(lcio::DataNotAvailableException zero) { }
		}

        LoadParam(_arbor_sel_params);
		Sequence = Arbor(inputABHit, _cepc_thresholds);
		ClusterBuilding( evtP, "EHBushes", inputHits, Trees, 0 );
		LinkVisulization(evtP, "Links_init", inputHits, InitLinks);
		LinkVisulization(evtP, "Links_iter_1", inputHits, IterLinks_1);
		LinkVisulization(evtP, "Links_iter_2", inputHits, IterLinks);
		LinkVisulization(evtP, "Links_init_Debug", inputHits, links_debug);

		for(unsigned int i2 = 0; i2 < IsoHitsIndex.size(); i2++)
		{
			CalorimeterHit* a_Isohit = inputHits[ IsoHitsIndex[i2] ];
			if(a_Isohit->getEnergy() > 0)	//Veto Trk End Hits
			{
				IsoHits.push_back(a_Isohit);
			}
		}

		MakeIsoHits(evtP, IsoHits, "AllIsolatedHits");
	}
}


void MarlinArbor::end()
{
    CloseROOT();
//    disABTree->Write();
//    disABFile->Close();
	std::cout<<"Arbor Ends. Good luck"<<std::endl;

}
