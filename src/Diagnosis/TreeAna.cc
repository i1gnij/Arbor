/*
 * used to readout the FD, FD_10First, FD_20Later, ... and other quantities of bush
 * objective: pattern tagging for ECAL as: penetrating mip, deep interaction CH, EM Cluster; 
 * 				  HCAL as: penetrating mip, early interaction H, later interaction H, etc.
 *	//Used as SP Analysis Code
 * */

#include <TreeAna.hh>
#include "ArborTool.hh"
#include "ArborToolLCIO.hh"
#include "DetectorPos.hh"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <values.h>
#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <EVENT/Cluster.h>
#include <EVENT/LCRelation.h>
#include <UTIL/CellIDDecoder.h>
#include <stdexcept>
#include <TFile.h> 
#include <TTree.h>
#include <Rtypes.h> 
#include <sstream>		
#include <TH1.h>
#include <TVector3.h>

using namespace std;

const string ECALCellIDDecoder = "M:3,S-1:3,I:9,J:9,K-1:6";

TreeAna aTreeAna ;
TreeAna::TreeAna()
	: Processor("TreeAna"),
	_output(0)
{
	_description = "Measure Bush Quantities" ;

	_treeFileName="TreeAna.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);

	std::vector<std::string> inputArborParticle;
        inputArborParticle.push_back(std::string("ArborTrkParticle_HQ_Barrel"));
        inputArborParticle.push_back(std::string("ArborTrkParticle_HQ_Endcap"));
        inputArborParticle.push_back(std::string("ArborTrkParticle_MQ_Barrel"));
	inputArborParticle.push_back(std::string("ArborTrkParticle_MQ_Endcap"));
	inputArborParticle.push_back(std::string("ArborTrkParticle_MQ_Vtx"));
	inputArborParticle.push_back(std::string("ArborTrkParticle_LQ"));
	inputArborParticle.push_back(std::string("ArborTrkParticle_NonMatch"));
	inputArborParticle.push_back(std::string("ArborNeutralParticle_ECAL"));
	inputArborParticle.push_back(std::string("ArborNeutralParticle_HCAL"));
	inputArborParticle.push_back(std::string("ArborNeutralParticle_LinkEHCAL"));
	registerProcessorParameter("InputArborParticle" ,
			"Name of Reconstructed Arbor Particle Collections" ,
			_inputArborParticle,
			inputArborParticle);


	std::vector<std::string> inputArborCollections;
	inputArborCollections.push_back("EcalTrees");	//Trees & Bushes
	inputArborCollections.push_back("EcalBushes");
	inputArborCollections.push_back("HcalTrees");
	inputArborCollections.push_back("HcalBushes");
	registerProcessorParameter("InputArborCollection" ,
			"Name of Reconstructed Arbor Collections" ,
			_inputArborCollections,
			inputArborCollections);


	_treeName="Arbor";
	registerProcessorParameter( "TreeName" , 
			"The name of the ROOT tree" ,
			_treeName ,
			_treeName);

	_overwrite=0;
	registerProcessorParameter( "OverwriteFile" , 
			"If zero an already existing file will not be overwritten." ,
			_overwrite ,
			_overwrite);

	_Num=0;

}

void TreeAna::init() {

	printParameters();

	TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));

	if (!tree_file->IsOpen()) {
		delete tree_file;
		tree_file=new TFile(_treeFileName.c_str(),"NEW");
	}

	_outputTree = new TTree(_treeName.c_str(),_treeName.c_str());
	_outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB

	_outputEvt = new TTree( "Evt", "Evt" );
	_outputEvt->SetAutoSave(32*1024*1024);

	_outputArbor = new TTree( "ArborPFO", "ArobrPFO");
	_outputArbor->SetAutoSave(32*1024*1024);

	_outputArbor->Branch("EventNr", &_eventNr, "EventNr/I");
        _outputArbor->Branch("Num", &_Num, "Num/I");
        _outputArbor->Branch("NTrack", &_NTrack, "NTrack/I");
	_outputArbor->Branch("En", &_En, "En/F");
	_outputArbor->Branch("P", _P, "P[3]/F");
	_outputArbor->Branch("ArborPID", &_ArborPID, "ArborPID/I");
	_outputArbor->Branch("MCPID", &_MCPID, "MCPID/I");
        _outputArbor->Branch("NClu", &_NClu, "NClu/I");
	_outputArbor->Branch("PosClu", _PosClu, "PosClu[3]/F");
	_outputArbor->Branch("E_Clu", &_E_Clu, "E_Clu/F");
	_outputArbor->Branch("EE_Clu", &_EE_Clu, "EE_Clu/F");
	_outputArbor->Branch("EH_Clu", &_EH_Clu, "EH_Clu/F");
	_outputArbor->Branch("NH_ECAL", &_NH_ECAL, "NH_ECAL/I");
	_outputArbor->Branch("NH_HCAL", &_NH_HCAL, "NH_HCAL/I");
	_outputArbor->Branch("FD_ECAL", &_FD_ECAL, "FD_ECAL/F");
        _outputArbor->Branch("FD_HCAL", &_FD_HCAL, "FD_HCAL/F");
	_outputArbor->Branch("FD_all", &_FD_all, "FD_all/F");
	_outputArbor->Branch("FD", _FD, "FD[8]/F");	//Divide into 6 parts; 
	_outputArbor->Branch("NH", _NH, "NH[8]/I");
	_outputArbor->Branch("NL", _NL, "NL[8]/I");
	_outputArbor->Branch("NLEcal", &_NLEcal, "NLEcal/I");
	_outputArbor->Branch("NLHcal", &_NLHcal, "NLHcal/I");
	_outputArbor->Branch("MaxDisHel", &_MaxDisHel, "MaxDisHel/F");	//Maximal Distance from Hit to Helix
	_outputArbor->Branch("MinDisHel", &_MinDisHel, "MinDisHel/F");
	_outputArbor->Branch("minDepth", &_minDepth, "minDepth/F");
	_outputArbor->Branch("maxDepth", &_maxDepth, "maxDepth/F");
	_outputArbor->Branch("nLFD01", &_nLFD01, "nLFD01/I");
        _outputArbor->Branch("nLNH20", &_nLNH20, "nLNH20/I");

	_outputEvt->Branch("EventNr", &_eventNr, "EventNr/I");
	_outputEvt->Branch("Num", &_Num, "Num/I");
	_outputEvt->Branch("NTrk", &_NTrk, "NTrk/I");
	_outputEvt->Branch("NChC", &_NChC, "NChC/I");
	_outputEvt->Branch("NNeC", &_NNeC, "NNeC/I");
	_outputEvt->Branch("MCPIDOri", &_MCPIDOri, "MCPIDOri/I");
	_outputEvt->Branch("MCPEn", &_MCPEn, "MCPEn/F");
	_outputEvt->Branch("MCPP", _MCPP, "MCPP[3]/F");

	_outputEvt->Branch("TotalEn", &_TotalEn, "TotalEn/F");
	_outputEvt->Branch("TotalP[3]", _TotalP, "TotalP[3]/F");
	_outputEvt->Branch("TotalClEn", &_TotalClEn, "TotalClEn/F");
	_outputEvt->Branch("ChEn", &_ChEn, "ChEn/F");
	_outputEvt->Branch("NeEn", &_NeEn, "NeEn/F");

	_outputEvt->Branch("LCEn", &_LCEn, "LCEn/F");   //Leading Cluster Quatities
	_outputEvt->Branch("LCSize", &_LCSize, "LCSize/I");
	_outputEvt->Branch("LCNH", _LCNH, "LCNH[4]/I");
	_outputEvt->Branch("LCFD", &_LCFD, "LCFD/F");
	_outputEvt->Branch("LCFD_E", &_LCFD_E, "LCFD_E/F");
	_outputEvt->Branch("LCFD_H", &_LCFD_H, "LCFD_H/F");
	_outputEvt->Branch("LCFD_E3", _LCFD_E3, "LCFD_E3[3]/F");
        _outputEvt->Branch("LCFD_H3", _LCFD_H3, "LCFD_H3[3]/F");

	_outputEvt->Branch("LCHCALEn", &_LCHCALEn, "LCHCALEn/F");
	_outputEvt->Branch("LCECALEn", &_LCECALEn, "LCECALEn/F");
	_outputEvt->Branch("LCHCALSize", &_LCHCALSize, "LCHCALSize/I");
        _outputEvt->Branch("LCECALSize", &_LCECALSize, "LCECALSize/I");
	_outputEvt->Branch("LCnLECAL", &_LCnLECAL, "LCnLECAL/I");
	_outputEvt->Branch("LCnLHCAL", &_LCnLHCAL, "LCnLHCAL/I");

	_outputTree->Branch("EventNr", &_eventNr, "EventNr/I");
	_outputTree->Branch("Num", &_Num, "Num/I");

	_outputTree->Branch("Type", &_type, "Type/I");
	_outputTree->Branch("Index", &_Index, "Index/I");

	_outputTree->Branch("CluEn", &_CluEn, "CluEn/F");
	_outputTree->Branch("CluSize", &_CluSize, "CluSize/I");
	_outputTree->Branch("CluFD", &_CluFD, "CluFD/F");
	_outputTree->Branch("CluPos", _Pos, "CluPos[3]/F");

	_outputTree->Branch("ECALEn", &_ECALEn, "ECALEn/F");	
	_outputTree->Branch("HCALEn", &_HCALEn, "HCALEn/F");
	_outputTree->Branch("ECALSize", &_ECALSize, "ECALSize/I");
	_outputTree->Branch("HCALSize", &_HCALSize, "HCALSize/I");
	_outputTree->Branch("ECALFD", &_ECALFD, "ECALFD/F");
	_outputTree->Branch("HCALFD", &_HCALFD, "HCALFD/F");

	_outputTree->Branch("nLHCAL", &_nLHCAL, "nLHCAL/I");
	_outputTree->Branch("nLECAL", &_nLECAL, "nLECAL/I");

//	_outputTree->Branch("nLFD01", &_nLFD01, "nLFD01/I");
//	_outputTree->Branch("nLNH20", &_nLNH20, "nLNH20/I");

	_outputTree->Branch("Depth", &_Depth, "Depth/F");
	_outputTree->Branch("LeadDepth", &_LeadDepth, "LeadDepth/F");
	_outputTree->Branch("MinDepth", &_MinDepth, "MinDepth/F");	

	//Trk Info

	/*
	   _outputTree->Branch("DirB", _DirB, "DirB[3]/F");
	   _outputTree->Branch("DirF", _DirF, "DirF[3]/F");
	   _outputTree->Branch("DirAngle1", &_DirAngle1, "DirAngle1/F");
	   _outputTree->Branch("DirAngle2", &_DirAngle2, "DirAngle2/F");
	   _outputTree->Branch("DirAngle3", &_DirAngle3, "DirAngle3/F");
	   */

}

void TreeAna::processEvent( LCEvent * evtP ) 
{		

	if (evtP) 								
	{	
		_eventNr=evtP->getEventNumber();
		if( _eventNr % 100 == 0)
			std::cout<<_eventNr<<" evts processed"<<std::endl;

		TVector3 CluPos, HitPos; //, TrkVTXPos, TrkEndPPos; 	//Define as the least deap position
		float currDepth = 0; 
		float HitEn = 0;
		float _Charge = 0;
		int NArbor = 0;
                float RecoEn = 0; 
		float CurrCEn = 0; 	//C En: Cluster Energy
		float MaxCEn = 0; 		

		float BushTime = 0;
		float BushDist[3] = {0, 0, 0};

		_NChC = 0; 
		_NNeC = 0; 

		_ChEn = 0; 
		_NeEn = 0;
		_TotalEn = 0; 
		_TotalClEn = 0; 	

		_LCEn = 0;
		_LCSize = 0;
		_LCFD = -1; 
		_LCFD_E = -1;
		_LCFD_H = -1; 

		Cluster* leadingCluster(0); 

		for(int k = 0; k < 4; k++)
		{
			if(k < 3)
			{
				_TotalP[k] = 0;
				_LCFD_E3[k] = -1;
				_LCFD_H3[k] = -1;
			}
			_LCNH[k] = 0;			//Number of HCAL Hit with bigger size
		}

		std::vector<CalorimeterHit*> Ecalhits; 
		std::vector<CalorimeterHit*> Hcalhits;

		std::vector<CalorimeterHit*> allhits; 

		std::vector<CalorimeterHit*> EH_1; 
		std::vector<CalorimeterHit*> EH_2; 
		std::vector<CalorimeterHit*> EH_3;			
		std::vector<CalorimeterHit*> HH_1; 
                std::vector<CalorimeterHit*> HH_2; 
                std::vector<CalorimeterHit*> HH_3;
		std::vector<CalorimeterHit*> HH_4; 
		std::vector<CalorimeterHit*> HH_5;

		try{

			LCCollection * MCP = evtP->getCollection("MCParticle");

			for(int i = 0; i < MCP->getNumberOfElements(); i++)
			{
				MCParticle * a_MCP = dynamic_cast<MCParticle*>(MCP->getElementAt(i));
				if( a_MCP->getParents().size() == 0 )
				{
					_MCPIDOri = a_MCP->getPDG();
					_MCPEn = a_MCP->getEnergy();
					_MCPP[0] = a_MCP->getMomentum()[0];
					_MCPP[1] = a_MCP->getMomentum()[1];
					_MCPP[2] = a_MCP->getMomentum()[2];
				}
			}

			LCCollection *Trks = evtP->getCollection("MarlinTrkTracks");
			_NTrk = Trks->getNumberOfElements();

			LCCollection *ChC = evtP->getCollection("ArborTrkChargedCluster");
			_NChC = ChC->getNumberOfElements();

			LCCollection *NeC = evtP->getCollection("ArborTrkNeutralCluster");
			_NNeC = NeC->getNumberOfElements();

			LCCollection *TrkMCLink = evtP->getCollection("MCTruthMarlinTrkTracksLink");
			int NLink = TrkMCLink->getNumberOfElements();

			CellIDDecoder<CalorimeterHit> idDecoder(ECALCellIDDecoder);

			for(unsigned int p0 = 0; p0 < _inputArborParticle.size(); p0++)
			{

				LCCollection * RecoCol = evtP->getCollection( _inputArborParticle[p0].c_str() );
				NArbor = RecoCol->getNumberOfElements();

				for(int i0 = 0; i0 < NArbor; i0++)
				{

					ReconstructedParticle* a_Arbor = dynamic_cast<ReconstructedParticle*>(RecoCol->getElementAt(i0));
					RecoEn = a_Arbor->getEnergy();

					_Charge = a_Arbor->getCharge();
					_TotalEn += RecoEn; 	


					_ArborPID = a_Arbor->getType();

					_MCPID = 0;
			
					_NTrack = a_Arbor->getTracks().size();
					_NClu = a_Arbor->getClusters().size();
					_En = a_Arbor->getEnergy();
					_P[0] = a_Arbor->getMomentum()[0];
					_P[1] = a_Arbor->getMomentum()[1];
					_P[2] = a_Arbor->getMomentum()[2];

					_E_Clu = 0; 
					_EE_Clu = 0; 
					_EH_Clu = 0; 
					_NH_ECAL = 0; 
					_NH_HCAL = 0; 

					_PosClu[0] = -1;
					_PosClu[1] = -1;
					_PosClu[2] = -1;

					allhits.clear();
					Ecalhits.clear();
                                        Hcalhits.clear();
                                        EH_1.clear();
                                        EH_2.clear();
                                        EH_3.clear();
                                        HH_1.clear();
                                        HH_2.clear();
                                        HH_3.clear();
                                        HH_4.clear();
                                        HH_5.clear();

					_maxDepth = -100;
					_minDepth = 1E6;
					_MaxDisHel = -1; 
					_MinDisHel = 1E10; 

					if(_NTrack == 1)
					{
						Track * a_trk = a_Arbor->getTracks()[0];
						HelixClass * currHelix = new HelixClass();
						currHelix->Initialize_Canonical(a_trk->getPhi(), a_trk -> getD0(), a_trk -> getZ0(), a_trk -> getOmega(), a_trk->getTanLambda(), 3.5);
						for(int j = 0; j < _NClu; j++)
						{

							Cluster * a_tmpClu = a_Arbor->getClusters()[j];
							_PosClu[0] = a_tmpClu->getPosition()[0];
							_PosClu[1] = a_tmpClu->getPosition()[1];
							_PosClu[2] = a_tmpClu->getPosition()[2];

							int NCluHits = a_tmpClu->getCalorimeterHits().size();

							for(int j1 = 0; j1 < NCluHits; j1++)
							{
								CalorimeterHit * a_hit = a_tmpClu->getCalorimeterHits()[j1];

								BushTime = currHelix->getDistanceToPoint((float*)a_hit->getPosition(), BushDist);

								if(BushTime > 0)
								{
									if(BushDist[2] > _MaxDisHel)
									{
										_MaxDisHel = BushDist[2];
									}
									if(BushDist[2] < _MinDisHel)
									{
										_MinDisHel = BushDist[2];
									}
								}
							}
						}

						for(int pp = 0; pp < NLink; pp++)
						{
							LCRelation *a_link = dynamic_cast<LCRelation*>(TrkMCLink->getElementAt(pp));
							if(a_link->getTo() == a_trk)
							{
								MCParticle * a_MCP = dynamic_cast<MCParticle*>(a_link->getFrom());
								_MCPID = a_MCP->getPDG(); 
								break; 
							}
						}
					}

					for(int s = 0; s < _NClu; s++)
					{

						Cluster * a_tmpClu = a_Arbor->getClusters()[s];

						_E_Clu += a_tmpClu->getEnergy();

						int NCluHits = a_tmpClu->getCalorimeterHits().size();

						for(int s1 = 0; s1 < NCluHits; s1++)
						{
							CalorimeterHit * a_hit = a_tmpClu->getCalorimeterHits()[s1];
							allhits.push_back(a_hit);

							int NLayer = idDecoder(a_hit)["K-1"];  

							HitPos = a_hit->getPosition();
							currDepth = DisSeedSurface(HitPos);
							if(currDepth > _maxDepth)
							{
								_maxDepth = currDepth;
							}
							if(currDepth < _minDepth)
							{
								_minDepth = currDepth;
							}

							if( fabs(a_hit->getEnergy() - DHCALCalibrationConstant) < 1.0E-6 )	//or other fancy judgements...
							{
								_NH_HCAL++;
								_EH_Clu += a_hit->getEnergy();
								Hcalhits.push_back(a_hit);
								if(NLayer < 10)
								{
									HH_1.push_back(a_hit);
								}
								else if(NLayer < 20)
								{
									HH_2.push_back(a_hit);
								}
								else if(NLayer < 30)
								{
									HH_3.push_back(a_hit);
								}
								else if(NLayer < 40)
								{
									HH_4.push_back(a_hit);
								}
								else
								{
									HH_5.push_back(a_hit);
								}

							}
							else
							{
								_NH_ECAL++; 
								_EE_Clu += a_hit->getEnergy();
								Ecalhits.push_back(a_hit);
								if(NLayer < 10)
								{
									EH_1.push_back(a_hit);
								}
								else if(NLayer < 20)
								{
									EH_2.push_back(a_hit);
								}
								else
								{
									EH_3.push_back(a_hit);
								}
							}
						}
					}

					_FD_all = FDV2(allhits, ECALCellIDDecoder);
					_FD_ECAL = FDV2(Ecalhits, ECALCellIDDecoder);
					_FD_HCAL = FDV2(Hcalhits, ECALCellIDDecoder);

					for(int p0 = 0; p0 < 8; p0++)
					{
						_NH[p0] = 0;
						_NL[p0] = 0;
						_FD[p0] = 0;
					}

					_FD[0] = FDV2(EH_1, ECALCellIDDecoder);
					_FD[1] = FDV2(EH_2, ECALCellIDDecoder);
					_FD[2] = FDV2(EH_3, ECALCellIDDecoder);
					_FD[3] = FDV2(HH_1, ECALCellIDDecoder);
					_FD[4] = FDV2(HH_2, ECALCellIDDecoder);
					_FD[5] = FDV2(HH_3, ECALCellIDDecoder);
					_FD[6] = FDV2(HH_4, ECALCellIDDecoder);
					_FD[7] = FDV2(HH_5, ECALCellIDDecoder);

					_NH[0] = EH_1.size();
					_NH[1] = EH_2.size();
					_NH[2] = EH_3.size();
					_NH[3] = HH_1.size();
					_NH[4] = HH_2.size();
					_NH[5] = HH_3.size();
					_NH[6] = HH_4.size();
					_NH[7] = HH_5.size();

					_NLEcal = ActiveLayers(Ecalhits, ECALCellIDDecoder);
					_NLHcal = ActiveLayers(Hcalhits, ECALCellIDDecoder);

					_NL[0] = ActiveLayers(EH_1, ECALCellIDDecoder);
					_NL[1] = ActiveLayers(EH_2, ECALCellIDDecoder);
					_NL[2] = ActiveLayers(EH_3, ECALCellIDDecoder);
					_NL[3] = ActiveLayers(HH_1, ECALCellIDDecoder);
					_NL[4] = ActiveLayers(HH_2, ECALCellIDDecoder);
					_NL[5] = ActiveLayers(HH_3, ECALCellIDDecoder);
					_NL[6] = ActiveLayers(HH_4, ECALCellIDDecoder);
					_NL[7] = ActiveLayers(HH_5, ECALCellIDDecoder);

					_nLFD01 = 0;
					_nLNH20 = 0;

					 for(int p1 = 0; p1 < 8; p1++)
                                        {
                                                if(_FD[p1] > 0.1)
                                                        _nLFD01 ++;
						if(_NH[p1] > 20)
							_nLNH20 ++;
                                        }

					allhits.clear();
					Ecalhits.clear();
					Hcalhits.clear();
					EH_1.clear();
					EH_2.clear();
					EH_3.clear();
					HH_1.clear();
					HH_2.clear();
					HH_3.clear();
					HH_4.clear();
					HH_5.clear();

					if(a_Arbor->getClusters().size() > 0)
					{
						CurrCEn = (a_Arbor->getClusters()[0])->getEnergy();
						_TotalClEn += CurrCEn; 						

						if(CurrCEn > MaxCEn)
						{
							MaxCEn = CurrCEn;
							leadingCluster = a_Arbor->getClusters()[0];
						}	
					}

					for(int s = 0; s < 3; s++)
					{
						_TotalP[s] += a_Arbor->getMomentum()[s];					
					}

					if(fabs(_Charge) > 0.01)
					{
						_ChEn += RecoEn;
					}
					else
					{
						_NeEn += RecoEn; 
					}


					_outputArbor->Fill();

				}

			}

			if(leadingCluster)
			{
				_LCEn = leadingCluster->getEnergy();
				_LCSize = leadingCluster->getCalorimeterHits().size();	
				_LCHCALSize = 0;
				_LCECALSize = 0; 
				_LCHCALEn = 0; 
				_LCECALEn = 0; 

				Ecalhits.clear();
				Hcalhits.clear();

				for(int i1 = 0; i1 < _LCSize; i1++)
				{
					CalorimeterHit * a_hit = leadingCluster->getCalorimeterHits()[i1];
					HitEn = a_hit->getEnergy();
					if( fabs(HitEn - DHCALCalibrationConstant) < 1.0E-6 )       //Should use some other fancy things...
					{
						_LCHCALSize ++;
						_LCHCALEn += HitEn;
						Hcalhits.push_back(a_hit);
					}
					else
					{
						_LCECALSize ++;
						_LCECALEn += HitEn;
						Ecalhits.push_back(a_hit);
					}
				}

				for(int s = 0; s < 4; s++)	//3, 5, 7, 9
				{
					_LCNH[s] = NHScaleV2(ECALCellIDDecoder, Hcalhits, 2*s + 3, 2*s + 3, 1);
				}

				_LCnLECAL = ActiveLayers(Ecalhits, ECALCellIDDecoder);
				_LCnLHCAL = ActiveLayers(Hcalhits, ECALCellIDDecoder);

				if(Ecalhits.size() > 0)
				{
					_LCFD_E = FDV2(Ecalhits, ECALCellIDDecoder);
					for(int p = 0; p < 3; p++)
					{
						_LCFD_E3[p] = FD_I(Ecalhits, ECALCellIDDecoder, p+3);
					}
				}
				if(Hcalhits.size() > 0)
				{
					_LCFD_H = FDV2(Hcalhits, ECALCellIDDecoder);
					for(int p = 0; p < 3; p++)
					{
						_LCFD_H3[p] = FD_I(Hcalhits, ECALCellIDDecoder, p+3);
					}
				}
			}	

			_outputEvt->Fill();


		}catch(lcio::DataNotAvailableException err){}


		for(unsigned int j = 0; j < _inputArborCollections.size(); j++ )
		{
			try
			{
				LCCollection *ClusterColl = evtP->getCollection( _inputArborCollections[j].c_str());

				_type = j; 				

				for(int k = 0; k < ClusterColl->getNumberOfElements(); k++)
				{
					Cluster * a_clu = dynamic_cast<Cluster*>(ClusterColl->getElementAt(k));
					_CluEn = a_clu->getEnergy();
					_CluSize = a_clu->getCalorimeterHits().size();
					_CluFD = -1; 
					if(_CluSize)
						_CluFD = FDV3(a_clu, ECALCellIDDecoder);
					CluPos = a_clu->getPosition();				

					_Pos[0] = CluPos.X();
					_Pos[1] = CluPos.Y();
					_Pos[2] = CluPos.Z();

					_Depth = DisSeedSurface(CluPos);

					_ECALSize = 0;
					_HCALSize = 0; 
					_HCALEn = 0; 
					_ECALEn = 0; 
					_nLECAL = 0; 
					_nLHCAL = 0; 

					Ecalhits.clear();
					Hcalhits.clear();

					_LeadDepth = 0;
					_MinDepth = 1.0E10;

					for(int i1 = 0; i1 < _CluSize; i1++)
					{
						CalorimeterHit * a_hit = a_clu->getCalorimeterHits()[i1];
						HitEn = a_hit->getEnergy();
						HitPos = a_hit->getPosition();
						currDepth = DisSeedSurface(HitPos);
						if(currDepth > _LeadDepth)
						{
							_LeadDepth = currDepth;
						}
						if(currDepth < _MinDepth)
						{
							_MinDepth = currDepth;
						}
						if( fabs(HitEn - DHCALCalibrationConstant) < 1.0E-6 )	//Should use some other fancy things...
						{
							_HCALSize ++;
							_HCALEn += HitEn; 
							Hcalhits.push_back(a_hit);
						}
						else 
						{
							_ECALSize ++;
							_ECALEn += HitEn; 
							Ecalhits.push_back(a_hit);
						}
					}

					_nLECAL = ActiveLayers(Ecalhits, ECALCellIDDecoder);
					_nLHCAL = ActiveLayers(Hcalhits, ECALCellIDDecoder);					

					_ECALFD = -1; 
					_HCALFD = -1;
					if(Ecalhits.size() > 0)
						_ECALFD = FDV2(Ecalhits, ECALCellIDDecoder);
					if(Hcalhits.size() > 0)
						_HCALFD = FDV2(Hcalhits, ECALCellIDDecoder);

					_outputTree->Fill();

				}


			}catch(lcio::DataNotAvailableException err){}
		}

		_Num++;

	}  	

}	

void TreeAna::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		delete tree_file;
	}

}


