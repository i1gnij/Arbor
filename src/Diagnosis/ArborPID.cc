/*
 * General PID
used to readout the FD, FD_10First, FD_20Later, ... and other quantities of bush
 * objective: pattern tagging for ECAL as: penetrating mip, deep interaction CH, EM Cluster; 
 * 				  HCAL as: penetrating mip, early interaction H, later interaction H, etc.
 *
 * */

#include <ArborPID.hh>
#include <DetectorPos.hh>
#include <ArborToolLCIO.hh>
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
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <EVENT/Cluster.h>
#include <EVENT/LCRelation.h>
#include <UTIL/CellIDDecoder.h>
#include <ArborTool.hh>
#include <stdexcept>
#include <TFile.h> 
#include <TTree.h>
#include <Rtypes.h> 
#include <sstream>		
#include <TH1.h>
#include <TVector3.h>


ArborPID aArborPID ;
ArborPID::ArborPID()
	: Processor("ArborPID"),
	_output(0)
{
	_description = "Measure Bush Quantities" ;

	_treeFileName="ArborPID.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);

	/*
	std::vector<std::string> BushCollections; 
	BushCollections.push_back(std::string("EcalBushes"));
	BushCollections.push_back(std::string("HcalBushes"));

	registerInputCollections( LCIO::CLUSTER,
			"BushCollections" ,
			"Bush to be analyzed" ,
			_BushCollections ,
			BushCollections);

	std::vector<std::string> BranchCollections;
        BranchCollections.push_back(std::string("EcalSortBranches"));
        BranchCollections.push_back(std::string("HcalSortBranches"));

        registerInputCollections( LCIO::CLUSTER,
                        "BranchCollections" ,
                        "Branch to be analyzed" ,
                        _BranchCollections ,
                        BranchCollections);	
	*/

	_treeName="Reco";
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

void ArborPID::init() {

	printParameters();

	TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));

	if (!tree_file->IsOpen()) {
		delete tree_file;
		tree_file=new TFile(_treeFileName.c_str(),"NEW");
	}

	_outputTree = new TTree(_treeName.c_str(),_treeName.c_str());
	_outputTree->SetAutoSave(32*1024*1024);  // autosave every 32MB

	_outputTree->Branch("EventNr", &_eventNr, "EventNr/I");
	_outputTree->Branch("Num", &_Num, "Num/I");

	_outputTree->Branch("MCPID", &_MCPID, "MCPID/I");
	_outputTree->Branch("MCE", &_MCEnergy, "MCE/F");
	_outputTree->Branch("MCP", _MCMomentum, "MCP[3]/F");
	_outputTree->Branch("RecoPID", &_RecoPID, "RecoPID/I");
	_outputTree->Branch("RecoE", &_RecoEnergy, "RecoE/F");
	_outputTree->Branch("RecoP", _RecoMomentum, "RecoP[3]/F");

	_outputTree->Branch("FDTotal", &_FDTotal, "FDTotal/F");
	_outputTree->Branch("FDEcal", &_FDEcal, "FDEcal/F");
	_outputTree->Branch("FDHcal", &_FDHcal, "FDHcal/F");
	_outputTree->Branch("FD", _FD, "FD[6]/F");
	
	_outputTree->Branch("NH", _NH, "NH[6]/I");
        _outputTree->Branch("NHPL", _NHPL, "NHPL[80]/I");       //ecal 30 + hcal 48 + n78: active ecal number of layer + n79, active hcal number of layer

	_outputTree->Branch("NBush",&_NBush,"NBush/I");
	_outputTree->Branch("Size", &_CluSize, "Size/I");	//Total Size
	_outputTree->Branch("CluTotalEn", &_CluTotalEn, "CluTotalEn/F");
	_outputTree->Branch("CluEn", _CluEn, "CluEn[6]/F");
	_outputTree->Branch("LCHitWeight", &_LCHitWeight, "LCHitWeight/F");
	_outputTree->Branch("LCEnWeight", &_LCEnWeight, "LCEnWeight/F");
	_outputTree->Branch("BushDepth", &_BushDepth, "BushDepth/F");

	_outputHits = new TTree("Hits", "Hits");
	_outputHits->Branch("EventNr", &_eventNr, "EventNr/I");
        _outputHits->Branch("Num", &_Num, "Num/I");
	_outputHits->Branch("CPid", &_CPid, "CPid/I");
	_outputHits->Branch("HitX", &_HitX, "HitX/F");
	_outputHits->Branch("HitY", &_HitY, "HitY/F");
	_outputHits->Branch("HitZ", &_HitZ, "HitZ/F");
	_outputHits->Branch("nLayer", &_nLayer, "nLayer/I");
	_outputHits->Branch("DepF", &_DepF, "DepF/I");

}


void ArborPID::processEvent( LCEvent * evtP ) 
{		

	if (evtP) 								
	{	
		try 	
		{    
			_eventNr=evtP->getEventNumber();
			if( _eventNr % 100 == 0)
				std::cout<<_eventNr<<" evts processed"<<std::endl;

			LCCollection * CaloHit = evtP->getCollection("ECALBarrel");	//just... to give the decoder
			CellIDDecoder<CalorimeterHit> idDecoder(CaloHit);
			int layernum = 0;
			TVector3 CurrHitP;

			LCCollection * RecoPMCPlink = evtP->getCollection("ArborCHPtoMCP");	//Only for chargedP...
			int NRel = RecoPMCPlink->getNumberOfElements();
			int LeadingCluSize = 0; 
			float LeadingCluEn = 0;
			float hitEn = 0;

			std::vector<Cluster *> attachedCluster; 

			for(int i0 = 0; i0 < NRel; i0++)
			{
				LCRelation *a_link = dynamic_cast<LCRelation*>( RecoPMCPlink->getElementAt( i0 ) );
				MCParticle *a_MCP = dynamic_cast<MCParticle*>(a_link->getFrom());
				ReconstructedParticle *a_RecoP = dynamic_cast<ReconstructedParticle*>(a_link->getTo());				
				_MCPID = a_MCP->getPDG();
				_MCEnergy = a_MCP->getEnergy();
				_MCMomentum[0] = a_MCP->getMomentum()[0];
				_MCMomentum[1] = a_MCP->getMomentum()[1];
				_MCMomentum[2] = a_MCP->getMomentum()[2];

				_RecoEnergy = a_RecoP->getEnergy();
				_RecoMomentum[0] = a_RecoP->getMomentum()[0];
				_RecoMomentum[1] = a_RecoP->getMomentum()[1];
				_RecoMomentum[2] = a_RecoP->getMomentum()[2];

				_NBush = a_RecoP->getClusters().size();

				attachedCluster.clear();

				EcalAll.clear();
				HcalAll.clear();
				EcalStart.clear();
				EcalMiddle.clear();
				EcalEnd.clear();
				HcalStart.clear();
				HcalMiddle.clear();
				HcalEnd.clear();

				_CluSize = 0;
				_CluTotalEn = 0;
				_FDTotal = 0;
				_FDEcal = 0;
				_FDHcal = 0;
				_LCHitWeight = 0;
				_LCEnWeight = 0;

				_CPid = i0;

				for(int j0 = 0; j0 < 80; j0++)
				{
					if(j0 < 6)
					{
						_FD[j0] = 0;
						_CluEn[j0] = 0;
						_NH[j0] = 0;
					}
					_NHPL[j0] = 0;
				}

				for(int i1 = 0; i1 < _NBush; i1++)
				{
					Cluster * a_Clu = a_RecoP->getClusters()[i1];
					int currCluSize = a_Clu->getCalorimeterHits().size();
					_CluSize += currCluSize; 
					_CluTotalEn += a_Clu->getEnergy();

					if(i1 == 0)
					{
						LeadingCluSize = currCluSize; 
						LeadingCluEn = a_Clu->getEnergy();
					}				
					attachedCluster.push_back(a_Clu);

					for(int i2 = 0; i2 < currCluSize; i2++)
					{
						CalorimeterHit * a_hit = a_Clu->getCalorimeterHits()[i2];
						CurrHitP = a_hit->getPosition();
						_HitX = a_hit->getPosition()[0];
						_HitY = a_hit->getPosition()[1];
						_HitZ = a_hit->getPosition()[2];
						_DepF = DepthFlag(CurrHitP);					
						_outputHits->Fill();					
	
						hitEn = a_hit->getEnergy();

						layernum = idDecoder(a_hit)["K-1"];
						_nLayer = layernum;

						if(DepthFlag(CurrHitP) == 0)	//EcalHit
						{
							_NHPL[layernum]++;
							EcalAll.push_back(a_hit);							

							if(layernum/10 == 0) 
							{
								EcalStart.push_back(a_hit);
							}
							else if(layernum/10 == 1) 
							{
								EcalMiddle.push_back(a_hit);	
							}
							else 
							{
								EcalEnd.push_back(a_hit);
							}
							_NH[ int(layernum/10) ]++;
							_CluEn[ int(layernum/10) ] += hitEn;
						}
						else if(DepthFlag(CurrHitP) == 1 || DepthFlag(CurrHitP) == 2)
						{
							_NHPL[layernum + 30]++;
							HcalAll.push_back(a_hit);

							if(layernum/16 == 0) 
							{	
								HcalStart.push_back(a_hit);
							}
							else if(layernum/16 == 1) 
							{
								HcalMiddle.push_back(a_hit);
							}			
							else
							{
								HcalEnd.push_back(a_hit);
							}
							_NH[ int(layernum/16) + 3 ]++;
							_CluEn[ int(layernum/16) + 3 ] += hitEn;

						}
						else
						{
							std::cout<<"DepthFlag Problematic"<<std::endl;
						}
					}

				}

				ClusterImpl * all_clu = NaiveMergeClu(attachedCluster);

				_FDTotal = FD(all_clu, CaloHit);

				/*
				if(_NH[0]) _FD[0] = FDV2( EcalStart, CaloHit );
				if(_NH[1]) _FD[1] = FDV2( EcalMiddle, CaloHit );
				if(_NH[2]) _FD[2] = FDV2( EcalEnd, CaloHit );

				if(_NH[3]) _FD[3] = FDV2( HcalStart, CaloHit );
				if(_NH[4]) _FD[4] = FDV2( HcalMiddle, CaloHit );
				if(_NH[5]) _FD[5] = FDV2( HcalEnd, CaloHit );

				if(EcalAll.size() != 0) 
					_FDEcal = FDV2(EcalAll, CaloHit);

				if(HcalAll.size() != 0)
					_FDHcal = FDV2(HcalAll, CaloHit);
				*/

				if(_CluSize && _CluTotalEn)
				{
					_LCHitWeight = float(LeadingCluSize)/_CluSize; 
					_LCEnWeight = LeadingCluEn/_CluTotalEn; 
				}
				_RecoPID = -1001; 

				_outputTree->Fill();
			}



		}		
		catch (lcio::DataNotAvailableException err) { }
	}  	

}	

void ArborPID::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		delete tree_file;
	}

}


