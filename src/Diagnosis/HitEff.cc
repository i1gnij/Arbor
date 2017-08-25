/*
 * General PID
used to readout the FD, FD_10First, FD_20Later, ... and other quantities of bush
 * objective: pattern tagging for ECAL as: penetrating mip, deep interaction CH, EM Cluster; 
 * 				  HCAL as: penetrating mip, early interaction H, later interaction H, etc.
 *
 * */

#include <HitEff.hh>
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


HitEff aHitEff ;
HitEff::HitEff()
	: Processor("HitEff"),
	_output(0)
{
	_description = "Measure Bush HitCollection Efficiency" ;

	_treeFileName="HitEff.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);

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

void HitEff::init() {

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
	_outputTree->Branch("MCPEn", &_MCPEn, "MCPEn/F");
	_outputTree->Branch("MCP", _MCP, "MCP[3]/F");

	_outputTree->Branch("NBush_E",&_NBush_E,"NBush_E/I");
	_outputTree->Branch("NBush_H",&_NBush_H,"NBush_H/I");

	_outputTree->Branch("NHit_E",&_NHit_E,"NHit_E/I");
	_outputTree->Branch("NHit_H",&_NHit_H,"NHit_H/I");
	_outputTree->Branch("NLimHit_E",&_NLimHit_E,"NLimHit_E/I");
        _outputTree->Branch("NLimHit_H",&_NLimHit_H,"NLimHit_H/I");
	_outputTree->Branch("NCluHit_E",&_NCluHit_E,"NCluHit_E/I");
	_outputTree->Branch("NCluHit_H",&_NCluHit_H,"NCluHit_H/I");
	_outputTree->Branch("NLCHit_E",&_NLCHit_E,"NLCHit_E/I");
	_outputTree->Branch("NLCHit_H",&_NLCHit_H,"NLCHit_H/I");

	_outputTree->Branch("En_E",&_En_E,"En_E/F");
        _outputTree->Branch("En_H",&_En_H,"En_H/F");
	_outputTree->Branch("LimEn_E",&_LimEn_E,"LimEn_E/F");
        _outputTree->Branch("LimEn_H",&_LimEn_H,"LimEn_H/F");
        _outputTree->Branch("CluEn_E",&_CluEn_E,"CluEn_E/F");
        _outputTree->Branch("CluEn_H",&_CluEn_H,"CluEn_H/F");
        _outputTree->Branch("LCEn_E",&_LCEn_E,"LCEn_E/F");
        _outputTree->Branch("LCEn_H",&_LCEn_H,"LCEn_H/F");

	_EcalInputHits.push_back("ECALBarrel");
	_EcalInputHits.push_back("ECALEndcap");
	_EcalInputHits.push_back("ECALOther");
	_EcalInputHits.push_back("LCAL");
	_EcalInputHits.push_back("LHCAL");

	_HcalInputHits.push_back("HCALBarrel");
        _HcalInputHits.push_back("HCALEndcap");
        _HcalInputHits.push_back("HCALOther");

	_Num = 0; 
}


void HitEff::processEvent( LCEvent * evtP ) 
{		

	if (evtP) 								
	{	
		try 	
		{    
			_eventNr=evtP->getEventNumber();

			// input Collections: MCParticle, EcalBush, HcalBush, EcalHits, HcalHits

			_MCPID = 0; _MCPEn = -1000; _MCP[0] = -1000; _MCP[1] = -1000; _MCP[2] = -1000;	
			_NHit_E = 0; _NHit_H = 0; _NLimHit_E = 0; _NLimHit_H = 0; _NCluHit_E = 0; _NCluHit_H = 0; _NLCHit_E = 0; _NLCHit_H = 0;
			_En_E = 0; _En_H = 0; _LimEn_E = 0; _LimEn_H = 0; _CluEn_E = 0; _CluEn_H = 0; _LCEn_E = 0; _LCEn_H = 0;		

			try{
			LCCollection * MCPart = evtP->getCollection("MCParticle");
                        int NMCP = MCPart->getNumberOfElements();
			LCCollection * EcalB = evtP->getCollection("EcalBushes");
			_NBush_E = EcalB->getNumberOfElements();
			LCCollection * HcalB = evtP->getCollection("HcalBushes");
                        _NBush_H = HcalB->getNumberOfElements();

			for(int i = 0; i < NMCP; i++)
			{
				MCParticle * a_MCP = dynamic_cast<MCParticle*>(MCPart->getElementAt(i));
				if(a_MCP->getParents().size() == 0)
				{
					_MCPID = a_MCP->getPDG();
					_MCPEn = a_MCP->getEnergy();
					_MCP[0] = a_MCP->getMomentum()[0];
					_MCP[1] = a_MCP->getMomentum()[1];
					_MCP[2] = a_MCP->getMomentum()[2];
					break; 
				}
			}

			int MaxBushSize = 0;
			int currBushSize = 0;

			for(int i1 = 0; i1 < _NBush_E; i1++)
			{
				Cluster * a_bush = dynamic_cast<Cluster*>(EcalB->getElementAt(i1));
				currBushSize = a_bush->getCalorimeterHits().size();
				_NCluHit_E += currBushSize; 
				_CluEn_E += a_bush->getEnergy();
				if(currBushSize > MaxBushSize)
				{
					_NLCHit_E = currBushSize;
					_LCEn_E = a_bush->getEnergy();
					MaxBushSize = currBushSize;
				}
			}

			MaxBushSize = 0;

			for(int i3 = 0; i3 < _NBush_H; i3++)
			{
				Cluster * a_bush = dynamic_cast<Cluster*>(HcalB->getElementAt(i3));
				currBushSize = a_bush->getCalorimeterHits().size();
				_NCluHit_H += currBushSize;
				_CluEn_H += a_bush->getEnergy();
				if(currBushSize > MaxBushSize)
				{
					_NLCHit_H = currBushSize;
					_LCEn_H = a_bush->getEnergy();
					MaxBushSize = currBushSize;
				}
			}
			}catch (lcio::DataNotAvailableException err) { }

			for(int j = 0; j < 5; j++)
			{
				try{
					LCCollection * A_ECALCol = evtP->getCollection(_EcalInputHits[j].c_str());
					_NHit_E += A_ECALCol->getNumberOfElements();
					for(int j1 = 0; j1 < int(A_ECALCol->getNumberOfElements()); j1++)
					{
						CalorimeterHit* a_hit = dynamic_cast<CalorimeterHit*>(A_ECALCol->getElementAt(j1));
						_En_E += a_hit->getEnergy();
					}
				}catch (lcio::DataNotAvailableException err) { }
			}

			for(int j2 = 0; j2 < 3; j2++)
			{
				try{
					LCCollection * A_HCALCol = evtP->getCollection(_HcalInputHits[j2].c_str());
					_NHit_H += A_HCALCol->getNumberOfElements();
					for(int j3 = 0; j3 < int(A_HCALCol->getNumberOfElements()); j3++)
					{
						CalorimeterHit* a_hit = dynamic_cast<CalorimeterHit*>(A_HCALCol->getElementAt(j3));
						_En_H += a_hit->getEnergy();
					}
				}catch (lcio::DataNotAvailableException err) { }
			}

			_outputTree->Fill();

			_Num++; 
		}		
		catch (lcio::DataNotAvailableException err) { }
	}  	

}	

void HitEff::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		delete tree_file;
	}

}


