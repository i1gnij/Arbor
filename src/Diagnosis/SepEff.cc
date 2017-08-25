/*
 * General PID
used to readout the FD, FD_10First, FD_20Later, ... and other quantities of bush
 * objective: pattern tagging for ECAL as: penetrating mip, deep interaction CH, EM Cluster; 
 * 				  HCAL as: penetrating mip, early interaction H, later interaction H, etc.
 *
 * */

#include <SepEff.hh>
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


SepEff aSepEff ;
SepEff::SepEff()
	: Processor("SepEff"),
	_output(0)
{
	_description = "Measure Bush HitCollection Efficiency" ;

	_treeFileName="SepEff.root";
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

void SepEff::init() {

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

	_outputTree->Branch("NPFO",&_NPFO,"NPFO/I");
        _outputTree->Branch("TotalRecoEn",&_TotalRecoEn,"TotalRecoEn/F");
	_outputTree->Branch("NCh",&_NCh,"NCh/I");
        _outputTree->Branch("NNe",&_NNe,"NNe/I");
	_outputTree->Branch("NChMC",&_NChMC,"NChMC/I");
        _outputTree->Branch("NNeMC",&_NNeMC,"NNeMC/I");
	_outputTree->Branch("sChEn", &_sChEn, "sChEn/F");
	_outputTree->Branch("sNeEn", &_sNeEn, "sNeEn/F");
	_outputTree->Branch("sChCluEn", &_sChCluEn, "sChCluEn/F");

	 _outputTree->Branch("sChHitEn", &_sChHitEn, "sChHitEn/F");
	 _outputTree->Branch("sNeHitEn", &_sNeHitEn, "sNeHitEn/F");

	_outputTree->Branch("LCNePFOEn", &_LCNePFOEn, "LCNePFOEn/F");
	_outputTree->Branch("LCChPFOEn", &_LCChPFOEn, "LCChPFOEn/F");
	_outputTree->Branch("LCCluEn", &_LCCluEn, "LCCluEn/F");
	_outputTree->Branch("Dis_SeedPro", &_Dis_SeedPro, "Dis_SeedPro/F");	//Seed Position between leading Charged/Neutral Cluster 
	_outputTree->Branch("PosCG_Ch", _PosCG_Ch, "PosCG_Ch[3]/F");
	_outputTree->Branch("PosCG_Ne", _PosCG_Ne, "PosCG_Ne[3]/F");
	_outputTree->Branch("PosCG_LNe", _PosCG_LNe, "PosCG_LNe[3]/F");

	_outputTree->Branch("HitCG_Ch", _HitCG_Ch, "HitCG_Ch[3]/F");
        _outputTree->Branch("HitCG_Ne", _HitCG_Ne, "HitCG_Ne[3]/F");

	_Num = 0; 
}


void SepEff::processEvent( LCEvent * evtP ) 
{		

	if (evtP) 								
	{	
		try 	
		{    
			_eventNr=evtP->getEventNumber();

			// input Collections: MCParticle, EcalBush, HcalBush, EcalHits, HcalHits
			_NNeMC = 0;
			_NChMC = 0;

			try{
				LCCollection * MCP = evtP->getCollection("MCParticle");
				int nMC = MCP->getNumberOfElements();
				TVector3 VtxPos, EndPPos;
				for(int i0 = 0; i0 < nMC; i0++)
				{
					MCParticle * a_MCP = dynamic_cast<MCParticle*>(MCP->getElementAt(i0));
					VtxPos = a_MCP->getVertex();
					EndPPos = a_MCP->getEndpoint();
					int PIDMCP=a_MCP->getPDG();
					int ChargeMCP = a_MCP->getCharge();
					if(((VtxPos.Perp() < 1750 && VtxPos.Perp() > 200 && fabs(VtxPos.Z()) < 2250) || (VtxPos.Perp() < 200 && fabs(VtxPos.Z()) < 1130)) && (EndPPos.Perp() > 1750 || fabs(EndPPos.Z()) > 2250 || (EndPPos.Perp() < 200 && fabs(EndPPos.Z()) > 1130)) && fabs(PIDMCP) != 12 &&  fabs(PIDMCP) != 14 &&  fabs(PIDMCP) != 16)
				{
					if(ChargeMCP == 0)
					{
						_NNeMC++;
					}
					else
					{
						_NChMC++;
					}
				}

				}
			}catch (lcio::DataNotAvailableException err) { }

			try{
				LCCollection * RecoP = evtP->getCollection("ArborPFOs");
				_NPFO = RecoP->getNumberOfElements();
				_TotalRecoEn = 0; 	
				_sChEn = 0; _sChCluEn = 0; _sNeEn = 0; _sChHitEn = 0; _sNeHitEn = 0; 
				_NCh = 0; _NNe = 0; _LCCluEn = 0; _LCNePFOEn = 0; _LCChPFOEn = 0; 

				float MaxNeEn = -1; 		
				float MaxChEn = -1; 

				TVector3 PosCh, PosNe, PosDiff, sumPosCh, sumPosNe, tmpPos; 
				PosCh.SetXYZ(0, 0, 0);
				PosNe.SetXYZ(0, 0, 0);
				sumPosCh.SetXYZ(0, 0, 0);
				sumPosNe.SetXYZ(0, 0, 0);
				Cluster* KK(0); 

				for(int p = 0; p < 3; p++)	
				{
					_HitCG_Ne[p] = 0;
					_HitCG_Ch[p] = 0;
				}
			
				for(int s = 0; s < _NPFO; s++)
				{
					ReconstructedParticle * a_RecoP = dynamic_cast<ReconstructedParticle*>(RecoP->getElementAt(s));
					float tmpEn = a_RecoP->getEnergy();

					for(unsigned int i0 = 0; i0 < a_RecoP->getClusters().size(); i0++)
					{
						Cluster * a_clu = a_RecoP->getClusters()[i0];

						for(unsigned int j0 = 0; j0 < a_clu->getCalorimeterHits().size(); j0++)
						{
							CalorimeterHit * a_hit = a_clu->getCalorimeterHits()[j0];
							if(a_RecoP->getCharge() != 0)
							{
								_sChHitEn += a_hit->getEnergy();
								_HitCG_Ch[0] += a_hit->getEnergy()*a_hit->getPosition()[0];
								_HitCG_Ch[1] += a_hit->getEnergy()*a_hit->getPosition()[1];
								_HitCG_Ch[2] += a_hit->getEnergy()*a_hit->getPosition()[2];
							}
							else
							{
								_sNeHitEn += a_hit->getEnergy();
								_HitCG_Ne[0] += a_hit->getEnergy()*a_hit->getPosition()[0];
                                                                _HitCG_Ne[1] += a_hit->getEnergy()*a_hit->getPosition()[1];
                                                                _HitCG_Ne[2] += a_hit->getEnergy()*a_hit->getPosition()[2];
							}
						}
					}

					if(a_RecoP->getClusters().size() > 0)
					{
						KK = a_RecoP->getClusters()[0];
					}
					_TotalRecoEn += tmpEn;
					
					if(a_RecoP->getType() == 501) continue;
					
					if(a_RecoP->getCharge() != 0)
					{
						_sChEn += tmpEn;
						if(a_RecoP->getClusters().size() > 0)
						{
							_sChCluEn += KK->getEnergy(); 
							tmpPos =  KK->getPosition();
							sumPosCh += KK->getEnergy() * tmpPos;
						}
						_NCh ++; 
						if(tmpEn > MaxChEn)
						{
							MaxChEn = tmpEn;
							_LCChPFOEn = tmpEn;
							if(a_RecoP->getClusters().size() > 0)
								PosCh = KK->getPosition();
						}
					}
					else
					{
						_sNeEn += tmpEn; 
						_NNe ++; 
						tmpPos = KK->getPosition();
						sumPosNe += tmpEn * tmpPos; 
						if(tmpEn > MaxNeEn)
						{
							MaxNeEn = tmpEn; 
							_LCNePFOEn = tmpEn; 
							if(a_RecoP->getClusters().size() > 0)
								PosNe = KK->getPosition();
						}
					}
				}

				if(_sNeEn > 0)
					sumPosNe = 1.0/_sNeEn*sumPosNe;
				if(_sChCluEn > 0)
					sumPosCh = 1.0/_sChCluEn*sumPosCh; 

				for(int p1 = 0; p1 < 3; p1++)
                                {
                                        _HitCG_Ch[p1] = _HitCG_Ch[p1] * 1.0/_sChHitEn;
					_HitCG_Ne[p1] = _HitCG_Ne[p1] * 1.0/_sNeHitEn;
                                }

				std::cout<<"Charged En: "<<_sChHitEn<<" : "<<_sChCluEn<<std::endl;
				std::cout<<"Neutral En: "<<_sNeHitEn<<" : "<<_sNeEn<<std::endl;

				_PosCG_Ch[0] = sumPosCh.X();
				_PosCG_Ch[1] = sumPosCh.Y();
				_PosCG_Ch[2] = sumPosCh.Z();

				_PosCG_Ne[0] = sumPosNe.X();
				_PosCG_Ne[1] = sumPosNe.Y();
				_PosCG_Ne[2] = sumPosNe.Z();

				_PosCG_LNe[0] = PosNe.X();
				_PosCG_LNe[1] = PosNe.Y();
				_PosCG_LNe[2] = PosNe.Z();

				PosDiff = PosNe - PosCh; 

				_Dis_SeedPro = PosDiff.Mag()*sin(PosDiff.Angle(PosCh));

				_outputTree->Fill();

			}catch (lcio::DataNotAvailableException err) { }

			_Num++; 
		}catch (lcio::DataNotAvailableException err) { }
	}  	

}	

void SepEff::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		delete tree_file;
	}

}


