#include <AnaTreeNew.hh>
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
#include <TVector3.h>

using namespace std;

AnaTreeNew aAnaTreeNew;
AnaTreeNew::AnaTreeNew()
	: Processor("AnaTreeNew"),
	_output(0)
{
	_description = "Arbor and Pandora Comparison for single charged particle" ;

	_treeFileName="AnaTreeNew.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);

	std::vector<std::string> inputMCParticle;
	inputMCParticle.push_back(std::string("MCParticle"));
	registerProcessorParameter("InputMCParticle" ,
			"Name of MC Particle Collections" ,
			_inputMCParticle,
			inputMCParticle);


	std::vector<std::string> inputArborCollections;
	inputArborCollections.push_back("ArborPFOs");
	registerProcessorParameter("InputArborCollection" ,
			"Name of Reconstructed Arbor Collections" ,
			_inputArborCollections,
			inputArborCollections);


	std::vector<std::string> inputPandoraCollections;
	inputPandoraCollections.push_back("PandoraPFOs");
	registerProcessorParameter("InputPandoraCollection" ,
			"Name of Reconstructed Pandora Collections" ,
			_inputPandoraCollections,
			inputPandoraCollections);


	_treeName="MCPart";
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


void AnaTreeNew::init() {

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
	_outputTree->Branch("nCHPFO_a", &nCHPFO_a, "nCHPFO_a/I");
	_outputTree->Branch("nCHPFO_p", &nCHPFO_p, "nCHPFO_p/I");
	_outputTree->Branch("nNEPFO_a", &nNEPFO_a, "nNEPFO_a/I");
        _outputTree->Branch("nNEPFO_p", &nNEPFO_p, "nNEPFO_p/I");
	_outputTree->Branch("nCHMCP", &nCHMCP, "nCHMCP/I");
	_outputTree->Branch("nNEMCP", &nNEMCP, "nNEMCP/I");
	
	_outputTree->Branch("LCEn",&_LCEn,"LCEn/F");
	_outputTree->Branch("TCEn",&_TCEn,"TCEn/F");
	
        _outputTree->Branch("EnPdP", &energyPdP, "EnPdP/F");
        _outputTree->Branch("PhiMCP", &phiMCP, "PhiMCP/F");
	_outputTree->Branch("EnMCP", &energyMCP, "EnMCP/F");
	_outputTree->Branch("ThetaMCP", &thetaMCP,"ThetaMCP/F");
	_outputTree->Branch("PMCP", &PMCP, "PMCP/F");
	_outputTree->Branch("MCPDGMother", &MCPDGMother, "MCPDGMother/I");

	_outputTree->Branch("ThetaPdP", &thetaPdP, "ThetaPdP/F");
        _outputTree->Branch("PhiPdP", &phiPdP, "PhiPdP/F");
	_outputTree->Branch("EnPdP", &energyPdP, "EnPdP/F");
	_outputTree->Branch("PdPID", &PdPID, "PdPID/I");
        _outputTree->Branch("nRecMu_p",&nRecMu_p,"nRecMu_p/I");
        _outputTree->Branch("nRecEle_p",&nRecEle_p,"nRecEle_p/I");
        _outputTree->Branch("nRecHad_p",&nRecHad_p,"nRecHad_p/I");
        _outputTree->Branch("nRecGamma_p",&nRecGamma_p,"nRecGamma_p/I");
        _outputTree->Branch("nRecOtherNeu_p",&nRecOtherNeu_p,"nRecOtherNeu_p/I");
        _outputTree->Branch("PdPIDMatrix",PdPIDMatrix,"PdPIDMatrix[5][5]/I");
	_outputTree->Branch("D0PdP", &PdD0, "D0PdP/F");
	_outputTree->Branch("Z0PdP", &PdZ0, "Z0PdP/F");
	_outputTree->Branch("dppPdP", &Pd_dpp, "dppPdP/F");
	_outputTree->Branch("LCEnPdP", &PdLCEn, "LCEnPdP/F");
	_outputTree->Branch("dEEPdP", &Pd_dEENeu, "dEEPdP/F");

	_outputTree->Branch("ThetaAbP", &thetaAbP, "ThetaAbP/F");
        _outputTree->Branch("PhiAbP", &phiAbP, "PhiAbP/F");
	_outputTree->Branch("EnAbP", &energyAbP, "EnAbP/F");
	_outputTree->Branch("AbPID", &AbPID, "AbPID/I");
	_outputTree->Branch("nRecMu_a",&nRecMu_a,"nRecMu_a/I");
        _outputTree->Branch("nRecEle_a",&nRecEle_a,"nRecEle_a/I");
        _outputTree->Branch("nRecHad_a",&nRecHad_a,"nRecHad_a/I");
        _outputTree->Branch("nRecGamma_a",&nRecGamma_a,"nRecGamma_a/I");
        _outputTree->Branch("nRecOtherNeu_a",&nRecOtherNeu_a,"nRecOtherNeu_a/I");
	_outputTree->Branch("AbPIDMatrix",AbPIDMatrix,"AbPIDMatrix[5][5]/I");
        _outputTree->Branch("D0AbP", &AbD0, "D0AbP/F");
        _outputTree->Branch("Z0AbP", &AbZ0, "Z0AbP/F");
        _outputTree->Branch("dppAbP", &Ab_dpp, "dppAbP/F");
        _outputTree->Branch("LCEnAbP", &AbLCEn, "LCEnAbP/F");
        _outputTree->Branch("dEEAbP", &Ab_dEENeu, "dEEAbP/F");

	_outputAPFO = new TTree("ArborPFOs","ArborPFOs");
	_outputAPFO->Branch("EventNr", &_eventNr, "EventNr/I");
        _outputAPFO->Branch("Num", &_Num, "Num/I");
	_outputAPFO->Branch("TypeA", &TypeA, "TypeA/I");	
	_outputAPFO->Branch("ChargeA", &ChargeA, "ChargeA/I");
	_outputAPFO->Branch("EnergyA", &EnergyA, "EnergyA/F");
	_outputAPFO->Branch("PA", PA, "PA[3]/F");
	_outputAPFO->Branch("CluEnA", &CluEnA, "CluEnA/F");
	_outputAPFO->Branch("CluEnComA", CluEnComA, "CluEnComA[2]/F");
	_outputAPFO->Branch("TrackHitA", &TrackHitA, "TrackHitA/I");
        _outputAPFO->Branch("StartPosA", StartPosA, "StartPosA[3]/F");
        _outputAPFO->Branch("EndPosA", EndPosA, "EndPosA[3]/F");
	

        _outputPPFO = new TTree("PandoraPFOs","PandoraPFOs");
	_outputPPFO->Branch("EventNr", &_eventNr, "EventNr/I");
        _outputPPFO->Branch("Num", &_Num, "Num/I");
	_outputPPFO->Branch("TypeP", &TypeP, "TypeP/I");	
	_outputPPFO->Branch("ChargeP", &ChargeP, "ChargeP/I");
	_outputPPFO->Branch("EnergyP", &EnergyP, "EnergyP/F");
	_outputPPFO->Branch("PP", PP, "P[3]/F");
	_outputPPFO->Branch("CluEnP", &CluEnP, "CluEnP/F");
	_outputPPFO->Branch("CluEnComP", CluEnComP, "CluEnComP[2]/F");
	_outputPPFO->Branch("TrackHitP", &TrackHitP, "TrackHitP/I");
        _outputPPFO->Branch("StartPosP", StartPosP, "StartPosP[3]/F");
        _outputPPFO->Branch("EndPosP", EndPosP, "EndPosP[3]/F");


	_Num = 0;
}

void AnaTreeNew::processEvent( LCEvent * evtP ) 
{
	if (evtP) 								
	{
		if(_Num%10 == 0) cout << "Event:" <<_Num<<endl;	

		_eventNr = evtP->getEventNumber();
		int index1 = -1;
		nCHPFO_a = 0; 
		nCHPFO_p = 0;
		nNEPFO_a = 0; 
		nNEPFO_p = 0;
		nCHMCP = 0;
		nNEMCP = 0;
		thetaMCP = -100;
		phiMCP = -100;
		energyMCP = -100;
		MCCharge = -2;
		thetaPdP = -100;
		phiPdP = -100;
		energyPdP = -100;
		thetaAbP = -100;
		phiAbP = -100;
		energyAbP = -100;
		MCPDGMother = 0;
		PdPID = 0;
		AbPID = 0;
		nRecMu_a = 0;
		nRecEle_a = 0;
		nRecHad_a = 0;
		nRecGamma_a = 0;
		nRecOtherNeu_a = 0;
		_LCEn = 0;
		_TCEn = 0;
		nRecMu_p = 0;
                nRecEle_p = 0;
                nRecHad_p = 0;
                nRecGamma_p = 0;
                nRecOtherNeu_p = 0;

		PdD0 = -10000;
		PdZ0 = -10000;
		PdLCEn = -100;
		Pd_dpp = -100;
		Pd_dEENeu = -100;

                AbD0 = -10000;
                AbZ0 = -10000;
                AbLCEn = -100;
                Ab_dpp = -100;
                Ab_dEENeu = -100;

		for(int i = 0;i<6;i++)
		{
			for(int k = 0; k < 5; k++)
			{
				AbPIDMatrix[i][k] = 0;
				PdPIDMatrix[i][k] = 0;
			}
		}


		try{
			LCCollection * MCP = evtP->getCollection("MCParticle");
			for(int s0 = 0; s0 < MCP->getNumberOfElements(); s0++)
			{
				
				MCParticle *a1_MCP = dynamic_cast<EVENT::MCParticle *>(MCP->getElementAt(s0));
                              	
				
				if(a1_MCP->getCharge() == 0) nNEMCP++;
				else nCHMCP++;

				if(s0 == 0)
				{
					TVector3 MCPM(a1_MCP->getMomentum()[0],a1_MCP->getMomentum()[1],a1_MCP->getMomentum()[2]);
					thetaMCP = MCPM.Theta();
					phiMCP = MCPM.Phi();
					PMCP = MCPM.Mag();
					energyMCP = a1_MCP->getEnergy();
					MCPDGMother = a1_MCP->getPDG();
					MCCharge = a1_MCP->getCharge();
				}
	
			}
			if(MCCharge != 0)
			{
				if(abs(MCPDGMother) == 13) index1 = 0;
				else if(abs(MCPDGMother) == 11) index1 = 1;
				else index1 = 2;
			}
			else
			{
				if(MCPDGMother == 22) index1 = 3;
				else index1 = 4;
			}
			

		}catch(lcio::DataNotAvailableException err) { }

		TVector3 MCDir(1,1,1);
		MCDir.SetTheta(thetaMCP);
		MCDir.SetPhi(phiMCP);

		try{

			LCCollection* col_RecoNeP = evtP->getCollection( "ArborPFOs" );
			float minAngA = 1.0E10;
			float minEnNeu = -1;
			for(int i1 = 0; i1 < col_RecoNeP->getNumberOfElements(); i1++)
			{
				for(int i0 = 0; i0 < 3; i0++)
				{
					PA[i0] = 0;
					StartPosA[i0] = 0;
					EndPosA[i0] = 0;
					if(i0 < 2)
					{
						CluEnComA[i0] = 0;
					}
				}

				ReconstructedParticle *a_RecoP = dynamic_cast<EVENT::ReconstructedParticle *>(col_RecoNeP->getElementAt(i1));	
				//if(a_RecoP->getClusters().size() == 0) continue;
				ChargeA = a_RecoP->getCharge();
				CluEnA = 0;
				if(ChargeA)
				{
					nCHPFO_a ++;
				}
				else
				{
					nNEPFO_a ++;
				}

				EnergyA = a_RecoP->getEnergy();
				PA[0] = a_RecoP->getMomentum()[0];
				PA[1] = a_RecoP->getMomentum()[1];
				PA[2] = a_RecoP->getMomentum()[2];

				TypeA = a_RecoP->getType();
				
                                if(abs(TypeA) == 13) nRecMu_a++; 
                                else if(abs(TypeA) == 11) nRecEle_a++;
                                else if(ChargeA != 0) nRecHad_a++;   
                                else if(abs(TypeA == 22)) nRecGamma_a++;
                                else nRecOtherNeu_a++;

                                if(ChargeA)
                                {
                                        Track * a_Trk = a_RecoP->getTracks()[0];
                                        TrackHitA = a_Trk->getTrackerHits().size();
                                        StartPosA[0] = (a_Trk->getTrackerHits()[0])->getPosition()[0];
                                        StartPosA[1] = (a_Trk->getTrackerHits()[0])->getPosition()[1];
                                        StartPosA[2] = (a_Trk->getTrackerHits()[0])->getPosition()[2];

                                        EndPosA[0] = (a_Trk->getTrackerHits()[TrackHitA - 1])->getPosition()[0];
                                        EndPosA[1] = (a_Trk->getTrackerHits()[TrackHitA - 1])->getPosition()[1];
                                        EndPosA[2] = (a_Trk->getTrackerHits()[TrackHitA - 1])->getPosition()[2];
                                }

				if(a_RecoP->getClusters().size() > 0)
				{
					Cluster * currClu = a_RecoP->getClusters()[0];
					CluEnA = currClu->getEnergy();
					
					_TCEn += CluEnA;
					if(CluEnA > _LCEn && ChargeA != 0) _LCEn = CluEnA;

					if(!ChargeA)
					{
						//CluEnComA[0] = currClu->getSubdetectorEnergies()[0];
						//CluEnComA[1] = currClu->getSubdetectorEnergies()[1];
					}
				}

				if(ChargeA && MCCharge)
				{
					TVector3 AbPM(a_RecoP->getMomentum()[0],a_RecoP->getMomentum()[1],a_RecoP->getMomentum()[2]);
					Track* a_trk = a_RecoP->getTracks()[0];
					float DifAngAMC = AbPM.Angle(MCDir);
					if(DifAngAMC < minAngA && a_RecoP->getClusters().size() > 0)
					{
						AbD0 = a_trk->getD0();
						AbZ0 = a_trk->getZ0();
						Ab_dpp = (PMCP-AbPM.Mag())/PMCP;
						AbPID = TypeA;
						thetaAbP = AbPM.Theta();
						phiAbP = AbPM.Phi();
						energyAbP = EnergyA;
					} 
				}

				if(ChargeA == 0 && MCCharge == 0)
				{
					TVector3 AbPM(a_RecoP->getMomentum()[0],a_RecoP->getMomentum()[1],a_RecoP->getMomentum()[2]);
					if(EnergyA > minEnNeu)
					{
						AbPID = TypeA;
						AbLCEn = EnergyA;
						Ab_dEENeu = (energyMCP - EnergyA)/energyMCP;
					}
				}				

				_outputAPFO->Fill();	
			}

		}catch(lcio::DataNotAvailableException err) { }

		try{
			LCCollection* col_RecoPandora = evtP->getCollection( "PandoraPFOs" );
			float minAngP = 1.0E10;
			float minEnNeu = -1;
			for(int i2 = 0; i2 < col_RecoPandora->getNumberOfElements(); i2++)
			{
				for(int i0 = 0; i0 < 3; i0++)
				{
					PP[i0] = 0;
					StartPosP[i0] = 0;
					EndPosP[i0] = 0;
					if(i0 < 2)
					{
						CluEnComP[i0] = 0;
					}
				}
				ReconstructedParticle *a_RecoP = dynamic_cast<EVENT::ReconstructedParticle *>(col_RecoPandora->getElementAt(i2));	
				if(a_RecoP->getClusters().size() == 0) continue;
				ChargeP = a_RecoP->getCharge();
				CluEnP = 0;
				if(ChargeP)
				{
					nCHPFO_p ++;
				}
				else
				{
					nNEPFO_p ++;
				}

				EnergyP = a_RecoP->getEnergy();
				PP[0] = a_RecoP->getMomentum()[0];
				PP[1] = a_RecoP->getMomentum()[1];
				PP[2] = a_RecoP->getMomentum()[2];

				TypeP = a_RecoP->getType();
				
				if(abs(TypeP) == 13) nRecMu_p++;
				else if(abs(TypeP) == 11) nRecEle_p++;
				else if(ChargeP != 0) nRecHad_p++;
				else if(abs(TypeP == 22)) nRecGamma_p++;
				else nRecOtherNeu_p++;		
		
                                if(ChargeP)
                                {
                                        Track * a_Trk = a_RecoP->getTracks()[0];
                                        TrackHitP = a_Trk->getTrackerHits().size();
                                        StartPosP[0] = (a_Trk->getTrackerHits()[0])->getPosition()[0];
                                        StartPosP[1] = (a_Trk->getTrackerHits()[0])->getPosition()[1];
                                        StartPosP[2] = (a_Trk->getTrackerHits()[0])->getPosition()[2];

                                        EndPosP[0] = (a_Trk->getTrackerHits()[TrackHitP - 1])->getPosition()[0];
                                        EndPosP[1] = (a_Trk->getTrackerHits()[TrackHitP - 1])->getPosition()[1];
                                        EndPosP[2] = (a_Trk->getTrackerHits()[TrackHitP - 1])->getPosition()[2];
                                }

				if(a_RecoP->getClusters().size() > 0)
				{
					Cluster * currClu = a_RecoP->getClusters()[0];
					CluEnP = currClu->getEnergy();

					if(!ChargeP)
					{
						CluEnComP[0] = currClu->getSubdetectorEnergies()[0];
						CluEnComP[1] = currClu->getSubdetectorEnergies()[1];
					}
				}
				if(ChargeP && MCCharge)
				{
					TVector3 PdPM(a_RecoP->getMomentum()[0],a_RecoP->getMomentum()[1],a_RecoP->getMomentum()[2]);
                                        Track* a_trk = a_RecoP->getTracks()[0];
					float DifAngAMC = PdPM.Angle(MCDir);
					if(DifAngAMC < minAngP && a_RecoP->getClusters().size() > 0)
					{
                                                PdD0 = a_trk->getD0();
                                                PdZ0 = a_trk->getZ0();
                                                Pd_dpp = (PMCP-PdPM.Mag())/PMCP;
						PdPID = TypeP;
						thetaPdP = PdPM.Theta();
						phiPdP = PdPM.Phi();
						energyPdP = EnergyP;
					} 
				}
				if(ChargeP == 0 && MCCharge == 0)
				{
                                        TVector3 PdPM(a_RecoP->getMomentum()[0],a_RecoP->getMomentum()[1],a_RecoP->getMomentum()[2]);
                                        if(EnergyP > minEnNeu)
                                        {
                                                PdPID = TypeP;
                                                PdLCEn = EnergyP;
                                                Pd_dEENeu = (energyMCP - EnergyP)/energyMCP;
                                        }

				}

				_outputPPFO->Fill();

				
			}   
		
		}catch(lcio::DataNotAvailableException err) { }
                
                if(index1 > -1)
		{	
			if(MCCharge)
 			{
				if(abs(AbPID) == 13) AbPIDMatrix[index1][0] = 1;
				else if(abs(AbPID) == 11) AbPIDMatrix[index1][1] = 1;
				else  AbPIDMatrix[index1][2] = 1;

                                if(abs(PdPID) == 13) PdPIDMatrix[index1][0] = 1;
                                else if(abs(PdPID) == 11) PdPIDMatrix[index1][1] = 1;
                                else  PdPIDMatrix[index1][2] = 1;
			}
			else
			{
				if(AbPID == 22) AbPIDMatrix[index1][3] = 1;
				else AbPIDMatrix[index1][4] = 1;

                                if(PdPID == 22) PdPIDMatrix[index1][3] = 1;
                                else PdPIDMatrix[index1][4] = 1;
			}
		}
	        	
		
		
		_outputTree->Fill();
		_Num++;



	}


}

void AnaTreeNew::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		delete tree_file;
	}

}
