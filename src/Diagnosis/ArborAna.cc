/*
 * used to readout the FD, FD_10First, FD_20Later, ... and other quantities of bush
 * objective: pattern tagging for ECAL as: penetrating mip, deep interaction CH, EM Cluster; 
 * 				  HCAL as: penetrating mip, early interaction H, later interaction H, etc.
 *
 * */

#include <ArborAna.hh>
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

ArborAna aArborAna ;
ArborAna::ArborAna()
	: Processor("ArborAna"),
	_output(0)
{
	_description = "Measure Bush Quantities" ;

	_treeFileName="ArborAna.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);

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

	std::vector<std::string> inputArborParticle;
        inputArborParticle.push_back(std::string("ArborTrkParticle_HQ_Barrel"));
        inputArborParticle.push_back(std::string("ArborTrkParticle_HQ_Endcap"));
        inputArborParticle.push_back(std::string("ArborTrkParticle_MQ_Barrel"));
	inputArborParticle.push_back(std::string("ArborTrkParticle_MQ_Endcap"));
        inputArborParticle.push_back(std::string("ArborTrkParticle_MQ_Vtx"));
        inputArborParticle.push_back(std::string("ArborTrkParticle_LQ"));
	inputArborParticle.push_back(std::string("ArborNeutralParticle"));
        registerProcessorParameter("InputArborCollection" ,
                        "Name of Reconstructed Arbor Collections" ,
                        _inputArborParticle,
                        inputArborParticle);

	_Num=0;

}

void ArborAna::init() {

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
	//MCPart
	_outputTree->Branch("Theta_Vis", &_ThetaVis, "Theta_Vis/F");    //Direction Tau Visible part
	_outputTree->Branch("Phi_Vis", &_PhiVis, "Phi_Vis/F");

	_outputTree->Branch("N_ChargeMCP", &_NChargeMCP, "N_ChargeMCP/I");
	_outputTree->Branch("N_NeutralMCP", &_NNeutralMCP, "N_NeutralMCP/I");
	_outputTree->Branch("N_Muon", &_NMuon, "NMuon/I");
	_outputTree->Branch("N_Pion", &_NPion, "NPion/I");
	_outputTree->Branch("N_EM", &_NEM, "NEM/I");
	_outputTree->Branch("N_Photon", &_NPhoton, "NPhoton/I");
	_outputTree->Branch("N_Neutrino", &_NNeutrino, "NNeutrino/I");
	_outputTree->Branch("E_ChargeMCP", &_EChargeMCP, "E_ChargeMCP/F");
	_outputTree->Branch("E_NeutralMCP", &_ENeutralMCP, "E_NeutralMCP/F");
	_outputTree->Branch("E_Neutrino", &_ENeutrino, "E_Neutrino/F");	
	_outputTree->Branch("EnSeg", _EnSeg, "EnSeg[10]/F");

	_outputTree->Branch("EnIsoFrag", _EnIsoFrag, "EnIsoFrag[2]/F");
	_outputTree->Branch("EnIsoHit", _EnIsoHit, "EnIsoHit[2]/F");		//EnCluster? 

	_outputTree->Branch("NE_sR", &_NE_sR, "NE_sR/F");
	_outputTree->Branch("CE_sR", &_CE_sR, "CE_sR/F");
	_outputTree->Branch("InvE_sR", &_InvE_sR, "InvE_sR/F");

	_outputTree->Branch("NE_lR", &_NE_lR, "NE_lR/F");
	_outputTree->Branch("CE_lR", &_CE_lR, "CE_lR/F");
	_outputTree->Branch("InvE_lR", &_InvE_lR, "InvE_lR/F");

	_outputTree->Branch("QuarkP", _QuarkP, "QuarkP[3]/F");
	_outputTree->Branch("CosTheta", &_CosTheta, "CosTheta/F");

	_outputTree->Branch("LeadingCMCP", _LeadingCMCP, "LeadingCMCP[3]/F");
	_outputTree->Branch("L_MCPID", &_L_MCPID, "L_MCPID/F");
	_outputTree->Branch("L_MCPEn", &_L_MCPEn, "L_MCPEn/F");
	_outputTree->Branch("L_ChargeArborP", _L_ChargeArborP, "L_ChargeArborP[3]/F");
	_outputTree->Branch("L_ChargeArborEn", &_L_ChargeArborEn, "L_ChargeArborEn/F");

	//Reco
	_outputTree->Branch("N_ChargeArbor", &_NChargeArbor, "N_ChargeArbor/I");
	_outputTree->Branch("N_NeutralArbor", &_NNeutralArbor, "N_NeutralArbor/I");
	_outputTree->Branch("N_EffChargeArbor", &_NEffChargeArbor, "N_EffChargeArbor/I");
	_outputTree->Branch("N_EffNeutralArbor", &_NEffNeutralArbor, "N_EffNeutralArbor/I");

	_outputTree->Branch("E_ChargeArbor", &_EChargeArbor, "E_ChargeArbor/F");
	_outputTree->Branch("E_NeutralArbor", &_ENeutralArbor, "E_NeutralArbor/F");
	_outputTree->Branch("TotalEVis", &_TotalEVis, "TotalEVis/F");
	_outputTree->Branch("E_noCluTrk", &_E_noCluTrk, "E_noCluTrk/F");
	_outputTree->Branch("HQ_Charged", &_HQ_Charged, "HQ_Charged/F");
	_outputTree->Branch("TotalPVis", _TotalPVis, "TotalPVis[3]/F");
	_outputTree->Branch("CoM", &_CoM, "CoM/F");

	_outputTree->Branch("InvMassVis", &_InvMassVis, "InvMassVis/F");
	_outputTree->Branch("Theta_Reco", &_ThetaReco, "Theta_Reco/F");    //Direction Tau
	_outputTree->Branch("Phi_Reco", &_PhiReco, "Phi_Reco/F");

	_outputTree->Branch("Theta_Calo", &_Theta_Calo, "Theta_Calo/F");
	_outputTree->Branch("Theta_Clu", &_Theta_Clu, "Theta_Clu/F");
	_outputTree->Branch("Phi_Calo", &_Phi_Calo, "Phi_Calo/F");
	_outputTree->Branch("Phi_Clu", &_Phi_Clu, "Phi_Clu/F");

	_outputTree->Branch("PandoraTotalEn", &_PandoraTotalEn, "PandoraTotalEn/F");
	_outputTree->Branch("PandoraChEn", &_PandoraChEn, "PandoraChEn/F");
	_outputTree->Branch("PandoraNeEn", &_PandoraNeEn, "PandoraNeEn/F");

	_outputPFO = new TTree("PFO", "PFO");
	_outputPFO->Branch("EventNr", &_eventNr, "EventNr/I");
	_outputPFO->Branch("Num", &_Num, "Num/I");
	_outputPFO->Branch("Index", &_Index, "Index/I");
	_outputPFO->Branch("Charge", &_Charge, "Charge/F");
	_outputPFO->Branch("TrkEn", &_TrkEn, "TrkEn/F");
	_outputPFO->Branch("CluEn", &_CluEn, "CluEn/F");
	_outputPFO->Branch("NClu", &_NClu, "NClu/I");
}

void ArborAna::processEvent( LCEvent * evtP ) 
{		

	if (evtP) 								
	{	
		try 	
		{    
			_eventNr=evtP->getEventNumber();
			if( _eventNr % 10 == 0)
				std::cout<<"ArborAna: "<<_eventNr<<" evts processed"<<std::endl;

			_Num++;

			_ThetaVis = 0; _PhiVis = 0;
			_NChargeMCP = 0; _NNeutralMCP = 0; _NMuon = 0; _NPion = 0; _NPion0 = 0;  _NEM = 0; _NPhoton = 0; _NNeutrino = 0;
			_EChargeMCP = 0; _ENeutralMCP = 0; _ENeutrino = 0;

			_NChargeArbor = 0; _NNeutralArbor = 0; 
			_NEffChargeArbor = 0; _NEffNeutralArbor = 0;
			_EChargeArbor = 0; _ENeutralArbor = 0; 
			_TotalEVis = 0; _InvMassVis = 0; 
			_ThetaReco = 0; _PhiReco = 0; 
			_Theta_Calo = 0; _Phi_Calo = 0; 
			_Theta_Clu = 0; _Phi_Clu = 0; 

			_NE_sR = 0; _NE_lR = 0; _CE_sR = 0; _CE_lR = 0; _InvE_sR = 0; _InvE_lR = 0;

			_L_MCPID = 0; _L_MCPEn = 0; _L_ChargeArborEn = 0; 

			int NArbor = 0;

			for(int p1 = 0; p1 < 10; p1++)
			{
				_EnSeg[p1] = 0;
				if(p1 < 3)
				{
					_LeadingCMCP[p1] = 0;
					_L_ChargeArborP[p1] = 0;
				}
			}

			LCCollection * IsoFrag_1 = evtP->getCollection("IsoFragEcalBushes");
			LCCollection * IsoFrag_2 = evtP->getCollection("IsoFragHcalBushes");
			// LCCollection * PandoraPFO = evtP->getCollection("PandoraPFOs");

			LCCollection * IsoHit_1 = evtP->getCollection("EcalIsoHits");
			LCCollection * IsoHit_2 = evtP->getCollection("HcalIsoHits");

			_EnIsoFrag[0] = 0;
			_EnIsoFrag[1] = 0;
			_EnIsoHit[0] = 0;
			_EnIsoHit[1] = 0;

			for(int i3 = 0; i3 < IsoFrag_1->getNumberOfElements(); i3++)
			{
				Cluster * a_clu = dynamic_cast<Cluster*>(IsoFrag_1 -> getElementAt(i3));
				_EnIsoFrag[0] += a_clu->getEnergy();
			}

			for(int i4 = 0; i4 < IsoFrag_2->getNumberOfElements(); i4++)
			{
				Cluster * a_clu = dynamic_cast<Cluster*>(IsoFrag_2 -> getElementAt(i4));
				_EnIsoFrag[1] += a_clu->getEnergy();
			}

			for(int j3 = 0; j3 < IsoHit_1->getNumberOfElements(); j3++)
			{
				CalorimeterHit * a_hit = dynamic_cast<CalorimeterHit*>(IsoHit_1->getElementAt(j3));
				_EnIsoHit[0] += a_hit->getEnergy();
			}

			for(int j4 = 0; j4 < IsoHit_2->getNumberOfElements(); j4++)
			{
				CalorimeterHit * a_hit = dynamic_cast<CalorimeterHit*>(IsoHit_2->getElementAt(j4));
				_EnIsoHit[1] += a_hit->getEnergy();
			}


			_PandoraTotalEn = 0;
			_PandoraChEn = 0;
			_PandoraNeEn = 0; 

			float RecoEn = 0; 
			
			/*
			for(int i5 = 0; i5 < PandoraPFO->getNumberOfElements(); i5++)
			{
				ReconstructedParticle* a_PFO = dynamic_cast<ReconstructedParticle*>(PandoraPFO->getElementAt(i5));
				RecoEn = a_PFO->getEnergy();
				_PandoraTotalEn += RecoEn;
				if( fabs(a_PFO->getCharge()) > 0.1 )
				{
					_PandoraChEn += RecoEn;
				}
				else
				{
					_PandoraNeEn += RecoEn;
				}
			}
			*/			

			for(unsigned int p0 = 0; p0 < _inputArborParticle.size(); p0++)
			{

				LCCollection * RecoCol = evtP->getCollection( _inputArborParticle[p0].c_str() );			
				NArbor = RecoCol->getNumberOfElements();

				_E_noCluTrk = 0;
				_HQ_Charged = 0;
				for(int k0 = 0; k0 < 3; k0++)
				{
					_TotalPVis[k0] = 0;
					_QuarkP[k0] = 0;
				}
				_CosTheta = -10; 
				_CoM = 0; 

				for(int i0 = 0; i0 < NArbor; i0++)
				{

					ReconstructedParticle* a_Arbor = dynamic_cast<ReconstructedParticle*>(RecoCol->getElementAt(i0));
					RecoEn = a_Arbor->getEnergy();

					_EnSeg[p0] += RecoEn; 

					_Charge = a_Arbor->getCharge();

					_NClu = a_Arbor->getClusters().size();
					_TotalPVis[0] += a_Arbor->getMomentum()[0];
					_TotalPVis[1] += a_Arbor->getMomentum()[1];
					_TotalPVis[2] += a_Arbor->getMomentum()[2];

					_CluEn = 0; 
					_Index = i0;

					// cout<<"Before Cluster Energy "<<endl; 

					for(int j0 = 0; j0 < _NClu; j0++)
					{
						Cluster * a_Clu = a_Arbor->getClusters()[j0];
						_CluEn += a_Clu->getEnergy();
					}

					// cout<<"After Cluster Energy "<<endl;

					if(fabs(_Charge) > 0.1)
					{

						if(RecoEn > _L_ChargeArborEn)
						{
							_L_ChargeArborEn = RecoEn; 
                        				_L_ChargeArborP[0] = a_Arbor->getMomentum()[0];
							_L_ChargeArborP[1] = a_Arbor->getMomentum()[1];
							_L_ChargeArborP[2] = a_Arbor->getMomentum()[2];
						}

						_EChargeArbor += RecoEn; 
						_NChargeArbor ++;
						_TrkEn = RecoEn; 
						if(_NClu == 0)
						{
							_E_noCluTrk += RecoEn; 
						}
						else
						{
							_HQ_Charged += RecoEn; 
						}

						if(RecoEn > 1.0)
						{
							_NEffChargeArbor++;
						}
					}       
					else    
					{       
						_ENeutralArbor += RecoEn;
						_NNeutralArbor ++;
						_TrkEn = 0; 
						//_CluEn = RecoEn; 
						if(RecoEn > 1.0)
						{
							_NEffNeutralArbor++;
						}
					}               

					// ArborMom += a_Arbor->getMomentum();
					_TotalEVis += RecoEn;

					_outputPFO->Fill();
				}

				_CoM = sqrt(_TotalEVis * _TotalEVis - _TotalPVis[0]*_TotalPVis[0] - _TotalPVis[1]*_TotalPVis[1] - _TotalPVis[2]*_TotalPVis[2] );

			}

			/*
			   std::vector<MCParticle*> FirstGeneDaughter;
			   FirstGeneDaughter.clear();
			   TVector3 TauMom, TauVisMom, ArborMom, CluMom, CluCaloMom; 
			   TauMom.SetXYZ(0, 0, 0);
			   TauVisMom.SetXYZ(0, 0, 0);
			   ArborMom.SetXYZ(0, 0, 0);
			   CluMom.SetXYZ(0, 0, 0);
			   CluCaloMom.SetXYZ(0, 0, 0);
			   float TauP = 0; 
			   */

			LCCollection * MCPCol = evtP->getCollection("MCParticle");
			int NMCP = MCPCol->getNumberOfElements();

			int PID = 0;
			float MCPEn = 0; 
			float MCPCharge = 0; 
			TVector3 VTX, EndP; 


			for(int i1 = 0; i1 < NMCP; i1++)
			{
				MCParticle *b_MCP = dynamic_cast<MCParticle*>(MCPCol->getElementAt(i1));
				VTX = b_MCP->getVertex();
				EndP = b_MCP->getEndpoint();
				MCPCharge = b_MCP->getCharge();
				PID = b_MCP->getPDG();
				MCPEn = b_MCP->getEnergy();

				if(abs(PID) < 6 && b_MCP->getParents().size() == 0)
				{
					_QuarkP[0] = b_MCP->getMomentum()[0];
					_QuarkP[1] = b_MCP->getMomentum()[1];
					_QuarkP[2] = b_MCP->getMomentum()[2];

					_CosTheta = _QuarkP[2]/MCPEn; 
				}

				if(b_MCP->getParents().size() == 0 && MCPEn > _L_MCPEn && abs(MCPCharge) > 0.001 )
				{
					_L_MCPEn = MCPEn;
					_L_MCPID = PID; 
					_LeadingCMCP[0] = b_MCP->getMomentum()[0];
					_LeadingCMCP[1] = b_MCP->getMomentum()[1];
					_LeadingCMCP[2] = b_MCP->getMomentum()[2];
				}

				if(VTX.Mag() < 1 && EndP.Mag() > 1)
				{
					if(fabs(MCPCharge) > 0.001)
					{
						_CE_sR += MCPEn; 
					}
					else if( abs(PID) == 12 || abs(PID) == 14 || abs(PID) == 16 )
					{
						_NNeutrino ++;
						_InvE_sR += MCPEn;
					}
					else
					{
						_NE_sR += MCPEn; 
					}
				}

				if(VTX.Mag() < 1800 && EndP.Mag() > 1800)
				{
					if(fabs(MCPCharge) > 0.001)
					{
						_CE_lR += MCPEn;
					}
					else if( abs(PID) == 12 || abs(PID) == 14 || abs(PID) == 16 )
					{
						_InvE_lR += MCPEn;
					}
					else
					{
						_NE_lR += MCPEn;
					}
				}
			}

			// cout<<"Before Tree Fill"<<endl;

			_outputTree->Fill();

		}		
		catch (lcio::DataNotAvailableException err) { }
	}  	

}	

void ArborAna::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		delete tree_file;
	}

}


