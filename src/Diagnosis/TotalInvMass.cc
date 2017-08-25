#include <TotalInvMass.hh>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/MCParticleImpl.h>
#include <values.h>
#include <string>
#include <iostream>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <stdexcept>
#include <TFile.h> 
#include <TTree.h>
#include <TVector3.h>
#include <TRandom.h>
#include <Rtypes.h> 
#include <sstream>		
#include <cmath>
#include <vector>
#include <TMath.h>
#include "TLorentzVector.h"

using namespace std;

const float sqrts = 250.0; 	//GeV

TotalInvMass a_TotalInvMass_instance;

TotalInvMass::TotalInvMass()
	: Processor("TotalInvMass"),
	_output(0)
{
	_description = "Print MC Truth" ;

	_treeFileName="MCTruth.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);

	/*
	_colName="MCParticle";
	registerProcessorParameter( "MCObjects" ,
			"The name of the PFOs" ,
			_colName ,
			_colName);
	*/
	
	_treeName="MCPart";
	registerProcessorParameter( "TreeName" , 
			"The name of the ROOT tree" ,
			_treeName ,
			_treeName);

	/*
	_leptonID = 13;
	registerProcessorParameter( "LeptonIDTag" ,
                        "Lepton ID that will be used in this analysis." ,
                        _leptonID ,
                        _leptonID);
	*/	

	_overwrite=0;
	registerProcessorParameter( "OverwriteFile" , 
			"If zero an already existing file will not be overwritten." ,
			_overwrite ,
			_overwrite);

}

void TotalInvMass::init() {

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
	
	_outputTree->Branch("OriQuarkID", &_OriQuarkID, "OriQuarkID/I");

	_outputTree->Branch("J1CosTheta", &_J1CosTheta, "J1CosTheta/F");
	_outputTree->Branch("J2CosTheta", &_J2CosTheta, "J2CosTheta/F");

	_outputTree->Branch("OQDir", &_OQDir, "OQDir/F");
	_outputTree->Branch("HDPID", &_HDPID, "HDPID/I");	//Higgs Daughter Type, Direction??
	_outputTree->Branch("HDir", &_HDir, "HDir/F");
	_outputTree->Branch("Mass_a", &_Mass_a, "Mass_a/F");
	_outputTree->Branch("Mass_p", &_Mass_p, "Mass_p/F");
	_outputTree->Branch("Mass_p_Pisr", &_Mass_p_Pisr, "Mass_p_Pisr");
	_outputTree->Branch("Mass_a_Pisr", &_Mass_a_Pisr, "Mass_a_Pisr");
	_outputTree->Branch("Mass_a_Plcal",&_Mass_a_Plcal, "Mass_a_Plcal");

	_outputTree->Branch("TotalP_a", TotalP_a, "TotalP_a[4]/F");
	_outputTree->Branch("ChP", ChP, "ChP[4]/F");
	_outputTree->Branch("PhP", PhP, "PhP[4]/F");	//Gamma
	_outputTree->Branch("NeP", NeP, "NeP[4]/F");
	_outputTree->Branch("FrP", FrP, "FrP[4]/F");

	_outputTree->Branch("FrPh", FrPh, "FrPh[4]/F");
	_outputTree->Branch("FrNe", FrNe, "FrNe[4]/F");
	
	_outputTree->Branch("KPF", KPF, "KPF[4]/F");

	_outputTree->Branch("UdP", UdP, "UdP[4]/F");

	_outputTree->Branch("TotalP_p", TotalP_p, "TotalP_p[4]/F");
	_outputTree->Branch("nCHPFO_a", &nCHPFO_a, "nCHPFO_a/I");
	_outputTree->Branch("nCHPFO_p", &nCHPFO_p, "nCHPFO_p/I");
	_outputTree->Branch("nNEPFO_a", &nNEPFO_a, "nNEPFO_a/I");
        _outputTree->Branch("nNEPFO_p", &nNEPFO_p, "nNEPFO_p/I");
	_outputTree->Branch("NeCaloE_a", NeCaloE_a, "NeCaloE_a[2]/F");
        _outputTree->Branch("NeCaloE_p", NeCaloE_p, "NeCaloE_p[2]/F");
	_outputTree->Branch("ElargeP", &ElargeP, "ElargeP[2]/F");
	_outputTree->Branch("TrkSumEn", &TrkSumEn, "TrkSumEn/F");
	_outputTree->Branch("EequP", &EequP, "EequP[2]/F");
	_outputTree->Branch("EsmallP", &EsmallP, "EsmallP[2]/F");

	_outputTree->Branch("EcalTotalE", &_EcalTotalE, "EcalTotalE/F");
	_outputTree->Branch("HcalTotalE", &_HcalTotalE, "HcalTotalE/F");
	_outputTree->Branch("EcalCluE", &_EcalCluE, "EcalCluE/F");
        _outputTree->Branch("HcalCluE", &_HcalCluE, "HcalCluE/F");

	_outputTree->Branch("EcalCluE_p", &_EcalCluE_p, "EcalCluE_p/F");
        _outputTree->Branch("HcalCluE_p", &_HcalCluE_p, "HcalCluE_p/F");

	_outputPFO = new TTree("PFO","PFO");
	_outputPFO->Branch("EventNr", &_eventNr, "EventNr/I");
        _outputPFO->Branch("Num", &_Num, "Num/I");
	_outputPFO->Branch("Type", &Type, "Type/I");	// 0 == Arbor & 1 == Pandora
	_outputPFO->Branch("Charge", &Charge, "Charge/I");
	_outputPFO->Branch("Energy", &Energy, "Energy/F");
	_outputPFO->Branch("P", P, "P[3]/F");
	_outputPFO->Branch("CluEn", &CluEn, "CluEn/F");
	_outputPFO->Branch("CluEnCom", CluEnCom, "CluEnCom[2]/F");

	_outputPFO->Branch("TrackHit", &TrackHit, "TrackHit/I");
        _outputPFO->Branch("StartPos", StartPos, "StartPos[3]/F");
        _outputPFO->Branch("EndPos", EndPos, "EndPos[3]/F");

	_Num = 0;
}

void TotalInvMass::processEvent( LCEvent * evtP ) 
{		

	if (evtP) 								
	{		
		TLorentzVector ArborTotalP(0, 0, 0, 0);
		TLorentzVector ArborChP(0, 0, 0, 0);
		TLorentzVector ArborPhP(0, 0, 0, 0);
		TLorentzVector ArborNeP(0, 0, 0, 0);
		TLorentzVector ArborFrP(0, 0, 0, 0);
		TLorentzVector ArborUdP(0, 0, 0, 0);
		TLorentzVector ArborFrPh(0, 0, 0, 0);
		TLorentzVector ArborFrNe(0, 0, 0, 0);
		TLorentzVector ArborKPF(0, 0, 0, 0);
		TLorentzVector PandoraTotalP(0, 0, 0, 0);
		TLorentzVector ArborLCAL(0, 0, 0, 0);
		TLorentzVector ArborISR(0, 0, 0, 0);
		TLorentzVector PandoraISR(0, 0, 0, 0);

		_eventNr = evtP->getEventNumber();

		for(int i0 = 0; i0 < 4; i0++)
		{
			TotalP_a[i0] = 0;
			TotalP_p[i0] = 0;
			ChP[i0] = 0;
			PhP[i0] = 0;
			NeP[i0] = 0;
			FrP[i0] = 0;
			UdP[i0] = 0;
			FrPh[i0] = 0;
			FrNe[i0] = 0;
			KPF[i0] = 0;
			if(i0 < 2)
			{
				NeCaloE_a[i0] = 0;
				NeCaloE_p[i0] = 0;
				CluEnCom[i0] = 0;
				ElargeP[i0] = 0;
				EequP[i0] = 0;
				EsmallP[i0] = 0;
			}
			if(i0 < 3)
			{
				P[i0] = 0;
			}
		}

		nCHPFO_a = 0; 
		nCHPFO_p = 0;
		nNEPFO_a = 0; 
		nNEPFO_p = 0;
		Type = -100;
		Charge = -100; 
		Energy = 0;
		CluEn = 0; 
		TrkSumEn = 0; 

		_EcalTotalE = 0; _HcalTotalE = 0; _EcalCluE = 0; _HcalCluE = 0; 
		_EcalCluE_p = 0; _HcalCluE_p = 0;

		_HDPID = 0; 
		_OriQuarkID = 0; 
		_OQDir = -10;
		_HDir = -10;

		_J1CosTheta = -2;
		_J2CosTheta = -2;

		std::vector<std::string> EcalHitColl;
		std::vector<std::string> HcalHitColl;
		EcalHitColl.push_back("ECALBarrel");
		EcalHitColl.push_back("ECALEndcap");
		EcalHitColl.push_back("ECALOther");
		EcalHitColl.push_back("LCAL");
		EcalHitColl.push_back("LHCAL");

		HcalHitColl.push_back("HCALBarrel");
                HcalHitColl.push_back("HCALEndcap");
                HcalHitColl.push_back("HCALOther");
	
		try{

			for(int t = 0; t< int(EcalHitColl.size()); t++)
			{
				LCCollection * aecalcoll = evtP->getCollection(EcalHitColl[t].c_str());
				for(int s = 0; s < aecalcoll->getNumberOfElements();s++)
				{
					CalorimeterHit * a_hit = dynamic_cast<EVENT::CalorimeterHit*>(aecalcoll->getElementAt(s));
					_EcalTotalE += a_hit->getEnergy();
				}
			}

			for(int t2 = 0; t2< int(HcalHitColl.size()); t2++)
                        {
                                LCCollection * ahcalcoll = evtP->getCollection(HcalHitColl[t2].c_str());
                                for(int s = 0; s < ahcalcoll->getNumberOfElements();s++)
                                {
                                        CalorimeterHit * a_hit = dynamic_cast<EVENT::CalorimeterHit*>(ahcalcoll->getElementAt(s));
                                        _HcalTotalE += a_hit->getEnergy();
                                }
                        }


		}catch(lcio::DataNotAvailableException err) { }




		try{
			LCCollection * MCP = evtP->getCollection("MCParticlesSkimmed");
			for(int s0 = 0; s0 < MCP->getNumberOfElements(); s0++)
			{
				MCParticle *a1_MCP = dynamic_cast<EVENT::MCParticle *>(MCP->getElementAt(s0));
				int tmpPID = a1_MCP->getPDG();
				int NParent = a1_MCP->getParents().size();
				int NDaughter = a1_MCP->getDaughters().size();

				if(tmpPID == 25 && NDaughter > 1 && NParent !=0 ) //Higgs
				{
					_HDPID = abs((a1_MCP->getDaughters()[0])->getPDG());
					_HDir = a1_MCP->getMomentum()[2]/a1_MCP->getEnergy();

					if(NDaughter == 2)
					{
						MCParticle *D1 = a1_MCP->getDaughters()[0];
						_J1CosTheta = D1->getMomentum()[2]/D1->getEnergy();
						MCParticle *D2 = a1_MCP->getDaughters()[1];
						_J2CosTheta = D2->getMomentum()[2]/D2->getEnergy();
					}
				}

				if(abs(tmpPID)<7 && NParent == 0)
				{
					_OriQuarkID = abs(tmpPID); 
					_OQDir = a1_MCP->getMomentum()[2]/a1_MCP->getEnergy();
				}	
			}
		}catch(lcio::DataNotAvailableException err) { }

		try{
			LCCollection* col_RecoNeP = evtP->getCollection( "ArborPFOs" );

			for(int i0 = 0; i0 < col_RecoNeP->getNumberOfElements(); i0++)
			{
				ReconstructedParticle *a_RecoP = dynamic_cast<EVENT::ReconstructedParticle *>(col_RecoNeP->getElementAt(i0));	
				TLorentzVector currP( a_RecoP->getMomentum()[0], a_RecoP->getMomentum()[1], a_RecoP->getMomentum()[2], a_RecoP->getEnergy());
				ArborTotalP += currP;	
				TVector3 currMom = a_RecoP->getMomentum();

				if(a_RecoP->getCharge() != 0)
				{
					ArborChP += currP;
				}
				else if(a_RecoP->getType() == 310)
				{
					ArborKPF += currP; 
				}
				else if(a_RecoP->getType() == 22)
				{
					if(a_RecoP->getEnergy() < 3.0)
						ArborFrPh += currP;
					else
						ArborPhP += currP;
				}
				else if(a_RecoP->getType() == 2112)
				{
					if(a_RecoP->getEnergy() < 3.0)
						ArborFrNe += currP;
					else
						ArborNeP += currP;
				}
				else if(a_RecoP->getEnergy() < 3.0)
				{
					ArborFrP += currP;
					cout<<"Undef "<<a_RecoP->getType()<<endl;  
				}
				else
				{
					ArborUdP += currP;
					cout<<"Undef "<<a_RecoP->getType() << "En "<<a_RecoP->getEnergy()<<endl;
				}

				if(a_RecoP->getClusters().size() > 0)
				{
					Cluster * a_clu = a_RecoP->getClusters()[0];
					TVector3 CluPos = a_clu->getPosition();
					if( CluPos.Perp() < 300 && abs(CluPos.Z()) < 1300 && a_RecoP->getCharge() == 0 )	// 1150-1300
						ArborLCAL += currP; 

					float MinAngleToCH = 1.0E6;
					float MinAngleToNE = 1.0E6;  

					if(a_RecoP->getEnergy() > 5 && a_RecoP->getCharge() == 0)
					{
						for(int i1 = 0; i1 < col_RecoNeP->getNumberOfElements(); i1++)
						{
							if(i1 != i0)
							{
								ReconstructedParticle *b_RecoP = dynamic_cast<EVENT::ReconstructedParticle *>(col_RecoNeP->getElementAt(i1));
								TVector3 tmpMom = b_RecoP->getMomentum();
								float tmpAngle = currMom.Angle(tmpMom);
								if( b_RecoP->getEnergy() > 3.0 )
								{
									if(b_RecoP->getCharge() != 0)
									{
										if(tmpAngle < MinAngleToCH)
										{
											MinAngleToCH = tmpAngle;
										}
									}
									else 
									{	
										if(tmpAngle < MinAngleToNE)
										{
											MinAngleToNE = tmpAngle;
										}
									}
								}
							}
						}

						if( MinAngleToNE > 0.5 || MinAngleToCH > 0.5 )
						{
							ArborISR += currP; 
						}
					}
				}

				TrackHit = 0; 
				for(int k = 0; k < 3; k++)	
				{
					StartPos[k] = 0;
					EndPos[k] = 0;
				}		

				Charge = int(a_RecoP->getCharge());
				CluEn = 0; 
				CluEnCom[0] = 0;
				CluEnCom[1] = 0;

				if(Charge)
				{
					Type = 1;
					nCHPFO_a ++;
				}
				else
				{
					Type = 0;
					nNEPFO_a ++;
				}

				Energy = a_RecoP->getEnergy();
				P[0] = a_RecoP->getMomentum()[0];
				P[1] = a_RecoP->getMomentum()[1];
				P[2] = a_RecoP->getMomentum()[2];

				if(a_RecoP->getClusters().size() > 0)
				{
					Cluster * currClu = a_RecoP->getClusters()[0];
					CluEn = currClu->getEnergy();

					if(!Charge)
					{
						CluEnCom[0] = currClu->getSubdetectorEnergies()[0];
						CluEnCom[1] = currClu->getSubdetectorEnergies()[1];	
						NeCaloE_a[0] += currClu->getSubdetectorEnergies()[0];
						NeCaloE_a[1] += currClu->getSubdetectorEnergies()[1];
					}

					_EcalCluE += currClu->getSubdetectorEnergies()[0];
					_HcalCluE += currClu->getSubdetectorEnergies()[1];

					if(Charge)
					{
						TrkSumEn += Energy; 

						Track * a_Trk = a_RecoP->getTracks()[0];
						TrackHit = a_Trk->getTrackerHits().size();
						StartPos[0] = (a_Trk->getTrackerHits()[0])->getPosition()[0];
						StartPos[1] = (a_Trk->getTrackerHits()[0])->getPosition()[1];
						StartPos[2] = (a_Trk->getTrackerHits()[0])->getPosition()[2];

						EndPos[0] = (a_Trk->getTrackerHits()[TrackHit - 1])->getPosition()[0];
						EndPos[1] = (a_Trk->getTrackerHits()[TrackHit - 1])->getPosition()[1];
						EndPos[2] = (a_Trk->getTrackerHits()[TrackHit - 1])->getPosition()[2];	

						if( Energy > CluEn + sqrt(Energy))
						{
							ElargeP[0] += Energy; 
							ElargeP[1] += CluEn; 
						}
						else if( fabs(Energy - CluEn) < sqrt(Energy) )
						{
							EequP[0] += Energy;
							EequP[1] += CluEn; 
						}
						else
						{
							EsmallP[0] += Energy;
							EsmallP[1] += CluEn; 
						}
					}
				}

				_outputPFO->Fill();

			}
		}catch (lcio::DataNotAvailableException err) { }


		try{
			LCCollection* col_RecoPandora = evtP->getCollection( "PandoraPFOs" );
			for(int i2 = 0; i2 < col_RecoPandora->getNumberOfElements(); i2++)
			{         
				ReconstructedParticle *a_RecoP = dynamic_cast<EVENT::ReconstructedParticle *>(col_RecoPandora->getElementAt(i2));
				TLorentzVector currP( a_RecoP->getMomentum()[0], a_RecoP->getMomentum()[1], a_RecoP->getMomentum()[2], a_RecoP->getEnergy());
				PandoraTotalP += currP;

				if(a_RecoP->getCharge())
				{
					nCHPFO_p++;
				}
				else
					nNEPFO_p++;

				TVector3 currMom = a_RecoP->getMomentum();

				if(a_RecoP->getClusters().size() > 0)
				{

					float MinAngleToCH = 1.0E6;
					float MinAngleToNE = 1.0E6;

					 Cluster * currClu = a_RecoP->getClusters()[0];
                                        CluEn = currClu->getEnergy();

 					_EcalCluE_p += CluEn;	
	                      //          _EcalCluE_p += currClu->getSubdetectorEnergies()[0];
                                        _HcalCluE_p += currClu->getSubdetectorEnergies()[1];

					if(a_RecoP->getEnergy() > 5 && a_RecoP->getCharge() == 0)
					{
					for(int i3 = 0; i3 < col_RecoPandora->getNumberOfElements(); i3++)
					{
						if(i3 != i2)
						{
							ReconstructedParticle *b_RecoP = dynamic_cast<EVENT::ReconstructedParticle *>(col_RecoPandora->getElementAt(i3));
							TVector3 tmpMom = b_RecoP->getMomentum();
							float tmpAngle = currMom.Angle(tmpMom);
							if( b_RecoP->getEnergy() > 3.0 )
							{
								if(b_RecoP->getCharge() != 0)
								{
									if(tmpAngle < MinAngleToCH)
									{
										MinAngleToCH = tmpAngle;
									}
								}
								else
								{
									if(tmpAngle < MinAngleToNE)
									{
										MinAngleToNE = tmpAngle;
									}
								}
							}
						}
					}

					if( MinAngleToNE > 0.5 || MinAngleToCH > 0.5 )
					{
						PandoraISR += currP;
					}
					}
				}



			}
		}catch (lcio::DataNotAvailableException err) { }

		_Mass_a = 0; 
		_Mass_p = 0; 
		_Mass_a_Pisr = 0;
		_Mass_a_Plcal = 0; 
		_Mass_p_Pisr = 0; 

		_Mass_a = ArborTotalP.M();
		_Mass_p = PandoraTotalP.M();

		_Mass_a_Pisr = (ArborTotalP - ArborISR).M();
		_Mass_a_Plcal = (ArborTotalP - ArborLCAL).M();
		_Mass_p_Pisr = (PandoraTotalP - PandoraISR).M();

		TotalP_a[0] = ArborTotalP.X();
		TotalP_a[1] = ArborTotalP.Y();
		TotalP_a[2] = ArborTotalP.Z();
		TotalP_a[3] = ArborTotalP.T();

		TotalP_p[0] = PandoraTotalP.X();
		TotalP_p[1] = PandoraTotalP.Y();
		TotalP_p[2] = PandoraTotalP.Z();
		TotalP_p[3] = PandoraTotalP.T();

		ChP[0] = ArborChP.X();
		ChP[1] = ArborChP.Y();
		ChP[2] = ArborChP.Z();
		ChP[3] = ArborChP.T();

		PhP[0] = ArborPhP.X();
		PhP[1] = ArborPhP.Y();
		PhP[2] = ArborPhP.Z();
		PhP[3] = ArborPhP.T();	

		NeP[0] = ArborNeP.X();
		NeP[1] = ArborNeP.Y();
		NeP[2] = ArborNeP.Z();
		NeP[3] = ArborNeP.T();

		FrP[0] = ArborFrP.X();
		FrP[1] = ArborFrP.Y();
		FrP[2] = ArborFrP.Z();
		FrP[3] = ArborFrP.T();

		UdP[0] = ArborUdP.X();
		UdP[1] = ArborUdP.Y();
		UdP[2] = ArborUdP.Z();
		UdP[3] = ArborUdP.T();

		FrPh[0] = ArborFrPh.X();
		FrPh[1] = ArborFrPh.Y();
		FrPh[2] = ArborFrPh.Z();
		FrPh[3] = ArborFrPh.T();

		FrNe[0] = ArborFrNe.X();
		FrNe[1] = ArborFrNe.Y();
		FrNe[2] = ArborFrNe.Z();
		FrNe[3] = ArborFrNe.T();

		KPF[0] = ArborKPF.X();
		KPF[1] = ArborKPF.Y();
		KPF[2] = ArborKPF.Z();
		KPF[3] = ArborKPF.T();

		cout<<_Mass_a<<" : "<<_Mass_p<<endl;		

		_outputTree->Fill();
		_Num++;

	}  	  

}	

void TotalInvMass::end()
{

	if (_outputTree) {

		TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		delete tree_file;
	}

}



