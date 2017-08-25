#include <PrintReco_MCP.hh>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <values.h>
#include <string>
#include <iostream>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <stdexcept>
#include <TFile.h> 
#include <TTree.h>
#include <Rtypes.h> 
#include <sstream>		
#include <cmath>
#include <TVector3.h>

PrintReco_MCP a_PrintReco_MCP_instance;

PrintReco_MCP::PrintReco_MCP()
	: Processor("PrintReco_MCP"),
	_output(0)
{
	_description = "Print MC Truth" ;



	_treeFileName="MCTruth.root";
	registerProcessorParameter( "TreeOutputFile" , 
			"The name of the file to which the ROOT tree will be written" ,
			_treeFileName ,
			_treeFileName);


	_colName="ArborPFOs";
	registerProcessorParameter( "ParticleFlowObject" ,
			"The name of the PFOs" ,
			_colName ,
			_colName);


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

}

void PrintReco_MCP::init() {

	printParameters();

	TFile *tree_file=new TFile(_treeFileName.c_str(),(_overwrite ? "RECREATE" : "UPDATE"));


	if (!tree_file->IsOpen()) {
		delete tree_file;
		tree_file=new TFile(_treeFileName.c_str(),"NEW");
	}

	_outputEvt = new TTree("Evt", "Evt");
        _outputEvt->SetAutoSave(32*1024*1024);  // autosave every 32MB
        _outputEvt->Branch("EventNr", &_eventNr, "EventNr/I");
	_outputEvt->Branch("NMCP",&_nMCP,"NMCP/I");
	_outputEvt->Branch("NPFOs",&_nPFOs,"NPFOs/I");
	_outputEvt->Branch("NChPFO", &_NChPFO, "NChPFO/I");
	_outputEvt->Branch("NNePFO", &_NNePFO, "NNePFO/I");

	_outputEvt->Branch("NSelMCP", &_nSelMCP, "NSelMCP/I");
	_outputEvt->Branch("TotalEn_MCP", &_TotalEnMCP, "TotalEnMCP/F");
	_outputEvt->Branch("P_MCP", &_P_MCP, "P_MCP[3]/F");
	_outputEvt->Branch("Mass_MCP", &_Mass_MCP, "Mass_MCP/F");

	_outputEvt->Branch("TotalEn_RecoP", &_TotalEnRecoP, "TotalEnRecoP/F");
        _outputEvt->Branch("P_RecoP", &_P_RecoP, "P_RecoP[3]/F");
        _outputEvt->Branch("Mass_RecoP", &_Mass_RecoP, "Mass_RecoP/F");


	_outputMCP = new TTree("MCP", "MCP");
	_outputMCP->SetAutoSave(32*1024*1024);  // autosave every 32MB
	_outputMCP->Branch("EventNr", &_eventNr, "EventNr/I");
	_outputMCP->Branch("NMCP",&_nMCP,"NMCP/I");
	_outputMCP->Branch("NSelMCP", &_nSelMCP, "NSelMCP/I");

	//_outputEvt

	_outputMCP->Branch("energy",&_energy,"energy/F");
	_outputMCP->Branch("mass",&_mass,"mass/F");
	_outputMCP->Branch("Charge",&_Charge,"Charge/F");
	_outputMCP->Branch("Px",&_Px,"Px/F");
	_outputMCP->Branch("Py",&_Py,"Py/F");
	_outputMCP->Branch("Pz",&_Pz,"Pz/F");
	_outputMCP->Branch("EndP_X",&_Posx,"EndP_X/F");
	_outputMCP->Branch("EndP_Y",&_Posy,"EndP_Y/F");
	_outputMCP->Branch("EndP_Z",&_Posz,"EndP_Z/F");
	_outputMCP->Branch("PID",&_PID,"PID/I");
	_outputMCP->Branch("Vex_X",&_Vex_x,"Vex_X/F");
	_outputMCP->Branch("Vex_Y",&_Vex_y,"Vex_Y/F");
	_outputMCP->Branch("Vex_Z",&_Vex_z,"Vex_Z/F");
	_outputMCP->Branch("GenStatus",&_GenStatus,"GenStatus/I");
	_outputMCP->Branch("SimStatus",&_SimStatus,"SimStatus/I");
	_outputMCP->Branch("NumofParent",&_parentnum,"NumofParent/I");
	_outputMCP->Branch("NumofDaughter",&_daughternum,"NumofDaughter/I");

	_outputPFO = new TTree("RecoP", "RecoP");
        _outputPFO->SetAutoSave(32*1024*1024);  // autosave every 32MB
        _outputPFO->Branch("EventNr", &_eventNr, "EventNr/I");

	_outputPFO->Branch("energy",&_energy,"energy/F");
	_outputPFO->Branch("CluEn", &_CluEn, "CluEn/F");
        _outputPFO->Branch("mass",&_mass,"mass/F");
        _outputPFO->Branch("Charge",&_Charge,"Charge/F");
        _outputPFO->Branch("Px",&_Px,"Px/F");
        _outputPFO->Branch("Py",&_Py,"Py/F");
        _outputPFO->Branch("Pz",&_Pz,"Pz/F");
	_outputPFO->Branch("Type",&_Type,"Type/I");
}

void PrintReco_MCP::processEvent( LCEvent * evtP ) 
{		

	if (evtP) 								
	{
		_Mass_RecoP = 0; 
		_Mass_MCP = 0; 
		
		try 	
		{    
			LCCollection* col_PFOs = evtP->getCollection( _colName ) ;

			_nPFOs=col_PFOs->getNumberOfElements();
			_NChPFO = 0; 
			_NNePFO = 0; 
			_eventNr=evtP->getEventNumber();
			_TotalEnRecoP = 0;
			_P_RecoP[0] = 0;
			_P_RecoP[1] = 0;
			_P_RecoP[2] = 0;

			for (int PFO_i=0; PFO_i < _nPFOs; PFO_i++)
			{		
				EVENT::ReconstructedParticle *a_PFO=dynamic_cast<EVENT::ReconstructedParticle *>(col_PFOs->getElementAt(PFO_i));

				_energy=a_PFO->getEnergy();
				_CluEn = 0;				

				for(unsigned int s = 0; s < a_PFO->getClusters().size(); s++)
				{
					_CluEn += (a_PFO->getClusters()[s])->getEnergy();
				}

				_mass=a_PFO->getMass();
				_Charge=a_PFO->getCharge();
				_Px=a_PFO->getMomentum()[0];
				_Py=a_PFO->getMomentum()[1];
				_Pz=a_PFO->getMomentum()[2];
				_Type = a_PFO->getType();

				if( fabs(a_PFO->getCharge()) < 0.01 )
				{
					_NNePFO++;
				}
				else
				{
					_NChPFO++;
				}

				_TotalEnRecoP += _energy;
				_P_RecoP[0] += _Px; 
				_P_RecoP[1] += _Py;
				_P_RecoP[2] += _Pz;	

				_outputPFO->Fill();
			}

			_Mass_RecoP = sqrt(_TotalEnRecoP*_TotalEnRecoP - _P_RecoP[0]*_P_RecoP[0] - _P_RecoP[1]*_P_RecoP[1] - _P_RecoP[2]*_P_RecoP[2]);
				
		}catch (lcio::DataNotAvailableException err) { }

		try{
			LCCollection* col_MCP = evtP->getCollection( "MCParticle" ) ;

			_nMCP=col_MCP->getNumberOfElements();
			_eventNr=evtP->getEventNumber();	
			TVector3 VtxPos, EndPPos; 
			_nSelMCP = 0;
			_TotalEnMCP = 0;
                        _P_MCP[0] = 0;
                        _P_MCP[1] = 0;
                        _P_MCP[2] = 0;

			for(int j = 0; j < _nMCP; j++)
			{
				MCParticle* a_MCP = dynamic_cast<MCParticle*>(col_MCP->getElementAt(j));

				_daughternum=a_MCP->getDaughters().size();
				_parentnum=a_MCP->getParents().size();
				VtxPos = a_MCP->getVertex();
				EndPPos = a_MCP->getEndpoint();
				_PID=a_MCP->getPDG();

				// if( _daughternum == 0 && _parentnum > 0 )	
				if( VtxPos.Mag() < 100 && EndPPos.Mag() > 100 && fabs(_PID) != 12 &&  fabs(_PID) != 14 &&  fabs(_PID) != 16 )
				{
					_nSelMCP++;
					_energy=a_MCP->getEnergy();
					_mass=a_MCP->getMass();
					_Charge=a_MCP->getCharge();
					_Px=a_MCP->getMomentum()[0];
					_Py=a_MCP->getMomentum()[1];
					_Pz=a_MCP->getMomentum()[2];

					_TotalEnMCP += _energy;
					_P_MCP[0] += _Px;
					_P_MCP[1] += _Py;
					_P_MCP[2] += _Pz;


					_Posx=a_MCP->getEndpoint()[0];
					_Posy=a_MCP->getEndpoint()[1];
					_Posz=a_MCP->getEndpoint()[2];
					_Vex_x=a_MCP->getVertex()[0];
					_Vex_y=a_MCP->getVertex()[1];
					_Vex_z=a_MCP->getVertex()[2];

					_GenStatus=a_MCP->getGeneratorStatus();
					_SimStatus=a_MCP->getSimulatorStatus();
					_outputMCP->Fill();
				}
			}

			_Mass_MCP = sqrt(_TotalEnMCP*_TotalEnMCP - _P_MCP[0]*_P_MCP[0] - _P_MCP[1]*_P_MCP[1] - _P_MCP[2]*_P_MCP[2] );

		}catch (lcio::DataNotAvailableException err) { }

		_outputEvt->Fill();
	}  	  
}	

void PrintReco_MCP::end()
{

	if (_outputMCP) {

		TFile *tree_file = _outputMCP->GetCurrentFile(); //just in case we switched to a new file
		tree_file->Write();
		delete tree_file;

	}

}



