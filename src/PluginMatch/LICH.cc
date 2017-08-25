#include <LICH.hh>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <EVENT/Cluster.h>
#include <EVENT/LCRelation.h>
#include <EVENT/ParticleID.h>

#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/ParticleIDImpl.h>
#include "UTIL/CellIDDecoder.h"
#include "HelixClass.hh"

#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <Rtypes.h>
#include <sstream>
#include <TH1.h>
#include <TVector3.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/Config.h"

#include "gear/BField.h"
#include "gear/CalorimeterParameters.h"
#include <marlin/Global.h>





using namespace std;

const float TightGeoThreshold = 30;
const string ECALCellIDDecoder  = "M:3,S-1:3,I:9,J:9,K-1:6";


LICH aLICH ;
LICH::LICH()
	: Processor("LICH"),
	_output(0)
{
	_description = "Produce single particle samples for training or produce the charged PID" ;
	_FileName="OUTPUT";
	registerProcessorParameter( "TreeOutputFile" ,
			        "The name of the file to which the pion ROOT tree will be written" ,
				_FileName ,
				_FileName);

	_sampleEn="SampleEnergy";
	registerProcessorParameter( "TrainingEn" ,
			        "The generated energy point choson for training" ,
				_sampleEn ,
				_sampleEn);

	_inputPFO="ArborPFOs";
	registerProcessorParameter("InputPFO",
			"Name of Inpout Reconstructed Particle Collection",
			_inputPFO,
			_inputPFO);

	_outputPFO="typedPFOs";
	registerProcessorParameter("OutputPFO",
			"Name of Output Reconstructed Particle Collection",
			_outputPFO,
			_outputPFO);

	_inputMCP="MCParticle";
	registerProcessorParameter("InputMCParticle",
			"Name of Inpout MC Particle Collection",
			_inputMCP,
			_inputMCP);

	_Training=0;
	registerProcessorParameter("TrainingFlag",
			"Produce sample for training when 0, otherwise proceed the PID",
			_Training,
			_Training);

	int en[]={1,2,3,5,7,10,20,30,40,50,70};
	float ang[]={0.0,0.3,0.55,0.75,1.0};
	NEn=sizeof(en)/sizeof(int);
	NPos=sizeof(ang)/sizeof(int);
	for(int i=0;i<NEn;i++){
		inputEn.push_back(en[i]);
	}
	for(int i=0;i<NPos;i++){
		inputPos.push_back(ang[i]);
	}

	registerProcessorParameter("InputEnergyPoints",
			"The interval of energy for PID",
			inputEn,
			inputEn);


	registerProcessorParameter("InputPositions",
			"The interval of angle for PID",
			inputPos,
			inputPos);


	inputDetectorModules.push_back("barrel1");
	inputDetectorModules.push_back("barrel2");
	inputDetectorModules.push_back("overlap");
	inputDetectorModules.push_back("endcap");
	registerProcessorParameter("InputDetectorModules" ,
			"DetectorModules" ,
			inputDetectorModules,
			inputDetectorModules);

	mvacut_e=0.5;
	mvacut_mu=0.5;
	mvacut_pi=0.5;

	registerProcessorParameter("mvacut_e",
			"mva value cut for electron",
			mvacut_e,
			mvacut_e);
	registerProcessorParameter("mvacut_pi",
			"mva value cut for electron",
			mvacut_pi,
			mvacut_pi);
	registerProcessorParameter("mvacut_mu",
			"mva value cut for electron",
			mvacut_mu,
			mvacut_mu);

	registerProcessorParameter("weightDir",
			"the weight file directory",
			weightDir,
			weightDir);

}

void LICH::gearPara(){

	const gear::CalorimeterParameters &ecalBarrelParameters = marlin::Global::GEAR->getEcalBarrelParameters();
	const gear::CalorimeterParameters &ecalEndCapParameters = marlin::Global::GEAR->getEcalEndcapParameters();
	const gear::CalorimeterParameters &hcalBarrelParameters = marlin::Global::GEAR->getHcalBarrelParameters();
	const gear::CalorimeterParameters &hcalEndCapParameters = marlin::Global::GEAR->getHcalEndcapParameters();

	 rminBarrelEcal=(ecalBarrelParameters.getExtent()[0]);
	 rmaxBarrelEcal=(ecalBarrelParameters.getExtent()[1]);
	 zminBarrelEcal=(ecalBarrelParameters.getExtent()[2]);
	 zmaxBarrelEcal=(ecalBarrelParameters.getExtent()[3]);

	 rminEndCapEcal=(ecalEndCapParameters.getExtent()[0]);
	 rmaxEndCapEcal=(ecalEndCapParameters.getExtent()[1]);
	 zminEndCapEcal=(ecalEndCapParameters.getExtent()[2]);
	 zmaxEndCapEcal=(ecalEndCapParameters.getExtent()[3]);

	 rminBarrelHcal=(hcalBarrelParameters.getExtent()[0]);
	 rmaxBarrelHcal=(hcalBarrelParameters.getExtent()[1]);
	 zminBarrelHcal=(hcalBarrelParameters.getExtent()[2]);
	 zmaxBarrelHcal=(hcalBarrelParameters.getExtent()[3]);

	 rminEndCapHcal=(hcalEndCapParameters.getExtent()[0]);
	 rmaxEndCapHcal=(hcalEndCapParameters.getExtent()[1]);
	 zminEndCapHcal=(hcalEndCapParameters.getExtent()[2]);
	 zmaxEndCapHcal=(hcalEndCapParameters.getExtent()[3]);

	TVector3 theTPCCenter(0.,0.,0.);
	TVector3 bfield = marlin::Global::GEAR->getBField().at(theTPCCenter);
	 BField=bfield.z();
}

int LICH::SubDeFlag(TVector3 inputPos){


	int FlagD(-1);
	if( fabs(inputPos[2]) > zmaxEndCapHcal || fabs(inputPos.Perp()) > rmaxBarrelHcal )
	{
		FlagD = 2;
	}
	else if( fabs(inputPos[2]) > zminEndCapHcal || fabs(inputPos.Perp()) > rminBarrelHcal)
	{
		FlagD = 1;          // Position outsider than DHCAL Region
	}
	else if( fabs(inputPos[2]) > zmaxEndCapEcal || fabs(inputPos.Perp()) > rmaxBarrelEcal)
	{
		FlagD = 10;
	}
	else if( fabs(inputPos[2]) > zminEndCapEcal || fabs(inputPos.Perp()) > rminBarrelEcal)
	{
		FlagD = 0;
	}

	else
	{
		FlagD = 11;         // Position inside Calo... Problematic for Seeds... But could be PreShower hits.
	}

	return FlagD;


}
int LICH::ActiveLayers(std::vector<CalorimeterHit*> clu, const std::string& encoder_str){

	std::vector<int> hitlayers;
	CellIDDecoder<CalorimeterHit> idDecoder(encoder_str);
	int NHits = clu.size();
	int tmpK = 0;
	int tmpS = 0;
	int tmpID = 0;

	for(int i = 0; i < NHits; i++)
	{
		CalorimeterHit *hit = dynamic_cast<CalorimeterHit*>( clu[i]);
		tmpK = idDecoder(hit)["K-1"]+1 ;
		tmpS = idDecoder(hit)["S-1"]+1 ;
		tmpID = tmpS * 50 + tmpK;
		if( std::find(hitlayers.begin(), hitlayers.end(), tmpID) == hitlayers.end() )
		{
			hitlayers.push_back(tmpID);
		}
	}
	return hitlayers.size();
}
int LICH::NHScaleV2( const std::string& encoder_str, std::vector<CalorimeterHit*> clu0, int RatioX, int RatioY, int RatioZ )
{

	int ReScaledNH = 0;
	int NumHit = clu0.size();
	int tmpI = 0;
	int tmpJ = 0;
	int tmpK = 0;
	float tmpEn = 0;
	int NewCellID0 = 0;

	CellIDDecoder<CalorimeterHit> idDecoder(encoder_str);      //Input Hits here refer to AllCleanHits collection

	std::map <double, float> testIDtoEnergy;

	for(int i = 0; i < NumHit; i++)
	{
		CalorimeterHit *hit = dynamic_cast<CalorimeterHit*>( clu0[i]);

		tmpI = idDecoder(hit)["I"]/RatioX;
		tmpJ = idDecoder(hit)["J"]/RatioY;
		tmpK = (idDecoder(hit)["K-1"]+1)/RatioZ;
		tmpEn = hit->getEnergy();

		NewCellID0 = (tmpK<<24) + (tmpJ<<12) + tmpI;

		if(testIDtoEnergy.find(NewCellID0) == testIDtoEnergy.end() )
		{
			testIDtoEnergy[NewCellID0] = tmpEn;
		}
		else
		{
			testIDtoEnergy[NewCellID0] += tmpEn;
		}
	}

	ReScaledNH = testIDtoEnergy.size();

	return ReScaledNH;
}

float LICH::FDV2(std::vector<CalorimeterHit*> clu, const std::string& encoder_str)
{
	float FractalDim = 0;
	int NReSizeHit[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	int Scale[10] = {20, 30, 40, 50, 60, 70, 80, 90, 100, 200};
	int OriNHit = clu.size();

	for(int j = 0; j < 10; j++)
	{
		NReSizeHit[j] = NHScaleV2(encoder_str, clu, Scale[j], Scale[j], 1);
		FractalDim += 0.1 * TMath::Log(float(OriNHit)/NReSizeHit[j])/TMath::Log(float(Scale[j]));
	}
	if(clu.size() == 0)
		FractalDim = -1;
	return FractalDim;
}
float LICH::DisSeedSurface( TVector3 SeedPos )
{

	float DisSS = 0;

	if( fabs(SeedPos[2]) > zmaxBarrelEcal )         //EcalEndcap hit start from 2350 + 100 = 2450
	{

		if( SeedPos.Perp() > rminBarrelEcal )
		{
			if( fabs(SeedPos[2])/SeedPos.Perp() > (zminEndCapEcal+3)/(rminBarrelEcal+ 100) )
			{
				DisSS = ( fabs(SeedPos[2]) - zminEndCapEcal - 3 ) * SeedPos.Mag()/fabs(SeedPos[2]);
			}
			else
			{
				DisSS = (SeedPos.Perp() - rminBarrelEcal - 100 )*SeedPos.Mag()/SeedPos.Perp();
			}
		}
		else
		{
			DisSS = fabs(SeedPos[2]) - zminEndCapEcal - 3;
		}

	}
	else if( SeedPos.Perp() > rminBarrelEcal + 400 )
	{
		DisSS = SeedPos.Perp() - rminBarrelEcal - 100;
	}
	else if( (SeedPos.Phi() > 0 && int(SeedPos.Phi() * 4/TMath::Pi() + 0.5) % 2 == 0 ) || (SeedPos.Phi() < 0 && int(SeedPos.Phi() * 4/TMath::Pi() + 8.5) % 2 == 0 ))
	{
		DisSS = min( fabs(fabs(SeedPos[0]) - rminBarrelEcal), fabs(fabs(SeedPos[1]) - rminBarrelEcal ) );
	}
	else
	{
		DisSS = min( fabs(fabs(SeedPos[0] + SeedPos[1])/TMath::Pi() -rminBarrelEcal), fabs(fabs(SeedPos[0] - SeedPos[1])/TMath::Pi() - rminBarrelEcal) );
	}

	return DisSS;
}
TVector3 LICH::ClusterCoG(Cluster * inputCluster)
{
	TVector3 CenterOfGravity; 

	int inputClusterSize = inputCluster->getCalorimeterHits().size();

	TVector3 tmphitPos; 
	float tmphitEnergy;
	float sumhitEnergy = 0; 

	for(int i = 0; i < inputClusterSize; i++)
	{
		CalorimeterHit * tmpHit = inputCluster->getCalorimeterHits()[i];
		tmphitPos = tmpHit->getPosition();
		tmphitEnergy = tmpHit->getEnergy();

		CenterOfGravity += tmphitPos*tmphitEnergy;
		sumhitEnergy += tmphitEnergy; 
	}

	CenterOfGravity = 1.0/sumhitEnergy * CenterOfGravity; 

	return CenterOfGravity; 
}
void LICH::init() {

	cout<<" _    _   ___  _   _ "<<endl;
	cout<<"| |  | | / _ \\| | | |"<<endl;
	cout<<"| |  | || / \\/| |_| |"<<endl;
	cout<<"| |  | || |   |     |"<<endl;
	cout<<"| |_ | || |   |  _  |"<<endl;
	cout<<"|   || || \\_/\\| | | | "<<endl;
	cout<<"|___||_| \\___/|_| |_|"<<endl;
	cout<<""<<endl;
	printParameters();


	Cluflag.setBit(LCIO::CHBIT_LONG);

	if(_Training)_treeFileName=_FileName+"_"+_sampleEn+".root";
	if(_Training==0)_treeFileName=_FileName+".root";

	tree_file=new TFile(_treeFileName.c_str(), "RECREATE");

	piTree = new TTree("Pion","Pion");

	piTree->Branch("EventNr", &eventNr, "EventNr/I");
	piTree->Branch("MCPDG", &MCPDG, "MCPDG/I");
	piTree->Branch("ClusterID", &ClusterID, "ClusterID/I");
	piTree->Branch("NPFO", &NPFO, "NPFO/I");
	piTree->Branch("EcalNHit",&EcalNHit,"EcalNHit/I");
	piTree->Branch("HcalNHit",&HcalNHit,"HcalNHit/I");
	piTree->Branch("CluNHit",&CluNHit,"CluNHit/I");
	piTree->Branch("NLEcal",&NLEcal,"NLEcal/I");
	piTree->Branch("NLHcal",&NLHcal,"NLHcal/I");
	piTree->Branch("maxDisHtoL",&maxDisHtoL,"maxDisHtoL/F");
	piTree->Branch("minDisHtoL",&minDisHtoL,"minDisHtoL/F");
	piTree->Branch("avDisHtoL",&avDisHtoL,"avDisHtoL/F");
	piTree->Branch("avEnDisHtoL",&avEnDisHtoL,"avEnDisHtoL/F");
	piTree->Branch("EcalEn",&EcalEn,"EcalEn/F");
	piTree->Branch("HcalEn",&HcalEn,"HcalEn/F");
	piTree->Branch("EClu",&EClu,"EClu/F");
	piTree->Branch("graDepth",&graDepth,"graDepth/F");
	piTree->Branch("cluDepth",&cluDepth,"cluDepth/F");
	piTree->Branch("graAbsDepth",&graAbsDepth,"graAbsDepth/F");
	piTree->Branch("maxDepth",&maxDepth,"maxDepth/F");
	piTree->Branch("minDepth",&minDepth,"minDepth/F");
	piTree->Branch("MaxDisHel",&MaxDisHel,"MaxDisHel/F");
	piTree->Branch("MinDisHel",&MinDisHel,"MinDisHel/F");
	piTree->Branch("FD_all",&FD_all,"FD_all/F");
	piTree->Branch("FD_ECAL",&FD_ECAL,"FD_ECAL/F");
	piTree->Branch("FD_HCAL",&FD_HCAL,"FD_HCAL/F");
	piTree->Branch("crdis",&crdis,"crdis/F");
	piTree->Branch("EEClu_L10",&EEClu_L10,"EEClu_L10/F");
	piTree->Branch("EEClu_R",&EEClu_R,"EEClu_R/F");
	piTree->Branch("EEClu_r",&EEClu_r,"EEClu_r/F");
	piTree->Branch("EEClu_p",&EEClu_p,"EEClu_p/F");
	piTree->Branch("rms_Ecal",&rms_Ecal,"rms_Ecal/F");
	piTree->Branch("rms_Hcal",&rms_Hcal,"rms_Hcal/F");
	piTree->Branch("rms_Ecal2",&rms_Ecal2,"rms_Ecal2/F");
	piTree->Branch("rms_Hcal2",&rms_Hcal2,"rms_Hcal2/F");
	piTree->Branch("av_NHE",&av_NHE,"av_NHE/F");
	piTree->Branch("av_NHH",&av_NHH,"av_NHH/F");
	piTree->Branch("AL_Ecal",&AL_Ecal,"AL_Ecal/I");
	piTree->Branch("AL_Hcal",&AL_Hcal,"AL_Hcal/I");
	piTree->Branch("FD_ECALF10",&FD_ECALF10,"FD_ECALF10/F");
	piTree->Branch("FD_ECALL20",&FD_ECALL20,"FD_ECALL20/F");
	piTree->Branch("NH_ECALF10",&NH_ECALF10,"NH_ECALF10/I");
	piTree->Branch("NH_ECALL20",&NH_ECALL20,"NH_ECALL20/I");
	piTree->Branch("dEdx",&dEdx,"dEdx/F");
	piTree->Branch("cosTheta",&cosTheta,"cosTheta/F");
	piTree->Branch("Phi",&Phi,"Phi/F");
	
	muTree = new TTree("Muon","Muon");

	muTree->Branch("EventNr", &eventNr, "EventNr/I");
	muTree->Branch("MCPDG", &MCPDG, "MCPDG/I");
	muTree->Branch("ClusterID", &ClusterID, "ClusterID/I");
	muTree->Branch("NPFO", &NPFO, "NPFO/I");
	muTree->Branch("EcalNHit",&EcalNHit,"EcalNHit/I");
	muTree->Branch("HcalNHit",&HcalNHit,"HcalNHit/I");
	muTree->Branch("CluNHit",&CluNHit,"CluNHit/I");
	muTree->Branch("NLEcal",&NLEcal,"NLEcal/I");
	muTree->Branch("NLHcal",&NLHcal,"NLHcal/I");
	muTree->Branch("maxDisHtoL",&maxDisHtoL,"maxDisHtoL/F");
	muTree->Branch("minDisHtoL",&minDisHtoL,"minDisHtoL/F");
	muTree->Branch("avDisHtoL",&avDisHtoL,"avDisHtoL/F");
	muTree->Branch("avEnDisHtoL",&avEnDisHtoL,"avEnDisHtoL/F");
	muTree->Branch("EcalEn",&EcalEn,"EcalEn/F");
	muTree->Branch("HcalEn",&HcalEn,"HcalEn/F");
	muTree->Branch("EClu",&EClu,"EClu/F");
	muTree->Branch("graDepth",&graDepth,"graDepth/F");
	muTree->Branch("cluDepth",&cluDepth,"cluDepth/F");
	muTree->Branch("graAbsDepth",&graAbsDepth,"graAbsDepth/F");
	muTree->Branch("maxDepth",&maxDepth,"maxDepth/F");
	muTree->Branch("minDepth",&minDepth,"minDepth/F");
	muTree->Branch("MaxDisHel",&MaxDisHel,"MaxDisHel/F");
	muTree->Branch("MinDisHel",&MinDisHel,"MinDisHel/F");
	muTree->Branch("FD_all",&FD_all,"FD_all/F");
	muTree->Branch("FD_ECAL",&FD_ECAL,"FD_ECAL/F");
	muTree->Branch("FD_HCAL",&FD_HCAL,"FD_HCAL/F");
	muTree->Branch("crdis",&crdis,"crdis/F");
	muTree->Branch("EEClu_L10",&EEClu_L10,"EEClu_L10/F");
	muTree->Branch("EEClu_R",&EEClu_R,"EEClu_R/F");
	muTree->Branch("EEClu_r",&EEClu_r,"EEClu_r/F");
	muTree->Branch("EEClu_p",&EEClu_p,"EEClu_p/F");
	muTree->Branch("rms_Ecal",&rms_Ecal,"rms_Ecal/F");
	muTree->Branch("rms_Hcal",&rms_Hcal,"rms_Hcal/F");
	muTree->Branch("rms_Ecal2",&rms_Ecal2,"rms_Ecal2/F");
	muTree->Branch("rms_Hcal2",&rms_Hcal2,"rms_Hcal2/F");
	muTree->Branch("av_NHE",&av_NHE,"av_NHE/F");
	muTree->Branch("av_NHH",&av_NHH,"av_NHH/F");
	muTree->Branch("AL_Ecal",&AL_Ecal,"AL_Ecal/I");
	muTree->Branch("AL_Hcal",&AL_Hcal,"AL_Hcal/I");
	muTree->Branch("FD_ECALF10",&FD_ECALF10,"FD_ECALF10/F");
	muTree->Branch("FD_ECALL20",&FD_ECALL20,"FD_ECALL20/F");
	muTree->Branch("NH_ECALF10",&NH_ECALF10,"NH_ECALF10/I");
	muTree->Branch("NH_ECALL20",&NH_ECALL20,"NH_ECALL20/I");
	muTree->Branch("dEdx",&dEdx,"dEdx/F");
	muTree->Branch("cosTheta",&cosTheta,"cosTheta/F");
	muTree->Branch("Phi",&Phi,"Phi/F");

	eTree = new TTree("Electron","Electron");

	eTree->Branch("EventNr", &eventNr, "EventNr/I");
	eTree->Branch("MCPDG", &MCPDG, "MCPDG/I");
	eTree->Branch("ClusterID", &ClusterID, "ClusterID/I");
	eTree->Branch("NPFO", &NPFO, "NPFO/I");
	eTree->Branch("EcalNHit",&EcalNHit,"EcalNHit/I");
	eTree->Branch("HcalNHit",&HcalNHit,"HcalNHit/I");
	eTree->Branch("CluNHit",&CluNHit,"CluNHit/I");
	eTree->Branch("NLEcal",&NLEcal,"NLEcal/I");
	eTree->Branch("NLHcal",&NLHcal,"NLHcal/I");
	eTree->Branch("maxDisHtoL",&maxDisHtoL,"maxDisHtoL/F");
	eTree->Branch("minDisHtoL",&minDisHtoL,"minDisHtoL/F");
	eTree->Branch("avDisHtoL",&avDisHtoL,"avDisHtoL/F");
	eTree->Branch("avEnDisHtoL",&avEnDisHtoL,"avEnDisHtoL/F");
	eTree->Branch("EcalEn",&EcalEn,"EcalEn/F");
	eTree->Branch("HcalEn",&HcalEn,"HcalEn/F");
	eTree->Branch("EClu",&EClu,"EClu/F");
	eTree->Branch("graDepth",&graDepth,"graDepth/F");
	eTree->Branch("cluDepth",&cluDepth,"cluDepth/F");
	eTree->Branch("graAbsDepth",&graAbsDepth,"graAbsDepth/F");
	eTree->Branch("maxDepth",&maxDepth,"maxDepth/F");
	eTree->Branch("minDepth",&minDepth,"minDepth/F");
	eTree->Branch("MaxDisHel",&MaxDisHel,"MaxDisHel/F");
	eTree->Branch("MinDisHel",&MinDisHel,"MinDisHel/F");
	eTree->Branch("FD_all",&FD_all,"FD_all/F");
	eTree->Branch("FD_ECAL",&FD_ECAL,"FD_ECAL/F");
	eTree->Branch("FD_HCAL",&FD_HCAL,"FD_HCAL/F");
	eTree->Branch("crdis",&crdis,"crdis/F");
	eTree->Branch("EEClu_L10",&EEClu_L10,"EEClu_L10/F");
	eTree->Branch("EEClu_R",&EEClu_R,"EEClu_R/F");
	eTree->Branch("EEClu_r",&EEClu_r,"EEClu_r/F");
	eTree->Branch("EEClu_p",&EEClu_p,"EEClu_p/F");
	eTree->Branch("rms_Ecal",&rms_Ecal,"rms_Ecal/F");
	eTree->Branch("rms_Hcal",&rms_Hcal,"rms_Hcal/F");
	eTree->Branch("rms_Ecal2",&rms_Ecal2,"rms_Ecal2/F");
	eTree->Branch("rms_Hcal2",&rms_Hcal2,"rms_Hcal2/F");
	eTree->Branch("av_NHE",&av_NHE,"av_NHE/F");
	eTree->Branch("av_NHH",&av_NHH,"av_NHH/F");
	eTree->Branch("AL_Ecal",&AL_Ecal,"AL_Ecal/I");
	eTree->Branch("AL_Hcal",&AL_Hcal,"AL_Hcal/I");
	eTree->Branch("FD_ECALF10",&FD_ECALF10,"FD_ECALF10/F");
	eTree->Branch("FD_ECALL20",&FD_ECALL20,"FD_ECALL20/F");
	eTree->Branch("NH_ECALF10",&NH_ECALF10,"NH_ECALF10/I");
	eTree->Branch("NH_ECALL20",&NH_ECALL20,"NH_ECALL20/I");
	eTree->Branch("dEdx",&dEdx,"dEdx/F");
	eTree->Branch("cosTheta",&cosTheta,"cosTheta/F");
	eTree->Branch("Phi",&Phi,"Phi/F");

	otherTree = new TTree("Other","Other");

	otherTree->Branch("EventNr", &eventNr, "EventNr/I");
	otherTree->Branch("MCPDG", &MCPDG, "MCPDG/I");
	otherTree->Branch("TrkEn", &TrkEn, "TrkEn/F");
	otherTree->Branch("ClusterID", &ClusterID, "ClusterID/I");
	otherTree->Branch("NPFO", &NPFO, "NPFO/I");
	otherTree->Branch("EcalNHit",&EcalNHit,"EcalNHit/I");
	otherTree->Branch("HcalNHit",&HcalNHit,"HcalNHit/I");
	otherTree->Branch("CluNHit",&CluNHit,"CluNHit/I");
	otherTree->Branch("NLEcal",&NLEcal,"NLEcal/I");
	otherTree->Branch("NLHcal",&NLHcal,"NLHcal/I");
	otherTree->Branch("maxDisHtoL",&maxDisHtoL,"maxDisHtoL/F");
	otherTree->Branch("minDisHtoL",&minDisHtoL,"minDisHtoL/F");
	otherTree->Branch("avDisHtoL",&avDisHtoL,"avDisHtoL/F");
	otherTree->Branch("avEnDisHtoL",&avEnDisHtoL,"avEnDisHtoL/F");
	otherTree->Branch("EcalEn",&EcalEn,"EcalEn/F");
	otherTree->Branch("HcalEn",&HcalEn,"HcalEn/F");
	otherTree->Branch("EClu",&EClu,"EClu/F");
	otherTree->Branch("graDepth",&graDepth,"graDepth/F");
	otherTree->Branch("cluDepth",&cluDepth,"cluDepth/F");
	otherTree->Branch("graAbsDepth",&graAbsDepth,"graAbsDepth/F");
	otherTree->Branch("maxDepth",&maxDepth,"maxDepth/F");
	otherTree->Branch("minDepth",&minDepth,"minDepth/F");
	otherTree->Branch("MaxDisHel",&MaxDisHel,"MaxDisHel/F");
	otherTree->Branch("MinDisHel",&MinDisHel,"MinDisHel/F");
	otherTree->Branch("FD_all",&FD_all,"FD_all/F");
	otherTree->Branch("FD_ECAL",&FD_ECAL,"FD_ECAL/F");
	otherTree->Branch("FD_HCAL",&FD_HCAL,"FD_HCAL/F");
	otherTree->Branch("crdis",&crdis,"crdis/F");
	otherTree->Branch("EEClu_L10",&EEClu_L10,"EEClu_L10/F");
	otherTree->Branch("EEClu_R",&EEClu_R,"EEClu_R/F");
	otherTree->Branch("EEClu_r",&EEClu_r,"EEClu_r/F");
	otherTree->Branch("EEClu_p",&EEClu_p,"EEClu_p/F");
	otherTree->Branch("rms_Ecal",&rms_Ecal,"rms_Ecal/F");
	otherTree->Branch("rms_Hcal",&rms_Hcal,"rms_Hcal/F");
	otherTree->Branch("rms_Ecal2",&rms_Ecal2,"rms_Ecal2/F");
	otherTree->Branch("rms_Hcal2",&rms_Hcal2,"rms_Hcal2/F");
	otherTree->Branch("av_NHE",&av_NHE,"av_NHE/F");
	otherTree->Branch("av_NHH",&av_NHH,"av_NHH/F");
	otherTree->Branch("AL_Ecal",&AL_Ecal,"AL_Ecal/I");
	otherTree->Branch("AL_Hcal",&AL_Hcal,"AL_Hcal/I");
	otherTree->Branch("FD_ECALF10",&FD_ECALF10,"FD_ECALF10/F");
	otherTree->Branch("FD_ECALL20",&FD_ECALL20,"FD_ECALL20/F");
	otherTree->Branch("NH_ECALF10",&NH_ECALF10,"NH_ECALF10/I");
	otherTree->Branch("NH_ECALL20",&NH_ECALL20,"NH_ECALL20/I");
	otherTree->Branch("dEdx",&dEdx,"dEdx/F");
	otherTree->Branch("cosTheta",&cosTheta,"cosTheta/F");
	otherTree->Branch("Phi",&Phi,"Phi/F");

	
	evtTree = new TTree("Evt","Evt");
	evtTree->Branch("EventNr", &eventNr, "EventNr/I");
	evtTree->Branch("NPFO", &NPFO, "NPFO/I");
	evtTree->Branch("NMuon", &NMuon, "NMuon/I");
	//evtTree->Branch("NPion", &NPion, "NPion/I");
	//evtTree->Branch("NElec", &NElec, "NElec/I");
	//evtTree->Branch("", &, "");
	//evtTree->Branch("", &, "");
	//evtTree->Branch("", &, "");

	if(_Training==0)
	{

		reader_all.clear();
		for(int i=0;i<11;i++){
			reader1.clear();
			for(int j=0;j<4;j++){
				TMVA::Reader *reader;
		      		TString wgtFname;
			
		   		TString slevel = Form("%i",inputEn[i]);
				TString spos=inputDetectorModules[j].c_str();
		   		TString fwgtFname = weightDir + "/weights/TMVAMulticlass_"+slevel+"GeV_"+spos+"_BDTG.weights.xml";
				if(!access(fwgtFname,0)){
					wgtFname= weightDir + "/weights/TMVAMulticlass_"+slevel+"GeV_"+spos+"_BDTG.weights.xml";
				}
				else{
					wgtFname=weightDir + "/weights/TMVAMulticlass_"+slevel+"GeV_all_BDTG.weights.xml";
				}
				TMVA::Tools::Instance();
		   		reader = new TMVA::Reader( "!Color:Silent" );
				reader->AddVariable("EcalNHit", &_EcalNHit);
				reader->AddVariable("HcalNHit", &_HcalNHit);
				reader->AddVariable("NLEcal", &_NLEcal);
				reader->AddVariable("NLHcal", &_NLHcal);
				reader->AddVariable("maxDisHtoL",&maxDisHtoL);
				reader->AddVariable("avDisHtoL",&avDisHtoL);
				reader->AddVariable("EE := EcalEn/EClu",&EE);
				reader->AddVariable("graDepth",&graDepth);
				reader->AddVariable("cluDepth",&cluDepth);
				reader->AddVariable("minDepth",&minDepth);
				reader->AddVariable("MaxDisHel",&MaxDisHel);
				reader->AddVariable("FD_all",&FD_all);
				reader->AddVariable("FD_ECAL",&FD_ECAL);
				reader->AddVariable("FD_HCAL",&FD_HCAL);
				reader->AddVariable("E_10 := EEClu_L10/EcalEn",&E_10);
				reader->AddVariable("E_R := EEClu_R/EcalEn",&E_R);
				reader->AddVariable("E_r := EEClu_r/EcalEn",&E_r);
				reader->AddVariable("rms_Hcal",&rms_Hcal);
				reader->AddVariable("av_NHH", &av_NHH);
				reader->AddVariable("AL_Ecal", &_AL_Ecal);
				reader->AddVariable("FD_ECALF10",&FD_ECALF10);
				reader->AddVariable("FD_ECALL20",&FD_ECALL20);
				reader->AddVariable("NH_ECALF10", &_NH_ECALF10);
				reader->AddVariable("dEdx",&dEdx);
			
				E_R=EEClu_R/EcalEn;
				E_r=EEClu_r/EcalEn;
				E_10=EEClu_L10/EcalEn;
				EE=EcalEn/EClu;
		   			reader->BookMVA("BDTG", wgtFname);
				reader1.push_back(reader);
			}
			reader_all.push_back(reader1);
     		} 
	}
}



void LICH::NewClusterFlag(Cluster* a_tree, Track* a_trk)
{

	int NH[16];
	
	dEdx=0;

	const string ECALCellIDDecoder = "M:3,S-1:3,I:9,J:9,K-1:6";
	CellIDDecoder<CalorimeterHit> idDecoder(ECALCellIDDecoder);
	const float mass = 0.139;       //Pion Mass

	HelixClass * TrkInit_Helix = new HelixClass();
	TrkInit_Helix->Initialize_Canonical(a_trk->getPhi(), a_trk -> getD0(), a_trk -> getZ0(), a_trk -> getOmega(), a_trk->getTanLambda(), BField);
	float TrackEn = mass*mass;


	for (int q3 = 0; q3 < 3; q3 ++)
	{
		TrackEn += (TrkInit_Helix->getMomentum()[q3])*(TrkInit_Helix->getMomentum()[q3]);
	}
	delete TrkInit_Helix;

	TrackEn = sqrt(TrackEn);
	int nSubTrk = a_trk->getTracks().size();
	int NHit=0;
	if(a_tree){
	NHit = a_tree->getCalorimeterHits().size();
	cout<<"a_tree"<<endl;
	}
	
	if ( (NHit > 4 && TrackEn > 1) || TrackEn <= 1 )
	{

		for(int t1 = 0; t1 < nSubTrk; t1++)
		{
			Track* a_SubTrk = a_trk->getTracks()[t1];
			std::vector<TrackerHit*> trhits = a_SubTrk->getTrackerHits();
                        int nhit = trhits.size();
                        int chose_high = int(0.3*nhit);
                        int chose_low = int(0.08*nhit);
                        float cen_high = 999.;
                        float cen_low = -999.;
                        int maxhit = 99999;
                        int minhit = 99999;
                        float totaled = 0.;
			for(int ihit = 0; ihit < chose_low;ihit++){

                        	minhit = ihit;
                                TrackerHit * hiti = dynamic_cast<TrackerHit*>(trhits[minhit]);
                                TrackerHit * hiti1 = dynamic_cast<TrackerHit*>(trhits[minhit+1]);
                                TVector3 hitposi = hiti->getPosition();
                                TVector3 lastposi = hiti1->getPosition();
                                float cudisi = (hitposi-lastposi).Mag();
                                float mindedx = hiti->getEDep()/cudisi;
                                for(int jhit = ihit+1; jhit<nhit-1; jhit++){

                                	TrackerHit * hitj = dynamic_cast<TrackerHit*>(trhits[jhit]);

                                        TrackerHit * hitj1 = dynamic_cast<TrackerHit*>(trhits[jhit+1]);
                                        TVector3 hitposj = hitj->getPosition();
                                        TVector3 lastposj = hitj1->getPosition();
                                        float cudisj = (hitposj-lastposj).Mag();
                                        if(hitj->getEDep()/cudisj<mindedx&&hitj->getEDep()/cudisj>cen_low){
                                        	minhit = jhit;
                                        	mindedx = hitj->getEDep()/cudisj;
                                        }
                                }
                                cen_low = mindedx;
                        }
                        for(int ihit = 0; ihit < chose_high;ihit++){

                        	maxhit = ihit;
                                TrackerHit * hiti = dynamic_cast<TrackerHit*>(trhits[maxhit]);
                                TrackerHit * hiti1 = dynamic_cast<TrackerHit*>(trhits[maxhit+1]);
                                TVector3 hitposi = hiti->getPosition();
                                TVector3 lastposi = hiti1->getPosition();
                                float cudisi = (hitposi-lastposi).Mag();
                                float maxdedx = hiti->getEDep()/cudisi;
                                for(int jhit = ihit+1; jhit<nhit-1; jhit++){

                                	TrackerHit * hitj = dynamic_cast<TrackerHit*>(trhits[jhit]);

                                        TrackerHit * hitj1 = dynamic_cast<TrackerHit*>(trhits[jhit+1]);
                                        TVector3 hitposj = hitj->getPosition();
                                        TVector3 lastposj = hitj1->getPosition();
                                        float cudisj = (hitposj-lastposj).Mag();
                                        if(hitj->getEDep()/cudisj>maxdedx&&hitj->getEDep()/cudisj<cen_high){
                                        	maxhit = jhit;
                                        	maxdedx = hitj->getEDep()/cudisj;
                                        }
                                }
                                cen_high = maxdedx;
                        }
                        int nhiteff = 0;


                        for(int ihit = 0; ihit < nhit-1; ihit++)
                        {
                             	TrackerHit * hit = dynamic_cast<TrackerHit*>(trhits[ihit]);

                                float mindis = 999.;
                                TrackerHit * last = dynamic_cast<TrackerHit*>(trhits[ihit+1]);
                                TVector3 hitpos = hit->getPosition();
                                TVector3 lastpos = last->getPosition();
                                float cudis = (hitpos-lastpos).Mag();
                                mindis= cudis;
                               	if(hit->getEDep()/mindis<cen_high&&hit->getEDep()/mindis>cen_low){
                                	float dedx = hit->getEDep()/mindis;
                                	totaled += dedx;
                                	nhiteff ++;
                        	}
                        }
                        dEdx = totaled/nhiteff;
		}
				
		if(a_tree){
		TVector3 CluPos;
		CluPos = a_tree->getPosition();
		cosTheta=CluPos[2]/CluPos.Mag();
		Phi=atan2(CluPos[1],CluPos[0]);
                TVector3 IntDir = ClusterCoG(a_tree)-CluPos;
		EClu = a_tree->getEnergy();
		EcalNHit = 0;
		HcalNHit = 0;
		CluNHit = 0;
		EcalEn = 0;
		HcalEn = 0;
                float currDepth = 0;
		maxDepth = -100;
		minDepth = 1E6;
                MaxDisHel = -1;   //maximal distance from Track to Helix
                MinDisHel = 1E10;

                EEClu_R = 0;
                EEClu_r = 0;
                EEClu_p = 0;
                EEClu_L10 = 0;


		std::vector<CalorimeterHit*> Ecalhits;
		std::vector<CalorimeterHit*> Hcalhits;
		std::vector<CalorimeterHit*> allhits;
		std::vector<CalorimeterHit*> EH_1;
		std::vector<CalorimeterHit*> EH_2;
		std::vector<CalorimeterHit*> EH_3;
		std::vector<CalorimeterHit*> EH_4;
                std::vector<CalorimeterHit*> EH_5;
                std::vector<CalorimeterHit*> EH_6;
		std::vector<CalorimeterHit*> HH_1;
		std::vector<CalorimeterHit*> HH_2;
		std::vector<CalorimeterHit*> HH_3;
		std::vector<CalorimeterHit*> HH_4;
		std::vector<CalorimeterHit*> HH_5;
		std::vector<CalorimeterHit*> HH_6;
                std::vector<CalorimeterHit*> HH_7;
                std::vector<CalorimeterHit*> HH_8;
                std::vector<CalorimeterHit*> HH_9;
                std::vector<CalorimeterHit*> HH_0;
                std::vector<CalorimeterHit*> Ecalf10hits;
                std::vector<CalorimeterHit*> Ecall20hits;



		allhits.clear();
		Ecalhits.clear();
		Hcalhits.clear();
		EH_1.clear();
		EH_2.clear();
		EH_3.clear();
		EH_4.clear();
                EH_5.clear();
                EH_6.clear();
		HH_1.clear();
		HH_2.clear();
		HH_3.clear();
		HH_4.clear();
		HH_5.clear();
		HH_6.clear();
                HH_7.clear();
                HH_8.clear();
                HH_9.clear();
                HH_0.clear();
                Ecalf10hits.clear();
                Ecall20hits.clear();


        	HelixClass * currHelix = new HelixClass();
                currHelix->Initialize_Canonical(a_trk->getPhi(), a_trk -> getD0(), a_trk -> getZ0(), a_trk -> getOmega(), a_trk->getTanLambda(), 3.5);  
                float BushDist[3] = {0, 0, 0};
                float BushTime = 0;
		
		std::vector<float> hitTheta;
		hitTheta.clear();
                                
                for(unsigned int j1 = 0; j1 < a_tree->getCalorimeterHits().size(); j1++)
                {
                          CalorimeterHit * a_hit = a_tree->getCalorimeterHits()[j1];
                          BushTime = currHelix->getDistanceToPoint((float*)a_hit->getPosition(), BushDist);
			  TVector3 tmpPos = a_hit->getPosition();
			  hitTheta.push_back(tmpPos.Theta());
                          if(BushTime > 0)
                          {
                                    if(BushDist[2] > MaxDisHel)
                                    {
                                           MaxDisHel = BushDist[2];
                                    }
                                    if(BushDist[2] < MinDisHel)
                                    {
                                           MinDisHel = BushDist[2];
                                    }
                          }
                }
                delete currHelix;

		float totTheta = 0;
		float avTheta = 0;
		float SDTheta;

                for(int t0 = 0; t0 < int(hitTheta.size()); t0++)
                {
                        float tmpTheta = hitTheta[t0];
                        totTheta += tmpTheta;
                }

                avTheta = totTheta/float(hitTheta.size());
                SDTheta = 0;

                for(int t1 = 0; t1 < int(hitTheta.size()); t1++)
                {
                        float tmpTheta = hitTheta[t1];
                        SDTheta += pow(tmpTheta-avTheta,2);
                }
                SDTheta = sqrt(SDTheta/float(hitTheta.size()));

		TVector3 HitPos;
		int currCluNHits = a_tree->getCalorimeterHits().size();
                CluNHit = currCluNHits;
		int index1 = 0, index2 = 0;

		for(int s1 = 0; s1 < currCluNHits; s1++)
		{
			CalorimeterHit * a_hit = a_tree->getCalorimeterHits()[s1];
			allhits.push_back(a_hit);
			int NLayer = idDecoder(a_hit)["K-1"];

			HitPos = a_hit->getPosition();

			currDepth = DisSeedSurface(HitPos);
			crdis = (CluPos-HitPos).Mag()*sin((CluPos-HitPos).Angle(IntDir));


			if(currDepth > maxDepth)
			{
				maxDepth = currDepth;
				index1 = s1;
			}
			if(currDepth < minDepth)
			{
				minDepth = currDepth;
				index2 = s1;
			}

			if(SubDeFlag(HitPos)==1)
			{
				HcalNHit++;
				HcalEn += a_hit->getEnergy();
				Hcalhits.push_back(a_hit);
				if(NLayer < 5)
				{
					HH_1.push_back(a_hit);
				}
				else if(NLayer < 10)
				{
					HH_2.push_back(a_hit);
				}
				else if(NLayer < 15)
                                {
                                        HH_3.push_back(a_hit);
                                }
                                else if(NLayer < 20)
                                {
                                        HH_4.push_back(a_hit);
                                }
                                else if(NLayer < 25)
                                {
                                        HH_5.push_back(a_hit);
                                }
                                else if(NLayer < 30)
                                {
                                        HH_6.push_back(a_hit);
                                }
				else if(NLayer < 35)
				{
					HH_7.push_back(a_hit);
				}
				else if(NLayer < 40)
				{
					HH_8.push_back(a_hit);
				}
				else if(NLayer < 45)
				{
					HH_9.push_back(a_hit);
				}
				else
				{
					HH_0.push_back(a_hit);
				}
			}
			else if(SubDeFlag(HitPos)==0)
			{
				EcalNHit++;
				EcalEn += a_hit->getEnergy();
				Ecalhits.push_back(a_hit);
                                if(NLayer< 10) Ecalf10hits.push_back(a_hit);
                                else Ecall20hits.push_back(a_hit);
				if(crdis < 22) EEClu_R += a_hit->getEnergy();
                                if(crdis < 11) EEClu_r += a_hit->getEnergy();
                                if(crdis < 6) EEClu_p += a_hit->getEnergy();
				if(NLayer < 5)
				{
					EH_1.push_back(a_hit);
					EEClu_L10 += a_hit->getEnergy();
				}
				else if(NLayer < 10)
				{
					EH_2.push_back(a_hit);
					EEClu_L10 += a_hit->getEnergy();
				}
				else if(NLayer < 15)
                                {
                                        EH_3.push_back(a_hit);
                                } 
				else if(NLayer < 20)
                                {
                                        EH_4.push_back(a_hit);
                                } 
				else if(NLayer < 25)
                                {
                                        EH_5.push_back(a_hit);
                                }
				else
				{
					EH_6.push_back(a_hit);
				}
			}
		}
                
		if(a_tree->getCalorimeterHits().size()>0){
			CalorimeterHit * maxdis_hit = a_tree->getCalorimeterHits()[index1];
			CalorimeterHit * mindis_hit = a_tree->getCalorimeterHits()[index2];
			TVector3 maxpos = maxdis_hit->getPosition();
			TVector3 minpos = mindis_hit->getPosition();
			TVector3 GraPos = ClusterCoG(a_tree);
			graDepth = DisSeedSurface(GraPos);
			cluDepth = (maxpos-minpos).Mag();
			graAbsDepth = (GraPos-minpos).Mag();



			maxDisHtoL = -100;
			minDisHtoL = 1E6;

			float totDisHtoL = 0;

			float totHitEn = 0;
			float totHitEnDis = 0;
			float HitEn;

			for(int s2 = 0; s2 < currCluNHits; s2++)
			{
				CalorimeterHit * a_hit2 = a_tree->getCalorimeterHits()[s2];
				HitPos = a_hit2->getPosition();
				HitEn  = a_hit2->getEnergy();
				TVector3 par1 = GraPos-minpos;
				TVector3 par2 = minpos-HitPos;
				TVector3 par3 = par1.Cross(par2);
				float disHtoL = par3.Mag()/par1.Mag();
				totDisHtoL+=disHtoL;
				totHitEn+=HitEn;
				totHitEnDis+=HitEn*disHtoL;
				if (disHtoL > maxDisHtoL) maxDisHtoL = disHtoL;
				if (disHtoL < minDisHtoL) minDisHtoL = disHtoL;

			}
			avDisHtoL = totDisHtoL/currCluNHits;
			avEnDisHtoL = totHitEnDis/totHitEn;
			FD_all = FDV2(allhits, ECALCellIDDecoder);
			FD_ECAL = FDV2(Ecalhits, ECALCellIDDecoder);
			FD_HCAL = FDV2(Hcalhits, ECALCellIDDecoder);
        	        FD_ECALF10 = FDV2(Ecalf10hits, ECALCellIDDecoder);
        	        NH_ECALF10 = Ecalf10hits.size();
        	        FD_ECALL20 = FDV2(Ecall20hits, ECALCellIDDecoder);
        	        NH_ECALL20 = Ecall20hits.size();

			NLEcal = 0;
			NLHcal = 0;
			for(int p0 = 0; p0 < 8; p0++)
        	        {
        	        	NH[p0] = 0;
        	        }

			NH[0] = EH_1.size();
			NH[1] = EH_2.size();
			NH[2] = EH_3.size();
			NH[3] = EH_4.size();
			NH[4] = EH_5.size();
			NH[5] = EH_6.size();
			NH[6] = HH_1.size();
			NH[7] = HH_2.size();
			NH[8] = HH_3.size();
			NH[9] = HH_4.size();
			NH[10] = HH_5.size();
			NH[11] = HH_6.size();
			NH[12] = HH_7.size();
			NH[13] = HH_8.size();
        	        NH[14] = HH_9.size();
        	        NH[15] = HH_0.size();

			NLEcal = ActiveLayers(Ecalhits, ECALCellIDDecoder);
			NLHcal = ActiveLayers(Hcalhits, ECALCellIDDecoder);

		

			float sum_NHE = 0, sum_NHH = 0;
        	        av_NHE = 0;
        	        av_NHH = 0;
        	        AL_Ecal = 0;
        	        AL_Hcal = 0;

        	        for(int r1 = 0; r1 < 16; r1++)
			{
				if(r1 < 6 && NH[r1]>0)
        	                {
        	                        sum_NHE += NH[r1];
        	                        AL_Ecal++;
        	                }
        	                if(r1 >= 6 && NH[r1]>0)
        	                {
        	                        sum_NHH += NH[r1];
        	                        AL_Hcal++;
        	                }
        	        }
        	        if(AL_Ecal > 0)
        	                av_NHE = sum_NHE/AL_Ecal;
        	        if(AL_Hcal > 0)
        	                av_NHH = sum_NHH/AL_Hcal;

        	        rms_Ecal = 0;
        	        rms_Hcal = 0;
			rms_Ecal2 = 0;
			rms_Hcal2 = 0;	
        	        for(int r0 = 0; r0 < 16; r0++)
        	        {
        	                if(r0 < 6)
        	                {
        	                        if(NH[r0] > 0)
        	                        {
        	                                rms_Ecal+=pow(NH[r0]-av_NHE,2);
        	                                rms_Ecal2 += pow(NH[r0],2);
        	                        }
        	                }
        	                else
        	                {
        	                        if(NH[r0] > 0)
        	                        {
        	                                rms_Hcal+=pow(NH[r0]-av_NHH,2);
        	                                rms_Hcal2 += pow(NH[r0],2);
        	                        }
        	                }
        	        }
        	        if(AL_Ecal > 0)
        	        {
        	                rms_Ecal2 = sqrt(rms_Ecal2/AL_Ecal);
        	                rms_Ecal = sqrt(rms_Ecal/AL_Ecal);
        	        }
        	        else
        	        {
        	                rms_Ecal2 = -1;
        	                rms_Ecal = -1;
        	        }
        	        if(AL_Hcal > 0)
        	        {
        	                rms_Hcal2 = sqrt(rms_Hcal2/AL_Ecal);
        	                rms_Hcal = sqrt(rms_Hcal/AL_Ecal);
        	        }
        	        else
        	        {
        	                rms_Hcal2 = -1;
        	                rms_Hcal = -1;
        	        }
		}
		}
		if(_Training==0){	
			_EcalNHit=(float)EcalNHit;
			_HcalNHit=(float)HcalNHit;
			_NLEcal=(float)NLEcal;
			_NLHcal=(float)NLHcal;
			_AL_Ecal=(float)AL_Ecal;
			_NH_ECALF10=(float)NH_ECALF10;
	
			E_R=EEClu_R/EcalEn;
			E_r=EEClu_r/EcalEn;
			E_10=EEClu_L10/EcalEn;
			EE=EcalEn/EClu;
			
			TMVA::Reader *reader(0);
			std::vector<TMVA::Reader*> reader_1;
			for(int i=0;i<NEn;i++){
				if(TrackEn>inputLevel[i]&&TrackEn<inputLevel[i+1]){
					reader_1=reader_all[i];
					for(int j=0;j<NPos-1;j++){
						if(abs(cosTheta)>inputPos[j]&&abs(cosTheta)<inputPos[j+1]){
							reader=reader_1[j];
						}
						else{
							reader=reader_1[2];
						}
					}
				}
			}
			
   			mvaVal_pi=0.;
   			mvaVal_mu=0.;
   			mvaVal_e=0.;

   			mvaVal_pi = (reader->EvaluateMulticlass("BDTG"))[0];
   			mvaVal_mu = (reader->EvaluateMulticlass("BDTG"))[1];
   			mvaVal_e = (reader->EvaluateMulticlass("BDTG"))[2];
			cout<<"e:mu:pi "<<mvaVal_e<<" "<<mvaVal_mu<<" "<<mvaVal_pi<<endl;
			if(mvaVal_e>mvacut_e&&mvaVal_mu<mvacut_mu&&mvaVal_pi<mvacut_pi){
				ClusterID = -11;
			}
			if(mvaVal_pi>mvacut_pi&&mvaVal_mu<mvacut_mu&&mvaVal_e<mvacut_e){
				ClusterID = 211;
			}
			if(mvaVal_mu>mvacut_mu&&mvaVal_pi<mvacut_pi&&mvaVal_e<mvacut_e){
				ClusterID = -13;
			}
		}	
	}
}

void LICH::SamplePro( LCEvent * evtP )
{
	if (evtP)
	{
		eventNr=evtP->getEventNumber();
		MCPDG=0;
		
		try{
			LCCollection * MCP = evtP->getCollection(_inputMCP.c_str());
			for(int i = 0; i < MCP->getNumberOfElements(); i++)
			{
				MCParticle * a_MCP = dynamic_cast<MCParticle*>(MCP->getElementAt(i));
				if( a_MCP->getParents().size() == 0 )
				{
					MCPDG=a_MCP->getPDG();
				}
			}
			LCCollection * RecoCol=evtP->getCollection(_inputPFO.c_str());
			NPFO=RecoCol->getNumberOfElements();

			if(NPFO==1&&_Training==1){
				ReconstructedParticle * a_pfo=dynamic_cast<ReconstructedParticle*>(RecoCol->getElementAt(0));
				Track * a_trk(0);
				if(a_pfo->getTracks().size())a_trk= a_pfo->getTracks()[0];

				Cluster * a_clu(0);
				if(a_pfo->getClusters().size())a_clu= a_pfo->getClusters()[0];
				if(a_pfo->getCharge()>0)NewClusterFlag(a_clu,a_trk);
				if(abs(MCPDG)==211){
					piTree->Fill();
				}
				else if(abs(MCPDG)==13){
					muTree->Fill();
				}
				else if(abs(MCPDG)==11){
					eTree->Fill();
				}
				else{
					otherTree->Fill();
				}

			}

			if(_Training==0){
				LCCollection *arborrecoparticle = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
				LCCollection * newMCP = evtP->getCollection("MCParticle");
				try{
				LCCollection * Trks = evtP->getCollection("ClupatraTracks");
				if(Trks->getNumberOfElements()>0){
				Track * a_cluTrk=dynamic_cast<Track*>(Trks->getElementAt(0));
				Cluster *a_empClu(0);
				if(NPFO==0){
					cout<<Trks->getNumberOfElements()<<"ClupatraTracks"<<endl;

					NewClusterFlag(a_empClu,a_cluTrk);
					otherTree->Fill();
				}
				}
				}catch(lcio::DataNotAvailableException err){}
				MCParticle * a_MCP = dynamic_cast<MCParticle*>(newMCP->getElementAt(0));
				TVector3 Mom=a_MCP->getMomentum();
				TrkEn=Mom.Mag();
				NMuon=0;
				cout<<eventNr<<"event have "<<NPFO<<" PFOs"<<endl;
				for(int i=0;i<NPFO;i++){

					ReconstructedParticleImpl * chargeparticle = new ReconstructedParticleImpl();
					ReconstructedParticle * b_pfo=dynamic_cast<ReconstructedParticle*>(RecoCol->getElementAt(i));
					ParticleIDImpl * pid=new ParticleIDImpl();
					Track * a_trk(0);
					if(b_pfo->getTracks().size())a_trk=b_pfo->getTracks()[0];
					Cluster * a_clu(0);
					if(b_pfo->getClusters().size())a_clu= b_pfo->getClusters()[0];
					if(b_pfo->getCharge()!=0){
						chargeparticle->addParticle(b_pfo);
						chargeparticle->setEnergy(b_pfo->getEnergy());
						chargeparticle->addTrack(b_pfo->getTracks()[0]);
						chargeparticle->setMomentum(b_pfo->getMomentum() );
						chargeparticle->setCharge(b_pfo->getCharge());
						chargeparticle->addCluster(b_pfo->getClusters()[0]);
						NewClusterFlag(a_clu,a_trk);
						otherTree->Fill();
						//cout<<ClusterID<<" "<<b_pfo->getType();
						//if(abs(b_pfo->getType())==13){NMuon++;}
						if(abs(ClusterID)==13){NMuon++;}
						pid->setPDG(int(ClusterID*b_pfo->getCharge()));
						chargeparticle->setType(int(ClusterID*b_pfo->getCharge()));
						chargeparticle->addParticleID(pid);
						arborrecoparticle->addElement(chargeparticle);
					}
					else{
						chargeparticle->addParticle(b_pfo);
						arborrecoparticle->addElement(chargeparticle);
					}
				}
						evtTree->Fill();
				cout<<NMuon<<endl;
				evtP->addCollection( arborrecoparticle, _outputPFO.c_str());

			}
		}catch(lcio::DataNotAvailableException err){}
	}
}


void LICH::processEvent( LCEvent * evtP )
{

	if (evtP)
	{
		int eventNr = evtP->getEventNumber();

		if(eventNr%1 == 0)
			cout<<"Nevts Processed: "<<eventNr<<endl;

		inputLevel.clear();
		NEn=inputEn.size();
		NPos=inputPos.size();
		for(int i=0;i<NEn;i++){
			if(i==0){inputLevel.push_back(0.0);}
			else{inputLevel.push_back(float(inputEn[i-1]+inputEn[i])/2.);}
		}
		inputLevel.push_back(9999.);


		LICH::gearPara();	
		LICH::SamplePro( evtP );	
	}
}

void LICH::end()
{
	std::cout<<"Bush Connection Finished, ArborObject Formed"<<std::endl;

		tree_file->Write();
		delete tree_file;


}

