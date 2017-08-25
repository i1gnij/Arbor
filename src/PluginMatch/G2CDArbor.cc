#include <G2CDArbor.hh>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include "UTIL/CellIDDecoder.h"
#include "DetectorPos.hh"

#include <values.h>
#include <string>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <Rtypes.h>
#include <TMath.h>
#include <TF1.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector3.h>
#include <TFile.h>

using namespace std;
using namespace lcio ;
using namespace marlin ;

static float digiEn1;
static float digiEn2;
static int eventID;
static float eSim ;
static float eDigi;
//static float hSim ;
static float hDigi;
static TVector3 posSim;
static TVector3 posDigi;
const double pi = acos(-1.0);

struct DigiHit {
    int digihitCellID0; 
    int digihitCellID1;
    float digihitCharge; 
    float digihitEnergyDepo; 
    int digihitNum1mmCell; 
    float PosX;
    float PosY;
    float PosZ; 
    float LeadChargeDepo; 
    float ChargeShare; 
    SimCalorimeterHit * LeadSimCaloHit; 
} ;

//std::map <int, std::pair<float, float> >WeightVector;

float ChanceOfKink = 0.1;
float KinkHitChargeBoost = 2.2; 

G2CDArbor aG2CDArbor ;
G2CDArbor::G2CDArbor()
    : Processor("G2CDArbor"),
    _output(0)
{
    _description = "Gaseous Calorimeter Digitizer: take 1mm Simulated Hits as input";

    std::vector<std::string> ecalCollections;
    ecalCollections.push_back(std::string("EcalBarrelSiliconCollection"));
    ecalCollections.push_back(std::string("EcalEndcapSiliconCollection"));
    ecalCollections.push_back(std::string("EcalEndcapRingCollection"));
    registerInputCollections( LCIO::SIMCALORIMETERHIT,
            "ECALCollections" ,
            "Input ECAL Hits Collection Names" ,
            _ecalCollections ,
            ecalCollections);

    std::vector<std::string> outputEcalCollections;
    outputEcalCollections.push_back(std::string("ECALBarrel"));
    outputEcalCollections.push_back(std::string("ECALEndcap"));
    outputEcalCollections.push_back(std::string("ECALOther"));
    registerProcessorParameter("DigiECALCollection" ,
            "Name of Digitized ECAL Hit Collections" ,
            _outputEcalCollections,
            outputEcalCollections);

    std::vector<std::string> EcalPreShowerCollections;
    EcalPreShowerCollections.push_back("EcalBarrelSiliconPreShowerCollection");
    EcalPreShowerCollections.push_back("EcalEndcapRingPreShowerCollection");
    EcalPreShowerCollections.push_back("EcalEndcapSiliconPreShowerCollection");
    registerInputCollections( LCIO::CALORIMETERHIT,      //adding a flag to Calo /SimCalo
            "HcalHitCollections" ,
            "Hit Collection Names" ,
            _EcalPreShowerCollections ,
            EcalPreShowerCollections);

    std::vector<float> calibrEcal;
    calibrEcal.push_back(40.91);
    calibrEcal.push_back(81.81);
    registerProcessorParameter("CalibrECAL" ,
            "Calibration coefficients for ECAL" ,
            _calibCoeffEcal,
            calibrEcal);

    registerProcessorParameter("NumThinEcalLayer" ,
            "Num of thiner Ecal layers" ,
            _NEcalThinLayer,
            20);

    registerProcessorParameter("ECALThreshold" ,
            "Threshold for ECAL Hits in GeV" ,
            _thresholdEcal,
            (float)5.0e-5);

    std::vector<std::string> hcalCollections;
    hcalCollections.push_back(std::string("HcalBarrelRegCollection"));
    hcalCollections.push_back(std::string("HcalEndcapsCollection"));
    hcalCollections.push_back(std::string("HcalEndcapRingsCollection"));
    registerInputCollections( LCIO::SIMCALORIMETERHIT,
            "HCALCollections" ,
            "HCAL Collection Names" ,
            _hcalCollections ,
            hcalCollections);

    std::vector<std::string> outputHcalCollections;
    outputHcalCollections.push_back(std::string("HCALBarrel"));
    outputHcalCollections.push_back(std::string("HCALEndcap"));
    outputHcalCollections.push_back(std::string("HCALOther"));
    registerProcessorParameter("DigiHCALCollection" ,
            "Name of Digitized HCAL Hit Collections" ,
            _outputHcalCollections,
            outputHcalCollections);

    std::vector<std::string> caloTruthLinkCollection;
    registerProcessorParameter("caloTruthLinkCollection",
            "caloTruthLinkCollection",
            _caloTruthLinkCollection,
            caloTruthLinkCollection);

    std::vector<float> ChargeSpatialDistri; 
    ChargeSpatialDistri.push_back(0.1);
    ChargeSpatialDistri.push_back(0.2);
    ChargeSpatialDistri.push_back(0.4);
    ChargeSpatialDistri.push_back(0.2);
    ChargeSpatialDistri.push_back(0.1);
    registerProcessorParameter("ChargeSpatialDistribution" ,
            "Spactial Distribution of MIP charge X*Y;" ,
            _ChargeSpatialDistri,
            ChargeSpatialDistri);

    std::vector<float> ShowerPositionShiftID; 
    ShowerPositionShiftID.push_back(0);
    ShowerPositionShiftID.push_back(0);
    ShowerPositionShiftID.push_back(0);
    registerProcessorParameter("PositionShiftID" ,
            "Global Position Shift For Overlay" ,
            _ShowerPositionShiftID,
            ShowerPositionShiftID);

    registerProcessorParameter( "DigiCellSize" ,
            "Size of Digitized Cell (in mm)" ,
            _DigiCellSize ,
            10);

    registerProcessorParameter( "ShiftInX" ,
            "Shift Distance in X directoin (in mm) NP only" ,
            _ShiftInX ,
            float(0.0));

    registerProcessorParameter( "UsingDefaultDetector",
            "Flag Parameter Setting (0 ~ self definition, 1 ~ MircoMegas, 2 ~ GRPC_PS, 3 ~ GRPC_SPS)",
            _UsingDefaultDetector ,
            0);

    registerProcessorParameter( "PolyaParaA" ,
            "Polya: x^A*exp(-b*x) + c" ,
            _PolyaParaA ,
            float(0.7));

    registerProcessorParameter( "PolyaParaB" ,
            "Polya: x^a*exp(-B*x) + c" ,
            _PolyaParaB ,
            float(0.045));

    registerProcessorParameter( "PolyaParaC" ,
            "Polya: x^a*exp(-b*x) + C" ,
            _PolyaParaC,
            float(0.03) );	

    registerProcessorParameter( "ChanceOfKink" ,
            "Chance of one boundary hit to create a multiple hit with boosted charge" ,
            _ChanceOfKink,
            float(0.0) );

    registerProcessorParameter( "KinkHitChargeBoost" ,
            "Scale factor of Charge on boosted multiple hits" ,
            _KinkHitChargeBoost,
            float(1.0) );

    std::string filename = "g2cd.root";
    registerProcessorParameter( "G2CD_ROOT" ,
            "G2CD_ROOT" ,
            _filename,
            filename);


}

void G2CDArbor::init() {

    printParameters();
    WeightVector.clear();
    CreateROOT(_filename);


    if( _outputHcalCollections.size() < _hcalCollections.size() )
    {
        cout<<"WARNING! Output Collections Array Smaller than Input"<<endl;
        exit(1);
    }
    else if( _ChargeSpatialDistri.size() % 2 == 0 )
    {
        cout<<"WARNING! Better Organize Charge Spatial distribution array in odd integer length"<<endl;
        exit(2);
    }

    if(_UsingDefaultDetector < 0 || _UsingDefaultDetector > 4)
    {
        cout<<"Parameter Flag Wrong. Need to be 0 (selfdefine), 1 (MircoMegas), 2 (GRPC_PS) or 3 (GRPC_SPS)"<<endl;
        exit(0);
    }
    else if(_UsingDefaultDetector == 1)		//MircoMegas Small Chamber, 2008 PS TB
    {
        _PolyaParaA = 0.7;
        _PolyaParaB = 0.045;
        _PolyaParaC = 0.03;
        _ChargeSpatialDistri.clear();
        _ChargeSpatialDistri.push_back(0.1);
        _ChargeSpatialDistri.push_back(8.0);
        _ChargeSpatialDistri.push_back(0.1);
        ChanceOfKink = 0.095;
        KinkHitChargeBoost = 2.2;
    }
    else if(_UsingDefaultDetector == 2)		//GRPC Cubic Meter, 2011 PS TB
    {
        _PolyaParaA = 4.0;
        _PolyaParaB = 4.5;
        _PolyaParaC = 0.0;
        _ChargeSpatialDistri.clear();
        _ChargeSpatialDistri.push_back(0.1);
        _ChargeSpatialDistri.push_back(0.1);
        _ChargeSpatialDistri.push_back(0.1);
        _ChargeSpatialDistri.push_back(0.1);
        _ChargeSpatialDistri.push_back(0.1);
        _ChargeSpatialDistri.push_back(0.1);
        _ChargeSpatialDistri.push_back(0.1);
        ChanceOfKink = 0.0;
        KinkHitChargeBoost = 1.0;
    }
    else if(_UsingDefaultDetector == 3)		//GRPC Bubic Meter, 2012 May SPS TB
    {
        _PolyaParaA = 1.1;
        _PolyaParaB = 1.0;
        _PolyaParaC = 0.0;
        _ChargeSpatialDistri.clear();
        _ChargeSpatialDistri.push_back(0.1);
        _ChargeSpatialDistri.push_back(0.2);
        _ChargeSpatialDistri.push_back(0.3);
        _ChargeSpatialDistri.push_back(0.6);
        _ChargeSpatialDistri.push_back(0.6);
        _ChargeSpatialDistri.push_back(0.6);
        _ChargeSpatialDistri.push_back(0.3);
        _ChargeSpatialDistri.push_back(0.2);
        _ChargeSpatialDistri.push_back(0.1);
        ChanceOfKink = 0.0;
        KinkHitChargeBoost = 1.0;
    }
    else
    {
        ChanceOfKink = _ChanceOfKink;
        KinkHitChargeBoost = _KinkHitChargeBoost; 
    }

    float PolyaDomain = 100;
    if(_PolyaParaB )
    {
        PolyaDomain = 20*_PolyaParaA/_PolyaParaB;
    }
    _QPolya = new TF1("QPolya", "x^[0]*exp(-1*x*[1]) + [2]", 0, PolyaDomain);
    _QPolya->SetParameters(_PolyaParaA, _PolyaParaB, _PolyaParaC);

    cout<<"Parameters After Para Review: "<<ChanceOfKink<<", "<<KinkHitChargeBoost <<", "<<_DigiCellSize<<endl;
    printParameters();

    int SizeCSD = _ChargeSpatialDistri.size();
    int CoverAreaLength = int((SizeCSD - 1)/2);
    if(CoverAreaLength > _DigiCellSize - 1)
    {
        cout<<"WARNING! CoverAreaLength is too large comparing to the Digitized Cell Size: to have maximally 4 multiple hit, CoverAreaLength must <= _DigiCellSize - 1"<<endl;
        exit(0);
    }
    int tmpIndex = 0;
    float NormalWeightFactor = 0;
    float NormalWeight[SizeCSD];
    int IndexA = 0; 
    int IndexB = 0;		//Used to denote charge distributions...
    float WeightA = 0; 
    float WeightB = 0; 

    for( int i0 = 0; i0 < SizeCSD; i0++ )
    {
        NormalWeightFactor += _ChargeSpatialDistri[i0];
    }

    for( int i1 = 0; i1 < SizeCSD; i1++ )
    {
        NormalWeight[i1] = _ChargeSpatialDistri[i1]/NormalWeightFactor;
    }

    int SignX(0);
    int SignY(0);
    pair<float, float> WMatrix;		//if Vec(A) != Vec(B)

    for(int i2 = 0; i2 < _DigiCellSize; i2++)
    {
        for(int j2 = 0; j2 < _DigiCellSize; j2++)
        {
            tmpIndex = _DigiCellSize*i2 + j2;
            IndexA = 0;
            IndexB = 0;
            WeightA = 0;
            WeightB = 0;

            if( i2 < CoverAreaLength ) 
            {
                IndexA = CoverAreaLength - i2;
                SignX = -1; 
            }
            else if( i2 > _DigiCellSize - CoverAreaLength - 1)
            {
                IndexA = CoverAreaLength + i2 - _DigiCellSize + 1;
                SignX = 1; 
            }

            if( j2 < CoverAreaLength ) 
            {
                IndexB = CoverAreaLength - j2;
                SignY = -1; 
            }
            else if( j2 > _DigiCellSize - CoverAreaLength - 1)
            {
                IndexB = CoverAreaLength + j2 - _DigiCellSize + 1;
                SignY = 1; 
            }

            for(int i3 = 0; i3 < CoverAreaLength; i3 ++)        
            {
                if(i3 < IndexA) WeightA += NormalWeight[i3];
                if(i3 < IndexB) WeightB += NormalWeight[i3];
            }

            WMatrix.first = SignX*WeightA; 
            WMatrix.second = SignY*WeightB; 
            WeightVector[tmpIndex] = WMatrix;	

            cout<<WMatrix.first<<"/"<<WMatrix.second<<", ";
        }
        cout<<endl;
    }
}

void G2CDArbor::processEvent( LCEvent * evtP )
{
    if (evtP)
    {
        try{

            _eventNr = evtP->getEventNumber();
            eventID = _eventNr;
            if(_eventNr % 100 == 0) cout<<"eventNr: "<<_eventNr<<endl;

            digiEn1 = 0.0;
            digiEn2 = 0.0;

            float HitEn = 0;
            float DigiHitEn = 0; 
            int LayerNum = 0;
            TVector3 HitPos; 
            int tmpM, tmpS, tmpI, tmpJ, tmpK;
            int CurrI = 0; 
            int CurrJ = 0; 
            int DHIndexI = 0;
            int DHIndexJ = 0; 
            float DHChargeWeight = 0;
            int DHCellID0 = 0; 
            float RndCharge = 0;
            int SingleMCPPID = 0; 
            //float SingleMCPPEn = 0; 
            float RefPosX = 0;
            float RefPosY = 0;
            float RefPosZ = 0;
            float DeltaPosI = 0;
            float DeltaPosJ = 0;
            std::vector<float> CurrWeightVec;
            float WeiI = 0;
            float WeiJ = 0;
            int DeltaI = 0;
            int DeltaJ = 0;
            int MapI = 0;
            int MapJ = 0;
            int MapIndex = 0;
            float FHitPos[3] = {0, 0, 0};
            float SimHitPos[3] = {0,0,0};
            float currTime = 0;
            float HitStepEn = 0;
            float EmaxStep = 0;


            LCCollectionVec *relcol = new LCCollectionVec(LCIO::LCRELATION);
            LCFlagImpl linkflag;
            linkflag.setBit(LCIO::CHBIT_LONG);
            relcol->setFlag(linkflag.getFlag());

            LCCollection * MCPCol = evtP->getCollection("MCParticle");
            for ( int s0(0); s0 < MCPCol->getNumberOfElements(); s0++)
            {
                MCParticle * a_MCP = dynamic_cast<MCParticle*> (MCPCol->getElementAt(s0));
                if( a_MCP->getParents().size() == 0 )
                {
                    SingleMCPPID = a_MCP->getPDG();
                    //SingleMCPPEn = a_MCP->getEnergy();
                    break;
                }
            }

            LCFlagImpl flag;
            flag.setBit(LCIO::CHBIT_LONG);                  //To set position & ID1
            flag.setBit(LCIO::CHBIT_ID1);
            flag.setBit(LCIO::RCHBIT_ENERGY_ERROR); //In order to use an additional FLOAT
            flag.setBit(LCIO::RCHBIT_TIME);

            for (unsigned int k0 = 0; k0 < _ecalCollections.size(); ++k0)
            {
                try{
                    LCCollection *Ecalcol = evtP->getCollection( _ecalCollections[k0].c_str() ) ;
                    CellIDDecoder<SimCalorimeterHit> idDecoder( Ecalcol );
                    int NumEcalhit = Ecalcol->getNumberOfElements();
                    LCCollectionVec *ecalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
                    ecalcol->setFlag(flag.getFlag());
                    string EcalinitString = Ecalcol->getParameters().getStringVal(LCIO::CellIDEncoding);
                    ecalcol->parameters().setValue(LCIO::CellIDEncoding, EcalinitString);

                    for(int k1 = 0; k1 < NumEcalhit; k1++)
                    {
                        SimCalorimeterHit * SimEcalhit = dynamic_cast<SimCalorimeterHit*>( Ecalcol->getElementAt( k1 ) ) ;
                        HitEn = SimEcalhit->getEnergy();
                        LayerNum = idDecoder(SimEcalhit)["K-1"];

                        currTime = 0;
                        EmaxStep = 0;
                        HitStepEn = 0;
                        for(int k=0; k<SimEcalhit->getNMCContributions(); k++)
                        {
                            HitStepEn = SimEcalhit->getEnergyCont(k);
                            if(HitStepEn > EmaxStep)
                            {
                                EmaxStep = HitStepEn;
                                currTime = SimEcalhit->getTimeCont(k);
                            }
                        }

                        // jing
                        if(LayerNum < _NEcalThinLayer)
                        {
                            DigiHitEn = HitEn * _calibCoeffEcal[0];
                            digiEn1 += HitEn;
                        }
                        else
                        {
                            DigiHitEn = HitEn * _calibCoeffEcal[1];
                            digiEn2 += HitEn;
                        }

                        // jing
                        //if (HitEn > _thresholdEcal) tecal->Fill();

                        if(HitEn > _thresholdEcal)
                        {
                            CalorimeterHitImpl * DigiEcalhit = new CalorimeterHitImpl();
                            FHitPos[0] = SimEcalhit->getPosition()[0] + _ShowerPositionShiftID[0]*10.0;
                            FHitPos[1] = SimEcalhit->getPosition()[1] + _ShowerPositionShiftID[1]*10.0;
                            FHitPos[2] = SimEcalhit->getPosition()[2] + _ShowerPositionShiftID[2]*10.0;

                            DigiEcalhit->setTime(currTime);
                            DigiEcalhit->setPosition( FHitPos );
                            DigiEcalhit->setCellID0(SimEcalhit->getCellID0());
                            DigiEcalhit->setCellID1(SimEcalhit->getCellID1());
                            DigiEcalhit->setEnergy(DigiHitEn);
                            ecalcol->addElement(DigiEcalhit);
                            LCRelationImpl *rel = new LCRelationImpl(DigiEcalhit, SimEcalhit, 1.0);    //only keep the leading contribution
                            relcol->addElement(rel);

                            eDigi = DigiHitEn;
                            eSim = HitEn;
                            tree->Fill();

                        }
                    }

                    evtP->addCollection(ecalcol,_outputEcalCollections[k0].c_str());

                }catch(lcio::DataNotAvailableException zero) { }
            }

            string k = evtP->getDetectorName();

            if(k=="CEPC_v1")	// 1cm cell simulation...
            {

                cout<<"Detector Tpype: "<<k<<endl; 

                for (unsigned int k3 = 0; k3 < _hcalCollections.size(); ++k3)
                {
                    try{
                        LCCollection *Hcalcol = evtP->getCollection( _hcalCollections[k3].c_str() ) ;
                        int NumHcalhit = Hcalcol->getNumberOfElements();
                        LCCollectionVec *hcalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
                        hcalcol->setFlag(flag.getFlag());
                        string HcalinitString = Hcalcol->getParameters().getStringVal(LCIO::CellIDEncoding);
                        hcalcol->parameters().setValue(LCIO::CellIDEncoding, HcalinitString);

                        for(int k4 = 0; k4 < NumHcalhit; k4++)
                        {
                            SimCalorimeterHit * SimHcalhit = dynamic_cast<SimCalorimeterHit*>( Hcalcol->getElementAt( k4 ) ) ;

                            currTime = 0;
                            EmaxStep = 0;
                            HitStepEn = 0;

                            for(int k=0; k<SimHcalhit->getNMCContributions(); k++)
                            {
                                HitStepEn = SimHcalhit->getEnergyCont(k);
                                if(HitStepEn > EmaxStep)
                                {
                                    EmaxStep = HitStepEn;
                                    currTime = SimHcalhit->getTimeCont(k);
                                }
                            }


                            if(SimHcalhit->getEnergy() > 1.e-9)	//some threshold can be added.
                            {
                                CalorimeterHitImpl * DigiHcalhit = new CalorimeterHitImpl();

                                FHitPos[0] = SimHcalhit->getPosition()[0] + _ShowerPositionShiftID[0]*10.0;
                                FHitPos[1] = SimHcalhit->getPosition()[1] + _ShowerPositionShiftID[1]*10.0;
                                FHitPos[2] = SimHcalhit->getPosition()[2] + _ShowerPositionShiftID[2]*10.0;


                                DigiHcalhit->setTime(currTime);
                                DigiHcalhit->setPosition( FHitPos );
                                DigiHcalhit->setCellID0(SimHcalhit->getCellID0());
                                DigiHcalhit->setCellID1(SimHcalhit->getCellID1());
                                DigiHcalhit->setEnergy(DHCALCalibrationConstant);
                                hcalcol->addElement(DigiHcalhit);
                                LCRelationImpl *rel = new LCRelationImpl(DigiHcalhit, SimHcalhit, 1.0);    //only keep the leading contribution
                                relcol->addElement(rel);
                                hDigi = SimHcalhit->getEnergy(); 
                                //Htree->Fill();
                            }
                        }

                        evtP->addCollection(hcalcol,_outputHcalCollections[k3].c_str());

                    }catch(lcio::DataNotAvailableException zero) { }
                }

                LCCollectionVec *ecalPScol = new LCCollectionVec(LCIO::CALORIMETERHIT);
                ecalPScol->setFlag(flag.getFlag());
                string EcalPSinitString = "M:3,S-1:3,I:9,J:9,K-1:6";
                ecalPScol->parameters().setValue(LCIO::CellIDEncoding, EcalPSinitString);

                //  parameter CellIDEncoding [string]: M:3,S-1:3,I:9,J:9,K-1:6, 

                for(unsigned int s0 = 0; s0 < _EcalPreShowerCollections.size(); s0++)   //I think should digitize and give a very small energy; even a new collection called PS Hits
                {
                    try
                    {
                        LCCollection * PSHitColl = evtP ->getCollection(_EcalPreShowerCollections[s0].c_str());
                        int NHitsCurrCol = PSHitColl->getNumberOfElements();
                        for(int s1 = 0; s1 < NHitsCurrCol; s1++)
                        {
                            SimCalorimeterHit * a_hit = dynamic_cast<SimCalorimeterHit*>(PSHitColl->getElementAt(s1));
                            currTime = 0;
                            EmaxStep = 0;
                            HitStepEn = 0;

                            for(int k=0; k<a_hit->getNMCContributions(); k++)
                            {
                                HitStepEn = a_hit->getEnergyCont(k);
                                if(HitStepEn > EmaxStep)
                                {
                                    EmaxStep = HitStepEn;
                                    currTime = a_hit->getTimeCont(k);
                                }
                            }

                            CalorimeterHitImpl * PShit = new CalorimeterHitImpl();

                            FHitPos[0] = a_hit->getPosition()[0] + _ShowerPositionShiftID[0]*10.0;
                            FHitPos[1] = a_hit->getPosition()[1] + _ShowerPositionShiftID[1]*10.0;
                            FHitPos[2] = a_hit->getPosition()[2] + _ShowerPositionShiftID[2]*10.0;

                            PShit->setTime(currTime);
                            PShit->setPosition( FHitPos );
                            PShit->setCellID0(a_hit->getCellID0());
                            PShit->setCellID1(a_hit->getCellID1());
                            PShit->setEnergy(a_hit->getEnergy());
                            ecalPScol->addElement(PShit);
                            LCRelationImpl *rel = new LCRelationImpl(PShit, a_hit, 1.0);
                            relcol->addElement(rel);
                            float tPShitEnergy = a_hit->getEnergy();
                        }
                    }catch(lcio::DataNotAvailableException zero){}
                }

                //evtP->addCollection(ecalPScol, "ECALPSHitCollection");
            }
            else
            {
                for (unsigned int i(0); i < _hcalCollections.size(); ++i) 
                {		//strictly follow barrel, endcap, ring order

                    std::map <int, DigiHit> IDtoDigiHit; 
                    IDtoDigiHit.clear();
                    std::map <int, TVector3> IDtoPos; 
                    IDtoPos.clear();
                    std::map <int, float> IDtoTime;
                    IDtoTime.clear();


                    try{
                        LCCollection * col = evtP->getCollection( _hcalCollections[i].c_str() ) ;

                        int numOrihits = col->getNumberOfElements();
                        CellIDDecoder<SimCalorimeterHit> idDecoder(col);
                        LCCollectionVec *hcalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
                        string initString = "M:3,S-1:3,I:9,J:9,K-1:6";		//Need to verify
                        hcalcol->parameters().setValue(LCIO::CellIDEncoding, initString);
                        hcalcol->setFlag(flag.getFlag());

                        for (int j(0); j < numOrihits; ++j) {
                            SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
                            HitEn = hit->getEnergy();

                            currTime = 0;
                            EmaxStep = 0;
                            HitStepEn = 0;

                            for(int k=0; k<hit->getNMCContributions(); k++)
                            {
                                HitStepEn = hit->getEnergyCont(k);
                                if(HitStepEn > EmaxStep)
                                {
                                    EmaxStep = HitStepEn;
                                    currTime = hit->getTimeCont(k);
                                }
                            }


                            tmpM = idDecoder(hit)["M"];
                            tmpS = idDecoder(hit)["S-1"];
                            tmpI = idDecoder(hit)["I"];
                            tmpJ = idDecoder(hit)["J"];
                            tmpK = idDecoder(hit)["K-1"];

                            RefPosX = hit->getPosition()[0];
                            RefPosY = hit->getPosition()[1];
                            RefPosZ = hit->getPosition()[2];

                            CurrI = tmpI/_DigiCellSize; 
                            CurrJ = tmpJ/_DigiCellSize;
                            MapI = tmpI % _DigiCellSize;
                            MapJ = tmpJ % _DigiCellSize;
                            MapIndex = MapI * _DigiCellSize + MapJ; 
                            WeiI = WeightVector[MapIndex].first; 
                            WeiJ = WeightVector[MapIndex].second;

                            DeltaPosI = (float(_DigiCellSize) - 1.0)/2 - MapI;
                            DeltaPosJ = (float(_DigiCellSize) - 1.0)/2 - MapJ;

                            // cout<<"DeltaI "<<DeltaPosI<<", "<<DeltaPosJ<<endl;

                            DeltaI = 0; 
                            DeltaJ = 0; 

                            if(fabs(WeiI) > 1E-9) 
                            {	
                                DeltaI = ((WeiI > 0) ? 1: -1);
                            }
                            if(fabs(WeiJ) > 1E-9)
                            {
                                DeltaJ = ((WeiJ > 0) ? 1: -1);
                            }

                            // cout<<"Weight... "<<WeiI<<"\t"<<DeltaI<<"\t"<<WeiJ<<"\t"<<DeltaJ<<endl;

                            RndCharge = _QPolya->GetRandom();

                            CurrWeightVec.clear();
                            CurrWeightVec.push_back((1 - fabs(WeiI))*(1 - fabs(WeiJ)));
                            if(DeltaI && !DeltaJ)
                            {
                                CurrWeightVec.push_back(fabs(WeiI));
                            }
                            else if(DeltaJ && !DeltaI)
                            {
                                CurrWeightVec.push_back(fabs(WeiJ));
                            }
                            else if(DeltaI && DeltaJ)
                            {
                                CurrWeightVec.push_back(fabs(WeiI) * (1 - fabs(WeiJ)));
                                CurrWeightVec.push_back(fabs(WeiJ) * (1 - fabs(WeiI)));
                                CurrWeightVec.push_back(fabs(WeiJ*WeiI));
                            }

                            for(int i0 = 0; i0 < int(CurrWeightVec.size()); i0++)
                            {
                                if(i0 == 0)
                                {
                                    DHIndexI = CurrI;
                                    DHIndexJ = CurrJ;
                                }
                                else if(i0 == 1)
                                {
                                    if(DeltaI)
                                    {
                                        DHIndexI = CurrI + DeltaI;
                                        DHIndexJ = CurrJ; 
                                    }
                                    if(DeltaJ)
                                    {
                                        DHIndexI = CurrI;
                                        DHIndexJ = CurrJ + DeltaJ;
                                    }
                                }
                                else if(i0 == 2)
                                {
                                    DHIndexI = CurrI + DeltaI;
                                    DHIndexJ = CurrJ;
                                }
                                else if(i0 == 3)
                                {
                                    DHIndexI = CurrI + DeltaI;
                                    DHIndexJ = CurrJ + DeltaJ;
                                }

                                DHChargeWeight = CurrWeightVec[i0];
                                DHCellID0 = (i<<30) + (tmpK<<24) + (DHIndexJ<<15) + (DHIndexI<<6) + (tmpS << 3) + tmpM;

                                if( IDtoDigiHit.find(DHCellID0) == IDtoDigiHit.end() && IDtoPos.find(DHCellID0) == IDtoPos.end() )
                                {
                                    IDtoDigiHit[DHCellID0].digihitCellID0 = DHCellID0;
                                    IDtoDigiHit[DHCellID0].digihitCharge = RndCharge * DHChargeWeight;
                                    IDtoDigiHit[DHCellID0].digihitEnergyDepo = HitEn * DHChargeWeight;	//Assumption...
                                    IDtoDigiHit[DHCellID0].digihitNum1mmCell = 1;
                                    IDtoDigiHit[DHCellID0].LeadChargeDepo = RndCharge * DHChargeWeight;
                                    IDtoDigiHit[DHCellID0].LeadSimCaloHit = hit;
                                    IDtoDigiHit[DHCellID0].ChargeShare = DHChargeWeight;

                                    if(i == 0)	//Barrel
                                    {
                                        IDtoDigiHit[DHCellID0].PosX = RefPosX + (DeltaPosI + int(DHIndexI - CurrI)*_DigiCellSize) * cos(tmpS*pi/4.0) + _ShowerPositionShiftID[0]*_DigiCellSize; //(mm in unit)
                                        IDtoDigiHit[DHCellID0].PosY = RefPosY + (DeltaPosI + int(DHIndexI - CurrI)*_DigiCellSize)* sin(tmpS*pi/4.0) + _ShowerPositionShiftID[1]*_DigiCellSize;
                                        IDtoDigiHit[DHCellID0].PosZ = RefPosZ + DeltaPosJ + int(DHIndexJ - CurrJ)*_DigiCellSize + _ShowerPositionShiftID[2]*_DigiCellSize;	//Rotation is needed, based on S; 
                                    }
                                    else	//endcap or ring;  
                                    {
                                        IDtoDigiHit[DHCellID0].PosX = RefPosX + DeltaPosI + (DHIndexI - CurrI)*_DigiCellSize + _ShowerPositionShiftID[0]*_DigiCellSize; //(mm in unit)
                                        IDtoDigiHit[DHCellID0].PosY = RefPosY + DeltaPosJ + (DHIndexJ - CurrJ)*_DigiCellSize + _ShowerPositionShiftID[1]*_DigiCellSize;
                                        IDtoDigiHit[DHCellID0].PosZ = RefPosZ;
                                    }

                                    HitPos.SetXYZ(IDtoDigiHit[DHCellID0].PosX, IDtoDigiHit[DHCellID0].PosY, IDtoDigiHit[DHCellID0].PosZ);

                                    IDtoPos[DHCellID0] = HitPos;
                                    IDtoTime[DHCellID0] = currTime;	
                                }
                                else
                                {
                                    IDtoDigiHit[DHCellID0].digihitCharge += RndCharge * DHChargeWeight;
                                    IDtoDigiHit[DHCellID0].digihitEnergyDepo += HitEn * DHChargeWeight;
                                    IDtoDigiHit[DHCellID0].digihitNum1mmCell ++;
                                    if(RndCharge * DHChargeWeight > IDtoDigiHit[DHCellID0].LeadChargeDepo)
                                    {
                                        IDtoDigiHit[DHCellID0].LeadChargeDepo = RndCharge * DHChargeWeight;
                                        IDtoDigiHit[DHCellID0].LeadSimCaloHit = hit; 
                                        IDtoDigiHit[DHCellID0].ChargeShare = DHChargeWeight;
                                    }

                                    if(IDtoDigiHit[DHCellID0].PosX != IDtoPos[DHCellID0].X() )
                                        cout<<"Position Changed.........."<<endl;
                                }

                                //cout<<"Com "<<IDtoPos.size()<<"\t"<<IDtoDigiHit.size()<<endl;
                            }
                        }

                        float DigiHitPos[3];

                        for(std::map <int, DigiHit>::iterator ff = IDtoDigiHit.begin(); ff!=IDtoDigiHit.end(); ff++)
                        {
                            CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
                            LCRelationImpl *rel = new LCRelationImpl(calhit, ff->second.LeadSimCaloHit, ff->second.ChargeShare);	//only keep the leading contribution
                            relcol->addElement(rel);

                            calhit->setCellID0( ff->first );	//Assume 100% efficiency
                            calhit->setEnergy( DHCALCalibrationConstant );		//Charge
                            calhit->setCellID1(SingleMCPPID);	//Use ID1 & Energy Error to denote the MCP info...
                            calhit->setEnergyError(ff->second.digihitCharge); // (SingleMCPPEn);
                            /*
                               DigiHitPos[0] = ff->second.PosX; 
                               DigiHitPos[1] = ff->second.PosY;
                               DigiHitPos[2] = ff->second.PosZ;
                               */
                            DigiHitPos[0] = IDtoPos[ff->first].X();
                            DigiHitPos[1] = IDtoPos[ff->first].Y();
                            DigiHitPos[2] = IDtoPos[ff->first].Z();

                            calhit->setTime(IDtoTime[ff->first]);
                            calhit->setPosition(DigiHitPos);		
                            hcalcol->addElement(calhit);
                        }

                        evtP->addCollection(hcalcol,_outputHcalCollections[i].c_str());
                        IDtoDigiHit.clear();
                    }
                    catch (lcio::DataNotAvailableException zero) { }
                }

            }

            evtP->addCollection(relcol, _caloTruthLinkCollection[0].c_str());
            tecal->Fill();
        }
        catch (lcio::DataNotAvailableException zero) { }
    }
}

void G2CDArbor::end()
{
    CloseROOT();
    std::cout<<"General Gas Digitizer FINISHED"<<std::endl;
}

void G2CDArbor::CreateROOT(std::string str){
    file = new TFile(str.c_str(),"RECREATE");
    tree = new TTree("tree","tree");
    tree->Branch("eSim",&eSim,"eSim/F");
    tree->Branch("eDigi",&eDigi,"eDigi/F");

    tecal = new TTree("tecal","tecal");
    tecal->Branch("digiEn1",&digiEn1,"digiEn1/F");
    tecal->Branch("digiEn2",&digiEn2,"digiEn2/F");
    tecal->Branch("eventID",&eventID,"eventID/I");
    //Htree = new TTree("Htree","Htree");
    //    Htree->Branch("hSim",&hSim,"hSim/F");
    //Htree->Branch("hDigi",&hDigi,"hDigi/F");

}

void G2CDArbor::CloseROOT(){
    file->cd();
    tecal->Write();
    tree->Write();
    //Htree->Write();
    std::cout << "ttree write here."<<std::endl;
    file->Close();
}
