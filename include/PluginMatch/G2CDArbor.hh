#ifndef _G2CDArbor_hh_
#define _G2CDArbor_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <EVENT/CalorimeterHit.h>
#include <TNtuple.h>
#include <TObject.h>

#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
class TTree;

// namespace CALICE {

class G2CDArbor  : public marlin::Processor
{
	public:

	Processor*  newProcessor() { return new G2CDArbor ; }

	G2CDArbor();

	~G2CDArbor() {};

	void init();

	void processEvent( LCEvent * evtP );

	void end();

    void CreateROOT(std::string str);
    void CloseROOT();
	protected:
    std::string _filename;
    TFile *file;
    TTree *tree;
    TTree *Htree;
    TTree *tecal;

	std::string _treeFileName;
	std::string _treeName;
	std::string _colName;
	std::vector<std::string> _caloTruthLinkCollection;
	std::vector<std::string> _hcalCollections;
	std::vector<std::string> _outputHcalCollections;
	std::vector<std::string> _ecalCollections;
        std::vector<std::string> _outputEcalCollections;
	std::vector<std::string> _EcalPreShowerCollections;
	std::vector<float> _ChargeSpatialDistri;
	std::vector<float> _thresholdHcal;
	std::vector<float> _calibCoeffEcal;
	std::vector<float> _ShowerPositionShiftID; //should be of the form deltaI, J, K
        std::map <int, std::pair<float, float> >WeightVector;

	float _thresholdEcal;
	int _NEcalThinLayer; 
	int _overwrite;
	int _DigiCellSize; 
	int _UsingDefaultDetector; 
	float _PolyaParaA, _PolyaParaB, _PolyaParaC; 
	float _ChanceOfKink, _KinkHitChargeBoost; 
	TTree *_outputTree;
	TH1F *_NH1stLayer, *_NH8thLayer; 
	TF1 * _QPolya; 

	int _Num;
	int _eventNr; 

	int _M, _S, _I, _J, _K, _Seg;
	float _PosX, _PosY, _PosZ; 
	float _EDepo, _Charge, _ShiftInX; 
	int _NHit1mm, _NHit1mmCenter, _NHit1mmCorner, _NHit1mmSide;
	int _TotalNHit1mm, _TotalNHit, _TotalMultiHit; 
	int _N1, _N2, _N3; 

	std::string _fileName;
	std::ostream *_output;
};

#endif


