#ifndef _Arbor_hh_
#define _Arbor_hh_

#include "ArborHit.hh"
#include <string>
#include <iostream>
#include <TVector3.h>
#include "ArborTool.hh"
#include <TTree.h>

static TFile *hits_file;
static TTree *tree;
void CreateROOT(std::string str);
void CloseROOT();

void LoadParam(const std::vector<Float_t>& vParam);
//Int_t GetIndex(Int_t layer, Int_t cell, Int_t idx);
//Float_t GetParaMean(Int_t layer, Int_t cell);
//Float_t GetParaRMS(Int_t layer, Int_t cell);

void init();

void HitsCleaning( std::vector<ArborHit> inputHits );

void HitsClassification( linkcoll inputLinks );

void BuildInitLink(std::vector<float> Thresholds);

void LinkIteration(int time);

void BranchBuilding(float SeedThreshold);

branchcoll Arbor( std::vector<ArborHit>, std::vector<float> Thresholds );

/*
 int NLayer_A, NStave_A, SubD_A; 
 int NLayer_B, NStave_B, SubD_B;
 float MagA, MagB, Depth_A, Depth_B, ECCorr, DisAB; 
 int FlagTrkPS, FlagEH; 
 int FlagPSEE, FlagHH;       
 int FlagStaveSame;
 int FlagStaveDiff;
 TVector3 PosA, PosB, PosDiffAB, PosDiffBA, linkDir;
*/

#endif


