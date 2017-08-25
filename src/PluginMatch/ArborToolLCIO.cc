#include "ArborToolLCIO.hh"
#include "ArborTool.hh"
#include "DetectorPos.hh"
#include <iostream>
#include <TMath.h>
#include "UTIL/CellIDDecoder.h"
#include "EVENT/CalorimeterHit.h"
#include "IMPL/ClusterImpl.h"
#include "IMPL/LCFlagImpl.h"
#include "HelixClass.hh"	

/*
//Fit Branch
#include "trajectory.h"
#include "fitting_root.h"
#include "geometry.h"
#include "point.h"
#include "segment3.h"
*/

using namespace std;

int NHScale( LCCollection *inputHit, Cluster * clu0, int RatioX, int RatioY, int RatioZ )
{

	int ReScaledNH = 0;
	int NumHit = clu0->getCalorimeterHits().size();
	int tmpI = 0;
	int tmpJ = 0;
	int tmpK = 0;
	float tmpEn = 0;
	int NewCellID0 = 0;
	int NewCellID1 = 0;

	CellIDDecoder<CalorimeterHit> idDecoder(inputHit);      //Input Hits here refer to AllCleanHits collection

	std::map <double, float> testIDtoEnergy;
	double testlongID = 0;

	for(int i = 0; i < NumHit; i++)
	{
		CalorimeterHit *hit = dynamic_cast<CalorimeterHit*>( clu0->getCalorimeterHits()[i]);

		tmpI = idDecoder(hit)["I"]/RatioX;
		tmpJ = idDecoder(hit)["J"]/RatioY;
		tmpK = (idDecoder(hit)["K-1"]+1)/RatioZ;
		tmpEn = hit->getEnergy();

		NewCellID0 = (tmpK<<24) + (tmpJ<<12) + tmpI;

		testlongID = NewCellID1*1073741824 + NewCellID0;
		if(testIDtoEnergy.find(testlongID) == testIDtoEnergy.end() )
		{
			testIDtoEnergy[testlongID] = tmpEn;
		}
		else
		{
			testIDtoEnergy[testlongID] += tmpEn;
		}
	}

	ReScaledNH = testIDtoEnergy.size();

	return ReScaledNH;

}

float FD( Cluster * clu, LCCollection *HitCollection )
{
	float FractalDim = 0;
	int NReSizeHit[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	int Scale[5] = {2, 3, 4, 5, 6};
	int OriNHit = clu->getCalorimeterHits().size();
	
	for(int j = 0; j < 5; j++)
	{
		NReSizeHit[j] = NHScale(HitCollection, clu, Scale[j], Scale[j], 1);
		FractalDim += 0.2 * TMath::Log(float(OriNHit)/NReSizeHit[j])/TMath::Log(float(Scale[j]));
	}

	return FractalDim;
}


bool HelixCluster_TightLink(float TrackEnergy, float currTrkTheta, float Dis)
{

/*	/besfs/groups/higgs/users/binsong/Software/Arbor/Deve/Arbor_v3/src/PluginMatch/FitSigMean.txt
0   -0.824153   6.13307   -7.27353   3.40849   -0.547119
0   0.883417   -1.56288   1.52334   -0.642984   0.100693
1   -1.03604   6.13911   -8.23187   3.99434   -0.63476
1   0.0644493   1.57968   -1.96033   0.921835   -0.145203
2   -0.750274   4.11594   -5.78318   2.83902   -0.450127
2   -0.155956   2.30511   -2.90032   1.36712   -0.214887
3   -0.447226   2.26029   -3.39447   1.68651   -0.265178
3   -0.228006   2.42051   -3.05811   1.44444   -0.227217
4   -0.208013   0.861233   -1.62842   0.851328   -0.133332
4   -0.0559016   1.52512   -1.92558   0.909111   -0.143227
5   -0.0156837   -0.16484   -0.327873   0.240926   -0.0381746
5   0.0553384   0.979369   -1.28817   0.625727   -0.100344
6   0.0428617   -0.540321   0.159683   0.00569509   -0.000808263
6   0.210229   0.284984   -0.445564   0.240462   -0.0412346
7   0.100135   -0.871974   0.590225   -0.199154   0.0314567
7   0.178681   0.323457   -0.430711   0.209505   -0.0337331
8   0.147677   -1.11072   0.869534   -0.324759   0.0504365
8   0.166636   0.279314   -0.32988   0.152814   -0.0240413
9   0.208789   -1.37252   1.16778   -0.460128   0.0719274
9   0.204023   0.127782   -0.131306   0.0551587   -0.00833517
*/	

	bool DisCut = false; 

	if(currTrkTheta < 3.142 && currTrkTheta > -0.001)
	{
		// float indexEn[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
		float p0m[10] = {-0.824153, -1.03604, -0.750274, -0.447226, -0.208013, -0.0156837, 0.0428617, 0.100135, 0.147677, 0.208789};
		float p1m[10] = {6.13307, 6.13911, 4.11594, 2.26029, 0.861233, -0.16484, -0.540321, -0.871974, -1.11072, -1.37252};
		float p2m[10] = {-7.27353, -8.23187, -5.78318, -3.39447, -1.62842, -0.327873, 0.159683, 0.590225, 0.869534, 1.16778};
		float p3m[10] = {3.40849, 3.99434, 2.83902, 1.68651, 0.851328, 0.240926, 0.00569509, -0.199154, -0.324759, -0.460128};
		float p4m[10] = {-0.547119, -0.63476, -0.450127, -0.265178, -0.133332, -0.0381746, -0.000808263, 0.0314567, 0.0504365, 0.0719274};

		float p0s[10] = {0.883417, 0.0644493, -0.155956, -0.228006, -0.0559016, 0.0553384, 0.210229, 0.178681, 0.166636, 0.204023};
		float p1s[10] = {-1.56288, 1.57968, 2.30511, 2.42051, 1.52512, 0.979369, 0.284984, 0.323457, 0.279314, 0.127782};
		float p2s[10] = {1.52334, -1.96033, -2.90032, -3.05811, -1.92558, -1.28817, -0.445564, -0.430711, -0.32988, -0.131306};
		float p3s[10] = {-0.642984, 0.921835, 1.36712, 1.44444, 0.909111, 0.625727, 0.240462, 0.209505, 0.152814, 0.0551587};
		float p4s[10] = {0.100693, -0.145203, -0.214887, -0.227217, -0.143227, -0.100344, -0.0412346, -0.0337331, -0.0240413, -0.00833517};

		float EnForCut = 0;
		if(TrackEnergy <= 1) EnForCut = 1.01;
		else if( TrackEnergy < 10) EnForCut = TrackEnergy;
		else EnForCut = 10.01;

		int EnIndex = int(EnForCut)-1;

		float meanEn1 = p0m[EnIndex] + p1m[EnIndex] * currTrkTheta + p2m[EnIndex] * pow(currTrkTheta,2) + p3m[EnIndex] * pow(currTrkTheta,3) + p4m[EnIndex] * pow(currTrkTheta,4);
		float meanEn2 = p0m[EnIndex+1] + p1m[EnIndex+1] * currTrkTheta + p2m[EnIndex+1] * pow(currTrkTheta,2) + p3m[EnIndex+1] * pow(currTrkTheta,3) + p4m[EnIndex+1] * pow(currTrkTheta,4);
		float Mean = (meanEn2-meanEn1)*(EnForCut-EnIndex)+meanEn1;

		float sigEn1 = p0s[EnIndex] + p1s[EnIndex] * currTrkTheta + p2s[EnIndex] * pow(currTrkTheta,2) + p3s[EnIndex] * pow(currTrkTheta,3) + p4s[EnIndex] * pow(currTrkTheta,4);
		float sigEn2 = p0s[EnIndex+1] + p1s[EnIndex+1] * currTrkTheta + p2s[EnIndex+1] * pow(currTrkTheta,2) + p3s[EnIndex+1] * pow(currTrkTheta,3) + p4s[EnIndex+1] * pow(currTrkTheta,4);
		float Sigma = (sigEn2-sigEn1)*(EnForCut-EnIndex)+sigEn1;

		DisCut = log10(Dis)-Mean < Sigma * 3;
	}

	return DisCut;
}

bool MIPFragFlag(Cluster * a_clu, float Dis, float TrkEn)
{
        float p0m[10] = {1.68984, 0.092209, -0.123757, 0.464582, -0.110746, 0.604014, -0.0155246, -0.107755, -0.267985, -0.0353545};
        float p1m[10] = {-0.00297416, 0.0053975, 0.00459457, 0.0018391, 0.00318535, 0.000772884, 0.00219315, 0.0020399, 0.00256288, 0.00157011};
        float p2m[10] = {1.1239e-06, -3.43177e-06, -2.80928e-06, -7.58474e-07, -1.6614e-06, -3.63653e-08, -9.30533e-07, -7.32234e-07, -1.23657e-06, -4.70653e-07};
        float sigma = 0.36;

        if(TrkEn < 1) TrkEn = 1;
        else if(TrkEn >= 10) TrkEn = 9.99;

        int EnIndex = int(TrkEn)-1;
        //bool flag = 0;
//       float CluFD = FDV3(a_clu,ECALCellIDDecoder);
        TVector3 SeedPos = a_clu->getPosition();
        float CluDepth = DisSeedSurface(SeedPos);
        //int CluSize = a_clu->getCalorimeterHits().size();
        float mean0 = p0m[EnIndex] + p1m[EnIndex]*CluDepth + p2m[EnIndex]*pow(CluDepth,2);
        float mean1 = p0m[EnIndex+1] + p1m[EnIndex+1]*CluDepth + p2m[EnIndex+1]*pow(CluDepth,2);
        float mean = mean0+(mean1-mean0)*(TrkEn-EnIndex);

        bool cutDepth = log10(Dis) - mean < 3 * sigma;

	// bool Zsuit; 
	// return Zsuit&&cutDepth

        return cutDepth;
}


int NHScaleV2( const std::string& encoder_str, std::vector<CalorimeterHit*> clu0, int RatioX, int RatioY, int RatioZ )
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

int ActiveLayers(  std::vector<CalorimeterHit*> clu, const std::string& encoder_str )
{
	std::vector<int> hitlayers; 
	CellIDDecoder<CalorimeterHit> idDecoder(encoder_str);

	int NHits = clu.size();	
	int tmpK = 0;	//Layer Number
	int tmpS = 0; 
	int tmpID = 0;

	for(int i = 0; i < NHits; i++)
	{
		CalorimeterHit *hit = dynamic_cast<CalorimeterHit*>( clu[i]);
		tmpK = idDecoder(hit)["K-1"]+1 ;
		tmpS = idDecoder(hit)["S-1"]+1 ;
		// cout<<"tmpK "<<tmpK<<endl; 
		tmpID = tmpS * 50 + tmpK;

		if( std::find(hitlayers.begin(), hitlayers.end(), tmpID) == hitlayers.end() )
		{
			hitlayers.push_back(tmpID);
		}
	}

	return hitlayers.size();

}

float FD_I( std::vector<CalorimeterHit*> clu, const std::string& encoder_str, int InitSize )	//Define initial size
{
	float FractalDim = 0;
	int NReSizeHit[5] = {0, 0, 0, 0, 0};
	int Scale[5] = {2, 3, 4, 5, 6};
	int OriNHit = NHScaleV2(encoder_str, clu, InitSize, InitSize, 1);;

	for(int j = 0; j < 5; j++)
	{
		NReSizeHit[j] = NHScaleV2(encoder_str, clu, Scale[j]*InitSize, Scale[j]*InitSize, 1);
		FractalDim += 0.2 * TMath::Log(float(OriNHit)/NReSizeHit[j])/TMath::Log(float(Scale[j]));
	}

	return FractalDim; 
}

float FDV2( std::vector<CalorimeterHit*> clu, const std::string& encoder_str )
{
	float FractalDim = 0;
	int NReSizeHit[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	int Scale[10] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 20};
	int OriNHit = clu.size();

	for(int j = 0; j < 10; j++)
	{
		NReSizeHit[j] = NHScaleV2(encoder_str, clu, Scale[j], Scale[j], 1);
		FractalDim += 0.1 * TMath::Log(float(OriNHit)/NReSizeHit[j])/TMath::Log(float(Scale[j]));
	}

	/*
	   if(FractalDim > 0.3)
	   {
	   cout<<FractalDim<<" at How many hits: "<<OriNHit<<endl; 
	   for(int t = 0; t < 10; t++)
	   {
	   cout<<"Scale & NH "<<Scale[t]<<" : "<<NReSizeHit[t]<<endl; 
	   }
	   }
	   */

	if(clu.size() == 0) 
		FractalDim = -1; 

	return FractalDim;
}

int NHScaleV3( const std::string& encoder_str, Cluster * clu0, int RatioX, int RatioY, int RatioZ )
{

	int ReScaledNH = 0;
	int NumHit = clu0->getCalorimeterHits().size();
	int tmpI = 0;
	int tmpJ = 0;
	int tmpK = 0;
	float tmpEn = 0;
	int NewCellID0 = 0;
	int NewCellID1 = 0;

	CellIDDecoder<CalorimeterHit> idDecoder(encoder_str);      //Input Hits here refer to AllCleanHits collection

	std::map <double, float> testIDtoEnergy;
	double testlongID = 0;

	for(int i = 0; i < NumHit; i++)
	{
		CalorimeterHit *hit = dynamic_cast<CalorimeterHit*>( clu0->getCalorimeterHits()[i]);

		tmpI = idDecoder(hit)["I"]/RatioX;
		tmpJ = idDecoder(hit)["J"]/RatioY;
		tmpK = (idDecoder(hit)["K-1"]+1)/RatioZ;
		tmpEn = hit->getEnergy();

		NewCellID0 = (tmpK<<24) + (tmpJ<<12) + tmpI;

		testlongID = NewCellID1*1073741824 + NewCellID0;
		if(testIDtoEnergy.find(testlongID) == testIDtoEnergy.end() )
		{
			testIDtoEnergy[testlongID] = tmpEn;
		}
		else
		{
			testIDtoEnergy[testlongID] += tmpEn;
		}
	}

	ReScaledNH = testIDtoEnergy.size();

	return ReScaledNH;

}

float FDV3( Cluster * clu, const std::string& encoder_str )
{
	float FractalDim = -1;
	int NReSizeHit[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	int Scale[5] = {2, 3, 4, 5, 6};
	int OriNHit = clu->getCalorimeterHits().size();

	if(clu->getCalorimeterHits().size() > 0)
	{
		FractalDim = 0.0;
		for(int j = 0; j < 5; j++)
		{
			NReSizeHit[j] = NHScaleV3(encoder_str, clu, Scale[j], Scale[j], 1);
			FractalDim += 0.2 * TMath::Log(float(OriNHit)/NReSizeHit[j])/TMath::Log(float(Scale[j]));
		}
	}

	return FractalDim;
}


float BushDis( Cluster *clu1, Cluster *clu2)
{
	float DisBetweenBush = 1.0E10; 

	int cluSize1 = clu1->getCalorimeterHits().size();
	int cluSize2 = clu2->getCalorimeterHits().size();

	TVector3 HitPos1, HitPos2; 
	TVector3 PosDiff; 
	// TVector3 XXXPos; 

	for(int i = 0; i < cluSize1; i++)
	{
		HitPos1 = (clu1->getCalorimeterHits()[i])->getPosition();
		for(int j = 0; j<cluSize2; j++)
		{
			HitPos2 = (clu2->getCalorimeterHits()[j])->getPosition();
			PosDiff = HitPos1 - HitPos2;

			if(PosDiff.Mag() < DisBetweenBush )
			{
				DisBetweenBush = PosDiff.Mag();
				//		XXXPos = PosDiff;
			}
		}
	}

	// cout<<" DisBetweenBush "<<XXXPos.X()<<" : "<<XXXPos.Y()<<" : "<<XXXPos.Z()<<" ~ "<<DisBetweenBush<<endl; 

	return DisBetweenBush; 
}

float DisPointToBush( TVector3 Pos1, Cluster * clu1)
{
	float Dis = 1.0E9; 
	float HitDis = 1.0E8;
	int clusize = clu1->getCalorimeterHits().size();

	TVector3 HitPos; 

	for(int s = 0; s < clusize; s++)
	{
		HitPos = (clu1->getCalorimeterHits()[s])->getPosition();
		HitDis = (HitPos - Pos1).Mag();
		if(HitDis < Dis) 
		{
			Dis = HitDis; 
		}
	}

	return Dis; 
}

TVector3 TrackOuterHit( Track* inputTrack )
{
	TVector3 OutHitPosition;

	int NTrkHits = inputTrack->getTrackerHits().size();

	float TrackHitMaxZ = -1;
	float tmpHitZ = -10;
	int OuterHitIndex = 0;

	for( int ihit = 0; ihit < NTrkHits; ihit++ )
	{
		TrackerHit * tmpHit = dynamic_cast<TrackerHit*> (inputTrack->getTrackerHits()[ihit]);
		tmpHitZ = fabs(tmpHit->getPosition()[2]);

		if(tmpHitZ > TrackHitMaxZ)
		{
			TrackHitMaxZ = tmpHitZ;
			OuterHitIndex = ihit;
		}
	}

	TrackerHit * OuterHit = dynamic_cast<TrackerHit*> (inputTrack->getTrackerHits()[OuterHitIndex]);

	OutHitPosition = OuterHit->getPosition();

	return OutHitPosition;
}


/*
   JSObj CaloTrackJS( Cluster *inputCluster )
   {
   JSObj calotrack; 

   int inputClusterSize = inputCluster->getCalorimeterHits().size();

   TVector3 tmphitPos; 
   TVector3 EndPointDiff; 

   calotrack.EndPos = inputCluster->getCalorimeterHits()[0]->getPosition();
   calotrack.BeginPos = inputCluster->getCalorimeterHits()[inputClusterSize -1]->getPosition();

   if(inputClusterSize == 1)
   {
   calotrack.BeginPos.SetXYZ(0, 0, 0);
   }

   if( calotrack.EndPos.Mag() < calotrack.BeginPos.Mag() )
   {
   tmphitPos = calotrack.EndPos;
   calotrack.EndPos = calotrack.BeginPos;
   calotrack.BeginPos = tmphitPos;
   }

   if(inputClusterSize < 300)
   {
   EndPointDiff = calotrack.EndPos - calotrack.BeginPos;
   EndPointDiff = float(1.0/(EndPointDiff.Mag() + 0.00001)) * EndPointDiff;
   calotrack.BeginMom = EndPointDiff;
   calotrack.EndMom = calotrack.BeginMom;
   }
   else		// Call Kalman Filter Or Anything Else
   {
   (calotrack.BeginMom).SetXYZ(0, 0, 0);
   (calotrack.EndMom).SetXYZ(0, 0, 0);
   }

   return calotrack; 
   }
   */

TVector3 ClusterCoG(Cluster * inputCluster)
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


LCCollection* ClusterVecColl( std::vector<Cluster*> inputClusters )
{
	LCCollection* vec_coll_Clusters = new LCCollectionVec(LCIO::CLUSTER);
	LCFlagImpl flag;
	flag.setBit(LCIO::CHBIT_LONG);
	vec_coll_Clusters->setFlag(flag.getFlag());

	int NClu = inputClusters.size();
	int CurrBranchSize = 0; 
	TVector3 SeedPos; 
	std::vector<float> CluEn; 
	std::vector<int> CluIndex; 

	for(int i0 = 0; i0 < NClu; i0++)
	{
		Cluster* a_clu = inputClusters[i0];
		CluEn.push_back(a_clu->getEnergy());
	}
	CluIndex = SortMeasure(CluEn, 1);

	for(int i1 = 0; i1 < NClu; i1++)
	{
		ClusterImpl * branchtmp = new ClusterImpl();
		Cluster* a_clu = inputClusters[CluIndex[i1]];

		CurrBranchSize = a_clu->getCalorimeterHits().size();

		for(int j1 = 0; j1 < CurrBranchSize; j1 ++)
		{
			CalorimeterHit * tmpHit = a_clu->getCalorimeterHits()[j1];
			branchtmp->addHit(tmpHit, float(1.0));
		}

		branchtmp->setPosition( (float*)a_clu->getPosition() );
		branchtmp->setEnergy(a_clu->getEnergy() );
		SeedPos = a_clu->getPosition();
		branchtmp->setITheta( SeedPos.Theta() );       //To be replaced, those worse than 1st order appro
		branchtmp->setIPhi( SeedPos.Phi()  );

		vec_coll_Clusters->addElement(branchtmp);
	}

	return vec_coll_Clusters;
}

std::vector<Cluster*> CollClusterVec(LCCollection * input_coll )
{
	std::vector<Cluster*> outputClusterVec; 

	bool isCluster = (input_coll->getTypeName() == LCIO::CLUSTER);

	if(! isCluster) 
	{	
		cout<<"Not Cluster Input Collection, Drop"<<endl;
		exit(5);
	}

	outputClusterVec.clear();

	for(int i = 0; i < input_coll->getNumberOfElements(); i++)	//We can have some sort here - according to depth/energy...
	{
		Cluster * a_clu = dynamic_cast<Cluster*>(input_coll->getElementAt(i));
		outputClusterVec.push_back(a_clu);
	}

	return outputClusterVec; 
}

std::vector<CalorimeterHit*> CollHitVec(LCCollection * input_coll, float EnergyThreshold)
{
	std::vector<CalorimeterHit*> outputHitVec;

	bool isCaloHit = (input_coll->getTypeName() == LCIO::CALORIMETERHIT);

	if(! isCaloHit)
	{
		cout<<"Not Cluster Input Collection, Drop"<<endl;
		exit(5);
	}

	outputHitVec.clear();

	for(int i = 0; i < input_coll->getNumberOfElements(); i++)      //We can have some sort here - according to depth/energy...
	{
		CalorimeterHit * a_hit = dynamic_cast<CalorimeterHit*>(input_coll->getElementAt(i));
		if(a_hit->getEnergy() > EnergyThreshold)
		{
			outputHitVec.push_back(a_hit);
		}
	}

	return outputHitVec;

}

std::vector<Cluster*> ClusterHitAbsorbtion( std::vector<Cluster*> MainClusters, std::vector<CalorimeterHit*> IsoHits, float DisThreshold )	// Projective Distance + Hit Depth correlation; 
{
	std::vector<Cluster*> outputClusterVec;

	int N_Core = MainClusters.size();
	int N_Hit = IsoHits.size();
	TVector3 HitPos, MBSeedPos;
	float currHitCoreDis = 0;  
	float MinHitCoreDis = 1.0E10; 
	int MinDisIndex = -1; 
	std::vector<std::pair<int, int> > Frag_Core_Links;
	std::pair<int, int> a_frag_core_link;
	Cluster* a_core; 

	for(int i0 = 0; i0 < N_Hit; i0++)
	{
		CalorimeterHit * a_hit = IsoHits[i0];
		HitPos = a_hit->getPosition();		
		MinHitCoreDis = 1.0E10;

		for(int j0 = 0; j0 < N_Core; j0++)
		{
			a_core = MainClusters[j0];
			currHitCoreDis = DisPointToBush(HitPos, a_core);
			if(currHitCoreDis < MinHitCoreDis)
			{
				MinHitCoreDis = currHitCoreDis; 
				MinDisIndex = j0;
			}
		}
		if(MinHitCoreDis < DisThreshold)
		{
			a_frag_core_link.first = i0;
			a_frag_core_link.second = MinDisIndex;
			Frag_Core_Links.push_back(a_frag_core_link);
		}
	}

	int N_frag_core_links = Frag_Core_Links.size();
	std::vector<CalorimeterHit*> tomerge_hits;
	float ClusterEn = 0; 

	for(int i2 = 0; i2 < N_Core; i2 ++)
	{
		a_core = MainClusters[i2];
		tomerge_hits.clear();

		for(int j4 = 0; j4 < N_frag_core_links; j4 ++)
		{
			a_frag_core_link = Frag_Core_Links[j4];
			if(a_frag_core_link.second == i2)
			{
				CalorimeterHit * a_frag = IsoHits[a_frag_core_link.first];
				tomerge_hits.push_back(a_frag);
			}
		}
		ClusterImpl *a_mergedfrag_core = new ClusterImpl();
		ClusterEn = 0; 

		for(unsigned int j2 = 0; j2 < a_core->getCalorimeterHits().size(); j2++)
		{
			CalorimeterHit * b_hit = a_core->getCalorimeterHits()[j2];
			a_mergedfrag_core->addHit(b_hit, float(1.0));
			ClusterEn += b_hit->getEnergy();
		}

		for(unsigned int j3 = 0; j3 < tomerge_hits.size(); j3++)
		{
			CalorimeterHit * c_hit = tomerge_hits[j3];
			a_mergedfrag_core->addHit(c_hit, float(1.0));
			ClusterEn += c_hit->getEnergy();
		}

		a_mergedfrag_core->setPosition( (float*)a_core->getPosition() );
		a_mergedfrag_core->setEnergy(ClusterEn);
		MBSeedPos = a_core->getPosition();
		a_mergedfrag_core->setITheta( MBSeedPos.Theta() );       //To be replaced, those worse than 1st order appro
		a_mergedfrag_core->setIPhi( MBSeedPos.Phi()  );

		outputClusterVec.push_back(a_mergedfrag_core);
	}

	return outputClusterVec; 
}


std::vector<Cluster*> ClusterAbsorbtion( std::vector<Cluster*> MainClusters, std::vector<Cluster*> FragClusters, float DisThreshold, float DepthSlope )	//ProjectiveDis
{
	/*
	   LCCollection * mergedcluster = new LCCollectionVec(LCIO::CLUSTER);

	   LCFlagImpl flag; 
	   flag.setBit(LCIO::CHBIT_LONG);
	   mergedcluster->setFlag(flag.getFlag());
	   */

	std::vector<Cluster*> outputClusterVec;

	int N_Core = MainClusters.size();
	int N_frag = FragClusters.size();

	//tag minimal distance

	Cluster* a_frag, *a_core; 
	float MinFragCoreDis = 1.0E10; 
	float CurrFragCoreDis = 0;
	int MinDisIndex = -1;  
	std::vector<std::pair<int, int> > Frag_Core_Links;
	std::pair<int, int> a_frag_core_link;
	std::map<int, int> TouchedFrag;
	TouchedFrag.clear();
	TVector3 fragPos; 

	for(int i0 = 0; i0 < N_frag; i0 ++)
	{
		a_frag = FragClusters[i0];
		fragPos = a_frag->getPosition();
		MinFragCoreDis = 1.0E10;
		for(int j0 = 0; j0 < N_Core; j0++)
		{
			a_core = MainClusters[j0];
			CurrFragCoreDis = BushDis(a_frag, a_core);
			if(CurrFragCoreDis < MinFragCoreDis)
			{
				MinFragCoreDis = CurrFragCoreDis;
				MinDisIndex = j0;
			}
		}

		// cout<<MinFragCoreDis<<endl; 

		if( MinFragCoreDis < DisThreshold + DepthSlope*DisSeedSurface(fragPos))
		{
			a_frag_core_link.first = i0;
			a_frag_core_link.second = MinDisIndex;
			Frag_Core_Links.push_back(a_frag_core_link);
		}
	}

	int N_frag_core_links = Frag_Core_Links.size();
	std::vector<Cluster*> tomerge_clu;

	for(int i4 = 0; i4 < N_Core; i4 ++)
	{
		a_core = MainClusters[i4];
		tomerge_clu.clear();
		tomerge_clu.push_back(a_core);

		for(int j4 = 0; j4 < N_frag_core_links; j4 ++)
		{
			a_frag_core_link = Frag_Core_Links[j4];
			if(a_frag_core_link.second == i4)
			{
				a_frag = FragClusters[a_frag_core_link.first];
				TouchedFrag[a_frag_core_link.first] = 1;
				tomerge_clu.push_back(a_frag);
			}
		}
		ClusterImpl *a_mergedfrag_core = NaiveMergeClu(tomerge_clu);
		// mergedcluster->addElement(a_mergedfrag_core);
		outputClusterVec.push_back(a_mergedfrag_core);
	}

	for(int i1 = 0; i1 < N_frag; i1++)
	{
		if(TouchedFrag.find(i1) == TouchedFrag.end())
		{
			a_frag = FragClusters[i1];
			ClusterImpl *a_isofrag = NaiveCluImpl(a_frag);
			// mergedcluster->addElement(a_isofrag);
			outputClusterVec.push_back(a_isofrag);
		}
	}

	// return mergedcluster; 
	return outputClusterVec;
}

LCCollection* ClusterVecMerge( std::vector<Cluster*> inputClusters, TMatrixF ConnectorMatrix  )
{
	LCCollection* mergedbranches = new LCCollectionVec(LCIO::CLUSTER);

	LCFlagImpl flag;
	flag.setBit(LCIO::CHBIT_LONG);
	mergedbranches->setFlag(flag.getFlag());

	int NinputClu = inputClusters.size();
	int Nrow = ConnectorMatrix.GetNrows();
	int Ncol = ConnectorMatrix.GetNcols();

	if(Ncol != NinputClu || Nrow != Ncol || Nrow != NinputClu)
	{
		cout<<"Size of Connector Matrix and inputClusterColl is not match"<<endl;
	}

	vector<Cluster*> branchToMerge;
	Cluster* Mergebranch_A;
	Cluster* Mergebranch_B;
	Cluster* tmpMergebranch;
	Cluster* Mainbranch (0);

	TVector3 tmpClusterSeedPos, MBSeedPos;	

	int CurrBranchSize = 0;
	float SeedPosMin = 1.0E10;
	float BranchEnergy = 0;

	int FlagBranchTouch[Nrow];

	for(int i0 = 0; i0 < Nrow; i0++)
	{
		FlagBranchTouch[i0] = 0;
	}

	for(int ibran = 0; ibran < Nrow ; ibran++)
	{
		if(FlagBranchTouch[ibran] == 0)
		{
			Mergebranch_A = inputClusters[ibran];
			branchToMerge.push_back(Mergebranch_A);
			FlagBranchTouch[ibran] = 1;
			BranchEnergy = 0;
			ClusterImpl * branchtmp = new ClusterImpl();

			for(int jbran = ibran + 1; jbran < Nrow; jbran++)
			{
				if(FlagBranchTouch[jbran] == 0)
				{
					Mergebranch_B = inputClusters[jbran];
					if( ConnectorMatrix(ibran, jbran) > 0.1 )
					{
						branchToMerge.push_back(Mergebranch_B);
						FlagBranchTouch[jbran] = 1;
					}
				}
			}

			SeedPosMin = 1.0E10;

			for(unsigned int i1 = 0; i1 < branchToMerge.size(); i1++)
			{
				tmpMergebranch = branchToMerge[i1];
				tmpClusterSeedPos = tmpMergebranch->getPosition();
				if( tmpClusterSeedPos.Mag() < SeedPosMin)
				{

					Mainbranch = tmpMergebranch;
					SeedPosMin = tmpClusterSeedPos.Mag();
				}

				CurrBranchSize = tmpMergebranch->getCalorimeterHits().size();
				BranchEnergy += tmpMergebranch->getEnergy();

				for(int j1 = 0; j1 < CurrBranchSize; j1 ++)
				{
					CalorimeterHit * tmpHit = tmpMergebranch->getCalorimeterHits()[j1];
					branchtmp->addHit(tmpHit, float(1.0));
				}
				/*
				   int NSubClu = tmpMergebranch->getClusters().size();
				   if(NSubClu == 0) branchtmp->addCluster(tmpMergebranch);
				   else
				   {
				   for(int aa = 0;aa < NSubClu; aa++)
				   {
				   branchtmp->addCluster(tmpMergebranch->getClusters()[aa]);
				   }
				   }
				   */
			}
			branchtmp->setPosition( (float*)Mainbranch->getPosition() );
			branchtmp->setEnergy(BranchEnergy);
			MBSeedPos = Mainbranch->getPosition();
			branchtmp->setITheta( MBSeedPos.Theta() );       //To be replaced, those worse than 1st order appro
			branchtmp->setIPhi( MBSeedPos.Phi()  );

			mergedbranches->addElement(branchtmp);
			branchToMerge.clear();
		}
	}

	return mergedbranches;

}

std::pair<TVector3, TVector3> ClosestPointPair(Cluster *a_clu, Cluster *b_clu)
{
	std::pair<TVector3, TVector3> a_hit_pair; 
	TVector3 Pos_Zero(0, 0, 0);
	a_hit_pair.first = Pos_Zero; 
	a_hit_pair.second = Pos_Zero; 

	float MinDis = 1.0E10; 	
	TVector3 Pos_A, Pos_B; 
	CalorimeterHit * a_hit, * b_hit; 	
	float CurrHitDis = 1.0E10;

	for(unsigned int i = 0; i < a_clu->getCalorimeterHits().size(); i++)
	{
		a_hit = a_clu->getCalorimeterHits()[i];
		Pos_A = a_hit->getPosition();
		for(unsigned int j = 0; j < b_clu->getCalorimeterHits().size(); j++)
		{
			b_hit = b_clu->getCalorimeterHits()[j];
			Pos_B = b_hit->getPosition();
			CurrHitDis = (Pos_B - Pos_A).Mag();			
			if(CurrHitDis < MinDis)
			{
				MinDis = CurrHitDis;
				a_hit_pair.first = Pos_A; 
				a_hit_pair.second = Pos_B; 
			}
		}
	}

	return a_hit_pair; 
}

LCCollection* ClusterMerge( LCCollection * inputClusterColl, TMatrixF ConnectorMatrix )
{
	LCCollection* mergedbranches = new LCCollectionVec(LCIO::CLUSTER);

	LCFlagImpl flag;
	flag.setBit(LCIO::CHBIT_LONG);
	// flag.setBit(LCIO::);
	mergedbranches->setFlag(flag.getFlag());

	int NinputClu = inputClusterColl->getNumberOfElements();
	int Nrow = ConnectorMatrix.GetNrows();
	int Ncol = ConnectorMatrix.GetNcols();

	if(Ncol != NinputClu || Nrow != Ncol || Nrow != NinputClu)
	{
		cout<<"Size of Connector Matrix and inputClusterColl is not match"<<endl;
	}

	vector<Cluster*> branchToMerge;
	Cluster* Mergebranch_A;
	Cluster* Mergebranch_B;
	Cluster* tmpMergebranch;
	Cluster* Mainbranch (0); 

	TVector3 tmpClusterSeedPos, MBSeedPos;

	//int MainBranchID = 0;
	int CurrBranchSize = 0;
	float SeedPosMin = 1.0E10; 
	float BranchEnergy = 0; 

	int FlagBranchTouch[Nrow];

	for(int i0 = 0; i0 < Nrow; i0++)
	{
		FlagBranchTouch[i0] = 0;
	}

	for(int ibran = 0; ibran < Nrow ; ibran++)
	{
		if(FlagBranchTouch[ibran] == 0)
		{
			Mergebranch_A = dynamic_cast<Cluster*>( inputClusterColl->getElementAt(ibran));
			branchToMerge.push_back(Mergebranch_A);
			FlagBranchTouch[ibran] = 1;
			BranchEnergy = 0;
			ClusterImpl * branchtmp = new ClusterImpl();

			for(int jbran = ibran + 1; jbran < Nrow; jbran++)
			{
				if(FlagBranchTouch[jbran] == 0)
				{
					Mergebranch_B = dynamic_cast<Cluster*>( inputClusterColl->getElementAt(jbran));
					if( ConnectorMatrix(ibran, jbran) > 0.1 )
					{
						branchToMerge.push_back(Mergebranch_B);
						FlagBranchTouch[jbran] = 1;
					}
				}
			}

			SeedPosMin = 1.0E10;

			//cout<<"Merge Size "<<branchToMerge.size()<<endl;

			float ECALTotalEn = 0;
			float HCALTotalEn = 0;

			for(unsigned int i1 = 0; i1 < branchToMerge.size(); i1++)
			{
				tmpMergebranch = branchToMerge[i1];
				tmpClusterSeedPos = tmpMergebranch->getPosition();
				if( tmpClusterSeedPos.Mag() < SeedPosMin)
				{
					//MainBranchID = i1;
					Mainbranch = tmpMergebranch;
					SeedPosMin = tmpClusterSeedPos.Mag();
				}

				CurrBranchSize = tmpMergebranch->getCalorimeterHits().size();
				BranchEnergy += tmpMergebranch->getEnergy();

				for(int j1 = 0; j1 < CurrBranchSize; j1 ++)
				{
					CalorimeterHit * tmpHit = tmpMergebranch->getCalorimeterHits()[j1];
					branchtmp->addHit(tmpHit, float(1.0));
					if(fabs(tmpHit->getEnergy() - DHCALCalibrationConstant) < 1.0E-6 )
					{
						HCALTotalEn += DHCALCalibrationConstant;
					}
					else
					{
						ECALTotalEn += tmpHit->getEnergy();
					}
				}
				/*
				   int NSubClu = tmpMergebranch->getClusters().size();
				   if(NSubClu == 0) branchtmp->addCluster(tmpMergebranch);
				   else
				   {
				   for(int aa = 0;aa < NSubClu; aa++)
				   {
				   branchtmp->addCluster(tmpMergebranch->getClusters()[aa]);
				   }
				   }
				   */
			}
			branchtmp->setPosition( (float*)Mainbranch->getPosition() );
			branchtmp->setEnergy(BranchEnergy);
			MBSeedPos = Mainbranch->getPosition();
			branchtmp->setITheta( MBSeedPos.Theta() );       //To be replaced, those worse than 1st order appro
			branchtmp->setIPhi( MBSeedPos.Phi()  );

			branchtmp->subdetectorEnergies().resize(6) ;
			branchtmp->subdetectorEnergies()[0] = ECALTotalEn ;
			branchtmp->subdetectorEnergies()[1] = HCALTotalEn ;
			branchtmp->subdetectorEnergies()[2] = 0 ;
			branchtmp->subdetectorEnergies()[3] = 0 ;
			branchtmp->subdetectorEnergies()[4] = 0 ;
			branchtmp->subdetectorEnergies()[5] = 0 ;

			mergedbranches->addElement(branchtmp);
			branchToMerge.clear();
		}
	}

	// cout<<"MBSize "<<mergedbranches->getNumberOfElements()<<endl;

	return mergedbranches;
}

/*
   std::vector<TVector3> BranchRefDir( Cluster* inputBranch )
   {
   std::vector<TVector3> RefDir; 
   TVector3 cluEndPos, cluBeginPos, cluDir_B, cluDir_E;

   int CurrPriBranchSize = inputBranch->getCalorimeterHits().size();

   double cone_crit = 20.; // Criteria for cone segmentation (in degrees)
   double chi2_crit = 50.;

   fit_type ft = HEL; // HEL (helix) sometimes is unstable for few points. Better use CIRC or SVD
   point axis = point(0,0,1);

   trajectory t(ft, axis, cone_crit, chi2_crit);
   vector <point> points;
   for(int j0 = 0; j0 < CurrPriBranchSize; j0++) // Now we add hits to t one by one
   {
   CalorimeterHit *a_hit = inputBranch->getCalorimeterHits()[j0];
   points.push_back( point( (double) a_hit->getPosition()[0],
   (double) a_hit->getPosition()[1],
   (double) a_hit->getPosition()[2]) );

   if(j0 == 0)
   cluEndPos = a_hit->getPosition();
   if(j0 == CurrPriBranchSize - 1)
   cluBeginPos = a_hit->getPosition();
   }

   unsigned inv;

   for(int j1 = 0; j1 < CurrPriBranchSize; j1++)
   {
   inv = points.size() - j1 - 1;
   t.add_point( points[ inv ] );
   }

   point RefD = t.dir_begin();
   point EndD = t.dir_end();

//Protection aganist Dir = (0, 0, +-1)

if( CurrPriBranchSize < 5 || abs(RefD.getZ()) > 0.999 || abs(EndD.getZ()) > 0.999)
{
if(CurrPriBranchSize > 5)
{
cout<<"Track Fit Anomally Tagged"<<endl; 
}

TVector3 Diff = cluEndPos - cluBeginPos;
Diff = 1.0/Diff.Mag()*Diff; 
cluDir_B.SetXYZ(Diff.X(), Diff.Y(), Diff.Z());
cluDir_E.SetXYZ(Diff.X(), Diff.Y(), Diff.Z());
}
else
{
cluDir_B.SetXYZ(RefD.getX(), RefD.getY(), RefD.getZ());
cluDir_E.SetXYZ(EndD.getX(), EndD.getY(), EndD.getZ());
}

cluDir_B = 1.0/cluDir_B.Mag() * cluDir_B;
cluDir_E = 1.0/cluDir_E.Mag() * cluDir_E;

RefDir.push_back(cluBeginPos);
RefDir.push_back(cluEndPos);
RefDir.push_back(cluDir_B);
RefDir.push_back(cluDir_E);

return RefDir; 
}
*/

ClusterImpl* NaiveMergeClu(std::vector<Cluster*> inputCluVec)
{
	ClusterImpl* MergedClu = new ClusterImpl();

	int NClu = inputCluVec.size();
	int CurrCluSize = 0; 
	float SeedDis = 1E9; 
	float MaxDis = 0; 
	float MergedCluEnergy = 0; 
	float HitEn = 0; 
	float SubDEn[6] = {0, 0, 0, 0, 0, 0};

	TVector3 CurrSeedPos, SeedPos, CurrHitPos, CluEndPos, CluRefDir;	//Seed Depth... CoG Comp...

	for(int i = 0; i < NClu; i++)
	{
		Cluster* a_Clu = inputCluVec[i];
		CurrSeedPos = a_Clu->getPosition();
		MergedCluEnergy += a_Clu->getEnergy();

		if(CurrSeedPos.Mag() < SeedDis)
		{
			SeedPos = CurrSeedPos; 
			SeedDis = CurrSeedPos.Mag();
		}

		CurrCluSize = a_Clu->getCalorimeterHits().size();
		/*	
			int NSubClu = a_Clu->getClusters().size();
			if(NSubClu == 0) MergedClu->addCluster(a_Clu);
			else
			{
			for(int aa = 0;aa < NSubClu; aa++)
			{
			MergedClu->addCluster(a_Clu->getClusters()[aa]);
			}
			}
			*/

		for(int j = 0; j < CurrCluSize; j++)
		{
			CalorimeterHit* a_hit = a_Clu->getCalorimeterHits()[j];
			MergedClu -> addHit(a_hit, 1.0);
			CurrHitPos = a_hit->getPosition();
			HitEn = a_hit->getEnergy();
			if(fabs(HitEn - DHCALCalibrationConstant) < 1.0E-6)	// ECAL, HCAL, Should use better criteria. 
			{
				SubDEn[1] += HitEn; 
			}
			else
			{
				SubDEn[0] += HitEn; 
			}

			if(CurrHitPos.Mag() > MaxDis)
			{
				MaxDis = CurrHitPos.Mag();
				CluEndPos = a_hit->getPosition();	
			}
		}
	}
	CluRefDir = (CluEndPos - SeedPos);

	float ClusterSeedPos[3] = {float(SeedPos.X()), float(SeedPos.Y()), float(SeedPos.Z())};
	MergedClu->setPosition(ClusterSeedPos);
	MergedClu->setEnergy(MergedCluEnergy);
	MergedClu->setITheta( CluRefDir.Theta() );
	MergedClu->setIPhi( CluRefDir.Phi() );

	MergedClu->subdetectorEnergies().resize(6) ;
	for(int i = 0; i < 6; i++)
	{
		MergedClu->subdetectorEnergies()[i] = SubDEn[i];
	}

	return MergedClu;
}

int JointsBetweenBush(Cluster* a_Clu, Cluster* b_Clu, float CellSize)
{
	int NJoint = 0; 
	int a_CluSize = a_Clu->getCalorimeterHits().size();
	int b_CluSize = b_Clu->getCalorimeterHits().size();
	TVector3 aHitPos, bHitPos, PosDiff, aCluPos, bCluPos; 	
	aCluPos = a_Clu->getPosition();
	bCluPos = b_Clu->getPosition();

	for(int i = 0; i < a_CluSize; i++)
	{
		CalorimeterHit * ahit = a_Clu->getCalorimeterHits()[i];
		aHitPos = ahit->getPosition();
		for(int j = 0; j < b_CluSize; j++)
		{
			CalorimeterHit * bhit = b_Clu->getCalorimeterHits()[j];
			bHitPos = bhit->getPosition();
			PosDiff = aHitPos - bHitPos; 
			if(PosDiff.Mag() < 1.5*CellSize)	//allow Diag connect... else use 1.2
			{
				// if((aCluPos - bHitPos).Mag() < 60 || (bCluPos - aHitPos).Mag() < 60)
				NJoint++;	//Change to NJoint Hit...
			}
		}
	}

	return NJoint; 
}


TVector3 TPCHitPos(MCParticle *inputMCP)
{

	TVector3 ExpectedHitPos, TmpHitPos1, TmpHitPos2;
	//Suppose all come from I.P with speed of light
	TVector3 MCPMomentum = inputMCP->getMomentum();
	//float theta = MCPMomentum.Theta();
	float TransMom = MCPMomentum.Pt();
	float phi = MCPMomentum.Phi();
	int charge = int(inputMCP->getCharge());

	float Radius = 0; //mm
	float eta = 0;
	float HitPosX = 0;
	float HitPosY = 0;
	float HitPosZ = 0;
	float HelixPhase = 0;

	if(charge == 0)
	{
		TmpHitPos1 = MCPMomentum*(ECALHalfZ/fabs(MCPMomentum[2]));
		TmpHitPos2 = MCPMomentum*(TPCRadius/TransMom);
	}
	else
	{
		Radius = MCPMomentum.Pt()*3333.33/4;   //mm C*B, assume unit charge

		eta = phi - charge*asin(TPCRadius/(2*Radius));
		TmpHitPos1.SetXYZ(cos(eta)*TPCRadius, sin(eta)*TPCRadius, MCPMomentum[2]/MCPMomentum.Pt()*TPCRadius);

		if(MCPMomentum[2] > 0) HitPosZ = ECALHalfZ;
		else HitPosZ = -1*ECALHalfZ;

		HelixPhase = charge*(MCPMomentum.Pt()*HitPosZ)/(MCPMomentum[2]*Radius);
		HitPosX = charge*Radius*(sin(phi) + sin(HelixPhase - phi ));
		HitPosY = charge*Radius*(cos(HelixPhase - phi ) - cos(phi));
		TmpHitPos2.SetXYZ(HitPosX, HitPosY, HitPosZ);
	}

	if(TmpHitPos1.Mag() < TmpHitPos2.Mag() )
	{
		ExpectedHitPos = TmpHitPos1;
	}
	else
	{
		ExpectedHitPos = TmpHitPos2;
	}

	return ExpectedHitPos;
}


TVector3 ECALHitPos(MCParticle *a_MCP, TVector3 &ExpHitDir)	//Problematic??
{
	TVector3 ExpHitPos, MCPVtxPos, MCPMomentum; 
	TVector3 NeutralPExpHitPos;

	MCPVtxPos = a_MCP->getVertex();
	MCPMomentum = a_MCP->getMomentum();
	float charge = a_MCP->getCharge();

	float Phi_v = MCPVtxPos.Phi();
	float Phi_P = MCPMomentum.Perp(); 	

	float v = MCPVtxPos.Perp(); 
	float r = MCPMomentum.Perp()*952.4;   //mm    10/3*3.5;
	float R = 2000.0;	//mm     out counter circle radius of ECAL
	float DisAV = 0;
	float tau1 = 0;		//Length parameter of the helix...
	float tau2 = 0; 
	float Zshift1 = 0;
	float Zshift2 = 0;
	float AxisPosX = 0; 
	float AxisPosY = 0;

	float BarrelX = 0;
	float BarrelY = 0;
	float BarrelZ = 0;
	float EndCapX = 0;
	float EndCapY = 0;
	float EndCapZ = 0;

	const double pi = acos(-1.0);

	if( fabs(charge) < 0.01 )
	{
		float alpha_N = Phi_P + asin(v/R*sin(Phi_P - Phi_v));

		BarrelX = R*cos(alpha_N);
		BarrelY = R*sin(alpha_N);
		float RefRadiusLength = sqrt((BarrelX-MCPVtxPos.X())*(BarrelX-MCPVtxPos.X()) + (BarrelY -MCPVtxPos.Y())*(BarrelY -MCPVtxPos.Y()));
		BarrelZ = MCPVtxPos.Z() + MCPMomentum.Z()/MCPMomentum.Perp()*RefRadiusLength;
		if(fabs(BarrelZ) < ECALHalfZ)	//EcalBarrelLength
		{
			ExpHitPos.SetXYZ(BarrelX, BarrelY, BarrelZ);
		}
		else
		{
			EndCapZ = ECALHalfZ*fabs(BarrelZ)/BarrelZ; //Sgn...
			EndCapX = MCPVtxPos.X() + (BarrelX - MCPVtxPos.X())*(EndCapZ - MCPVtxPos.Z())/(BarrelZ - MCPVtxPos.Z());
			EndCapY = MCPVtxPos.Y() + (BarrelY - MCPVtxPos.Y())*(EndCapZ - MCPVtxPos.Z())/(BarrelZ - MCPVtxPos.Z());
			ExpHitPos.SetXYZ(EndCapX, EndCapY, EndCapZ);
		}
		ExpHitDir = MCPMomentum;
	}
	else if( fabs(charge) > 0.9 )
	{
		//Charged Case: alpha: Hit Position Azimuth Angle; satisfy C*Cos(Phi_A) - B*Sin(Phi_A) = A

		AxisPosX = MCPVtxPos.X() + r*charge*sin(Phi_P);
		AxisPosY = MCPVtxPos.Y() - r*charge*cos(Phi_P);

		if( 2*r + v - R > 0 )
		{
			float Coff_A = v*v + R*R + 2*r*charge*v*sin(Phi_P - Phi_v);
			float Coff_C = 2*R*(sin(Phi_P)*r*charge + cos(Phi_v)*v);
			float Coff_BCM = 2*R*sqrt(r*r + v*v + 2*r*charge*v*sin(Phi_P - Phi_v));

			float beta = asin(Coff_C/Coff_BCM);
			float alpha = beta - asin(Coff_A/Coff_BCM);	
			BarrelX = R*cos(alpha);		//Charge already into account.
			BarrelY = R*sin(alpha);
			DisAV = sqrt((BarrelX - MCPVtxPos.X())*(BarrelX - MCPVtxPos.X()) + (BarrelY - MCPVtxPos.Y())*(BarrelY - MCPVtxPos.Y()) );
			tau1 = 2*atan(0.5*DisAV/r);	
			Zshift1 = tau1*r*MCPMomentum.Z()/MCPMomentum.Perp();
			BarrelZ = MCPVtxPos.Z() + Zshift1; 

			if(fabs(BarrelZ) < ECALHalfZ)
			{
				ExpHitPos.SetXYZ(BarrelX, BarrelY, BarrelZ);
				tau2 = tau1; 
			}
			else
			{
				EndCapZ =ECALHalfZ*fabs(BarrelZ)/BarrelZ; 
				Zshift2 = EndCapZ - MCPVtxPos.Z();	
				tau2 = tau1*Zshift2/Zshift1; 
				EndCapX = AxisPosX - r*sin(Phi_P - tau2);
				EndCapY = AxisPosY + r*cos(Phi_P - tau2);
				ExpHitPos.SetXYZ(EndCapX, EndCapY, EndCapZ);
			}
			ExpHitDir.SetXYZ(cos(Phi_P - tau2), sin(Phi_P - tau2), MCPMomentum.Z()/MCPMomentum.Perp());
		}
		else
		{
			EndCapZ = ECALHalfZ*MCPMomentum.Z()/fabs(MCPMomentum.Z());
			Zshift2 = EndCapZ - MCPVtxPos.Z();
			tau2 = Zshift2*MCPMomentum.Perp()/(MCPMomentum.Z()*2*pi*r);
			EndCapX = AxisPosX - r*sin(Phi_P - tau2);
			EndCapY = AxisPosY + r*cos(Phi_P - tau2);
			ExpHitPos.SetXYZ(EndCapX, EndCapY, EndCapZ);
			ExpHitDir.SetXYZ(cos(Phi_P - tau2), sin(Phi_P - tau2), MCPMomentum.Z()/MCPMomentum.Perp());
		}
	}

	ExpHitDir = 1.0/ExpHitDir.Mag()*ExpHitDir;

	//cout<<"HitExpPos "<<ExpHitPos.X()<<", "<<ExpHitPos.Y()<<",  "<<ExpHitPos.Z()<<", Radius "<<ExpHitPos.Perp()<<endl;

	return ExpHitPos; 
}

float DisHelixCluster(Cluster* a_clu, HelixClass* a_helix, float &Time)
{
	float BushDist[3] = {0, 0, 0}; 
	float MinDis = 1E9;
	int Nhits = a_clu->getCalorimeterHits().size();
	
	for(int s = 0; s < Nhits; s++)
	{
		CalorimeterHit * a_hit = a_clu->getCalorimeterHits()[s];
		Time = a_helix->getDistanceToPoint((float*)a_hit->getPosition(), BushDist);		
		if(Time > -100 && BushDist[2] < MinDis)
		{
			MinDis = BushDist[2];
		}

	}

	return MinDis; 
}


ClusterImpl* NaiveCluImpl(Cluster* a0_clu)
{
	ClusterImpl * b0_clu = new ClusterImpl();
	b0_clu->setPosition((float*)a0_clu->getPosition());
	b0_clu->setEnergy(a0_clu->getEnergy());
	int NCaloHit = a0_clu->getCalorimeterHits().size();
	float HitEn = 0; 
	float SubDEn[6] = {0, 0, 0, 0, 0, 0};

	for(int t0 = 0; t0 < NCaloHit; t0++)
	{
		CalorimeterHit * a0_hit = a0_clu->getCalorimeterHits()[t0];
		b0_clu->addHit(a0_hit, 1.0);
		HitEn = a0_hit->getEnergy();
		if(fabs(HitEn - DHCALCalibrationConstant) < 1.0E-6)
		{
			SubDEn[1] += HitEn;
		}
		else
		{	
			SubDEn[0] += HitEn;
		}
	}

	b0_clu->subdetectorEnergies().resize(6) ;
	for(int i = 0; i < 6; i++)
	{
		b0_clu->subdetectorEnergies()[i] = SubDEn[i];
	}
	/*
	   int NSubClu = a0_clu->getClusters().size();
	   if(NSubClu == 0) b0_clu->addCluster(a0_clu);
	   else
	   {
	   for(int aa = 0;aa < NSubClu; aa++)
	   {
	   b0_clu->addCluster(a0_clu->getClusters()[aa]);
	   }
	   }
	   */
	return b0_clu; 
}

TVector3 CluEndP(Cluster* a_clu)	//add algorithm to calculate the MIP end
{
	TVector3 PosEndP, CurrHitP; 
	float MaxDis = 0;
	int NCaloHit = a_clu->getCalorimeterHits().size();

	for(int i = 0; i < NCaloHit; i++)
	{
		CalorimeterHit * a_hit = a_clu->getCalorimeterHits()[i];
		CurrHitP = a_hit->getPosition();
		if(CurrHitP.Mag() > MaxDis)
		{
			MaxDis = CurrHitP.Mag();
			PosEndP = CurrHitP;
		}
	}

	return PosEndP;
}

std::vector<TVector3> CluRefDir(Cluster* a_clu)	// a general analysis code as well...
{
	std::vector<TVector3> RefDirections; 

	TVector3 Pos, EndP, tmpDirection, MIPexpoPoint, HitPos;

	Pos = a_clu->getPosition();
	EndP = CluEndP(a_clu);
	tmpDirection = EndP - Pos; 
	tmpDirection = 1.0/tmpDirection.Mag()*tmpDirection;  
	RefDirections.push_back(tmpDirection);

	int clusize = a_clu->getCalorimeterHits().size();
	int NH[50];
	float XMean[50];
	float YMean[50];
	float ZMean[50];
	int LayerNum = 0; 
	int LastHitLayer = 0;
	int CellID0 = 0; 	//direct decoding!!!!!!!!!

	for(int i0 = 0; i0 < 50; i0++)
	{
		NH[i0] = 0;
		XMean[i0] = 0;
		YMean[i0] = 0;
		ZMean[i0] = 0;
	}

	for(int i1 = 0; i1 < clusize; i1++)	//cannot be used for E-H combined. or, define position.
	{
		CalorimeterHit *a_hit = a_clu->getCalorimeterHits()[i1];
		CellID0 = a_hit->getCellID0();
		LayerNum = (CellID0 & 0x3F000000)>>24;	
		HitPos = a_hit->getPosition();
		if(LayerNum < 50 && LayerNum > -1)
		{
			NH[LayerNum] ++;	
			XMean[LayerNum] += HitPos.X();
			YMean[LayerNum] += HitPos.Y();
			ZMean[LayerNum] += HitPos.Z();
		}
	}

	for(int i2 = 0; i2 < 50; i2++)
	{
		if( NH[i2] )
		{
			XMean[i2] /= NH[i2];
			YMean[i2] /= NH[i2];
			ZMean[i2] /= NH[i2];
			LastHitLayer = i2; 
		}
	}

	//tag the interaction layer: for MIP-exposition case; 
	int NActiveLayer = 0; 

	for(int i3 =  1; i3 < 50; i3++)
	{
		if(NH[i3 - 1] < 3 && NH[i3] > 2 )	//only valid for mips without multiplicity...	
			NActiveLayer = i3;
		break;
	}

	if(NActiveLayer == 0)
		NActiveLayer = LastHitLayer;

	MIPexpoPoint.SetXYZ(XMean[NActiveLayer], YMean[NActiveLayer], ZMean[NActiveLayer]);

	TVector3 tmpDir = MIPexpoPoint - Pos;
	tmpDir = 1.0/tmpDir.Mag()*tmpDir;

	RefDirections.push_back(tmpDirection);

	return RefDirections;
}

float CluProjectiveDis(TVector3 Pos, TVector3 Dir, Cluster* a_clu)
{
	TVector3 currhitPos;
	Dir = 1.0/Dir.Mag()*Dir; 	//Normalization
	int Size = a_clu->getCalorimeterHits().size();
	float currhitDis = 0;
	float MinimalDis = 1E9; 

	for(int i = 0; i < Size; i++)
	{
		CalorimeterHit * a_hit = a_clu->getCalorimeterHits()[i];
		currhitPos = a_hit->getPosition();
		if( (currhitPos-Pos).Mag() < 150.0 )	//distance threshold...
		{
			currhitDis = ((currhitPos - Pos).Cross(Dir)).Mag();
			if(currhitDis < MinimalDis)
			{
				MinimalDis = currhitDis; 
			}
		}
	}	
	return MinimalDis; 
}

void ClusterCollBuild( LCEvent * evtPP, std::string Name, std::vector<Cluster*> inputClusters )	//Operator Overload
{
	LCCollection *currclustercoll = new LCCollectionVec(LCIO::CLUSTER);
	LCFlagImpl flag;
	flag.setBit(LCIO::CHBIT_LONG);
	currclustercoll->setFlag(flag.getFlag());

	int NCluster = inputClusters.size();

	for(int i = 0; i < NCluster; i++)
	{
		Cluster * a_clu = inputClusters[i];
		ClusterImpl *a_bush = NaiveCluImpl(a_clu);
		currclustercoll->addElement(a_bush);
	}

	evtPP->addCollection( currclustercoll, Name );
}

void ClusterBuilding( LCEvent * evtPP, std::string Name, std::vector<CalorimeterHit*> Hits, std::vector< std::vector<int> > BranchOrder, int DHCALFlag )
{
	LCCollection *currbranchcoll = new LCCollectionVec(LCIO::CLUSTER);
	LCFlagImpl flag;
	flag.setBit(LCIO::CHBIT_LONG);
	currbranchcoll->setFlag(flag.getFlag());

	int NBranch = BranchOrder.size();
	int BranchSize = 0;
	float currBranchEnergy = 0;
	TVector3 SeedPos, currPos;
	float MinMag = 1E9;
	float currMag = 0; 
	float ECALTotalEn = 0; 
	float HCALTotalEn = 0;

	for(int i0 = 0; i0 < NBranch; i0++)
	{
		ClusterImpl* a_branch = new ClusterImpl();
		std::vector<int> currbranchorder = BranchOrder[i0];
		BranchSize = currbranchorder.size();
		currBranchEnergy = 0;
		ECALTotalEn = 0;
		HCALTotalEn = 0;
		// CalorimeterHit *Seedhit = Hits[currbranchorder[BranchSize - 1]];
		MinMag = 1E9;

		for(int j = 0; j < BranchSize; j++)
		{
			CalorimeterHit * a_hit = Hits[currbranchorder[j]];
			currPos = a_hit->getPosition();
			// currMag = DisSeedSurface(currPos);
			currMag = currPos.Mag();

			if( currMag < MinMag )
			{
				MinMag = currMag;
				SeedPos = currPos;
			}

			if(fabs(a_hit->getEnergy() - DHCALCalibrationConstant) < 1.0E-6 )
			{
				HCALTotalEn += DHCALCalibrationConstant;
			}
			else
			{
				ECALTotalEn += a_hit->getEnergy();
			}

			currBranchEnergy +=  a_hit->getEnergy();

			a_branch->addHit(a_hit, (float)1.0);
		}

		a_branch->setEnergy(currBranchEnergy);
		float ArraySeedPos[3] = { float(SeedPos.X()), float(SeedPos.Y()), float(SeedPos.Z()) };
		a_branch->setPosition( ArraySeedPos );

		a_branch->subdetectorEnergies().resize(6) ;
		a_branch->subdetectorEnergies()[0] = ECALTotalEn ;
		a_branch->subdetectorEnergies()[1] = HCALTotalEn ;
		a_branch->subdetectorEnergies()[2] = 0 ;
		a_branch->subdetectorEnergies()[3] = 0 ;
		a_branch->subdetectorEnergies()[4] = 0 ;
		a_branch->subdetectorEnergies()[5] = 0 ;

		currbranchcoll -> addElement(a_branch);
	}

	evtPP->addCollection(currbranchcoll, Name);
}
/*
float *SimpleDisTrackClu(Track * a_trk, Cluster * a_clu)
{
	static float Distance[3] = {1.0E9,1.0E9,1.0E9}; 
	float minDis = 1.0E9;
	float BushDist[3] = {0, 0 ,0};
	HelixClass * a_Helix = new HelixClass();
	//float refPoint[3] = a_Helix->getReferencePoint();
	a_Helix->Initialize_Canonical(a_trk->getPhi(), a_trk -> getD0(), a_trk -> getZ0(), a_trk -> getOmega(), a_trk->getTanLambda(), BField);
	int NCaloHits = a_clu->getCalorimeterHits().size();


	for(int i1 = 0; i1 < NCaloHits; i1++)
	{
		CalorimeterHit *a_hit = a_clu->getCalorimeterHits()[i1];

		float BushTime = a_Helix->getDistanceToPoint((float*)a_hit->getPosition(), BushDist);
		if(BushTime > 0 && BushDist[2] < minDis )
		{
			minDis = BushDist[2];
			Distance[0] = BushDist[0];
			Distance[1] = BushDist[1];
			Distance[2] = BushDist[2];
		}

		//float phi = atan2(a_hit->getPosition()[1]-_yCentre,a_hit->getPosition()[0]-_xCentre);
        	//float phi0 = atan2(_refPoint[1]-_yCentre,_refPoint[0]-_xCentre);


		//if(fabs(a_Helix->getMomentum()[2]) > 0.000001)
		//	cout<<BushTime<<" : "<<BushDist[1]/a_Helix->getMomentum()[2]<<endl; 
	}

	return Distance;
}
*/

float* SimpleDisTrackClu(Track * a_trk, Cluster * a_clu)
{
	float* Distance = new float[3];
       	Distance[0]	= 1.0E9;
       	Distance[1]	= 1.0E9;
       	Distance[2]	= 1.0E9;
	float minDis = 1.0E9;
	float BushDist[3] = {0, 0 ,0};
	HelixClass * a_Helix = new HelixClass();
	//float refPoint[3] = a_Helix->getReferencePoint();
	a_Helix->Initialize_Canonical(a_trk->getPhi(), a_trk -> getD0(), a_trk -> getZ0(), a_trk -> getOmega(), a_trk->getTanLambda(), BField);
	int NCaloHits = a_clu->getCalorimeterHits().size();


	for(int i1 = 0; i1 < NCaloHits; i1++)
	{
		CalorimeterHit *a_hit = a_clu->getCalorimeterHits()[i1];

		float BushTime = a_Helix->getDistanceToPoint((float*)a_hit->getPosition(), BushDist);
		if(BushTime > 0 && BushDist[2] < minDis )
		{
			minDis = BushDist[2];
			Distance[0] = BushDist[0];
			Distance[1] = BushDist[1];
			Distance[2] = BushDist[2];
		}
	}
	delete a_Helix;

	return Distance;
}
float SimpleBushTimeTrackClu(Track * a_trk, Cluster * a_clu)
{
        float Distance = 1.0E9;
        float Time = 0;
        float BushDist[3] = {0, 0 ,0};
        HelixClass * a_Helix = new HelixClass();
        a_Helix->Initialize_Canonical(a_trk->getPhi(), a_trk -> getD0(), a_trk -> getZ0(), a_trk -> getOmega(), a_trk->getTanLambda(), BField);
        int NCaloHits = a_clu->getCalorimeterHits().size();

        for(int i1 = 0; i1 < NCaloHits; i1++)
        {
                CalorimeterHit *a_hit = a_clu->getCalorimeterHits()[i1];

                float BushTime = a_Helix->getDistanceToPoint((float*)a_hit->getPosition(), BushDist);
                if(BushTime > 0 && BushDist[2] < Distance )
                {
                        Time = BushTime;
			Distance = BushDist[2];
                }
        }
	delete a_Helix;
        return Time;
}

int SimpleBushNC(Track * a_trk, Cluster * a_clu)
{
	float Distance = 1.0E9;
	//float Time = 0; 
	int NC = 0;
	float BushDist[3] = {0, 0 ,0};
	HelixClass * a_Helix = new HelixClass();
	a_Helix->Initialize_Canonical(a_trk->getPhi(), a_trk -> getD0(), a_trk -> getZ0(), a_trk -> getOmega(), a_trk->getTanLambda(), BField);
	int NCaloHits = a_clu->getCalorimeterHits().size();

	for(int i1 = 0; i1 < NCaloHits; i1++)
	{
		CalorimeterHit *a_hit = a_clu->getCalorimeterHits()[i1];

		float BushTime = a_Helix->getDistanceToPoint((float*)a_hit->getPosition(), BushDist);
		int nCircles = a_Helix->getNCircle((float*)a_hit->getPosition());
		if(BushTime > 0 && BushDist[2] < Distance )
		{
			//Time = BushTime;
			NC = nCircles;
			Distance = BushDist[2];
		}
	}

	delete a_Helix;
	return NC;
}

float *DisTrackClu(Track * a_trk, Cluster * a_clu)
{
	static float Distance[2] = {1.0E10, 1.0E10};
	float BushDist[3] = {0, 0, 0};
	float BushTime = -100;
	float BushTimeHit = -100;  
	float MinTrackCluHitDis = 1.0E10; 

	HelixClass * a_Helix = new HelixClass();
	a_Helix->Initialize_Canonical(a_trk->getPhi(), a_trk -> getD0(), a_trk -> getZ0(), a_trk -> getOmega(), a_trk->getTanLambda(), BField);

	/*
	   float CluEn = a_clu->getEnergy();	
	   float TrkEn = 0;
	   for (int q3 = 0; q3 < 3; q3 ++)
	   {
	   TrkEn += (a_Helix->getMomentum()[q3])*(a_Helix->getMomentum()[q3]);
	   }
	   TrkEn = sqrt(TrkEn);
	   */

	BushTime = a_Helix->getDistanceToPoint((float*)a_clu->getPosition(), BushDist);
	if(BushTime > 0)
	{
		Distance[0] = BushDist[2];
	}

	int NCaloHits = a_clu->getCalorimeterHits().size();
	for(int i1 = 0; i1 < NCaloHits; i1++)
	{
		CalorimeterHit *a_hit = a_clu->getCalorimeterHits()[i1];
		BushTimeHit = a_Helix->getDistanceToPoint((float*)a_hit->getPosition(), BushDist);
		if(BushTimeHit > 0 && BushDist[2] < MinTrackCluHitDis)
		{
			MinTrackCluHitDis = BushDist[2];
		}
	}	

	if(MinTrackCluHitDis < 1.0E6)
	{
		Distance[1] = MinTrackCluHitDis;
	}
	delete a_Helix;

	return Distance;
}



int ClusterFlag(Cluster* a_tree, Track* a_trk, LCCollection *col_TPCTrk)
{
	// give each charged core cluster a flag
	//  Fragmentation:       999
	//  MIP: matched         130
	//       non-matched     131
	//  EM:  matched         110
	//       non-matched     111
	//  HAD: matched         211
	//       non-matched     212

	int CluIDFlag = 999;
	int ClusterID = 211;
	int EcalNHit, HcalNHit, NLEcal, NLHcal, NH[16];
	float avEnDisHtoL;
	float EcalEn, HcalEn, EClu, cluDepth, maxDepth, minDepth, MaxDisHel, MinDisHel, FD_all, FD_ECAL, FD_HCAL;
	float crdis, EEClu_L10, EEClu_R, EEClu_r, EEClu_p, rms_Ecal, rms_Hcal, rms_Ecal2, rms_Hcal2, av_NHE, av_NHH;
	int AL_Ecal, AL_Hcal;
	float FD_ECALF10;
	float dEdx = 0;

	const string ECALCellIDDecoder = "M:3,S-1:3,I:9,J:9,K-1:6";
	CellIDDecoder<CalorimeterHit> idDecoder(ECALCellIDDecoder);
	const float mass = 0.139;       //Pion Mass

	HelixClass * TrkInit_Helix = new HelixClass();
	TrkInit_Helix->Initialize_Canonical(a_trk->getPhi(), a_trk -> getD0(), a_trk -> getZ0(), a_trk -> getOmega(), a_trk->getTanLambda(), BField);
	float TrackEn = mass*mass;

	int ntpctrk = col_TPCTrk->getNumberOfElements();

	for (int q3 = 0; q3 < 3; q3 ++)
	{
		TrackEn += (TrkInit_Helix->getMomentum()[q3])*(TrkInit_Helix->getMomentum()[q3]);
	}
	delete TrkInit_Helix;

	TrackEn = sqrt(TrackEn);
	int nSubTrk = a_trk->getTracks().size();

	//int NHit = a_tree->getCalorimeterHits().size();
	if(1 == 1) //if ( (NHit > 4 && TrackEn > 1) || TrackEn <= 1 )
	{

		for(int t0 = 0; t0 < ntpctrk; t0++)
		{
			Track *a_tpctrk = dynamic_cast<EVENT::Track *>(col_TPCTrk->getElementAt(t0));
			float TPCD0 = a_tpctrk->getD0();
			float TPCZ0 = a_tpctrk->getZ0();
			float TPCPhi = a_tpctrk->getPhi();
			float TPCTL = a_tpctrk->getTanLambda();
			for(int t1 = 0; t1 < nSubTrk; t1++)
			{
				Track* a_SubTrk = a_trk->getTracks()[t1];
				float SubD0 = a_SubTrk->getD0();
				float SubZ0 = a_SubTrk->getZ0();
				float SubPhi = a_SubTrk->getPhi();
				float SubTL = a_SubTrk->getTanLambda();
				if (fabs(SubTL-TPCTL) < 1.0E-6 && fabs(SubD0-TPCD0) < 1.0E-6 && fabs(SubZ0-TPCZ0) < 1.0E-6 && fabs(SubPhi-TPCPhi) < 1.0E-6)
				{
					std::vector<TrackerHit*> trhits = a_tpctrk->getTrackerHits();
					int nhit = trhits.size();
					int chose_high = int(0.3*nhit);
					float cen_high = 999.;
					int maxhit = 99999;
					float totaled = 0.;

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
						if(hit->getEDep()/mindis<cen_high){
							float dedx = hit->getEDep()/mindis;
							totaled += dedx;
							nhiteff ++;
						}
					}
					dEdx = totaled/nhiteff;
				}
			}

		}

		TVector3 CluPos;
		CluPos = a_tree->getPosition();
		TVector3 IntDir = ClusterCoG(a_tree)-CluPos;
		EClu = a_tree->getEnergy();
		EcalNHit = 0;
		HcalNHit = 0;
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
		if(currCluNHits == 0) return CluIDFlag;
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

			if( fabs(a_hit->getEnergy() - DHCALCalibrationConstant) < 1.0E-6 )      //or other fancy judgements...^M
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
			else
			{
				EcalNHit++;
				EcalEn += a_hit->getEnergy();
				Ecalhits.push_back(a_hit);
				if(NLayer< 10) Ecalf10hits.push_back(a_hit);
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

		CalorimeterHit * maxdis_hit = a_tree->getCalorimeterHits()[index1];
		CalorimeterHit * mindis_hit = a_tree->getCalorimeterHits()[index2];
		TVector3 maxpos = maxdis_hit->getPosition();
		TVector3 minpos = mindis_hit->getPosition();
		TVector3 GraPos = ClusterCoG(a_tree);
		cluDepth = (maxpos-minpos).Mag();

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
			totHitEn+=HitEn;
			totHitEnDis+=HitEn*disHtoL;
		}
		avEnDisHtoL = totHitEnDis/totHitEn;
		FD_all = FDV2(allhits, ECALCellIDDecoder);
		FD_ECAL = FDV2(Ecalhits, ECALCellIDDecoder);
		FD_HCAL = FDV2(Hcalhits, ECALCellIDDecoder);
		FD_ECALF10 = FDV2(Ecalf10hits, ECALCellIDDecoder);

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




		bool cutmu1 = EcalNHit+2.*HcalNHit < 500;
		bool cutmu2;

		if(cluDepth < 1100) cutmu2 = 0;
		else if(cluDepth < 1400) cutmu2 = FD_all < 0.5/pow(400,2)*pow((cluDepth-1000),2);
		else cutmu2 = 1;

		bool cutmu3;
		if(TrackEn > 70) cutmu3 = FD_all < (600./avEnDisHtoL + 20)/100.;
		else if (TrackEn > 10) cutmu3 = FD_all < (600./avEnDisHtoL -10 + 0.5*TrackEn)/100.;
		else if (TrackEn > 7.5) cutmu3 = FD_all < (600./avEnDisHtoL -10)/100.;
		else cutmu3 = 1;
		bool cutmu3b;
		if (TrackEn > 10) cutmu3b = avEnDisHtoL < 25;
		else if (TrackEn > 4.5) cutmu3b = avEnDisHtoL < 25 + 10 * (10 - TrackEn);
		else cutmu3b = 1;
		bool cutmu10en = cutmu1 && cutmu2 && cutmu3 && cutmu3b;


		bool cutmu4 = FD_HCAL >= 0;
		bool cutmu5 = cluDepth > 750 - 20/TrackEn;
		bool cutmu6 = cluDepth > 1200;
		bool cutmu7 = FD_all < 0.3/sqrt(400.)*sqrt(cluDepth-1200.);
		bool cutmu8 = rms_Hcal < 10 && FD_HCAL < -0.25/600.*MaxDisHel+0.25;
		bool cutmu9 = FD_all < 0.35/sqrt(400)*sqrt(400-MaxDisHel) && MaxDisHel < 400 && EcalNHit+2.*HcalNHit > 85;
		bool cutmu;

		if(TrackEn > 9.5) cutmu = cutmu10en;
		else if(TrackEn < 1.5)
		{
			if(cutmu5) cutmu = 1;
			else cutmu = 0;
		}
		else if(TrackEn < 3.5)
		{
			if(cutmu4 && cutmu6 && cutmu7) cutmu = 1;
			else cutmu = 0;
		}
		else
		{
			if(cutmu3b && cutmu4 && cutmu6 && cutmu7 && cutmu8 && cutmu9) cutmu = 1;
			else cutmu = 0;
		}

		//bool cute1 = EcalEn/(EcalEn+HcalEn) > 0.9;
		bool cute2 = (dEdx > 0.17e-6 && dEdx < 0.3e-6)||dEdx==0;
		bool cute3 = FD_all > 0.9*sin(cluDepth*1.57/800.) && cluDepth < 800;
		bool cute4 = FD_ECALF10 > 0.9*sin((cluDepth-200)*1.57/600.);
		bool cute5 = EEClu_r/TrackEn > 0.8 * sin((avEnDisHtoL-15)*1.57/20.) && avEnDisHtoL < 35;
		bool cute6 = FD_ECAL > 0.2 * log10(EcalNHit) && log10(EcalNHit) > 1.5;
		bool cute7 = FD_all >= 0.6/1.5*(log10(EcalNHit+2.*HcalNHit)-1.2-0.4*TrackEn/100.);
		bool cute8 = SDTheta < 0.012/1.5*(log10(EcalNHit+2.*HcalNHit)-1);
		bool cute;
		if(TrackEn < 1.5) cute = cute2;
		else cute = cute2 && cute3 && cute4 && cute5 && cute6 && cute7 && cute8;


		if(cutmu) ClusterID = 13;
		else if (cute)  ClusterID = 11;

		if(ClusterID == 13 )
		{
			if(NLEcal+NLHcal < 30) CluIDFlag = 131;
			else CluIDFlag = 130;
		}

		if(ClusterID == 11 )
		{
			if(EClu < 0.8*TrackEn) CluIDFlag = 111;
			else CluIDFlag = 110;
		}

		if(ClusterID == 211 )
		{
			if(EClu < 0.8*TrackEn) CluIDFlag = 212;
			else CluIDFlag = 211;
		}

	}

	return CluIDFlag;

}


int PhotonTag(Cluster *a_clu)
{
	const string ECALCellIDDecoder = "M:3,S-1:3,I:9,J:9,K-1:6";

	float CluEn, CluFD, d_angle, ECALEn, MaxLength, Depth;
	TVector3 cogdir;
	int CluSize, nLECAL, minlay, lay;
	CluEn = 0.;
	CluFD = -1;
	d_angle = 0.;
	ECALEn = 0.;
	MaxLength = 0.;
	nLECAL =0;
	minlay = 30;
	lay = 0;
	TVector3 CluPos;
	std::vector<CalorimeterHit*> Ecalhits;

	CluEn = a_clu->getEnergy();
	CluSize = a_clu->getCalorimeterHits().size();
	if(CluSize)
		CluFD = FDV3(a_clu, ECALCellIDDecoder);
	CluPos = a_clu->getPosition();
	Depth = DisSeedSurface(CluPos);
	Ecalhits.clear();
	float cluE_R = 0;
	float cluE_r = 0;
	float cluE_p = 0;


	CellIDDecoder<CalorimeterHit> idDecoder(ECALCellIDDecoder);
	for(int i1 = 0; i1 < CluSize; i1++)
	{
		TVector3 HitPos;
		//		float crdis, cluE_R, cluE_r, cluE_p;
		CalorimeterHit * a_hit = a_clu->getCalorimeterHits()[i1];
		float HitEn = a_hit->getEnergy();
		HitPos = a_hit->getPosition();
		cogdir = ClusterCoG(a_clu)-CluPos;
		float crdis = (CluPos-HitPos).Mag() * sin((CluPos-HitPos).Angle(cogdir));
		if( fabs(HitEn - DHCALCalibrationConstant) > 1.0E-6 ){
			if(crdis<22){
				cluE_R += a_hit->getEnergy();
			}
			if(crdis<11){
				cluE_r += a_hit->getEnergy();
			}
			if(crdis<6){
				cluE_p += a_hit->getEnergy();
			}

			if(crdis > MaxLength)
			{
				MaxLength = crdis;
			}
			int thislay = idDecoder(a_hit)["K-1"];
			if(thislay<minlay){
				minlay=thislay;
			}
			ECALEn += HitEn;
			Ecalhits.push_back(a_hit);
		}


	}
	d_angle = cogdir.Angle(CluPos);

	lay=minlay;
	nLECAL = ActiveLayers(Ecalhits, ECALCellIDDecoder);
	if((CluEn<5&&cluE_p/ECALEn<0.764+0.24*exp(CluEn*-0.54)&&Depth<80&&MaxLength>-44.7+48.1*pow(CluEn+0.9,0.15)&&d_angle<0.75&&CluFD>-12.2972+12.1928/pow(CluEn+1.67389,-0.0221877)&&CluFD<-1.576+2.289/pow(CluEn+0.0407,-0.02993)&&nLECAL>-53.54+1.662/pow(CluEn+38.18,-0.9672)&&lay<12)||(CluEn>5&&cluE_p/ECALEn<0.68+0.106*exp(CluEn*-0.019)&&Depth<80&&MaxLength>-44.7+48.1*pow(CluEn+0.9,0.15)&&d_angle<0.75&&CluFD>-7.273+7.434/pow(CluEn-1.367,-0.0176)&&CluFD<-1.54+2.289/pow(CluEn+0.0407,-0.0287)&&nLECAL>min(-28.55+3.499/pow(CluEn+86.68,-0.5276),28.-lay)&&lay<12))
	{ 
		return 1;
	}
	else return 0;
}

float EnUltraHotCorr(float En, Cluster *a_tree)
{
	float EnUltraHot = 0;
	for(int i = 0; i < int( a_tree->getCalorimeterHits().size() ); i++)
	{
		CalorimeterHit* a_hit = a_tree->getCalorimeterHits()[i];
		if(a_hit->getEnergy() > 1)
			EnUltraHot += a_hit->getEnergy();	
	}
	return En - EnUltraHot; 
}

int ClusterFlag1st(Cluster* a_tree)
{

	int ClusterID = 211;
	int EcalNHit, HcalNHit, NH_ECALF10;
	float avEnDisHtoL;
	float EcalEn, HcalEn, EClu, cluDepth, maxDepth, minDepth, FD_all, FD_ECAL, FD_HCAL;
	float FD_ECALF10;


	const string ECALCellIDDecoder = "M:3,S-1:3,I:9,J:9,K-1:6";
	CellIDDecoder<CalorimeterHit> idDecoder(ECALCellIDDecoder);

		TVector3 CluPos;
		CluPos = a_tree->getPosition();
                
		EClu = a_tree->getEnergy();
		EcalNHit = 0;
		HcalNHit = 0;
		EcalEn = 0;
		HcalEn = 0;
                float currDepth = 0;
		maxDepth = -100;
		minDepth = 1E6;
		
		std::vector<CalorimeterHit*> allhits;

		std::vector<CalorimeterHit*> Ecalhits;
		std::vector<CalorimeterHit*> Hcalhits;

                std::vector<CalorimeterHit*> Ecalf10hits;



		allhits.clear();
		Ecalhits.clear();
		Hcalhits.clear();

                Ecalf10hits.clear();



		std::vector<float> hitTheta;
		hitTheta.clear();
                                
                for(unsigned int j1 = 0; j1 < a_tree->getCalorimeterHits().size(); j1++)
                {
                          CalorimeterHit * a_hit = a_tree->getCalorimeterHits()[j1];
                          
			  TVector3 tmpPos = a_hit->getPosition();
			  hitTheta.push_back(tmpPos.Theta());
                       
                }

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
		if(currCluNHits == 0) return 1;
		int index1 = 0, index2 = 0;
		int HitDepth50 = 1.0e4;

		for(int s1 = 0; s1 < currCluNHits; s1++)
		{
			CalorimeterHit * a_hit = a_tree->getCalorimeterHits()[s1];
			allhits.push_back(a_hit);
			int NLayer = idDecoder(a_hit)["K-1"];
			float tmpHitEn = a_hit->getEnergy();
			HitPos = a_hit->getPosition();

			currDepth = DisSeedSurface(HitPos);
			if(tmpHitEn > 0.05 && currDepth < HitDepth50)
			{
				HitDepth50 = currDepth;
			}


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

			if( fabs(a_hit->getEnergy() - DHCALCalibrationConstant) < 1.0E-6 )      //or other fancy judgements...^M
			{
				HcalNHit++;
				HcalEn += a_hit->getEnergy();
				Hcalhits.push_back(a_hit);
			}
			else
			{
				EcalNHit++;
				EcalEn += a_hit->getEnergy();
				Ecalhits.push_back(a_hit);
                                if(NLayer< 10) Ecalf10hits.push_back(a_hit);
			}
		}
                
		CalorimeterHit * maxdis_hit = a_tree->getCalorimeterHits()[index1];
		CalorimeterHit * mindis_hit = a_tree->getCalorimeterHits()[index2];
		TVector3 maxpos = maxdis_hit->getPosition();
		TVector3 minpos = mindis_hit->getPosition();
		TVector3 GraPos = ClusterCoG(a_tree);
		cluDepth = (maxpos-minpos).Mag();

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
			totHitEn+=HitEn;
			totHitEnDis+=HitEn*disHtoL;
		}
		avEnDisHtoL = totHitEnDis/totHitEn;
		FD_all = FDV2(allhits, ECALCellIDDecoder);
		FD_ECAL = FDV2(Ecalhits, ECALCellIDDecoder);
		FD_HCAL = FDV2(Hcalhits, ECALCellIDDecoder);
                FD_ECALF10 = FDV2(Ecalf10hits, ECALCellIDDecoder);
		NH_ECALF10 = Ecalf10hits.size();


                bool cute1 = log10(FD_all/log10(EClu)) > (log10(cluDepth/log10(EClu))-2.9);
		bool cute2 = FD_ECAL > 0.1 && FD_ECAL > 0.2*log10(avEnDisHtoL-3)+0.05*sqrt(avEnDisHtoL-3);
		float x;
		if (EcalNHit == 0) x = 0;
		else x = float(NH_ECALF10)/float(EcalNHit);
                bool cute3;
		if (x < 0.9)  cute3 = FD_ECALF10 > 0.3/sqrt(sqrt(0.9))*sqrt(sqrt(0.9-x));
		else cute3 = 1;
		bool cute4 = HcalNHit/EClu < 0.3;                
                bool cute;
		cute = cute1 && cute2 && cute3 && cute4;


		bool cutmu1 = (cluDepth > 1300 || EClu < 4/1000.*cluDepth+0.9);
                bool cutmu2 = EClu < 15;
                bool cutmu3 = avEnDisHtoL*SDTheta/EClu <0.02 && FD_all < 0.4/0.025*(0.025-avEnDisHtoL*SDTheta/EClu);
                
                
                bool cutmu;

                cutmu = cutmu1 && cutmu2 && cutmu3;
                

		if(currCluNHits <= 4) ClusterID = 1;
                else if(cute) ClusterID = 11;
                else if (cutmu)  ClusterID = 13;
		else 
		{
			bool cutef1 = FD_HCAL == -1;
			bool cutef2 = minDepth < 50;
			bool cutef3 = FD_ECAL > 0.2;
			bool cutef3b = 1;
			if (FD_ECAL < 0.6) cutef3b = FD_ECAL > 6/3.5*(1.75-log10((EcalNHit+HcalNHit)/EClu));
			bool cutef4;
			if (avEnDisHtoL <= 8)
				cutef4 = FD_ECALF10 > 0.2;
			else
				cutef4 = FD_ECALF10 > 0.2+0.8*sqrt(log10(avEnDisHtoL/8));
			bool cutef5 = FD_ECAL/EClu > (log10(EcalNHit)-1.6)*(log10(EcalNHit)-1.6)*4+0.25 && FD_ECALF10 > (avEnDisHtoL-9)*(avEnDisHtoL-9)*0.006+0.2;
			bool cutef;
			
			cutef = (cutef1 && cutef2 && cutef3 && cutef3b && cutef4) || cutef5;
			if(cutef) ClusterID = 11;  
		}

		if(ClusterID == 211)
		{
			bool cutmuf1 = FD_all == 0 && log10((EcalNHit+2*HcalNHit)*EClu) > 1.2;
			bool cutmuf2 = log10((EcalNHit+2*HcalNHit)*EClu) > 1.9 && log10((EcalNHit+2*HcalNHit)*EClu) < 3.1;
			bool cutmuf3 = log10(FD_all/cluDepth) < -223.074+431.315*log10((EcalNHit+2*HcalNHit)*EClu)-331.72*pow(log10((EcalNHit+2*HcalNHit)*EClu),2)+124.588*pow(log10((EcalNHit+2*HcalNHit)*EClu),3)-22.8786*pow(log10((EcalNHit+2*HcalNHit)*EClu),4)+1.64577*pow(log10((EcalNHit+2*HcalNHit)*EClu),5);
			bool cutmuf = cutmuf1 || (cutmuf2 && cutmuf3);
			if(cutmuf) ClusterID = 13;
			else if(minDepth < 0.77+0.23*EClu) ClusterID = 211;
			else ClusterID = 2; 
		}

	return ClusterID;
}

float EMClusterEE( Cluster *inputCluster )
{
	const string ECALCellIDDecoder = "M:3,S-1:3,I:9,J:9,K-1:6";
	float aaa = -50.2288, ab = 219.398,ac =0.17679,ad =  0.00241144;
	float ba =  -56.6164,bbb = 162.647,bc = 0.679974,bd  = 0.00423267,be  = -0.324786;
	float ff1a = -21.878,ff1b = 1146.3, ff1c =  0.0267898, ff1d  =  0.000712903;
	float ff2a = -85.8291,ff2b = 2466.59,ff2c = 0.55722, ff2d  =  0.00159572;
	float ArEhito10=0, ArEhite10=0, ArEhito20=0, ArEhite20m=0, ArEhite20p=0;
	int   ArNhito10=0, ArNhite10=0, ArNhito20=0, ArNhite20m=0, ArNhite20p=0;
	float diff=1.0;
	float x = 0, y= 0, f1 = 0, f2 = 0;
	float _costheta = 0;
	float Ethetacorr = 1;
	float Ephicorr = 1;
	float EMC = inputCluster->getEnergy(); 
	TVector3 CluPos = inputCluster->getPosition(); 
	float CluTheta = CluPos.Theta();
	float CluPhi = CluPos.Phi()*5.72957795130823229e+01;	


	int NCluHits = inputCluster->getCalorimeterHits().size();
	CellIDDecoder<CalorimeterHit> idDecoder(ECALCellIDDecoder);
	for(int s1 = 0; s1 < NCluHits; s1++)
	{
		CalorimeterHit * a_hit = inputCluster->getCalorimeterHits()[s1];
		int NLayer = idDecoder(a_hit)["K-1"];
		if(NLayer > 20){
			if (NLayer%2 ==0){
				ArNhite10 ++;
				ArEhite10 +=a_hit->getEnergy();
			}
			else if (NLayer%2 ==1){
				ArNhito10 ++;
				ArEhito10 +=a_hit->getEnergy();
			}

		}
		else if(NLayer <= 20){
			if (NLayer%2 ==0){
				if(NLayer ==0){
					ArNhite20m ++;
					ArEhite20m +=a_hit->getEnergy();
				}
				if(NLayer >0){
					ArNhite20p ++;
					ArEhite20p +=a_hit->getEnergy();
				}
			}
			else if (NLayer%2 ==1){
				ArNhito20 ++;
				ArEhito20 +=a_hit->getEnergy();
			}
		}
	}
	while(diff>0.01 && EMC > 0)
	{
		float temp_a=sqrt(EMC*EMC);
		x=exp((-EMC+aaa)/ab)+ac/sqrt(EMC)+ad*EMC;
		y=1/(exp((-EMC+ba)/bbb)+bc/sqrt(EMC)+bd*EMC)+be;
		f1=1/(exp((-EMC+ff1a)/ff1b)+ff1c/sqrt(EMC)+ff1d*EMC);
		f2=1/(exp((-EMC+ff2a)/ff2b)+ff2c/sqrt(EMC)+ff2d*EMC);
		EMC=x*ArEhito20+f1*ArEhite20p+y*ArEhito10+f2*ArEhite10;
		float temp_b=sqrt(EMC*EMC);
		diff=fabs(temp_a-temp_b);

		
	}
	_costheta=cos(CluTheta);
	/*
	Ethetacorr=(50/(0.294596*fabs(_costheta)*fabs(_costheta)-1.58336*fabs(_costheta)+49.9219));
	for(int n=0;n<8;n++){
		if(CluPhi>=(20.5+n*45) && CluPhi<(62.5+n*45) ) Ephicorr=(50/(1/(exp((-(CluPhi-n*45)-15.1254)/0.918483)+0.00839675/sqrt(CluPhi-n*45)+0.000369287*logf(CluPhi-n*45)+0.0174033)));
		if(CluPhi>=(17.5+n*45) && CluPhi<(20.5+n*45) ) Ephicorr=(50/(0.4255992*(CluPhi-n*45)*(CluPhi-n*45)-0.4748*(CluPhi-n*45)+63));
		if(CluPhi>=(20.5+n*45) && CluPhi<(22.5+n*45) ) Ephicorr=(50/(0.486258*(CluPhi-n*45)*(CluPhi-n*45)-0.413742*(CluPhi-n*45)+12.6491));
	}
	*/

	if(EMC/inputCluster->getEnergy() < 0.5) EMC = inputCluster->getEnergy();

	 Ethetacorr=(50/(0.294596*fabs(_costheta)*fabs(_costheta)-1.58336*fabs(_costheta)+49.9219));

        float ModuleCluPhi = 0;
        if(CluPhi > 17.5)
        {
                ModuleCluPhi = CluPhi - int((CluPhi -17.5)/45)*45;
        }
        else
        {
                ModuleCluPhi = CluPhi - int((CluPhi -17.5)/45)*45 + 45;
        }

        if(ModuleCluPhi>=17.5 && ModuleCluPhi<20.5 ) Ephicorr=(50/(-0.122*ModuleCluPhi*ModuleCluPhi+2.76*ModuleCluPhi+38.04));
        if(ModuleCluPhi>=20.5 && ModuleCluPhi<22.5 ) Ephicorr=(50/(0.22*ModuleCluPhi*ModuleCluPhi-6.43*ModuleCluPhi+81.69));
        if(ModuleCluPhi>=22.5 && ModuleCluPhi<62.5 ) Ephicorr=50*(0.00839675/sqrt(ModuleCluPhi)+0.000369287*logf(ModuleCluPhi)+0.0174033);

	EMC=EMC*Ethetacorr*Ephicorr;

	return EMC;
}


float ClusterEE(Cluster* inputCluster)
{
	float ClusterEnergy = 0;
	float tmpCluEn = 0;
	int CluType = 211;
	if(ClusterFlag1st(inputCluster) == 11) CluType = 22;

	int NCluHit = inputCluster->getCalorimeterHits().size();
	float hitEn = 0;
	float EnCorrector = 1;

	if(CluType == 22)
	{
		ClusterEnergy = EMClusterEE(inputCluster);
	}
	else
	{
		for(int i = 0; i < NCluHit; i++)
		{
			CalorimeterHit * a_hit = inputCluster->getCalorimeterHits()[i];
			hitEn = a_hit->getEnergy();

			if(CluType == 211)
			{
				if( hitEn < 1.5 )       // Or some function of the initCluEn
				{
					tmpCluEn += hitEn;
				}
				else
				{
					cout<<"Ultra Hot Hit"<<endl;	// Use ShuZhen's Hot Hit Finding function...
					tmpCluEn += 0.05;       // MIP Hit Energy, Value to be determined
				}
			}
			else    // For EM & MIP: should veto accordingly
			{
				tmpCluEn += hitEn;
			}
		}

		EnCorrector = 1;
		if(tmpCluEn > 1.5 && tmpCluEn < 22 && 1 == 0)
		{
			EnCorrector = 0.6*(1 + 1.0/log10(tmpCluEn)); 
		}

		//ClusterEnergy = tmpCluEn*EnCorrector;
		ClusterEnergy = HADClusterEE(tmpCluEn,inputCluster);
	}

	return ClusterEnergy;
}

float HADClusterEE(float En, Cluster *inputCluster)
{
	float outputEn;

	float factorTheta = 1.0;
	float factorEn = 1.0;
	TVector3 CluPos = inputCluster->getPosition();
	float tmpTheta = CluPos.Theta();
	if(fabs(cos(tmpTheta)) < 0.988)
	{ 
		factorTheta = 0.501563+2.27856*tmpTheta-6.47055*pow(tmpTheta,2)+7.82752*pow(tmpTheta,3)-4.54267*pow(tmpTheta,4)+1.25052*pow(tmpTheta,5)-0.131147*pow(tmpTheta,6);
	}
	//else
	//{
	//	factorTheta = 0.2+0.554/0.2*(1.0-fabs(cos(tmpTheta)));
	//}
	if(En < 100)
	{
		factorEn = 0.690245+0.0600232*En-0.00478908*pow(En,2)+0.000204681*pow(En,3)-4.8215e-06*pow(En,4)+6.23579e-08*pow(En,5)-4.14229e-10*pow(En,6)+1.10348e-12*pow(En,7);
	}
	else
	{
		factorEn = 0.690245+0.0600232*100-0.00478908*pow(100,2)+0.000204681*pow(100,3)-4.8215e-06*pow(100,4)+6.23579e-08*pow(100,5)-4.14229e-10*pow(100,6)+1.10348e-12*pow(100,7);
	}
	
	//outputEn = En/(factorTheta*factorEn);
	outputEn = En/factorEn;
	cout << "f1:f2 " << factorTheta <<":"<<factorEn<<endl;
	return outputEn;
	
}




