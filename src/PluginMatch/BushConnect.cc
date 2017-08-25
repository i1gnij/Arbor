#include <BushConnect.hh>
#include <ArborTool.hh>
#include <ArborToolLCIO.hh>
#include <DetectorPos.hh>
#include <LICH.hh>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <EVENT/Cluster.h>
#include <EVENT/LCRelation.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/ClusterImpl.h>
#include "UTIL/CellIDDecoder.h"
#include "HelixClass.hh"		//in Marlin Util
#include <values.h>
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
//#include <ArborTrack.hh>

using namespace std;

const float TightGeoThreshold = 30;
const string ECALCellIDDecoder  = "M:3,S-1:3,I:9,J:9,K-1:6";
int Diag = 0; 

BushConnect aBushConnect ;
BushConnect::BushConnect()
	: Processor("BushConnect"),
        _output(0)
{
	_description = "Track Cluster Linking. Track info represented by MCTruth at this moment" ;
}

void BushConnect::init() {
	printParameters();
	Cluflag.setBit(LCIO::CHBIT_LONG);
}


void BushConnect::Clean(){

	Track_Energy.clear();
	Track_Type.clear();
	Track_Phi.clear();
	Track_Theta.clear();
	SortedTracks.clear();
	SortedSMBushes.clear();
	ClusterType_1stID.clear();
	ChCoreID.clear();

	ecalnecore_EM.clear();
	ecalnecore_NonEM.clear();
	ecalfrag.clear();
	ecalundef.clear();
	ecalfrag_TBM_CH.clear();
	ecalfrag_TBM_NE_EM.clear();
	ecalfrag_TBM_NE_NonEM.clear();
	ecalpotentialbackscattering.clear();
	ecalundef_iso.clear();
	trkendposition.clear();
	chargedclustercore.clear();
	non_chargedclustercore.clear();		//all clusters besides charged cluster core
	pem_neutral_core.clear();
	non_charged_pem_neutral_core.clear();

	TrackEndPoint.clear();
	TrackStartPoint.clear();
	CluFD.clear();
	CluEnergy.clear();
}

void BushConnect::TrackSort(LCEvent* evtPP) //, &std::map<Track*, int>Track_Tpye, &std::map<Track*, float> Track_Energy)
{

	LCCollection * TrackColl = evtPP->getCollection("MarlinTrkTracks");

	int NTrk = TrackColl->getNumberOfElements();
	float D0 = 0;
	float Z0 = 0;
	int NTrkHit = 0;
	const float mass = 0.139;	//Pion Mass
	TVector3 EndPointPos, StartPointPos; 
	int TrackType = 0; 

	std::vector<Track*> tracks_HQ_Barrel; 
	std::vector<Track*> tracks_HQ_Endcap;
	std::vector<Track*> tracks_HQ_Shoulder;
	std::vector<Track*> tracks_HQ_Forward; 
	std::vector<Track*> tracks_MQ_Barrel;
	std::vector<Track*> tracks_MQ_Endcap;
	std::vector<Track*> tracks_MQ_Shoulder;
	std::vector<Track*> tracks_MQ_Forward;
	std::vector<Track*> tracks_Vtx; 
	std::vector<Track*> tracks_LQ; 
	std::vector<Track*> tracks_LE; 
	std::vector<Track*> curr_tracks;

	trkendposition.clear();

	tracks_HQ_Barrel.clear();
	tracks_HQ_Endcap.clear();
	tracks_HQ_Shoulder.clear();
	tracks_HQ_Forward.clear();
	tracks_MQ_Barrel.clear();
	tracks_MQ_Endcap.clear();
	tracks_MQ_Shoulder.clear();
	tracks_MQ_Forward.clear();
	tracks_Vtx.clear();
	tracks_LQ.clear();
	tracks_LE.clear();

	std::vector<Track*> tracks_ILL;
	tracks_ILL.clear();
	std::vector<Track*> tracks_preInteraction;
	tracks_preInteraction.clear();	//Used to denote pion and electron interaction inside TPC/Tracker. Simply vetoed for avoid double counting... but muon may still be problematic. Better way of treating would be find the cascade photons & tracks - clusters, and veto all the daughters instead of mother. Similar can done for Kshort...
	// Condition, tracks_head to others tail. head position far from boundary. and, track energy >= sum of cascade

	std::vector<int> TrackOrder; 
	TrackOrder.clear();	
	std::map<Track*, int> Track_Index; 
	Track_Index.clear();
	Track_Energy.clear();
	Track_Type.clear();
	Track_P3.clear();

	TrackEndPoint.clear();
	TrackStartPoint.clear();

	for(int i0 = 0; i0 < NTrk; i0++)
	{
		Track* a_Trk = dynamic_cast<Track*>( TrackColl->getElementAt( i0 ) );
		NTrkHit = a_Trk->getTrackerHits().size();		
		EndPointPos = (a_Trk->getTrackerHits()[NTrkHit - 1])->getPosition();	
		StartPointPos = (a_Trk->getTrackerHits()[0])->getPosition();
		TrackEndPoint[a_Trk] = EndPointPos;
		TrackStartPoint[a_Trk] = StartPointPos;

		HelixClass * TrkInit_Helix = new HelixClass();
		TrkInit_Helix->Initialize_Canonical(a_Trk->getPhi(), a_Trk -> getD0(), a_Trk -> getZ0(), a_Trk -> getOmega(), a_Trk->getTanLambda(), BField);
		float TrackEn = mass*mass;

		for (int q3 = 0; q3 < 3; q3 ++)
		{
			TrackEn += (TrkInit_Helix->getMomentum()[q3])*(TrkInit_Helix->getMomentum()[q3]);
		}
		TVector3 TrkMom(TrkInit_Helix->getMomentum()[0],TrkInit_Helix->getMomentum()[1],TrkInit_Helix->getMomentum()[2]);
		
		TrackEn = sqrt(TrackEn);
		Track_Energy[a_Trk] = TrackEn;
		Track_Theta[a_Trk] = TrkMom.Theta();
		Track_Phi[a_Trk] = TrkMom.Phi();
		Track_P3[a_Trk] = TrkMom;		
		
		delete TrkInit_Helix;
	}

	TVector3 currEp, currSp;
	float currMotherEn = 0;
	float sumDauEn = 0; 

	for(int i1 = 0; i1 < NTrk; i1++)
	{
		Track* a_Trk = dynamic_cast<Track*>( TrackColl->getElementAt( i1 ) );		
		currEp = TrackEndPoint[a_Trk];

		if( currEp.Perp() < 1600 && currEp.Perp() > 400 && abs(currEp.Z()) < 2000 )	//Only check 
		{
			currMotherEn = Track_Energy[a_Trk];
			sumDauEn = 0;	
			for(int i2 = 0; i2 < NTrk; i2++)
			{
				Track* b_Trk = dynamic_cast<Track*>( TrackColl->getElementAt( i2 ) );
				if(i2 != i1)
				{
					currSp = TrackStartPoint[b_Trk];
					if( (currEp - currSp).Mag() < 40  )
						sumDauEn += Track_Energy[b_Trk];
				}
			}
			if(currMotherEn + 0.1 > 0.9*sumDauEn && currMotherEn > 3 && sumDauEn > 0 )	//Some protection is always needed...
			{
				tracks_preInteraction.push_back(a_Trk);
				// cout<<"PPPPPPPPPPPPPPPPPPPRRRRRRRRRRRRRRRRRRRRRRinteraction Track found "<<currMotherEn<<" to "<<sumDauEn<<endl; 
			}
		}
	}

	for(int t0 = 0; t0 < NTrk; t0++)
	{
		Track* a_Trk = dynamic_cast<Track*>( TrackColl->getElementAt( t0 ) );
		D0 = a_Trk->getD0();
		Z0 = a_Trk->getZ0();
		NTrkHit = a_Trk->getTrackerHits().size();
		TrackerHit * last_hit = a_Trk->getTrackerHits()[NTrkHit - 1];
		EndPointPos = last_hit->getPosition();
		trkendposition[a_Trk] = EndPointPos;
		StartPointPos = (a_Trk->getTrackerHits()[0])->getPosition();
		Track_Index[a_Trk] = t0;

		if( NTrkHit > 9 || (fabs(EndPointPos.Z()) > LStar - 500 && EndPointPos.Perp() < TPCInnerRadius ) || fabs(EndPointPos.Z()) > ECALHalfZ - 200  )		// Min requirement for track quality
		{	// LStar - 500, suppose to be the last Disk Position

			if( find(tracks_preInteraction.begin(), tracks_preInteraction.end(), a_Trk ) != tracks_preInteraction.end() )
			{
				cout<<"So We Drop it! "<<Track_Energy[a_Trk]<<endl; 
				continue; 
			}

			TrackType = 0;
			if((Track_Energy[a_Trk] < 1.0 && fabs(Track_Theta[a_Trk]-1.57)< 0.4) || (fabs(Track_Theta[a_Trk]-1.57) >= 0.4 && log10(Track_Energy[a_Trk]) < -(fabs(Track_Theta[a_Trk]-1.57)-0.4)*0.2/0.3 ))
			{
				TrackType = 100;
			}
			else if( fabs(EndPointPos.Z()) > ECALHalfZ - 500 && EndPointPos.Perp() > TPCOuterRadius - 300  )	//Shoulder
			{
				TrackType = 30;
			}
			else if( fabs(EndPointPos.Z()) > LStar - 500 && EndPointPos.Perp() < TPCInnerRadius )		//Forward
			{
				TrackType = 40;
			}
			else if( EndPointPos.Perp() > TPCOuterRadius - 100 )		//Barrel
			{
				TrackType = 10;
			}
			else if( fabs(EndPointPos.Z()) > ECALHalfZ - 200 )		//Endcap
			{
				TrackType = 20; 
			}

			if( fabs(D0) < 1 && fabs(Z0) < 1 )
			{
				TrackType += 1;
			}

			Track_Type[a_Trk] = TrackType; 

			if(TrackType == 11)
				tracks_HQ_Barrel.push_back(a_Trk);
			else if(TrackType == 21)
				tracks_HQ_Endcap.push_back(a_Trk);
			else if(TrackType == 31)
				tracks_HQ_Shoulder.push_back(a_Trk);
			else if(TrackType == 41)
				tracks_HQ_Forward.push_back(a_Trk);
			else if(TrackType == 10)
				tracks_MQ_Barrel.push_back(a_Trk);
			else if(TrackType == 20)
				tracks_MQ_Endcap.push_back(a_Trk);
			else if(TrackType == 30)
				tracks_MQ_Shoulder.push_back(a_Trk);
			else if(TrackType == 40)
				tracks_MQ_Forward.push_back(a_Trk);
			else if(TrackType == 1)
				tracks_Vtx.push_back(a_Trk);
			else if(TrackType == 101)
				tracks_LE.push_back(a_Trk);
			else if( (StartPointPos.Mag() > 50 && EndPointPos.Mag() < 1000 && NTrkHit < 50) || TrackType == 100  )
				tracks_ILL.push_back(a_Trk);
			else
				tracks_LQ.push_back(a_Trk);
		}
	}

	std::vector<float > currTrkMomentum;
	std::vector<int> currTrkIndex;

	for(int t1 = 0; t1 < 11; t1++)
	{
		currTrkMomentum.clear();
		currTrkIndex.clear();
		curr_tracks.clear();
		if(t1 == 0)
			curr_tracks = tracks_HQ_Endcap;
		else if(t1 == 1)
			curr_tracks = tracks_HQ_Barrel;
		else if(t1 == 2)
			curr_tracks = tracks_MQ_Endcap;
		else if(t1 == 3)
			curr_tracks = tracks_MQ_Barrel;
		else if(t1 == 4)
			curr_tracks = tracks_HQ_Shoulder;
		else if(t1 == 5)
			curr_tracks = tracks_MQ_Shoulder;
		else if(t1 == 6)
			curr_tracks = tracks_HQ_Forward;
		else if(t1 == 7)
			curr_tracks = tracks_MQ_Forward;
		else if(t1 == 8)
			curr_tracks = tracks_Vtx;
		else if(t1 == 9)			
			curr_tracks = tracks_LQ; 
		else if(t1 == 10)			
			curr_tracks = tracks_LE; 


		int N_currTrack = curr_tracks.size();

		for(int t2 = 0; t2 < N_currTrack; t2++)
		{
			Track* tmpTrk = curr_tracks[t2];
			currTrkMomentum.push_back(Track_Energy[tmpTrk]);
		}

		currTrkIndex = SortMeasure(currTrkMomentum, 1);

		for(int t3 = 0; t3 < N_currTrack; t3++)
		{
			Track* b_tmpTrk = curr_tracks[currTrkIndex[t3]];
			if(t1 < 9 || Track_Energy[b_tmpTrk] < 10)
				TrackOrder.push_back(Track_Index[b_tmpTrk]);
		}
	}

	for(unsigned int t4 = 0; t4 < TrackOrder.size(); t4++)
	{
		Track* b_Trk = dynamic_cast<Track*>( TrackColl->getElementAt( TrackOrder[t4] ) );
		SortedTracks.push_back(b_Trk);
	}
}

void CheckConsistence(LCEvent *evtQQ, std::vector<std::string> SubClass)
{
	float inputEn = 0;
	int inputNCl = 0;
	float TotalEn = 0;
	int TotalNCl = 0;

	for(int c1 = 0; c1 < int(SubClass.size()); c1++)
	{

		LCCollection * LocalCol = evtQQ->getCollection( SubClass[c1].c_str() );

		int CurrCollSize = LocalCol->getNumberOfElements();
		float CurrEnergy = 0;

		for(int c2 = 0; c2 < CurrCollSize; c2++)
		{
			Cluster * a_clu = dynamic_cast<Cluster*>(LocalCol->getElementAt(c2));
			CurrEnergy += a_clu->getEnergy();
		}

		if(c1 < 2)
		{
			inputEn += CurrEnergy;
			inputNCl += CurrCollSize;
		}
		else
		{
			TotalNCl += CurrCollSize;
			TotalEn += CurrEnergy;
		}

		cout<<"CurrColl "<<SubClass[c1].c_str()<<" has "<<CurrCollSize<<" Clusters with total En = "<<CurrEnergy<<" GeV"<<endl;
	}

	cout<<endl<<"Comparison Energy "<<inputEn<<" =? "<<TotalEn<<endl;
	cout<<"Comparison Cluster Size "<<inputNCl<<" =? "<<TotalNCl<<endl<<endl;
}

void BushConnect::BushSelfMerge(LCEvent * evtPP)
{
	LCCollection * CaloClu = evtPP->getCollection("EHBushes");      //A sort here should be helpful
	int NClu = CaloClu->getNumberOfElements();

	std::vector<Cluster* > Core_1st; 
	std::vector<Cluster* > Frag_1st;
	std::vector<Cluster* > UnId_1st; 
	Core_1st.clear();
	Frag_1st.clear();
	UnId_1st.clear();

	float CluDepth = 0; 
	float CluEn = 0;
	int CluSize = 0; 
	TVector3 PosCluSeed, PosSeedDiff, PosSeedA, PosSeedB; 

	// Maybe first ID? Should be helpful

	float TotalCluEn = 0;
	float TotalCluEn_1stAB = 0;

	int NJoints = 0; 	
	int SmallCluSize = 0; 
	float Depth_A = 0; 
	float Depth_B = 0;
	int Size_A = 0; 
	int Size_B = 0; 

	TMatrixF FlagMerge(NClu, NClu);

	for(int i0 = 0; i0 < NClu; i0++)
	{
		Cluster * a_clu = dynamic_cast<Cluster*>(CaloClu->getElementAt(i0));
		float currCluFD = FDV3(a_clu, ECALCellIDDecoder);
		CluFD[a_clu] = currCluFD;
	}

	for(int s0 = 0; s0 < NClu; s0++)
	{
		Cluster * a_clu = dynamic_cast<Cluster*>(CaloClu->getElementAt(s0));
		PosSeedA = a_clu->getPosition();
		Depth_A = DisSeedSurface(PosSeedA);
		Size_A = a_clu->getCalorimeterHits().size();

		for(int s1 = s0 + 1; s1 < NClu; s1++)
		{
			Cluster *b_clu = dynamic_cast<Cluster*>(CaloClu->getElementAt(s1));
			NJoints = JointsBetweenBush(a_clu, b_clu, 4);
			PosSeedB = b_clu->getPosition();
			Depth_B = DisSeedSurface(PosSeedB);
			PosSeedDiff = PosSeedA - PosSeedB;
			Size_B = a_clu->getCalorimeterHits().size();
			float DeeperDepth = std::max(Depth_A, Depth_B);
			if(NJoints && PosSeedDiff.Perp() < 120 + 0.05*DeeperDepth )	//And depth...
			{
				SmallCluSize = std::min( Size_A, Size_B );

				if( ( ( NJoints > 4 || (NJoints > 1 && SmallCluSize < 10) ) && DeeperDepth > 30 ) || NJoints > 8 )
				{	
					FlagMerge[s0][s1] = 1.0;
					FlagMerge[s1][s0] = 1.0;
				}
			}
			//Head Tail Connection. Could be more sophsticate && should be very strict.
			if( PosSeedA.Angle(PosSeedB) < 0.1 && PosSeedDiff.Mag() < 1000 && PosSeedDiff.Mag()*PosSeedA.Angle(PosSeedB) < 60 + 0.02*DeeperDepth && ((CluFD[a_clu] < 0.2 && Size_A > 6) || (CluFD[b_clu] < 0.2 && Size_B > 6)) )
			{
				if( (PosSeedA.Mag() > PosSeedB.Mag() && PosSeedA.Angle(PosSeedB - PosSeedA) < 0.2) || (PosSeedB.Mag() > PosSeedA.Mag() && PosSeedA.Angle(PosSeedA - PosSeedB) < 0.2) )
				{
					FlagMerge[s0][s1] = 2.0;
					FlagMerge[s1][s0] = 2.0;
					// cout<<"tail found"<<endl; 
				}
			}

			/*			
			//Check if Head-Tail Connection
			if( CluFD[a_clu] < 0.2 || CluFD[b_clu] < 0.2 )	// Applied 
			{
				std::pair<TVector3, TVector3> PointPair = ClosestPointPair(a_clu, b_clu);
				// cout<<"CLOSESSSSSSSSSSSSSSSSSSSS Pair"<< (PointPair.first - PointPair.second).Mag() <<endl; 
				// cout<<PointPair.first.X()<<" : "<<PointPair.first.Y()<<" : "<<PointPair.first.Z()<<endl;
				// cout<<PointPair.second.X()<<" : "<<PointPair.second.Y()<<" : "<<PointPair.second.Z()<<endl;

				TVector3 CoGDeep, CoGShallow;
				if( PointPair.first.Mag() < PointPair.second.Mag() )
				{
					CoGDeep = ClusterCoG(b_clu);
					CoGShallow = ClusterCoG(a_clu);
					// cout<<"AB " <<(CoGDeep - PointPair.second).Angle(PointPair.second - PointPair.first)<<endl; 
				}
				else
				{
					CoGDeep = ClusterCoG(a_clu);
                                        CoGShallow = ClusterCoG(b_clu);
                                        // cout<<"BA " <<(CoGDeep - PointPair.first).Angle(PointPair.first - PointPair.second)<<endl;
				}
			}
			*/
		}
		TotalCluEn += a_clu->getEnergy();
	}

	std::vector<Cluster*> OriInputEHBushes = CollClusterVec(CaloClu);
	TMatrixF MergeSYM = MatrixSummarize(FlagMerge);
	LCCollection* CloseMergedCaloClu = ClusterVecMerge( OriInputEHBushes, MergeSYM);

	// 1st iteration absorbtion: large energy threshold, and small merge region, keep purity

	std::map<Cluster*,float> MinDisSeedToBush;
	MinDisSeedToBush.clear();
	for(int i0 = 0; i0 < CloseMergedCaloClu->getNumberOfElements(); i0++)
	{
		Cluster * a_clu = dynamic_cast<Cluster*>(CloseMergedCaloClu->getElementAt(i0));
		PosCluSeed = a_clu->getPosition();
		float tmpmindis = 1e10;
		for(int i1 = 0; i1 < CloseMergedCaloClu->getNumberOfElements(); i1++)
		{
			Cluster * b_clu = dynamic_cast<Cluster*>(CloseMergedCaloClu->getElementAt(i1));
			if(i1 != i0)
			{
				if(DisPointToBush(PosCluSeed,b_clu) < tmpmindis) tmpmindis = DisPointToBush(PosCluSeed,b_clu);  
			}
		}
		MinDisSeedToBush[a_clu] = tmpmindis;
	}

	for(int i0 = 0; i0 < CloseMergedCaloClu->getNumberOfElements(); i0++)
	{
		Cluster * a_clu = dynamic_cast<Cluster*>(CloseMergedCaloClu->getElementAt(i0));
		PosCluSeed = a_clu->getPosition();
		CluDepth = DisSeedSurface(PosCluSeed);
		CluEn = a_clu->getEnergy();
		CluSize = a_clu->getCalorimeterHits().size();

		if( CluEn > 2.0 + 0.002*CluDepth || (CluSize > 10 && CluDepth < 20) || CluEn > 5.0)
		{
			Core_1st.push_back(a_clu);
		}
		else if( (CluSize > 10 || CluEn > 0.2) && MinDisSeedToBush[a_clu] > 20 && (PhotonTag(a_clu) == 1 || ClusterFlag1st(a_clu) == 11  ))
		{
			Core_1st.push_back(a_clu);
		}
		else if( CluSize < 5 && CluEn < 0.3 && CluDepth > 40 )
		{
			Frag_1st.push_back(a_clu);
		}
		else
		{
			UnId_1st.push_back(a_clu);
		}
	}

	std::vector<Cluster* > UndefFrag_1stAB = ClusterAbsorbtion(UnId_1st, Frag_1st, 50, 0.02);
	std::vector<Cluster* > CoreFrag_1stAB = ClusterAbsorbtion(Core_1st, UndefFrag_1stAB, 50, 0.02);	


	float MaxCluEn = 0;

	for(int s = 0; s < int(CoreFrag_1stAB.size()); s++)
	{
		Cluster * a_clu = CoreFrag_1stAB[s];
		TotalCluEn_1stAB += a_clu->getEnergy();
		if(MaxCluEn < a_clu->getEnergy())
			MaxCluEn = a_clu->getEnergy();
	}
	

	// 2nd iteration, considering the cores

	std::vector<Cluster* > Core_2nd;
	std::vector<Cluster* > Frag_2nd;
	std::vector<Cluster* > UnId_2nd;
	Core_2nd.clear();
	Frag_2nd.clear();
	UnId_2nd.clear();

	std::map<Cluster*,float> MinDisSeedToBush2;
	MinDisSeedToBush2.clear();
	for(int i0 = 0; i0 < int(CoreFrag_1stAB.size()); i0++)
	{
		Cluster * a_clu = CoreFrag_1stAB[i0];
		PosCluSeed = a_clu->getPosition();
		float tmpmindis = 1e10;
		for(int i1 = 0; i1 < int(CoreFrag_1stAB.size()); i1++)
		{
			Cluster * b_clu = CoreFrag_1stAB[i1];
			if(i1 != i0)
			{
				if(DisPointToBush(PosCluSeed,b_clu) < tmpmindis) tmpmindis = DisPointToBush(PosCluSeed,b_clu);
			}
		}
		MinDisSeedToBush2[a_clu] = tmpmindis;
	}

	for(int i2 = 0; i2 < int(CoreFrag_1stAB.size()); i2++)
	{
		Cluster * a_clu = CoreFrag_1stAB[i2];
		PosCluSeed = a_clu->getPosition();
		CluDepth = DisSeedSurface(PosCluSeed);
		CluEn = a_clu->getEnergy();
		CluSize = a_clu->getCalorimeterHits().size();
		if( CluEn > 1.5 + 0.003*CluDepth || ( (PhotonTag(a_clu) == 1 || ClusterFlag1st(a_clu) == 11) && (CluSize > 10 || CluEn > 0.2) && MinDisSeedToBush2[a_clu] > 20 && CluDepth < 25) )
		{
			Core_2nd.push_back(a_clu);
		}
		else
		{
			UnId_2nd.push_back(a_clu);
		}
	}	

	float MinDisToCore = 1.0E10; 
	float tmpUnIdCoreDis = 1.0E10; 

	for(int i3 = 0; i3 < int(UnId_2nd.size()); i3++)
	{
		Cluster * a_unId = UnId_2nd[i3];
		CluEn = a_unId->getEnergy();
		PosCluSeed = a_unId->getPosition();
		CluDepth = DisSeedSurface(PosCluSeed);
		CluSize = a_unId->getCalorimeterHits().size();
		MinDisToCore = 1.0E10;
		tmpUnIdCoreDis = 1.0E10;

		for(int j3 = 0; j3 < int(Core_2nd.size()); j3++)
		{
			Cluster *a_core = Core_2nd[j3];
			tmpUnIdCoreDis = BushDis(a_core, a_unId);	//Maybe Giveback the vector is even better. Surely!
			if( tmpUnIdCoreDis < MinDisToCore )
			{
				MinDisToCore = tmpUnIdCoreDis;
			}			
		}
		if( (1 - 0.0006*CluDepth) * CluEn * MinDisToCore > 30 ) //|| (CluEn > 1.5 && CluDepth < 100) )	//function and para need to be optimized, Projective distance should make more sense
		{	
			Core_2nd.push_back(a_unId);
		}
		else
		{
			Frag_2nd.push_back(a_unId);	//Should be merged to the minimal index one...
		}
	}

	std::vector<Cluster* > CoreFrag_2ndAB = ClusterAbsorbtion(Core_2nd, Frag_2nd, 50, 0.01);

	LCCollection *SMBush_2nd = ClusterVecColl(CoreFrag_2ndAB);
	evtPP->addCollection(SMBush_2nd, "SMBush_2ndIt");

	//Enable the Shoulder Merge. 
	
	std::vector<Cluster* > CluCat1;
        std::vector<Cluster* > CluCat2;
        std::vector<Cluster* > CluCat3;
	std::vector<Cluster* > CluCat4;
        CluCat1.clear();
        CluCat2.clear();
        CluCat3.clear();
	CluCat4.clear();

	int nClu2edIt = CoreFrag_2ndAB.size();
	for (int i3 = 0; i3 < nClu2edIt;i3++)
	{
		Cluster *a_clu = CoreFrag_2ndAB[i3];
		TVector3 CluPos = a_clu->getPosition();
		if(fabs(CluPos[2]) > 2645 && CluPos.Perp() > 2000 ) 
		{
			CluCat1.push_back(a_clu);
		}
		else if(fabs(CluPos[2]) > 2400 && CluPos.Perp() > 2000 ) 
		{
			CluCat2.push_back(a_clu);
		}
		else if(fabs(CluPos[2]) > 2000 && fabs(CluPos[2]) < 2400) 
		{
			CluCat3.push_back(a_clu);
		}
		else
		{
			CluCat4.push_back(a_clu);
		}
	}

	std::vector<Cluster* > CluCat12 = ClusterAbsorbtion(CluCat2, CluCat1, 30, 0);
	std::vector<Cluster* > CluCat123 = ClusterAbsorbtion(CluCat3, CluCat12, 200, 0);
	std::vector<Cluster* > CluShoulderMerged = ClusterAbsorbtion(CluCat4, CluCat123, 0, 0);

	LCCollection *SMBush_3rd = ClusterVecColl(CluShoulderMerged);
	evtPP->addCollection(SMBush_3rd, "SMBush_3rdIt");

	SortedSMBushes = CollClusterVec(SMBush_3rd);
	//Set ID...		
}

void BushConnect::TagCore(LCEvent * evtPP) 
{
	LCCollection * ECALClu = evtPP->getCollection("SMBush_3rdIt");	//A sort here should be helpful
	int NClu_Ecal = ECALClu->getNumberOfElements();
	int NTrk = SortedTracks.size();
	std::vector<Cluster* > ecalbushes = CollClusterVec(ECALClu);
	TVector3 CluPos;

	std::map<Cluster*, int> BushTouchFlag; 
	std::map<Track*, int> TrkTouchFlag; 
	std::map<Track*, Cluster*> FCMap_Track_ABSCluster; 	//Careful about 1 - 1 on to...
	std::map<Track*, float> Map_Track_ClusterEnergy; 
	std::map<Track*, Cluster*> FCMap_Track_CHCore;
	BushTouchFlag.clear();
	TrkTouchFlag.clear();
	FCMap_Track_ABSCluster.clear();
	FCMap_Track_CHCore.clear();
	Map_Track_ClusterEnergy.clear();

	float currTrkEn = 0; 
	int currTrackType = 0;
	//float MinimalTrkBushDis_E = 1.0E9;
	float DisMatrix_Track_Clu_E[NTrk][NClu_Ecal];
	float TimeMatrix_Track_Clu_E[NTrk][NClu_Ecal];
	float CluDepth = 0; 

	LCCollection *chcorecluster = new LCCollectionVec(LCIO::CLUSTER);
	chcorecluster->setFlag(Cluflag.getFlag());

	for(int s0 = 0; s0 < NTrk; s0++)
	{
		for(int s1 = 0; s1 < NClu_Ecal; s1++)
		{
			DisMatrix_Track_Clu_E[s0][s1] = 1.0E10;
			TimeMatrix_Track_Clu_E[s0][s1] = 1.0E10; 
		}
	}

	std::vector<Cluster*> TightLinkedCluster; 

	//float EcalCoreEnergy = 0; 
	float CoreMergeDistanceDepthCorrector = 0; 
	float MinDisToEntrance = 1.0E10; 

	TVector3 TrkEndPoint(0, 0, 0);
	static float TrkEndP[3] = {0,0,0};
	static float TrkEndM[3] = {0,0,0};
	TVector3 RefECALEntrance(0, 0, 0);

//~~~~~~~ find the closest cluster first...

	std::map<int, int> Closest_Trk_Clu_Map;
	Closest_Trk_Clu_Map.clear();

	for(int g0 = 0; g0 < NTrk; g0++)
	{
		Track* a_trk = SortedTracks[g0];
		float ClosestDis = 1.0E9;  
		int ClosestCluIndex = -1; 
		int ClosestNC = 1E9;
		float ThetaDiff = 0;
		TrkEndPoint = trkendposition[a_trk];
		TrkEndP[0]=TrkEndPoint(0);
		TrkEndP[1]=TrkEndPoint(1);
		TrkEndP[2]=TrkEndPoint(2);

		HelixClass * trkHelix = new HelixClass();
		trkHelix->Initialize_Canonical(a_trk->getPhi(), a_trk -> getD0(), a_trk -> getZ0(), a_trk -> getOmega(), a_trk->getTanLambda(), BField);
		trkHelix->getExtrapolatedMomentum(TrkEndP,TrkEndM);
		TVector3 TrkEndMomentum(TrkEndM);
		
		currTrackType = Track_Type[a_trk];
		TVector3 TrkP3 = Track_P3[a_trk];
		
		for(int g1 = 0; g1 < NClu_Ecal; g1++)
		{
			Cluster *fccand_bush = ecalbushes[g1];
			float* Dis = SimpleDisTrackClu(a_trk, fccand_bush);
			float Time = SimpleBushTimeTrackClu(a_trk, fccand_bush);
			int NC = SimpleBushNC(a_trk, fccand_bush);
			TVector3 CluPos = fccand_bush->getPosition();
			//ThetaDiff = abs(TrkEndMomentum.Theta() - (CluPos - TrkEndPoint).Theta());
			ThetaDiff=TrkEndMomentum.Angle(CluPos - TrkEndPoint);
			if(Dis[2] > -0.1)
			{
				DisMatrix_Track_Clu_E[g0][g1] = Dis[2];
				TimeMatrix_Track_Clu_E[g0][g1] = Time;
				if( Dis[2] < ClosestDis && ThetaDiff < 0.05)
				{
					ClosestDis = Dis[2]; 
					ClosestCluIndex = g1;
					ClosestNC = NC;
				}
			}
		}

		// if( (ClosestNC < 3 || ClosestDis < 5) && abs(TrkP3.Theta() - 1.57) < 0.01 )
		// if(ClosestDis < 10 && ClosestCluIndex > -0.1 && log10(ClosestTime) < 3.5 && ClosestTime > 0)	// to be replace by Binsong's poly-4 curve
		if( ClosestDis < 10 && ClosestCluIndex > -0.1 && (ClosestNC < 3 || abs(TrkP3.Theta() - 1.57) < 0.01 ) ) 
		{
			Cluster * candiclu = ecalbushes[ClosestCluIndex];
			TVector3 CluPos = candiclu->getPosition();
			float TrackEndPDis = (TrkEndPoint - CluPos).Mag();
			float AngDiff = TrkEndPoint.Angle(CluPos-TrkEndPoint);
			if(TrackEndPDis < 400 && AngDiff < 0.4 && (fabs(Track_Energy[a_trk] - candiclu->getEnergy() ) < 3.0*sqrt(Track_Energy[a_trk]) + 1.0 ||  candiclu->getEnergy() < 8) )
			{
				Closest_Trk_Clu_Map[g0] = ClosestCluIndex;
				BushTouchFlag[candiclu] = 0;
			
			}
		}
	}

	//~~~~~~~ end of finding close cluster

	for(int i0 = 0; i0 < NTrk; i0++)  //Dropped Size can exist
	{
		Track* a_trk = SortedTracks[i0];
		currTrackType = Track_Type[a_trk];
		currTrkEn = Track_Energy[a_trk];

		// EcalCoreEnergy = 0;
		// MinimalTrkBushDis_E = 1.0E9;
		// CluIndex = -1;

		TrkEndPoint = trkendposition[a_trk];
		MinDisToEntrance = 1.0E10; 
		TightLinkedCluster.clear();
		float fccanden = 0; 
		if( Closest_Trk_Clu_Map.find(i0) != Closest_Trk_Clu_Map.end() )
		{
			Cluster * closeClu = ecalbushes[Closest_Trk_Clu_Map[i0]];
			TightLinkedCluster.push_back(closeClu);
		}

		for(int j0 = 0; j0 < NClu_Ecal; j0++)
		{
			Cluster *fccand_bush = ecalbushes[j0];			
			float Dis = DisMatrix_Track_Clu_E[i0][j0]; //SimpleDisTrackClu(a_trk, fccand_bush);
			float BushTime = TimeMatrix_Track_Clu_E[i0][j0];
			CluPos = fccand_bush->getPosition();
			float DisToEntrance = (CluPos - RefECALEntrance).Mag();			
			CluDepth = DisSeedSurface(CluPos);
			int currCluType = ClusterFlag1st(fccand_bush);

			//if(currCluType == 11 && Dis > 5) continue;
			
			CoreMergeDistanceDepthCorrector = 0;
			if(CluDepth > 20)
				CoreMergeDistanceDepthCorrector = 20;
			else if(CluDepth > 10)
				CoreMergeDistanceDepthCorrector = 10;

			float TrackEndPDis = (TrkEndPoint - CluPos).Mag();
			//cout<<BushTime<<" "<<currTrackType<<" "<<Dis<<" "<<CoreMergeDistanceDepthCorrector<<" "<<currTrkEn<<" "<<fccand_bush->getEnergy()<<" "<<TrackEndPDis<<endl;
			if(log10(BushTime) < 3.5 && BushTime > 0 && currTrackType != 101 && Dis < 7 + CoreMergeDistanceDepthCorrector && Dis > -0.1 && BushTouchFlag.find(fccand_bush) == BushTouchFlag.end() && (currTrkEn > 3 || fccand_bush->getEnergy() < 5 || currCluType == 13 ) && (TrackEndPDis < 400 || currTrackType != 11))
			{
				TightLinkedCluster.push_back(fccand_bush);
				BushTouchFlag[fccand_bush] = currTrackType;
				fccanden += fccand_bush->getEnergy();
			}
			if(DisToEntrance < MinDisToEntrance)
			{
				MinDisToEntrance = DisToEntrance;
			}
		}
		if( TightLinkedCluster.size() > 0 ) // && EcalCoreEnergy + HcalCoreEnergy < 2.0*currTrkEn )...
		{
			ClusterImpl * chcorecluster_eh =  NaiveMergeClu(TightLinkedCluster);
			chcorecluster->addElement(chcorecluster_eh);
			FCMap_Track_CHCore[a_trk] = chcorecluster_eh;
			Map_Track_ClusterEnergy[a_trk] = chcorecluster_eh->getEnergy();
			chargedclustercore.push_back(chcorecluster_eh);
		}
	}

	evtPP->addCollection(chcorecluster, "CHARGEDCORE_BeforeABS");

	float DisToClosestTrack = 1.0E10;
	//int CluSize = 0; 
	int ClosetTrackIndex = -1;
	float CluEnergy = 0; 
	float DepthShiftTBM = 0; 
	float DisToCore = 0; 
	int FlagBackScattering = 0; 

	for(int i1 = 0; i1 < NClu_Ecal; i1++)	// Could be replaced by CluEnergy, CluType, CluSize Maps. 
	{
		Cluster *a_bush = ecalbushes[i1];

		if( BushTouchFlag.find(a_bush) == BushTouchFlag.end() )	//Might be a neutral core
		{
			DisToClosestTrack = 1.0E10;
			DisToCore = 1.0E10; 
			CluPos = a_bush->getPosition();
			CluDepth = DisSeedSurface(CluPos);
			CluEnergy = a_bush->getEnergy();
			//CluSize = a_bush->getCalorimeterHits().size();
			ClosetTrackIndex = -1;
			FlagBackScattering = 0;	

			for(int j1 = 0; j1 < NTrk; j1++)
			{
				if( DisMatrix_Track_Clu_E[j1][i1] < DisToClosestTrack )
				{
					DisToClosestTrack = DisMatrix_Track_Clu_E[j1][i1];
					ClosetTrackIndex = j1; 
				}
			}

			if(DisToClosestTrack < 100 + CluDepth)
			{
				Cluster *currcore = FCMap_Track_CHCore[ SortedTracks[ClosetTrackIndex] ];
				if(currcore)
				{
					int Currcoresize = currcore->getCalorimeterHits().size();
					DisToCore = BushDis(currcore, a_bush);
					if(DisToCore < 7 && CluDepth < 10 && Currcoresize > 20 && PhotonTag(a_bush) == 0 && ClusterFlag1st(a_bush) != 11 )	//Cell Size, close attached
					{
						FlagBackScattering = 1; 	//Create a dedicated collection; then
						ecalpotentialbackscattering.push_back(a_bush);
					}	
				}
			}

			if(!FlagBackScattering)	//Only Potential
			{
				if(CluEnergy < 0.05 && CluDepth > 30)
				{
					ecalfrag.push_back(a_bush);
				}
				else if( (DisToClosestTrack > 15 && DisToCore > 10) || CluEnergy > 2.0 )     //FD ~ see if ... should enlarge this distances...
				{
					if( ClusterFlag1st(a_bush) == 11 || PhotonTag(a_bush) == 1 )
					{
						ecalnecore_EM.push_back(a_bush);
					}
					else if(CluEnergy > 0.5)
					{
						ecalnecore_NonEM.push_back(a_bush);
					}
					else
					{
						ecalundef.push_back(a_bush);
					}
				}
				else if(DisToClosestTrack < 60 + CluDepth)	//T.B.M	
				{
					if(CluDepth > 40) 
						DepthShiftTBM = 20;
					else 
						DepthShiftTBM = 0.5*CluDepth;					

					if( DisToCore < 7 + DepthShiftTBM )	//Neighbours?... ~ same level
					{
						ecalfrag_TBM_CH.push_back(a_bush);		//Direct Merging
					}
					else	// Potentially loose Core
					{
						ecalundef.push_back(a_bush);
					}
				}
				else
				{
					ecalundef.push_back(a_bush);
				}
			}		
		}
	}

	int NUndef = ecalundef.size();
	int NNeCoreEM = ecalnecore_EM.size();
	int NNeCoreNonEM = ecalnecore_NonEM.size();
	int NChCore = chcorecluster->getNumberOfElements();
	//int NCore = NChCore + NNeCoreEM + NNeCoreNonEM;
	//float DisMatrix[NUndef][NCore];

	float tmpDis = 0; 
	float MinDistanceToChCore = 1.0E10; 
	float MinDistanceToNeCore_EM = 1.0E10; 
	float MinDistanceToNeCore_NonEM = 1.0E10; 

	for(int i2 = 0; i2 < NUndef; i2++)
	{
		Cluster *a_Undef = ecalundef[i2];
		MinDistanceToChCore = 1.0E10;
		MinDistanceToNeCore_EM = 1.0E10;
		MinDistanceToNeCore_NonEM = 1.0E10; 

		CluPos = a_Undef->getPosition();
		CluDepth = DisSeedSurface(CluPos);

		for(int j2 = 0; j2 < NChCore; j2++)
		{
			Cluster *a_chcore = dynamic_cast<Cluster*>(chcorecluster->getElementAt(j2));
			tmpDis = BushDis(a_chcore, a_Undef);
			//DisMatrix[i2][j2] = BushDis(a_chcore, a_Undef);
			if(tmpDis < MinDistanceToChCore)
			{
				MinDistanceToChCore = tmpDis;
			}
		}
		for(int j3 = 0; j3 < NNeCoreEM; j3++)
		{
			Cluster *a_necore = ecalnecore_EM[j3];
			tmpDis = BushDis(a_necore, a_Undef);
			//DisMatrix[i2][j3 + NChCore] = BushDis(a_necore, a_Undef);
			if(tmpDis < MinDistanceToNeCore_EM)
			{
				MinDistanceToNeCore_EM = tmpDis;
			}
		}
		for(int j4 = 0; j4 < NNeCoreNonEM; j4++)
		{
			Cluster *a_necore = ecalnecore_NonEM[j4];
			tmpDis = BushDis(a_necore, a_Undef);
			//DisMatrix[i2][j4 + NChCore + NNeCoreEM] = BushDis(a_necore, a_Undef);
			if(tmpDis < MinDistanceToNeCore_NonEM)
			{
				MinDistanceToNeCore_NonEM = tmpDis;
			}
		}

		if(MinDistanceToChCore < 20 && MinDistanceToChCore < MinDistanceToNeCore_EM && MinDistanceToChCore < MinDistanceToNeCore_NonEM  )
		{
			ecalfrag_TBM_CH.push_back(a_Undef);
		}
		else if(MinDistanceToNeCore_EM < 40 && MinDistanceToNeCore_EM < MinDistanceToChCore && MinDistanceToNeCore_EM < MinDistanceToNeCore_NonEM)	//Distance should be depends on Energy & Depth As well
		{
			ecalfrag_TBM_NE_EM.push_back(a_Undef);
		}
		else if(MinDistanceToNeCore_NonEM < 40 && MinDistanceToNeCore_NonEM < MinDistanceToChCore && MinDistanceToNeCore_NonEM < MinDistanceToNeCore_EM)
		{
			ecalfrag_TBM_NE_NonEM.push_back(a_Undef);
		}
		else
		{
			ecalundef_iso.push_back(a_Undef);
		}
	}

	chargedclustercore_abs = ClusterAbsorbtion(chargedclustercore, ecalfrag_TBM_CH, 1000, 0);

	LCCollection *absorbedChCore = ClusterVecColl(chargedclustercore_abs);
	evtPP->addCollection(absorbedChCore, "ChargedCore");
	
	int N_ChartgedCore = chargedclustercore_abs.size();

	std::vector<Cluster* > em_core = ClusterAbsorbtion(ecalnecore_EM, ecalfrag_TBM_NE_EM, 1000, 0);
	std::vector<Cluster* > nonem_core = ClusterAbsorbtion(ecalnecore_NonEM, ecalfrag_TBM_NE_NonEM, 1000, 0);

	LCCollection *absorbedNeCore_EM = ClusterVecColl(em_core);
	LCCollection *absorbedNeCore_NonEM = ClusterVecColl(nonem_core);

	evtPP->addCollection(absorbedNeCore_EM, "EMCore");
	evtPP->addCollection(absorbedNeCore_NonEM, "NonEMCore");

	//Re-create the Map between Tracks and Charged Core: which is stupid...

	LCCollection *arborrecoparticle_ch = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	LCCollection *arborrecoparticle_ne = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	LCCollection *arborchargedcorecluster = new LCCollectionVec(LCIO::CLUSTER);
	LCCollection *arborneutralcorecluster = new LCCollectionVec(LCIO::CLUSTER);
	arborchargedcorecluster->setFlag(Cluflag.getFlag());
	arborneutralcorecluster->setFlag(Cluflag.getFlag());

	//Non_Core_Neutral_Cluster	Push...
	//Core_Neutral_Cluster		Push...
	//
	//	Tight Absorb Coll
	//	Loose Absorb Coargedclustercore_absl

	LCCollection * col_TPCTrk =  evtPP->getCollection("ClupatraTracks");

	for(int j5 = 0; j5 < NTrk; j5++)
	{
		Track* a_trk = SortedTracks[j5];

		int Track_Core_ID = 0;

		if(Map_Track_ClusterEnergy[a_trk] > 0)		// If already a ... trk	
		{
			currTrackType = Track_Type[a_trk];
			currTrkEn = Track_Energy[a_trk];
			float MinTrkCluDis = 1.0E10;
			int CloseCluIndex = -1;
			//float CoreClusterEnergy = 0; 

			for(int k5 = 0; k5 < N_ChartgedCore; k5++)
			{
				Cluster * a_chargedcore_abs = chargedclustercore_abs[k5];
				float* Dis=new float[3];
				Dis[0]=0;
				Dis[1]=0;
				Dis[2]=0;

				
				Dis = SimpleDisTrackClu(a_trk, a_chargedcore_abs);	
				if(Dis[2] < MinTrkCluDis && Dis[2] > -0.1 && Map_Track_ClusterEnergy[a_trk] <= a_chargedcore_abs->getEnergy() )   //  ) && a_chargedcore_abs->getEnergy() < currTrkEn + 2*sqrt(currTrkEn) + 2)
				{
					MinTrkCluDis = Dis[2]; 
					CloseCluIndex = k5; 
				}
				//delete Dis;
			}
			if(CloseCluIndex > -1)
			{
				Cluster *a_close_tree = chargedclustercore_abs[CloseCluIndex];

				Track_Core_ID = ClusterFlag(a_close_tree, a_trk, col_TPCTrk);
				if(BushTouchFlag.find(a_close_tree) == BushTouchFlag.end())
				{
					BushTouchFlag[a_close_tree] = currTrackType;
					FCMap_Track_ABSCluster[a_trk] = a_close_tree;
					//CoreClusterEnergy = a_close_tree->getEnergy();
				}
				
			}
		}
		
		//Shall I veto some crazy track...

		ReconstructedParticleImpl * chargeparticle = new ReconstructedParticleImpl();
		chargeparticle->setEnergy( Track_Energy[a_trk] );
		chargeparticle->setCharge(a_trk -> getOmega()/fabs(a_trk -> getOmega()));
		HelixClass * currHelix = new HelixClass();
		currHelix->Initialize_Canonical(a_trk->getPhi(), a_trk -> getD0(), a_trk -> getZ0(), a_trk -> getOmega(), a_trk->getTanLambda(), BField);
		chargeparticle->setMomentum( currHelix->getMomentum() );
		chargeparticle->addTrack( a_trk );
		chargeparticle->setType(Track_Core_ID);
		if( FCMap_Track_ABSCluster[a_trk] )
		{
			ClusterImpl * chargedarborcluster =  NaiveCluImpl(FCMap_Track_ABSCluster[a_trk]);
			chargeparticle->addCluster(chargedarborcluster);
			arborchargedcorecluster->addElement(chargedarborcluster);
		}
		arborrecoparticle_ch->addElement(chargeparticle);
		ChCoreID[chargeparticle] = Track_Core_ID;

	}

	evtPP->addCollection( arborchargedcorecluster, "ClusterChargedCore" );
	evtPP->addCollection( arborrecoparticle_ch, "ArborChargedCore" );

	//~~~~~~~~~~~~~~~~~~~~~
	//Reabsorbtion also defined over here...

	//Cores:
	int NEMCore = em_core.size();
	int NNonEMCore = nonem_core.size();
	int NNeutralCore = NEMCore + NNonEMCore;
	float NAMom[3] = {0, 0, 0};
	TVector3 BushSeedPos; 
	int tmp_ClusterTypeID = 0; 
	int NBKS = ecalpotentialbackscattering.size();
	int NEcalIso = ecalundef_iso.size();
	int NEcalFrag = ecalfrag.size();
	int NFrags = NEcalIso + NEcalFrag + NBKS;

	//Frags:

	for(int j7 = 0; j7 < NFrags; j7++)
	{
		Cluster * a_clu(0);

		if(j7 < NEcalFrag )
		{
			a_clu = ecalfrag[j7];
			tmp_ClusterTypeID = 11;
		}
		else if(j7 < NEcalFrag + NEcalIso)
		{
			a_clu = ecalundef_iso[j7 - NEcalFrag];
			tmp_ClusterTypeID = 12;
		}
		else
		{
			a_clu = ecalpotentialbackscattering[j7 - NEcalFrag - NEcalIso];
			tmp_ClusterTypeID = 14;
		}
		non_chargedclustercore.push_back(a_clu);
		ClusterType_1stID[a_clu] = tmp_ClusterTypeID;
	}

	for(int j6 = 0; j6< NNeutralCore; j6++)	
	{
		Cluster * a_clu(0); 

		if(j6 < NEMCore)
		{
			a_clu = em_core[j6];
			tmp_ClusterTypeID = 1;
		}
		else if(j6 < NEMCore + NNonEMCore)
		{
			a_clu = nonem_core[j6 - NEMCore];
			tmp_ClusterTypeID = 2;
			TVector3 Pos=a_clu->getPosition();

			if( a_clu->getEnergy() > 2.5 && DisSeedSurface(Pos) < 20 )	// PUT BY HAND...
			{
				tmp_ClusterTypeID = 1;
			}
		}

		non_chargedclustercore.push_back(a_clu);
		ClusterType_1stID[a_clu] = tmp_ClusterTypeID; 

		BushSeedPos = a_clu->getPosition();

		ReconstructedParticleImpl * neutralparticle = new ReconstructedParticleImpl();
		neutralparticle->setEnergy( a_clu->getEnergy() );
		neutralparticle->setMass( 0.0 );
		neutralparticle->setCharge( 0.0 );
		NAMom[0] = CluEnergy*1.0/BushSeedPos.Mag()*BushSeedPos.X();
		NAMom[1] = CluEnergy*1.0/BushSeedPos.Mag()*BushSeedPos.Y();
		NAMom[2] = CluEnergy*1.0/BushSeedPos.Mag()*BushSeedPos.Z();
		neutralparticle->setMomentum( NAMom );
		ClusterImpl * neutralarborcluster =  NaiveCluImpl(a_clu);
		neutralparticle->addCluster(neutralarborcluster);
		arborneutralcorecluster->addElement(neutralarborcluster);
		arborrecoparticle_ne->addElement(neutralparticle);
	}

	evtPP->addCollection( arborneutralcorecluster, "ClusterNeutralCore");
	evtPP->addCollection( arborrecoparticle_ne, "ArborNeutralCore");
}

void BushConnect::ParticleReco( LCEvent * evtPP )
{
	LCCollection *arborrecoparticle = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	LCCollection *mergedclu_ch = new LCCollectionVec(LCIO::CLUSTER);
	LCCollection *mergedclu_ne = new LCCollectionVec(LCIO::CLUSTER);
	mergedclu_ch->setFlag(Cluflag.getFlag());
	mergedclu_ne->setFlag(Cluflag.getFlag());

	LCCollection * ChargedCore = evtPP->getCollection("ArborChargedCore");
	LCCollection * col_TPCTrk = evtPP->getCollection("ClupatraTracks");
	LCCollection * col_IsoHit = evtPP->getCollection("AllIsolatedHits");
	std::vector<CalorimeterHit*> IsoHits = CollHitVec(col_IsoHit, 0);

	int NChargedObj = ChargedCore->getNumberOfElements();
	int NNeutralCluster = non_chargedclustercore.size();
	double DisMatrix_Core_Neutral[NChargedObj][NNeutralCluster][2];		//Define different types of distances; 

	float TotalChEn = 0;
	Track * a_chargedTrk(0); 
	Track * a_neighbourTrk(0);
	Cluster * a_chargedClu(0), *a_NeCandiClu(0); 
	float CluDepth = 0;
	std::map<Cluster*, double> CluDepthMap; 
	CluDepthMap.clear();
	int currChargeCoreType = 0;  
	TVector3 CluPos; 

	// Per Track usage...
	std::vector<Cluster*> loosecandicluster; 
	std::vector<Cluster*> tightcandicluster;		//Muon potential candi?
	std::vector<Cluster*> mergedcluster; 			//tmp for each charged P
	std::vector<Cluster*> chargedclustercore_merged; 	//overall

	chargedclustercore_merged.clear();

	std::vector<double> reftightdis; 
	std::vector<double> refloosedis; 

	std::map<Cluster*, int> NNCTouchFlag; 
	std::vector<Track*> SecondIterTracks;
	SecondIterTracks.clear();

	TVector3 currTrkEnd, neighbourTrkEnd, LeadP; 

	for(int i = 0; i < NChargedObj; i++)
	{
		ReconstructedParticle * a_recoP_ch = dynamic_cast<ReconstructedParticle*>(ChargedCore->getElementAt(i));

		loosecandicluster.clear();
		tightcandicluster.clear();
		mergedcluster.clear();
		reftightdis.clear();
		refloosedis.clear();
		a_chargedTrk = a_recoP_ch->getTracks()[0];
		currTrkEnd = trkendposition[a_chargedTrk];
		currChargeCoreType = ChCoreID[a_recoP_ch];
		int currTrkType = Track_Type[a_chargedTrk];

		float CurrClusterEnergy = 0;
		float CurrTrackEnergy = Track_Energy[a_chargedTrk];
		if(a_recoP_ch->getClusters().size() != 0)
		{
			a_chargedClu = a_recoP_ch->getClusters()[0];
			CurrClusterEnergy = a_chargedClu->getEnergy();
			mergedcluster.push_back(a_chargedClu);		//Actually can use this chance to question if previous energy are balance...
		}

		float MinDisToNoClusterTrk = 1.0E10; 
		float MinDisToOtherTrack = 1.0E10;

		for( int is = 0; is < NChargedObj; is++ )
		{
			if(is != i)
			{
				ReconstructedParticle * b_recoP_ch = dynamic_cast<ReconstructedParticle*>(ChargedCore->getElementAt(is));
				a_neighbourTrk = b_recoP_ch->getTracks()[0];
				neighbourTrkEnd = trkendposition[a_neighbourTrk];
				float currDD = (neighbourTrkEnd - currTrkEnd).Mag();
				if( currDD < MinDisToOtherTrack )
				{
					MinDisToOtherTrack = currDD;
				}
			}
		}

		for(int j = 0; j < NNeutralCluster; j++)
		{
			a_NeCandiClu = non_chargedclustercore[j];
			float NeCandEn = a_NeCandiClu->getEnergy(); 
			CluPos = a_NeCandiClu->getPosition();
			CluDepth = DisSeedSurface(CluPos);
			CluDepthMap[a_NeCandiClu] = CluDepth; 	

			if( ClusterType_1stID[a_NeCandiClu] == 1 )   continue; 

			for(int k = 0; k < 2; k++)
			{
				DisMatrix_Core_Neutral[i][j][k] = 1.0E9;
			}

			if(CurrClusterEnergy > 1E-6)	//put by hand...
			{
				DisMatrix_Core_Neutral[i][j][0] = BushDis(a_chargedClu, a_NeCandiClu);
			}
			float* Dis = SimpleDisTrackClu(a_chargedTrk, a_NeCandiClu);
			DisMatrix_Core_Neutral[i][j][1] = Dis[2];

			if( NNCTouchFlag.find(a_NeCandiClu) == NNCTouchFlag.end() && ( currChargeCoreType == 0 || DisMatrix_Core_Neutral[i][j][0] < 1000 ) && currTrkType != 101)
			{			
				if( currChargeCoreType == 130 )			//Matched Muon, should ignore
				{
					if( DisMatrix_Core_Neutral[i][j][1] < 0.2*CluDepth && CluDepth > 200  )	//&& FD?
					{
						tightcandicluster.push_back(a_NeCandiClu);
						reftightdis.push_back( DisMatrix_Core_Neutral[i][j][1] );	//dependence on Cluster Flag & Clu Depth. use some more fancy sort algorithm...
					}
				}
				else if( currChargeCoreType == 131 )
				{
					if( DisMatrix_Core_Neutral[i][j][1] < 0.3*CluDepth && CluDepth > 150 )
					{
						tightcandicluster.push_back(a_NeCandiClu);
						reftightdis.push_back( DisMatrix_Core_Neutral[i][j][1] );	
					}
					else if( DisMatrix_Core_Neutral[i][j][1] < 0.5*CluDepth && CluDepth > 100 )
					{
						loosecandicluster.push_back(a_NeCandiClu);
						refloosedis.push_back( DisMatrix_Core_Neutral[i][j][1] );
					}
				}	
				else if( currChargeCoreType == 110  )		// Electron
				{
					if( DisMatrix_Core_Neutral[i][j][0] < 0.15*CluDepth + 15 )
					{
						tightcandicluster.push_back(a_NeCandiClu);
						reftightdis.push_back( DisMatrix_Core_Neutral[i][j][0] );
					}
				}			
				else if( currChargeCoreType == 111 )		// look behind... might be pion...
				{
					if( DisMatrix_Core_Neutral[i][j][0] < 0.1*CluDepth + 15 && DisMatrix_Core_Neutral[i][j][1] < 0.1*CluDepth + 10 )	//Define Brems Photon region for correct
					{
						tightcandicluster.push_back(a_NeCandiClu);
						if(DisMatrix_Core_Neutral[i][j][0] < DisMatrix_Core_Neutral[i][j][1])	// not fully adequate.
						{
							reftightdis.push_back( DisMatrix_Core_Neutral[i][j][0] );
						}
						else
						{
							reftightdis.push_back( DisMatrix_Core_Neutral[i][j][1] );
						}
					}
					else if( DisMatrix_Core_Neutral[i][j][0] < 0.2*CluDepth + 15 || DisMatrix_Core_Neutral[i][j][1] < 0.2*CluDepth + 15  )
					{	
						loosecandicluster.push_back(a_NeCandiClu);

						if(DisMatrix_Core_Neutral[i][j][0] < DisMatrix_Core_Neutral[i][j][1])   // not fully adequate.
						{
							refloosedis.push_back( DisMatrix_Core_Neutral[i][j][0] );
						}
						else
						{
							refloosedis.push_back( DisMatrix_Core_Neutral[i][j][1] );
						}
					}
				}
				else if( currChargeCoreType == 211 )	//Main Cluster distance oriented
				{
					if(DisMatrix_Core_Neutral[i][j][0] < 0.2*CluDepth)
					{
						tightcandicluster.push_back(a_NeCandiClu);
						reftightdis.push_back( DisMatrix_Core_Neutral[i][j][0] );
					}
				}
				else if( currChargeCoreType == 212 )	//Non_Matched
				{
					if(DisMatrix_Core_Neutral[i][j][0] < 10 + 0.5*CluDepth )	//Energy Dependence...
					{
						tightcandicluster.push_back(a_NeCandiClu);
						reftightdis.push_back( DisMatrix_Core_Neutral[i][j][0] );
					}
					else if(DisMatrix_Core_Neutral[i][j][0] < 10 + 0.4*CluDepth || DisMatrix_Core_Neutral[i][j][1] < 20 + 0.5*CluDepth )
					{
						loosecandicluster.push_back(a_NeCandiClu);
						refloosedis.push_back( DisMatrix_Core_Neutral[i][j][0] );
					}
				}
				else if( currChargeCoreType == 0 ) // && a_recoP_ch->getEnergy() < 3 ) // && !FlagLowPtPz )	//
				{
					if(CluDepth < 20)
					{
						if(DisMatrix_Core_Neutral[i][j][1] < MinDisToNoClusterTrk)	//Tag minimal distance cluster... and see if it can be potentially linked.
						{
							MinDisToNoClusterTrk = DisMatrix_Core_Neutral[i][j][1];
						}
						if( MinDisToNoClusterTrk < 300 && abs(a_recoP_ch->getEnergy() - NeCandEn) < 1.5*a_recoP_ch->getEnergy() )	//some hard cut
						{
							tightcandicluster.push_back(a_NeCandiClu);
							reftightdis.push_back( DisMatrix_Core_Neutral[i][j][1] );
						}
					}
				}
				else
				{
					cout<<"Over balanced/Un matched/defined case: "<<a_recoP_ch->getEnergy()<<" ??? "<<currChargeCoreType<<endl; 
				}
			}
		}

		float totaltightcandiEn = 0; 
		float totalloosecandiEn = 0; 
		for(unsigned int s = 0; s < tightcandicluster.size(); s++)
		{
			totaltightcandiEn += tightcandicluster[s]->getEnergy();
		}

		for(unsigned int s = 0; s < loosecandicluster.size(); s++)
		{
			totalloosecandiEn += loosecandicluster[s]->getEnergy();
		}

		if( currChargeCoreType == 130 )
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				Cluster* a_clu = tightcandicluster[i1];
				if( ClusterType_1stID[a_clu] >= 10 ) //  && CurrClusterEnergy < CurrTrackEnergy + 1.0*sqrt(CurrTrackEnergy) )	//Frags...
				{
					mergedcluster.push_back( a_clu );		
					CurrClusterEnergy += a_clu->getEnergy();
				}
				else if( ClusterType_1stID[a_clu] < 10 && (CurrClusterEnergy + a_clu->getEnergy() < CurrTrackEnergy + 1.0*sqrt(CurrTrackEnergy) ))
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}
		}
		else if( currChargeCoreType == 131 )
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				Cluster* a_clu = tightcandicluster[i1];
				if( ClusterType_1stID[a_clu] >= 10 || (ClusterType_1stID[a_clu] < 10 && (CurrClusterEnergy + a_clu->getEnergy() < CurrTrackEnergy + 1.0*sqrt(CurrTrackEnergy)))  ) 
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}

			//Maybe Some ID over here?...	//layers & numbers...	//BS Case ID

			for(unsigned int i2 = 0; i2 < loosecandicluster.size(); i2++)
			{
				Cluster* a_clu = loosecandicluster[i2];
				if( ClusterType_1stID[a_clu] >= 10 || CurrClusterEnergy + a_clu->getEnergy() < CurrTrackEnergy + 1.0*sqrt(CurrTrackEnergy))       //Frags...Or some minmal hit cut
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}
		}
		else if( currChargeCoreType == 110 )
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				Cluster* a_clu = tightcandicluster[i1];
				if( ClusterType_1stID[a_clu] >= 10 || (ClusterType_1stID[a_clu] < 10 && CurrClusterEnergy + a_clu->getEnergy() < CurrTrackEnergy + 0.5*sqrt(CurrTrackEnergy))  )
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}
		}
		else if( currChargeCoreType == 111 )
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				Cluster* a_clu = tightcandicluster[i1];
				if( ClusterType_1stID[a_clu] >= 10 || (ClusterType_1stID[a_clu] < 10 && CurrClusterEnergy + a_clu->getEnergy() < CurrTrackEnergy + 0.5*sqrt(CurrTrackEnergy))  )
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}

			for(unsigned int i2 = 0; i2 < loosecandicluster.size(); i2++)
			{
				Cluster* a_clu = loosecandicluster[i2];
				if( ClusterType_1stID[a_clu] >= 10 || fabs(CurrClusterEnergy + a_clu->getEnergy() - CurrTrackEnergy) < fabs(CurrClusterEnergy - CurrTrackEnergy) )       //Frags...Or some minmal hit cut
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}	
		}
		else if( currChargeCoreType == 211 )	// Matched
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				Cluster* a_clu = tightcandicluster[i1];
				if( ClusterType_1stID[a_clu] >= 10 || (ClusterType_1stID[a_clu] < 10 && CurrClusterEnergy + a_clu->getEnergy() < CurrTrackEnergy + 1.5*sqrt(CurrTrackEnergy))  )
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}
		}
		else if( currChargeCoreType == 212)
		{
			for(unsigned int i1 = 0; i1 < tightcandicluster.size(); i1++)
			{
				Cluster* a_clu = tightcandicluster[i1];
				if( ClusterType_1stID[a_clu] >= 10 || (ClusterType_1stID[a_clu] < 10 && CurrClusterEnergy + a_clu->getEnergy() < CurrTrackEnergy + 1.5*sqrt(CurrTrackEnergy))  )
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}

			for(unsigned int i2 = 0; i2 < loosecandicluster.size(); i2++)
			{
				Cluster* a_clu = loosecandicluster[i2];
				if( ClusterType_1stID[a_clu] >= 10 || fabs(CurrClusterEnergy + a_clu->getEnergy() - CurrTrackEnergy) < fabs(CurrClusterEnergy - CurrTrackEnergy) )       //Frags...Or some minmal hit cut
				{
					mergedcluster.push_back( a_clu );
					CurrClusterEnergy += a_clu->getEnergy();
				}
			}
		}
		else if( currChargeCoreType == 0 && reftightdis.size() > 0)
		{
			float mindis = 1.0E10;
			int minindex = 0; 

			for(unsigned int i1 = 0; i1 < reftightdis.size(); i1 ++)
			{
				if(reftightdis[i1] < mindis)
				{
					mindis = reftightdis[i1];
					minindex = i1; 
				}
			}

			Cluster* a_clu = tightcandicluster[minindex];	// Only 1? ...

			mergedcluster.push_back( a_clu );
		}
		else
		{
			cout<<"No_match"<<endl; 
		}

		float CHCluEnergy = 0;

		for(int is = 0; is < int(mergedcluster.size()); is++)
		{       
			Cluster* a_TBM_clu = mergedcluster[is]; 
			CHCluEnergy += EnUltraHotCorr(a_TBM_clu->getEnergy(), a_TBM_clu);
		}
		if( !( CHCluEnergy < 1 && CurrTrackEnergy > 5 ) || (MinDisToOtherTrack < 100) )	// Need to check if exist nearby larger charged cluster: maybe absorbed together and only left tiny MIP tail in the ECAL //* bool closeTonearByEnergeticChargedShower = 0  // MIP like; should also protect against energies
		{
			for(int i2 = 0; i2 < int(mergedcluster.size()); i2++)
			{
				Cluster* a_TBM_clu = mergedcluster[i2];
				NNCTouchFlag[a_TBM_clu]	= 2; 		// can make use of this intereting flag...
			}

			float charge = a_chargedTrk -> getOmega()/fabs(a_chargedTrk -> getOmega());

			ReconstructedParticleImpl * chargeparticle = new ReconstructedParticleImpl();
			chargeparticle->setEnergy( CurrTrackEnergy );
			chargeparticle->setCharge(charge);
			HelixClass * currHelix = new HelixClass();
			currHelix->Initialize_Canonical(a_chargedTrk->getPhi(), a_chargedTrk -> getD0(), a_chargedTrk -> getZ0(), a_chargedTrk -> getOmega(), a_chargedTrk->getTanLambda(), BField);
			chargeparticle->setMomentum( currHelix->getMomentum() );
			chargeparticle->addTrack( a_chargedTrk );

			TotalChEn += CurrTrackEnergy;
			int flagEnergyFlow = 0;

			//Clustermerging
			ClusterImpl * chclustermerged =  NaiveMergeClu(mergedcluster);
			mergedclu_ch->addElement(chclustermerged);
			chargeparticle->addCluster(chclustermerged);
			chargedclustercore_merged.push_back(chclustermerged);

			int currChargeCoreType2;
			
			currChargeCoreType2 = ClusterFlag(chclustermerged, a_chargedTrk, col_TPCTrk);

			if(currChargeCoreType2 == 130 || currChargeCoreType2 == 131)
			{
				chargeparticle->setType( int(-13*charge) );
			}
			else if(currChargeCoreType2 == 110 || currChargeCoreType2 == 111)
			{
				chargeparticle->setType( int(-11*charge) );
				if(CHCluEnergy > CurrTrackEnergy + 0.5*sqrt(CurrTrackEnergy) + 1)
				{
					flagEnergyFlow = 1; 
				}
			}
			else
			{
				chargeparticle->setType( int(211*charge) );
				if(CHCluEnergy > CurrTrackEnergy + 1.2*sqrt(CurrTrackEnergy) + 1)
				{
					flagEnergyFlow = 2;
				}
			}

			//Energy Flow Procedure

			if( flagEnergyFlow )
			{
				ReconstructedParticleImpl * a_Ef_Ne_particle = new ReconstructedParticleImpl();
				a_Ef_Ne_particle->setEnergy( CHCluEnergy - CurrTrackEnergy );
				TVector3 corePos = chclustermerged->getPosition();
				float WFactor = (CHCluEnergy - CurrTrackEnergy)/corePos.Mag(); 
				float PFNEMom[3] = {WFactor*float(corePos.X()), WFactor*float(corePos.Y()), WFactor*float(corePos.Z())};
				a_Ef_Ne_particle->setMomentum(PFNEMom);
				a_Ef_Ne_particle->setMass( 0.0 );
				a_Ef_Ne_particle->setCharge( 0.0 );
				a_Ef_Ne_particle->setType(501);
				arborrecoparticle->addElement(a_Ef_Ne_particle);

				cout<<"Energy Flow Neutral Tagged "<<CHCluEnergy - CurrTrackEnergy<<endl; 
			}

			arborrecoparticle->addElement(chargeparticle);
		}
		else	// push non valid tracks, etc to second iteration, as those for PreInteracting ones
		{
			SecondIterTracks.push_back(a_chargedTrk);
			cout<<"Second Iter Track Found"<<endl; 
		}	
	}
	evtPP->addCollection(mergedclu_ch, "ArborCharged");

	std::vector<Cluster*> Ab_or_veto_clu;
	std::vector<Cluster*> BBCore; 
	Ab_or_veto_clu.clear();
	BBCore.clear();

	for(int p6 = 0; p6 < NNeutralCluster; p6 ++)
	{
		Cluster * c_clu = non_chargedclustercore[p6];
		if( NNCTouchFlag.find(c_clu) == NNCTouchFlag.end() )
		{
			if( ClusterType_1stID[c_clu] < 10 || c_clu->getEnergy() > 0.05 + 0.001*CluDepthMap[c_clu] )	//Cores
			{
				BBCore.push_back(c_clu);
			}
		} 
	}

	float NAMom[3] = {0, 0, 0};

	//Final Re-absorption
	std::vector<Cluster*> NBBNeutral; 
	NBBNeutral.clear();
	//End 

	for(int s = 0; s < int (BBCore.size()); s++)
	{
		Cluster * a_clu = BBCore[s];
		TVector3 PosClu = a_clu->getPosition();
		float Depth = DisSeedSurface(PosClu);
		float CoreEnCorr = ClusterEE(a_clu);

		if(ClusterFlag1st(a_clu) == 11)
		{
			cout<<"PPPP "<<CoreEnCorr<<" : "<<Depth<<endl; 
			TVector3 BushSeedPos = a_clu->getPosition();
			ReconstructedParticleImpl * neutralparticle = new ReconstructedParticleImpl();
			neutralparticle->setType(22);
			TVector3 PP = ClusterCoG(a_clu);
			NAMom[0] = CoreEnCorr*1.0/PP.Mag()*PP.X();
			NAMom[1] = CoreEnCorr*1.0/PP.Mag()*PP.Y();
			NAMom[2] = CoreEnCorr*1.0/PP.Mag()*PP.Z();
			neutralparticle->setEnergy( CoreEnCorr );
			neutralparticle->setMass( 0.0 );
			neutralparticle->setCharge( 0.0 );
			neutralparticle->setMomentum( NAMom );
			ClusterImpl * a_neclu = NaiveCluImpl(a_clu);
			a_neclu->setEnergy( CoreEnCorr );	//Reset...
			neutralparticle->addCluster(a_neclu);
			mergedclu_ne->addElement(a_neclu);
			arborrecoparticle->addElement(neutralparticle);
		}
		else	// Distance to Charged Core > sth;
		{
			float MinDisToChCore = 1.0E9;
			float currDis = 0; 
			int NChCore = mergedclu_ch->getNumberOfElements();
			float closestChCluEn = 0; 			
			for(int t = 0; t < NChCore; t++)
			{
				Cluster * a_chclu = dynamic_cast<Cluster*>(mergedclu_ch->getElementAt(t));
				currDis = BushDis(a_chclu, a_clu);
				if(currDis < MinDisToChCore)
				{
					MinDisToChCore = currDis;
					closestChCluEn = a_chclu->getEnergy();	// Or the Trk En??
				}
			}
			if( MinDisToChCore > 0.4*(15 + closestChCluEn + Depth*0.01) || a_clu->getEnergy() > 2.0 )	//Joint Depth??
			{
				NBBNeutral.push_back(a_clu);
			}
		}
	}

	// Add: Neural Core Remerge & Energy Scale Recalculate
	// IsoHit Abso

	std::vector<Cluster*> NBBAbs = ClusterHitAbsorbtion(NBBNeutral, IsoHits, 100);	// Huge??
	
	std::vector<float> BBAbsEn; 
	BBAbsEn.clear();

	for(unsigned s1 = 0; s1 < NBBAbs.size(); s1++)
	{
		BBAbsEn.push_back(NBBAbs[s1]->getEnergy());
	}

	std::vector<int> BBAbsIndex = SortMeasure(BBAbsEn, 1);

	std::vector<Cluster *> NeutronCore;
	std::vector<Cluster *> NeutronFlag;
	NeutronCore.clear();
	NeutronFlag.clear();	

	for(unsigned int s2 = 0; s2 < NBBAbs.size(); s2++)	//Sort it; the first one must be a neutral core?
	{
		Cluster * p_clu = NBBAbs[BBAbsIndex[s2]];
		float currCluEn = p_clu->getEnergy();
		cout<<p_clu->getEnergy()<<" vs "<<NBBNeutral[s2]->getEnergy()<<endl; 
		if( currCluEn > 1.0 || (currCluEn > 0.5 && s2 < 2) )
		{
			NeutronCore.push_back(p_clu);
		}
		else
		{
			NeutronFlag.push_back(p_clu);
		}
	}

	std::vector<Cluster *> Neutrons = ClusterAbsorbtion(NeutronCore, NeutronFlag, 200, 0.01);
	// std::cout<<"BBB "<<Neutrons.size()<<" < "<< NBBAbs.size() <<" =?= "<<NeutronCore.size()<<" + "<<NeutronFlag.size()<<std::endl; 

	for(unsigned int s3 = 0; s3 < Neutrons.size(); s3++)
	{
		Cluster * a_clu = Neutrons[s3];
		float CoreEnCorr = ClusterEE(a_clu);
		TVector3 SeedPos = a_clu->getPosition();
		float Depth = DisSeedSurface(SeedPos);
	
		cout<<"SSS "<<CoreEnCorr<<" : "<<Depth<<endl; 

		if( CoreEnCorr > 0.1 + 0.003*Depth || a_clu->getCalorimeterHits().size() > 20)
		{
			if( ClusterFlag1st(a_clu) == 11 )	// Photon ID
				cout<<"WARNING... Photons after neutron merge merged"<<endl; 
			TVector3 PosClu = a_clu->getPosition();
			ReconstructedParticleImpl * neutralparticle = new ReconstructedParticleImpl();
			neutralparticle->setType(2112);
			TVector3 PP = ClusterCoG(a_clu);
			NAMom[0] = CoreEnCorr*1.0/PP.Mag()*PP.X();
			NAMom[1] = CoreEnCorr*1.0/PP.Mag()*PP.Y();
			NAMom[2] = CoreEnCorr*1.0/PP.Mag()*PP.Z();
			neutralparticle->setEnergy( CoreEnCorr );
			neutralparticle->setMass( 0.0 );
			neutralparticle->setCharge( 0.0 );
			neutralparticle->setMomentum( NAMom );
			ClusterImpl * a_neclu = NaiveCluImpl(a_clu);
			a_neclu->setEnergy( CoreEnCorr );       //Reset...
			neutralparticle->addCluster(a_neclu);
			mergedclu_ne->addElement(a_neclu);
			arborrecoparticle->addElement(neutralparticle);
		}
	}

	evtPP->addCollection( arborrecoparticle, "ArborPFOs");
	evtPP->addCollection( mergedclu_ne, "ArborNeutral" );
}


void BushConnect::processEvent( LCEvent * evtP )
{
	if (evtP)
	{

		_eventNr = evtP->getEventNumber();

		if(_eventNr%1 == 0)
			cout<<"Nevts Processed: "<<_eventNr<<endl;

		BushConnect::Clean();	
		BushConnect::TrackSort( evtP ); 
		BushConnect::BushSelfMerge( evtP ); 	
		BushConnect::TagCore( evtP );		//ECAL
		BushConnect::ParticleReco( evtP );
	}
}

void BushConnect::end()
{
	std::cout<<"Bush Connection Finished, ArborObject Formed"<<std::endl;	
}

