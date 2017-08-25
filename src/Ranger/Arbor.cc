#include <ArborHit.hh>
#include <Arbor.hh>
#include <TTree.h>
#include <TH1F.h>
#include <TFile.h>
#include <algorithm>
#include <TMath.h>
#include "KDTreeLinkerAlgoT.h"
#include <unordered_map>
#include "Rtypes.h"

//#define DEBUG

using namespace std; 

typedef std::unordered_multimap<unsigned,unsigned> HitLinkMap;

std::vector<ArborHit> cleanedHits;

std::vector<int> LeafHitsIndex; 
std::vector<int> JointHitsIndex; 
std::vector<int> StarJointHitsIndex; 
std::vector<int> IsoHitsIndex;
std::vector<int> SimpleSeedHitsIndex;
std::vector<int> StarSeedHitsIndex;
std::map<int, int> HitsFlag; 	// Tell Hits type; 

//common block
static Float_t  sDisAB;
static Int_t    sSubD_A;
static Int_t    sSubD_B;
static Float_t  sDepth_A;
static Float_t  sDepth_B;
static Int_t    sFlagEE_J;
static Int_t    sFlagPSEE;
static Float_t  sEnergyA;
static Float_t  sEnergyB;
static Float_t  sthetaA;
static Float_t  sthetaB;
static TVector3 sPosA;
static TVector3 sPosB;


static Float_t fParaMeanofDisABL;
static Float_t fParaRMSofDisABL;
static Float_t fParaMeanofDABL;
static Float_t fParaRMSofDABL;
static Float_t fParaMeanofDisABR;
static Float_t fParaRMSofDisABR;
static Float_t fParaMeanofDABR;
static Float_t fParaRMSofDABR;
static std::vector<Float_t> sParamTable;
void LoadParam(const std::vector<Float_t>& vParam){
    fParaMeanofDisABL=vParam[0];
    fParaRMSofDisABL=vParam[1];
    fParaMeanofDABL=vParam[2];
    fParaRMSofDABL=vParam[3];

    fParaMeanofDisABR=vParam[4];
    fParaRMSofDisABR=vParam[5];
    fParaMeanofDABR=vParam[6];
    fParaRMSofDABR=vParam[7];
}
/*Int_t GetIndex(Int_t layer, Int_t cell, Int_t idx)
  {
  const Int_t NELEMENT_PER_ROW = 8;
  Int_t row(0), col(0);
  if(layer<=16) row=0;
  else if(layer>16&&layer<=20) row=1;
  else if(layer>20&&layer<=26) row=2;
  else if(layer>26&&layer<=30) row=3;
  if(cell<=1) col=0;
  else if(cell>1&&cell<=2)  col=1;
  else if(cell>2&&cell<=5)  col=2;
  else if(cell>5&&cell<=10) col=3;
  return row*NELEMENT_PER_ROW + col + idx;
  }
  Float_t GetParaMean(Int_t layer, Int_t cell, Int_t idx)
  {
  Int_t loc = GetIndex(layer, cell, idx);
  return sParamTable[loc];
  }
  Float_t GetParaRMS(Int_t layer, Int_t cell, Int_t idx)
  {
  Int_t loc = GetIndex(layer, cell, idx);
  return sParamTable[loc];
  }*/
//
/*
   HitLinkMap BackLinksMap;
   HitLinkMap alliterBackLinksMap;
   */
HitLinkMap IterBackLinks;

linkcoll Links;
linkcoll InitLinks;
linkcoll IterInputLinks; 
linkcoll alliterlinks;
linkcoll links_debug; 
linkcoll IterLinks; 
linkcoll IterLinks_1; 
branchcoll LengthSortBranchCollection;
branchcoll Trees; 

int NHits = 0; 

void init() {
    cleanedHits.clear();
    LeafHitsIndex.clear();
    JointHitsIndex.clear();
    StarJointHitsIndex.clear();
    IsoHitsIndex.clear();
    SimpleSeedHitsIndex.clear();
    StarSeedHitsIndex.clear();
    HitsFlag.clear();

    Links.clear();
    InitLinks.clear();
    IterLinks_1.clear();
    IterLinks.clear();
    IterInputLinks.clear();
    IterBackLinks.clear();
    LengthSortBranchCollection.clear();
    Trees.clear();
    alliterlinks.clear();
    links_debug.clear();
}

void HitsCleaning( std::vector<ArborHit> inputHits )
{
    cleanedHits = inputHits;        //Cannot Really do much things here. Mapping before calling
    NHits = cleanedHits.size();
}

void HitsClassification( linkcoll inputLinks )
{
    int NLinks =  inputLinks.size();

    LeafHitsIndex.clear();
    JointHitsIndex.clear();
    StarJointHitsIndex.clear();
    IsoHitsIndex.clear();
    SimpleSeedHitsIndex.clear();
    StarSeedHitsIndex.clear();

    std::pair<int, int> a_link;

    int BeginIndex[NHits];
    int EndIndex[NHits];

    for(int i0 = 0; i0 < NHits; i0++)
    {
        BeginIndex[i0] = 0;
        EndIndex[i0] = 0;
    }

    for(int j0 = 0; j0 < NLinks; j0++)
    {
        BeginIndex[ (inputLinks[j0].first) ] ++;
        EndIndex[ (inputLinks[j0].second) ] ++;
    }

    for(int i1 = 0; i1 < NHits; i1++)
    {
        if(BeginIndex[i1] == 0 && EndIndex[i1] == 1)
        {
            LeafHitsIndex.push_back(i1);
            HitsFlag[i1] = 5;
        }
        else if(BeginIndex[i1] == 1 && EndIndex[i1] == 1)
        {
            JointHitsIndex.push_back(i1);
            HitsFlag[i1] = 4; 
        }
        else if(BeginIndex[i1] > 1 && EndIndex[i1] == 1)
        {
            StarJointHitsIndex.push_back(i1);
            HitsFlag[i1] = 3;
        }
        else if(BeginIndex[i1] == 1 && EndIndex[i1] == 0)
        {
            SimpleSeedHitsIndex.push_back(i1);
            HitsFlag[i1] = 2;
        }
        else if(BeginIndex[i1] > 1 && EndIndex[i1] == 0)
        {
            StarSeedHitsIndex.push_back(i1);
            HitsFlag[i1] = 1;
        }
        else if(BeginIndex[i1] == 0 && EndIndex[i1] == 0)
        {
            IsoHitsIndex.push_back(i1);
            HitsFlag[i1] = 0;
        }
        else
        {
            cout<<"WARNING: UNCLASSIFIED HITS, Begin Index: "<<BeginIndex[i1]<<",  End Index:  "<<EndIndex[i1]<<endl; 
        }
    }

#ifdef DEBUG
    cout<<"Verification of Hits Classification: "<<endl;
    cout<<"Seed - Simple/Star: "<<SimpleSeedHitsIndex.size()<<" : "<<StarSeedHitsIndex.size()<<endl;
    cout<<"Joint - Simple/Star: "<<JointHitsIndex.size()<<" : "<<StarJointHitsIndex.size()<<endl;
    cout<<"Leaves: "<<LeafHitsIndex.size()<<endl;
    cout<<"IsoHits: "<<IsoHitsIndex.size()<<endl; 
    cout<<"TotalHits: "<<NHits<<endl; 
#endif
}

linkcoll LinkClean( std::vector<ArborHit> allhits, linkcoll alllinks )
{
    linkcoll cleanedlinks; 

    int NLinks = alllinks.size();
    int Ncurrhitlinks = 0;
    int MinAngleIndex = -1;
    float MinAngle = 1E6;
    float tmpOrder = 0;
    float DirAngle = 0;

    std::pair<int, int> SelectedPair;

    TVector3 PosA, PosB, PosDiffAB;

    std::vector< std::vector<int> > LinkHits;
    LinkHits.clear();
    for(int s1 = 0; s1 < NHits; s1++)
    {
        std::vector<int> hitlink;
        for(int t1 = 0; t1 < NLinks; t1++)
        {
            if(alllinks[t1].second == s1)
            {
                hitlink.push_back(alllinks[t1].first);
            }
        }
        LinkHits.push_back(hitlink);
    }

    for(int i1 = 0; i1 < NHits; i1++)
    {
        PosB = cleanedHits[i1].GetPosition();
        MinAngleIndex = -10;
        MinAngle = 1E6;

        std::vector<int> currhitlink = LinkHits[i1];

        Ncurrhitlinks = currhitlink.size();

        for(int k1 = 0; k1 < Ncurrhitlinks; k1++)
        {
            PosA = cleanedHits[ currhitlink[k1] ].GetPosition();
            DirAngle = (PosA + PosB).Angle(PosB - PosA);		//Replace PosA + PosB with other order parameter ~ reference direction
            tmpOrder = (PosB - PosA).Mag() * (DirAngle + 0.1);
            if( tmpOrder < MinAngle ) // && DirAngle < 2.5 )
            {
                MinAngleIndex = currhitlink[k1];
                MinAngle = tmpOrder;
            }
        }

        if(MinAngleIndex > -0.5)
        {
            SelectedPair.first = MinAngleIndex;
            SelectedPair.second = i1;
            cleanedlinks.push_back(SelectedPair);
        }
    }

#ifdef DEBUG
    cout<<"NStat "<<NHits<<" : "<<NLinks<<" InitLinks "<<InitLinks.size()<<endl;
#endif
    return cleanedlinks;
}

void BuildInitLink( std::vector<float>Thresholds )
{
    KDTreeLinkerAlgo<unsigned,3> kdtree;
    typedef KDTreeNodeInfoT<unsigned,3> KDTreeNodeInfo;
    std::array<float,3> minpos{ {0.0f,0.0f,0.0f} }, maxpos{ {0.0f,0.0f,0.0f} };
    std::vector<KDTreeNodeInfo> nodes, found;

    for(int i0 = 0; i0 < NHits; ++i0 )
    {
        const auto& hit = cleanedHits[i0].GetPosition();
        nodes.emplace_back(i0,(float)hit.X(),(float)hit.Y(),(float)hit.Z());
        if( i0 == 0 )
        {
            minpos[0] = hit.X(); minpos[1] = hit.Y(); minpos[2] = hit.Z();
            maxpos[0] = hit.X(); maxpos[1] = hit.Y(); maxpos[2] = hit.Z();
        }
        else
        {
            minpos[0] = std::min((float)hit.X(),minpos[0]);
            minpos[1] = std::min((float)hit.Y(),minpos[1]);
            minpos[2] = std::min((float)hit.Z(),minpos[2]);
            maxpos[0] = std::max((float)hit.X(),maxpos[0]);
            maxpos[1] = std::max((float)hit.Y(),maxpos[1]);
            maxpos[2] = std::max((float)hit.Z(),maxpos[2]);
        }
    }
    KDTreeCube kdXYZCube(minpos[0],maxpos[0],
            minpos[1],maxpos[1],
            minpos[2],maxpos[2]);
    kdtree.build(nodes,kdXYZCube);
    nodes.clear();

    Links.clear();	//all tmp links

    // jing

    TVector3 PosA, PosB, PosDiffAB, PosDiffBA;
    int NLayer_A = 0;
    int NLayer_B = 0;
    int NStave_A = 0;
    int NStave_B = 0;
    int SubD_A = 0;
    int SubD_B = 0;
    float Depth_A = 0;
    float Depth_B = 0;
    float DisAB = 0;
    float MagA = 0;
    float MagB = 0;
    float ECCorr = 0;
    int FlagTrkPS = 0;
    int FlagPSEE = 0;
    int FlagEH = 0;
    int FlagHH = 0;
    int FlagStaveSame = 0;
    int FlagStaveDiff = 0;
    int FlagEE_J = 0;
    float EnergyA = 0.0;
    float EnergyB = 0.0;
    float thetaA = 0.0;
    float thetaB = 0.0;
    //TFile *hits_root_file = new TFile("hits_th.root","update");
    //TTree *tree = (TTree*)hits_root_file->Get("tt");
    /*
       TH1F *ps = new TH1F("TrkPS","SubD=0,3",20,0,15);
       TH1F *ee_sub_dis = new TH1F("ee_sub_dis","SubD=1,3 && DisAB<EE+0.05Depths",20,0,15);
       TH1F *ee_sub = new TH1F("ee_sub","SubD=1,3",20,0,15);
       TH1F *eh_sub_dis = new TH1F("eh_sub_dis","SubD=1,2 && DisAB<EH+ECCorr",20,0,500);
       TH1F *eh_sub = new TH1F("eh_sub","SubD=1,2",20,0,500);
       TH1F *hh_sub_dis = new TH1F("hh_sub_dis","SubD=2,2 && DisAB<HH+0.03Depths",20,0,500);
       TH1F *hh_sub = new TH1F("hh_sub","SubD=2,2",20,0,500);
       */
    for(int i0 = 0; i0 < NHits; i0++)
    {

        found.clear();

        PosA = cleanedHits[i0].GetPosition();
        NLayer_A = cleanedHits[i0].GetLayer();
        NStave_A = cleanedHits[i0].GetStave();
        SubD_A = cleanedHits[i0].GetSubD();	//SubD_Index, 0 = track endpoint, 1 = Ecal, 2 = Hcal, 3 = EcalPS
        Depth_A = cleanedHits[i0].GetDepth();
        EnergyA = cleanedHits[i0].GetEnergy();
        thetaA = PosA.Theta();

        const float side = 200;	//could also be sub detector dependent
        const float xplus(PosA.X() + side), xminus(PosA.X() - side);
        const float yplus(PosA.Y() + side), yminus(PosA.Y() - side);
        const float zplus(PosA.Z() + side), zminus(PosA.Z() - side);
        const float xmin(std::min(xplus,xminus)), xmax(std::max(xplus,xminus));
        const float ymin(std::min(yplus,yminus)), ymax(std::max(yplus,yminus));
        const float zmin(std::min(zplus,zminus)), zmax(std::max(zplus,zminus));
        KDTreeCube searchcube( xmin, xmax,
                ymin, ymax,
                zmin, zmax );
        kdtree.search(searchcube,found);

        for(unsigned int j0 = 0; j0 < found.size(); j0++)
        {
            if( found[j0].data <= (unsigned)i0 ) continue;
            PosB = cleanedHits[found[j0].data].GetPosition();
            NLayer_B = cleanedHits[found[j0].data].GetLayer();
            NStave_B = cleanedHits[found[j0].data].GetStave();
            SubD_B = cleanedHits[found[j0].data].GetSubD();
            Depth_B = cleanedHits[found[j0].data].GetDepth();
            EnergyB = cleanedHits[j0].GetEnergy();
            thetaB = PosB.Theta();

            PosDiffAB = PosA - PosB;
            PosDiffBA = PosB - PosA;
            DisAB = PosDiffAB.Mag();
            ECCorr = 0;
            if( (fabs(PosA.Z()) - 2600 )*(fabs(PosB.Z()) - 2600) < 0 )
                ECCorr = 40;

            //std::cout <<DisAB<<endl;
            // For the XX seed, using triangle method...
            FlagTrkPS = 0; FlagPSEE = 0; FlagEH = 0; FlagHH = 0; FlagStaveSame = 0; FlagStaveDiff = 0;FlagEE_J=0;

            if( SubD_A*SubD_B == 0 && SubD_A + SubD_B == 3 && DisAB < 180 && (PosDiffAB.Angle(PosA) < 0.1 || PosDiffBA.Angle(PosA) < 0.1) )
                FlagTrkPS = 1;
            else if( (SubD_A == 1 || SubD_A == 3) && (SubD_B == 1 || SubD_B == 3) )
            {
                FlagEE_J=1;
                if ( DisAB < Thresholds[0] + 0.05*(Depth_A + Depth_B) ) 
                   FlagPSEE = 1;
            //    Float_t normDis, normDAB;
            //    if(0.5*(Depth_A + Depth_B)< 100){
            //        normDis = (DisAB-fParaMeanofDisABL)/fParaRMSofDisABL;
            //        normDAB = (0.5*(Depth_A + Depth_B)-fParaMeanofDABL)/fParaRMSofDABL;
            //    }else{
            //        normDis = (DisAB-fParaMeanofDisABR)/fParaRMSofDisABR;
            //        normDAB = (0.5*(Depth_A + Depth_B)-fParaMeanofDABR)/fParaRMSofDABR;
            //    }
            //    if(normDis*normDis + normDAB*normDAB< 50)//3 sigma deviation
            //        FlagPSEE =1 ;
            }
            else if( (SubD_A * SubD_B == 2 && DisAB < Thresholds[1] + ECCorr) )
                FlagEH = 1;
            else if( SubD_A == 2 && SubD_B == 2 && DisAB < Thresholds[2] + 0.03*(Depth_A + Depth_B) )
                FlagHH = 1;
            //if (tree->Fill()) std::cout <<"tree fill"<<endl;
            sDisAB=DisAB;
            sSubD_A=SubD_A;
            sSubD_B=SubD_B;
            sDepth_A=Depth_A;
            sDepth_B=Depth_B;
            sFlagEE_J=FlagEE_J;
            sFlagPSEE=FlagPSEE;
            sEnergyA=EnergyA;
            sEnergyB=EnergyB;
            sthetaA=thetaA;
            sthetaB=thetaB;
            sPosA=PosA;
            sPosB=PosB;
            tree->Fill();
            /*
               if( SubD_A*SubD_B == 0 && SubD_A + SubD_B == 3 && DisAB < 180 && (PosDiffAB.Angle(PosA) < 0.1 || PosDiffBA.Angle(PosA) < 0.1) )
               {
               FlagTrkPS = 1;
               ps->Fill(DisAB);
               }
               else if( (SubD_A == 1 || SubD_A == 3) && (SubD_B == 1 || SubD_B == 3) )
               {
               ee_sub->Fill(DisAB);
               if (DisAB < Thresholds[0] + 0.05*(Depth_A + Depth_B) )
               {
               FlagPSEE = 1;
               ee_sub_dis->Fill(DisAB);
               }
               }
               else if( SubD_A * SubD_B == 2 )
               {
               eh_sub->Fill(DisAB);
               if ( DisAB < Thresholds[1] + ECCorr)
               {
               FlagEH = 1;
               eh_sub_dis->Fill(DisAB);
               }
               }
               else if( SubD_A == 2 && SubD_B == 2  )
               {
               hh_sub->Fill(DisAB);
               if (DisAB < Thresholds[2] + 0.03*(Depth_A + Depth_B) )
               {
               FlagHH = 1;
               hh_sub_dis->Fill(DisAB);
               }
               }
               dis_ab->Fill(DisAB);
               */
            if( FlagTrkPS || FlagPSEE || FlagEH || FlagHH )
            {
                std::pair<int, int> a_Link;
                MagA = PosA.Mag();
                MagB = PosB.Mag();

                if( NStave_A != NStave_B || ( NLayer_A == 0 && NLayer_B != 0 ) || ( NLayer_B == 0 && NLayer_A != 0 ) )
                    FlagStaveDiff = 1;
                else if( NLayer_A != NLayer_B && NStave_A == NStave_B )
                    FlagStaveSame = 1;

                if( FlagStaveDiff || FlagStaveSame || FlagTrkPS || FlagEH )
                {
                    if( MagA > MagB )
                    {
                        a_Link.first = found[j0].data;
                        a_Link.second = i0;
                    }
                    else
                    {
                        a_Link.first = i0;
                        a_Link.second = found[j0].data;
                    }
                    Links.push_back(a_Link);
                }
            }
        }
    }
    links_debug = Links;
    std::cout<<"Links Size: "<<Links.size()<<endl;
    /*
       dis_ab->Write();
       ps->Write();
       ee_sub->Write();
       ee_sub_dis->Write();
       eh_sub->Write();
       eh_sub_dis->Write();
       hh_sub->Write();
       hh_sub_dis->Write();
       */
    //tree->Write("",TObject::kOverwrite);
    //hits_root_file->Close();
}

void LinkIteration( int time )	//Energy corrections, semi-local correction
{

    KDTreeLinkerAlgo<unsigned,3> kdtree;
    typedef KDTreeNodeInfoT<unsigned,3> KDTreeNodeInfo;
    std::array<float,3> minpos{ {0.0f,0.0f,0.0f} }, maxpos{ {0.0f,0.0f,0.0f} };
    std::vector<KDTreeNodeInfo> nodes, found;

    TVector3 RefDir[NHits];
    int Nin_hit[NHits];
    int Nout_hit[NHits];

    for(int i0 = 0; i0 < NHits; i0++)
    {
        const auto& hit = cleanedHits[i0].GetPosition();
        RefDir[i0] = 1.0/hit.Mag() * hit;
        Nin_hit[i0] = 0;
        Nout_hit[i0] = 0;

        nodes.emplace_back(i0,(float)hit.X(),(float)hit.Y(),(float)hit.Z());
        if( i0 == 0 )
        {
            minpos[0] = hit.X(); minpos[1] = hit.Y(); minpos[2] = hit.Z();
            maxpos[0] = hit.X(); maxpos[1] = hit.Y(); maxpos[2] = hit.Z();
        }
        else
        {
            minpos[0] = std::min((float)hit.X(),minpos[0]);
            minpos[1] = std::min((float)hit.Y(),minpos[1]);
            minpos[2] = std::min((float)hit.Z(),minpos[2]);
            maxpos[0] = std::max((float)hit.X(),maxpos[0]);
            maxpos[1] = std::max((float)hit.Y(),maxpos[1]);
            maxpos[2] = std::max((float)hit.Z(),maxpos[2]);
        }
    }

    KDTreeCube kdXYZCube(minpos[0],maxpos[0],
            minpos[1],maxpos[1],
            minpos[2],maxpos[2]);
    kdtree.build(nodes,kdXYZCube);
    nodes.clear();

    IterBackLinks.clear();

    if(time == 1)
    {
        IterLinks_1.clear();
        alliterlinks = InitLinks;
    }
    else
    {
        IterLinks.clear();
        alliterlinks = IterLinks_1;
    }
    int NcurrLinks = alliterlinks.size();

    TVector3 hitPos, PosA, PosB, PosDiffAB, PosDiffBA, linkDir;
    int NLayer_A = 0;
    int NLayer_B = 0;
    int NStave_A = 0;
    int NStave_B = 0;
    int SubD_A = 0;
    int SubD_B = 0;
    int FlagEE = 0;
    int FlagHH = 0;
    int AngleAccIndex = 0;
    float DisAB = 0;
    float MagA = 0;
    float MagB = 0;
    std::pair<int, int> currlink;
    std::pair<int, int> a_Link;
    std::pair<int, int> a_tmpLink, b_tmpLink;
    int FlagNoJoint = 0;
    int FlagNoIso = 0;
    //	int FlagNoExisting = 0;

    for(int i = 0; i < NHits; i++)
    {
        hitPos = cleanedHits[i].GetPosition();
        RefDir[i] = 1.0/hitPos.Mag() * hitPos;
        Nin_hit[i] = 0;
        Nout_hit[i] = 0;
    }

    // std::vector<int> ---> all the hits linked to this hit
    std::vector< std::vector<int> > hitLinksArray;
    hitLinksArray.resize(NHits);
    for(int j = 0; j < NcurrLinks; j++)
    {
        currlink = alliterlinks[j];
        PosA = cleanedHits[ currlink.first ].GetPosition();
        PosB = cleanedHits[ currlink.second ].GetPosition();
        linkDir = (PosA - PosB);		//Links are always from first point to second - verify
        linkDir *= 1.0/linkDir.Mag();
        RefDir[currlink.first] += 4*linkDir; 	//Weights... might be optimized...
        RefDir[currlink.second] += 6*linkDir;
        Nin_hit[currlink.first] ++;
        Nout_hit[currlink.second] ++;
        hitLinksArray[currlink.first].push_back(currlink.second);
    }

    //Classification of cases: Branch to IsoHit, Branch, Seed
    for(int i1 = 0; i1 < NHits; i1++)
    {
        found.clear();
        RefDir[i1] *= 1.0/RefDir[i1].Mag();
        PosA = cleanedHits[i1].GetPosition();
        NLayer_A = cleanedHits[i1].GetLayer();
        NStave_A = cleanedHits[i1].GetStave();
        SubD_A = cleanedHits[i1].GetSubD();

        const float side = 200 + 100*time;
        const float xplus(PosA.X() + side), xminus(PosA.X() - side);
        const float yplus(PosA.Y() + side), yminus(PosA.Y() - side);
        const float zplus(PosA.Z() + side), zminus(PosA.Z() - side);
        const float xmin(std::min(xplus,xminus)), xmax(std::max(xplus,xminus));
        const float ymin(std::min(yplus,yminus)), ymax(std::max(yplus,yminus));
        const float zmin(std::min(zplus,zminus)), zmax(std::max(zplus,zminus));
        KDTreeCube searchcube( xmin, xmax,
                ymin, ymax,
                zmin, zmax );
        kdtree.search(searchcube,found);

        for(unsigned int j1 = 0; j1 < found.size(); j1++)
        {
            if( found[j1].data <= (unsigned)i1 ) continue;

            FlagNoJoint = Nout_hit[j1] * Nin_hit[i1] * Nout_hit[i1] * Nin_hit[j1];
            FlagNoIso = Nout_hit[j1] + Nin_hit[i1] + Nout_hit[i1] + Nin_hit[j1];
            a_tmpLink.first = i1;
            a_tmpLink.second = j1;
            b_tmpLink.first = j1;
            b_tmpLink.second = i1;

            bool isConnected = false;

            std::vector<int>& linksOfHitI = hitLinksArray[i1];
            std::vector<int>& linksOfHitJ = hitLinksArray[j1];


            for(auto it=linksOfHitI.begin(); it!=linksOfHitI.end(); ++it) {
                int hitConnected = *it;

                if(hitConnected==(int)j1) {
                    isConnected = true;
                    break;
                }
            }

            // a piece of code copied from upper lines :(
            if(!isConnected) {
                for(auto it=linksOfHitJ.begin(); it!=linksOfHitJ.end(); ++it) {
                    int hitConnected = *it;

                    if(hitConnected==(int)i1) {
                        isConnected = true;
                        break;
                    }
                }
            }


            if( FlagNoJoint == 0 && FlagNoIso != 0 )
            {
                if(!isConnected)
                {
                    PosB = cleanedHits[found[j1].data].GetPosition();
                    NLayer_B = cleanedHits[found[j1].data].GetLayer();
                    NStave_B = cleanedHits[found[j1].data].GetStave();
                    SubD_B = cleanedHits[found[j1].data].GetSubD();
                    PosDiffAB = PosB - PosA;
                    PosDiffBA = PosA - PosB;
                    DisAB = PosDiffAB.Mag();

                    FlagEE = 0;
                    FlagHH = 0;
                    AngleAccIndex = 0;

                    if( PosDiffAB.Angle(RefDir[i1]) < 0.6/time )
                        AngleAccIndex = 1;
                    else if( PosDiffAB.Angle(RefDir[j1]) < 0.6/time )
                        AngleAccIndex = 2;
                    else if( PosDiffBA.Angle(RefDir[i1]) < 0.6/time )
                        AngleAccIndex = 3;
                    else if( PosDiffBA.Angle(RefDir[j1]) < 0.6/time )
                        AngleAccIndex = 4;

                    if(AngleAccIndex)
                    {
                        if(SubD_A == 1 && SubD_B == 1 && DisAB < 15*(time+1) ) //  && DisAB > 15*time)
                            FlagEE = 1;
                        if(SubD_A > 2 && SubD_B == 2 && DisAB < 50*(time+1) ) // && DisAB > 50*time )
                            FlagHH = 1;
                    }

                    if(FlagEE || FlagHH)
                    {
                        MagA = PosA.Mag();
                        MagB = PosB.Mag();

                        //	if(NLayer_A != NLayer_B)
                        //	{
                        if( NStave_A != NStave_B || ( NLayer_A == 0 && NLayer_B != 0 ) || ( NLayer_B == 0 && NLayer_A != 0 ) )
                        {
                            if( MagA > MagB && AngleAccIndex < 3 )
                            {
                                a_Link.first = found[j1].data;
                                a_Link.second = i1;
                            }
                            else if( MagA < MagB && AngleAccIndex > 2 )
                            {
                                a_Link.first = i1;
                                a_Link.second = found[j1].data;
                            }
                            alliterlinks.push_back(a_Link);
                            hitLinksArray[a_Link.first].push_back(a_Link.second);
                        }
                        else if( NLayer_A != NLayer_B && NStave_A == NStave_B)
                        {
                            if( MagA > MagB && AngleAccIndex < 3 )
                            {
                                a_Link.first = found[j1].data;
                                a_Link.second = i1;
                            }
                            else if( MagA < MagB && AngleAccIndex > 2 )
                            {
                                a_Link.first = i1;
                                a_Link.second = found[j1].data;
                            }
                            alliterlinks.push_back(a_Link);
                            hitLinksArray[a_Link.first].push_back(a_Link.second);
                        }
                        //	}
                    }
                }
            }
        }
    }

    //Reusage of link iteration codes?

    int NLinks = alliterlinks.size();
    int MinAngleIndex = -10;
    int Ncurrhitlinks = 0; 
    float MinAngle = 1E6; 
    float tmpOrder = 0;
    float DirAngle = 0; 
    std::pair<int, int> SelectedPair; 

    std::vector< std::vector<int> > LinkHits;
    LinkHits.clear();
    for(int s1 = 0; s1 < NHits; s1++)
    {
        std::vector<int> hitlink;
        for(int t1 = 0; t1 < NLinks; t1++)
        {
            if(alliterlinks[t1].second == s1)
            {
                hitlink.push_back(alliterlinks[t1].first);
            }
        }
        LinkHits.push_back(hitlink);
    }

    for(int i2 = 0; i2 < NHits; i2++)
    {
        PosB = cleanedHits[i2].GetPosition();
        MinAngleIndex = -10;
        MinAngle = 1E6;

        std::vector<int> currhitlink = LinkHits[i2];

        Ncurrhitlinks = currhitlink.size();

        for(int j2 = 0; j2 < Ncurrhitlinks; j2++)
        {
            PosA = cleanedHits[ currhitlink[j2] ].GetPosition();
            DirAngle = (RefDir[i2]).Angle(PosA - PosB);
            tmpOrder = (PosB - PosA).Mag() * (DirAngle + 0.6);
            if(tmpOrder < MinAngle) //  && DirAngle < 1.0)
            {
                MinAngleIndex = currhitlink[j2];
                MinAngle = tmpOrder;
            }
        }

        if(MinAngleIndex > -0.5)
        {
            SelectedPair.first = MinAngleIndex;
            SelectedPair.second = i2;
            if(SelectedPair.first == SelectedPair.second)
            {
                cout<<"WTTTTTTTTFFFFFFFFFFFFFF"<<endl;
                continue; 
            }
            if(time == 1)
            {
                IterLinks_1.push_back(SelectedPair);
            }
            else
            {
                IterLinks.push_back(SelectedPair);
                IterBackLinks.emplace(SelectedPair.second,SelectedPair.first);
            }
        }
    }	

#ifdef DEBUG
    cout<<"Init-Iter Size "<<InitLinks.size()<<" : "<<IterLinks.size()<<endl;
#endif
}

void BranchBuilding(float SeedThreshold)
{

    int NLinks = IterLinks.size();
    int NBranches = 0;
    std::map <int, int> HitBeginIndex;
    std::map <int, int> HitEndIndex;
    std::vector< std::vector<int> > InitBranchCollection;
    std::vector< std::vector<int> > PrunedInitBranchCollection;
    std::vector< std::vector<int> > TmpBranchCollection;
    TVector3 PosA, PosB;

    for(int i1 = 0; i1 < NHits; i1++)
    {
        HitBeginIndex[i1] = 0;
        HitEndIndex[i1] = 0;
    }

    for(int j1 = 0; j1 < NLinks; j1++)
    {
        HitBeginIndex[ (IterLinks[j1].first) ] ++;
        HitEndIndex[ (IterLinks[j1].second) ] ++;
    }

    int iterhitindex = 0;
    int LL = 0;     //Uplimt to be set to twice the total layer thickness...

    for(int i2 = 0; i2 < NHits; i2++)
    {
        if(HitEndIndex[i2] > 1)
            cout<<"WARNING OF INTERNAL LOOP with more than 1 link stopped at the same Hit"<<endl;

        if(HitBeginIndex[i2] == 0 && HitEndIndex[i2] == 1)        //EndPoint
        {
            NBranches ++;
            std::vector<int> currBranchhits;      //array of indexes 

            iterhitindex = i2;
            currBranchhits.push_back(i2);
            LL = 0;

            while(HitEndIndex[iterhitindex] != 0 && LL < 300)      // 100 put by hand
            {
                /*
                   for(int j2 = 0; j2 < NLinks; j2++)
                   {
                   std::pair<int, int> PairIterator = IterLinks[j2];
                   if(  PairIterator.second == iterhitindex && std::find(currBranchhits.begin(), currBranchhits.end(), PairIterator.first) == currBranchhits.end() )
                   {
                   currBranchhits.push_back(PairIterator.first);
                   iterhitindex = PairIterator.first;
                   break;
                   }
                   }  
                   */
                auto iterlink_range = IterBackLinks.equal_range(iterhitindex);
                assert( std::distance(iterlink_range.first,iterlink_range.second) == 1 );
                iterhitindex = iterlink_range.first->second;
                currBranchhits.push_back(iterhitindex);
                LL++;
            }
            InitBranchCollection.push_back(std::move(currBranchhits) );
        }
    }

    PrunedInitBranchCollection.resize(InitBranchCollection.size());

    std::vector<float> BranchSize;
    std::vector<float> cutBranchSize;
    std::vector<int> SortedBranchIndex;
    std::vector<int> SortedcutBranchIndex;
    std::vector<int> currBranch;
    std::vector<int> iterBranch;
    std::vector<int> touchedHits;
    std::vector<bool> touchedHitsMap(NHits,false);
    std::vector<bool> seedHitsMap(NHits,false);
    std::vector<unsigned> seedHits;
    std::vector<int> leadingbranch;

    KDTreeLinkerAlgo<unsigned,3> kdtree;
    typedef KDTreeNodeInfoT<unsigned,3> KDTreeNodeInfo;
    std::array<float,3> minpos{ {0.0f,0.0f,0.0f} }, maxpos{ {0.0f,0.0f,0.0f} };
    std::vector<KDTreeNodeInfo> nodes, found;
    bool needInitPosMinMax = true;

    HitLinkMap seedToBranchesMap;
    std::unordered_map<int,int> branchToSeed;
    std::map<branch, int> SortedBranchToOriginal;
    SortedBranchToOriginal.clear();
    branchToSeed.clear();

    int currBranchSize = 0;
    int currHit = 0;

    for(int i3 = 0; i3 < NBranches; i3++)
    {
        currBranch = InitBranchCollection[i3];
        BranchSize.push_back( float(currBranch.size()) );
    }

    SortedBranchIndex = std::move( SortMeasure(BranchSize, 1) );

    for(int i4 = 0; i4 < NBranches; i4++)
    {
        currBranch = InitBranchCollection[SortedBranchIndex[i4]];
        currBranchSize = currBranch.size();
        iterBranch.clear();

        for(int j4 = 0; j4 < currBranchSize; j4++)
        {
            currHit = currBranch[j4];

            if( !touchedHitsMap[currHit] )
            {
                iterBranch.push_back(currHit);
                touchedHitsMap[currHit] = true;
            }
        }
        const auto theseed = currBranch[currBranchSize - 1];
        SortedBranchToOriginal[iterBranch] = theseed;     //Map to seed...
        branchToSeed.emplace(SortedBranchIndex[i4],theseed);
        seedToBranchesMap.emplace(theseed,SortedBranchIndex[i4]); // map seed to branches
        if( !seedHitsMap[theseed] )
        {
            seedHitsMap[theseed] = true;
            const auto& hit = cleanedHits[theseed].GetPosition();
            nodes.emplace_back(theseed,(float)hit.X(),(float)hit.Y(),(float)hit.Z());
            if( needInitPosMinMax ) {
                needInitPosMinMax = false;
                minpos[0] = hit.X(); minpos[1] = hit.Y(); minpos[2] = hit.Z();
                maxpos[0] = hit.X(); maxpos[1] = hit.Y(); maxpos[2] = hit.Z();
            } else {
                minpos[0] = std::min((float)hit.X(),minpos[0]);
                minpos[1] = std::min((float)hit.Y(),minpos[1]);
                minpos[2] = std::min((float)hit.Z(),minpos[2]);
                maxpos[0] = std::max((float)hit.X(),maxpos[0]);
                maxpos[1] = std::max((float)hit.Y(),maxpos[1]);
                maxpos[2] = std::max((float)hit.Z(),maxpos[2]);
            }
        }

        TmpBranchCollection.push_back(iterBranch);
        cutBranchSize.push_back( float(iterBranch.size()) );
        PrunedInitBranchCollection[SortedBranchIndex[i4]] = std::move(iterBranch);
    }

    SortedcutBranchIndex = std::move( SortMeasure(cutBranchSize, 1) );

    for(int i6 = 0; i6 < NBranches; i6++)
    {
        currBranch.clear();
        currBranch = TmpBranchCollection[ SortedcutBranchIndex[i6]];
        LengthSortBranchCollection.push_back(currBranch);;
    }

    std::vector<bool> link_helper(NBranches*NBranches,false);
    TVector3 DisSeed;
    KDTreeCube kdXYZCube(minpos[0],maxpos[0],
            minpos[1],maxpos[1],
            minpos[2],maxpos[2]);
    kdtree.build(nodes,kdXYZCube);
    nodes.clear();

    QuickUnion qu(NBranches);

    for(int i7 = 0; i7 < NBranches; i7++)
    {
        auto SeedIndex_A = branchToSeed[i7];
        auto shared_branches = seedToBranchesMap.equal_range(SeedIndex_A);

        for( auto itr = shared_branches.first; itr != shared_branches.second; ++itr ) {
            const auto foundSortedIdx = itr->second;
            if( foundSortedIdx <= (unsigned)i7 ) continue;
            if( link_helper[NBranches*i7 + foundSortedIdx] ||
                    link_helper[NBranches*foundSortedIdx + i7]    ) continue;
            if( !qu.connected(i7,foundSortedIdx) ) {
                qu.unite(i7,foundSortedIdx);
            }
            link_helper[NBranches*i7 + foundSortedIdx] = true;
            link_helper[NBranches*foundSortedIdx + i7] = true;
        }

        const auto& seedpos = cleanedHits[SeedIndex_A].GetPosition();
        const float side = 10*SeedThreshold;
        const float xplus(seedpos.X() + side), xminus(seedpos.X() - side);
        const float yplus(seedpos.Y() + side), yminus(seedpos.Y() - side);
        const float zplus(seedpos.Z() + side), zminus(seedpos.Z() - side);
        const float xmin(std::min(xplus,xminus)), xmax(std::max(xplus,xminus));
        const float ymin(std::min(yplus,yminus)), ymax(std::max(yplus,yminus));
        const float zmin(std::min(zplus,zminus)), zmax(std::max(zplus,zminus));
        KDTreeCube searchcube( xmin, xmax,
                ymin, ymax,
                zmin, zmax );
        found.clear();
        kdtree.search(searchcube,found);
        for(unsigned j7 = 0; j7 < found.size(); j7++) {
            DisSeed = seedpos - cleanedHits[ found[j7].data ].GetPosition();

            if( DisSeed.Mag() < SeedThreshold )
            {
                auto seed_branches = seedToBranchesMap.equal_range(found[j7].data);
                for( auto itr = seed_branches.first; itr != seed_branches.second; ++itr ){
                    const auto foundSortedIdx = itr->second;
                    if( foundSortedIdx <= (unsigned)i7 ) continue;
                    if( link_helper[NBranches*i7 + foundSortedIdx] ||
                            link_helper[NBranches*foundSortedIdx + i7]    ) continue;
                    if( !qu.connected(i7,foundSortedIdx) ) {
                        qu.unite(i7,foundSortedIdx);
                    }
                    link_helper[NBranches*i7 + foundSortedIdx] = true;
                    link_helper[NBranches*foundSortedIdx + i7] = true;
                }
            }
        }
    }

    Trees.clear();
    std::unordered_map<unsigned,branch> merged_branches(qu.count());
    Trees.reserve(qu.count());

    for( unsigned i = 0; i < (unsigned)NBranches; ++i ) {
        unsigned root = qu.find(i);
        const auto& branch = PrunedInitBranchCollection[i];
        auto& merged_branch = merged_branches[root];
        merged_branch.insert(merged_branch.end(),branch.begin(),branch.end());
    }

    unsigned total_hits = 0;
    for( auto& final_branch : merged_branches ) {
        total_hits += final_branch.second.size();
        Trees.push_back(std::move(final_branch.second));
    }
}

void BushMerging()
{
    cout<<"Merging branch"<<endl;
}

//HitThresholds: EEThreshold, EHThreshold, HHThreshold, EE_Seed_Threshold;
std::vector< std::vector<int> > Arbor( std::vector<ArborHit> inputHits, std::vector<float>Thresholds)
{

    if(Thresholds.size() != 4)
    {
        cout<<"Threshold Set Wrong, Threshold Size is "<<Thresholds.size()<<" Should be 4"<<endl; 
        exit(2); 
    }

    init();	
    HitsCleaning(inputHits);
    BuildInitLink(Thresholds);
    InitLinks = LinkClean( cleanedHits, Links );

    LinkIteration(1);
    LinkIteration(2);
    HitsClassification(IterLinks);	
    BranchBuilding(Thresholds[3]);	

    return LengthSortBranchCollection;
}

static int root_count=0;
void CreateROOT(std::string str){
    hits_file = new TFile(str.c_str(),"RECREATE");
    tree = new TTree("DisAB","DisAB");
    root_count++;
    std::cout<<"root_count++"<<endl;
    tree->Branch("DisAB",&sDisAB,"DisAB/F");
    tree->Branch("SubD_A",&sSubD_A,"SubD_A/I");
    tree->Branch("SubD_B",&sSubD_B,"SubD_B/I");
    tree->Branch("Depth_A",&sDepth_A,"Depth_A/F");
    tree->Branch("Depth_B",&sDepth_B,"Depth_B/F");
    tree->Branch("FlagEE_J",&sFlagEE_J,"FlagEE_J/I");
    tree->Branch("FlagPSEE",&sFlagPSEE,"FlagPSEE/I");
    tree->Branch("EnergyA",&sEnergyA,"EnergyA/F");
    tree->Branch("EnergyB",&sEnergyB,"EnergyB/F");
    tree->Branch("thetaA",&sthetaA,"thetaA/F");
    tree->Branch("thetaB",&sthetaB,"thetaB/F");
    tree->Branch("PosA",&sPosA[0],"PosA[3]/D");
    tree->Branch("PosB",&sPosB[0],"PosB[3]/D");
}

void CloseROOT(){
    std::cout <<"tree entries:"<< tree->GetEntries()<<endl;
    hits_file->cd();
    tree->Write("",TObject::kOverwrite);
    hits_file->Close();
    std::cout<<"root_count="<<root_count<<endl;
}
