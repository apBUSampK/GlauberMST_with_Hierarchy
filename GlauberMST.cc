#include "TGlauber/TGlauberMC.hh"
#include "TGlauber/TGlauNucleon.hh"
#include "TGlauber/TGlauNucleus.hh"
#include "include/GMSTManager.hh"
#include "include/GMSTClustering.hh"
#include "TVector3.h"
#include "TObjArray.h"
#include "TObject.h"
 #include "TFile.h"
 #include "TH1D.h"
 #include "TH2D.h"
 #include "TTree.h"
 #include "TRandom.h" 

int main(){

int NpartA = 0;
int NpartB = 0;
int Ncoll = 0;
int Ncollnn = 0;
int Ncollpp = 0;
int Ncollpn = 0;

double b = 0;
std::vector<int> A;
std::vector<int> Z;
std::vector<int> A_t;
std::vector<int> Z_t;


GMSTManager* manager = new GMSTManager();
manager->GetTree(false)->Branch("NpartA", &NpartA, "NpartA/I");
manager->GetTree(false)->Branch("NpartB", &NpartB, "NpartB/I");
manager->GetTree(false)->Branch("Ncoll", &Ncoll, "Ncoll/I");
manager->GetTree(false)->Branch("b", &b, "b/D");
manager->GetTree(false)->Branch("A", "std::vector" ,&A);
manager->GetTree(false)->Branch("Z", "std::vector" ,&Z);
manager->GetTree(true)->Branch("NpartA", &NpartA, "NpartA/I");
manager->GetTree(true)->Branch("NpartB", &NpartB, "NpartB/I");
manager->GetTree(true)->Branch("Ncoll", &Ncoll, "Ncoll/I");
manager->GetTree(true)->Branch("b", &b, "b/D");
manager->GetTree(true)->Branch("A", "std::vector" ,&A_t);
manager->GetTree(true)->Branch("Z", "std::vector" ,&Z_t);

GMSTClustering* clust_manager = new GMSTClustering();
clust_manager->SetCD(manager->GetCriticalDistance());

TRandom* rand = new TRandom(0);
TGlauberMC *mcg=new TGlauberMC(manager->GetNucleusA(),manager->GetNucleusB(),manager->GetSigmaNN(),-1, rand->GetSeed());
  mcg->SetMinDistance(0);
  mcg->SetNodeDistance(0);
  mcg->SetCalcLength(0);
  mcg->SetCalcArea(0);
  mcg->SetCalcCore(0);
  mcg->SetDetail(99);
 // mcg->SetBmin(0);
 // mcg->SetBmax(18);

for(int count = 0; count < manager->GetIterations()+1; count++){
    mcg->Run(1);
    NpartA = mcg->GetNpartA();
    NpartB = mcg->GetNpartB();
    Ncoll  = mcg->GetNcoll();
    b = mcg->GetB();

    TObjArray* nucleons=mcg->GetNucleons();
    clust_manager->SetUp(nucleons);
    GMSTClusterVector clusters_output;
    clusters_output = clust_manager->GetClusters();

    for(auto & i : clusters_output) {
        A.push_back(i.GetA());
        Z.push_back(i.GetZ());
    }


    clusters_output = clust_manager->GetClusters_HSilhouette();

    for(auto & i : clusters_output) {
        A_t.push_back(i.GetA());
        Z_t.push_back(i.GetZ());
    }


    if(!(count%100)){std::cout<< count <<" events were calculated!    \r" << std::flush;}
    manager->GetTree(true)->Fill();
    manager->GetTree(false)->Fill();
    A.clear();
    Z.clear();
    A_t.clear();
    Z_t.clear();
}

delete manager;
delete clust_manager;
return 0;
}
