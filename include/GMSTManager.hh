
// Igor Pshenichnov 06.08.2008

#ifndef GMSTManager_h
#define GMSTManager_h 1
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

/*class TFile;
class TH1D;
class TH2D;
class TTree;*/

class GMSTManager
{
  public:

  GMSTManager();
   ~GMSTManager();

  public:

  //TH1D* GetHisto(G4int id) {return histo[id];};
  //TH2D* GetHisto2(G4int id) {return histo2[id];};

       
  void Book();
  void Save();
  
  inline TTree* GetTree(bool t) {return t ? tree_t : tree;};
  inline const char* GetNucleusA() {return NucA.c_str();}; 
  inline const char* GetNucleusB() {return NucB.c_str();}; 
  inline double GetSigmaNN() {return sigmaNN;};
  inline double GetCriticalDistance() {return CritDist;}
  inline double GetSingleSilhouette() {return single_silh;}
  inline double GetVariation() {return variation;}
  inline int GetZ() {return sourceZ;};
  inline int GetA() {return sourceA;};
  inline int GetIterations()  {return iterations;};


  
  private:

  
    TFile* fFile;
    TFile* fFile_t;
    TH1D*  histo[20];
    TH2D*  histo2[10];
    TTree* tree;
    TTree* tree_t;

    int sourceZ;
    int sourceA;
    int Z;
    int A;
    int iterations;  
    int MForNot; 
    double sigmaNN;
    double CritDist;
    double single_silh;
    double variation;
    std::string NucA;
    std::string NucB;

    std::string fileName;
    std::string fileType;
    std::string fileOpenPath;

    std::string fileFullName;
    int    compressionFactor;

    int    binsExEn;
    int    eventsPerBin;      

};

#endif
