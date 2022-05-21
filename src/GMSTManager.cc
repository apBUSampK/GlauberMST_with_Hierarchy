

#include "../include/GMSTManager.hh"

 #include "TFile.h"
 #include "TH1D.h"
 #include "TH2D.h" 


GMSTManager::GMSTManager()
  : NucA("Pb"), NucB("Pb"), sigmaNN(-1), CritDist(-1), iterations(-1), single_silh(0), variation(-1)
{ 
  std::cout << "------------------------------------------------\n";

  std::cout << "|######### GlauberMC + MST clustering #########|" << std::endl;
  
  std::cout << "------------------------------------------------\n";

  std::cout << "Please enter colliding nucleus name (side A). U, Pb, Pbpnrw(with neutron skin) Au, Xe, Ag, Br, Cu, Al, O, C is available : ";
  std::cin >> NucA;
 
  std::cout<<"\n";

  std::cout << "Please enter colliding nucleus name (side B). U, Pb, Pbpnrw(with neutron skin) Au, Xe, Ag, Br, Cu, Al, O, C is available : ";
  std::cin >> NucB;

  std::cout<<"\n";

  while(sigmaNN < 0) { 
   std::cout<<"Please enter inelastic cross section (in mb) : ";
   std::cin >> sigmaNN;
  } 

  std::cout<<"\n";

  while(CritDist < 0) { 
   std::cout<<"Please enter critical distance (in fm) : ";
   std::cin >> CritDist;
  } 

  std::cout<<"\n";

  while ( iterations < 1 || iterations>10000000 ) {
    std::cout<<"Please enter number of iterations: ";
    std::cin >> iterations;

  }

  std::cout<<"\n";

  while (single_silh < -1 || single_silh > 1) {
      std::cout << "Please enter silhouette of a single nucleon: ";
      std::cin >> single_silh;
  
  std::cout<<"\n";
  }

  while (variation < 0) {
      std::cout << "Please enter the acceptable variation of critical distance (in percent): ";
      std::cin >> variation;
      variation /= 100;
  }

  std::cout<<"\n";

  std::cout << "Please enter the file name to write histograms (.root will be supplied): ";
  std::cin >> fileName;

  std::cout<<"\n";

  Book();
}


GMSTManager::~GMSTManager()
{
  Save();
}


void GMSTManager::Book()
{

 if ( fileName == "") fileName = "GlauberMST";
 fileType = "root";
 fileFullName = fileName+"."+fileType;
 compressionFactor = 8;
 fFile = new TFile(fileFullName.c_str(), "RECREATE", fileName.c_str(), 8);
 std::cout << "TTree will be written to " << fileFullName << std::endl;
 std::cout<<"\n";

 tree = new TTree("GMST", "TTree to store output of a toy model");
 tree_t = new TTree("GMST_test", "TTree to store output of a test model");

}


void GMSTManager::Save()
{
  fFile->Write();
  std::cout << "\n----> TTree was written into the file " << fileFullName << std::endl;
  delete fFile;
}

