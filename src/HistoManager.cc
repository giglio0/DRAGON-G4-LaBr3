//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <CLHEP/Units/SystemOfUnits.h>

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
:rootFile(0)
{
      
  // Histograms
  for (G4int k=0; k<MaxHisto; k++) histo[k] = 0;
  for (G4int k=0; k<Max2Histo; k++) histo2[k] = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  if ( rootFile ) delete rootFile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{ 
// Creating a tree container to handle histograms.
// This tree is associated with an output file.
//
 G4String fileName = "LaBr3CeResTestG4.root";
 rootFile = new TFile(fileName,"RECREATE");
 if(!rootFile) {
   G4cout << " HistoManager::book :" 
          << " problem creating the ROOT TFile "
          << G4endl;
   return;
 }
// 1-D Histos  
 histo[1] = new TH1D("Energy 1", "LaBr3:Ce Detector 1 Spectrum", 1000, 0*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[1]) G4cout << "\n can't create histo 1" << G4endl;

histo[1]->GetXaxis()->SetTitle("Energy (MeV)");
histo[1]->GetYaxis()->SetTitle("Entries");

 histo[2] = new TH1D("Energy Resolution 1", "LaBr3:Ce Detector 1 Spectrum with Energy Resolution", 1000, 0*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[2]) G4cout << "\n can't create histo 2" << G4endl;
 
histo[2]->GetXaxis()->SetTitle("Energy (MeV)");
histo[2]->GetYaxis()->SetTitle("Entries");

 histo[3] = new TH1D("Energy 2", "LaBr3:Ce Detector 2 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[3]) G4cout << "\n can't create histo 3" << G4endl;

histo[3]->GetXaxis()->SetTitle("Energy (MeV)");
histo[3]->GetYaxis()->SetTitle("Entries");

 histo[4] = new TH1D("Energy Resolution 2", "LaBr3:Ce Detector 2 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[4]) G4cout << "\n can't create histo 4" << G4endl;
 
histo[4]->GetXaxis()->SetTitle("Energy (MeV)");
histo[4]->GetYaxis()->SetTitle("Entries");

 histo[5] = new TH1D("Energy 3", "LaBr3:Ce Detector 3 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[5]) G4cout << "\n can't create histo 5" << G4endl;

histo[5]->GetXaxis()->SetTitle("Energy (MeV)");
histo[5]->GetYaxis()->SetTitle("Entries");

 histo[6] = new TH1D("Energy Resolution 3", "LaBr3:Ce Detector 3 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[6]) G4cout << "\n can't create histo 6" << G4endl;
 
histo[6]->GetXaxis()->SetTitle("Energy (MeV)");
histo[6]->GetYaxis()->SetTitle("Entries");

 histo[7] = new TH1D("Energy 4", "LaBr3:Ce Detector 4 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[7]) G4cout << "\n can't create histo 7" << G4endl;

histo[7]->GetXaxis()->SetTitle("Energy (MeV)");
histo[7]->GetYaxis()->SetTitle("Entries");

 histo[8] = new TH1D("Energy Resolution 4", "LaBr3:Ce Detector 4 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[8]) G4cout << "\n can't create histo 8" << G4endl;
 
histo[8]->GetXaxis()->SetTitle("Energy (MeV)");
histo[8]->GetYaxis()->SetTitle("Entries");

 histo[9] = new TH1D("Energy 5", "LaBr3:Ce Detector 5 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[9]) G4cout << "\n can't create histo 9" << G4endl;

histo[9]->GetXaxis()->SetTitle("Energy (MeV)");
histo[9]->GetYaxis()->SetTitle("Entries");

 histo[10] = new TH1D("Energy Resolution 5", "LaBr3:Ce Detector 5 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[10]) G4cout << "\n can't create histo 10" << G4endl;
 
histo[10]->GetXaxis()->SetTitle("Energy (MeV)");
histo[10]->GetYaxis()->SetTitle("Entries");

 histo[11] = new TH1D("Energy 6", "LaBr3:Ce Detector 6 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[11]) G4cout << "\n can't create histo 11" << G4endl;

histo[11]->GetXaxis()->SetTitle("Energy (MeV)");
histo[11]->GetYaxis()->SetTitle("Entries");

 histo[12] = new TH1D("Energy Resolution 6", "LaBr3:Ce Detector 6 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[12]) G4cout << "\n can't create histo 12" << G4endl;
 
histo[12]->GetXaxis()->SetTitle("Energy (MeV)");
histo[12]->GetYaxis()->SetTitle("Entries");

 histo[13] = new TH1D("Energy 7", "LaBr3:Ce Detector 7 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[13]) G4cout << "\n can't create histo 13" << G4endl;

histo[13]->GetXaxis()->SetTitle("Energy (MeV)");
histo[13]->GetYaxis()->SetTitle("Entries");

 histo[14] = new TH1D("Energy Resolution 7", "LaBr3:Ce Detector 7 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[14]) G4cout << "\n can't create histo 14" << G4endl;
 
histo[14]->GetXaxis()->SetTitle("Energy (MeV)");
histo[14]->GetYaxis()->SetTitle("Entries");

 histo[15] = new TH1D("Energy 8", "LaBr3:Ce Detector 8 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[15]) G4cout << "\n can't create histo 15" << G4endl;

histo[15]->GetXaxis()->SetTitle("Energy (MeV)");
histo[15]->GetYaxis()->SetTitle("Entries");

 histo[16] = new TH1D("Energy Resolution 8", "LaBr3:Ce Detector 8 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[16]) G4cout << "\n can't create histo 16" << G4endl;
 
histo[16]->GetXaxis()->SetTitle("Energy (MeV)");
histo[16]->GetYaxis()->SetTitle("Entries");

 histo[17] = new TH1D("Energy 9", "LaBr3:Ce Detector 9 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[17]) G4cout << "\n can't create histo 17" << G4endl;

histo[17]->GetXaxis()->SetTitle("Energy (MeV)");
histo[17]->GetYaxis()->SetTitle("Entries");

 histo[18] = new TH1D("Energy Resolution 9", "LaBr3:Ce Detector 9 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[18]) G4cout << "\n can't create histo 18" << G4endl;
 
histo[18]->GetXaxis()->SetTitle("Energy (MeV)");
histo[18]->GetYaxis()->SetTitle("Entries");

 histo[19] = new TH1D("Energy 10", "LaBr3:Ce Detector 10 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[19]) G4cout << "\n can't create histo 19" << G4endl;

histo[19]->GetXaxis()->SetTitle("Energy (MeV)");
histo[19]->GetYaxis()->SetTitle("Entries");

 histo[20] = new TH1D("Energy Resolution 10", "LaBr3:Ce Detector 10 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[20]) G4cout << "\n can't create histo 20" << G4endl;
 
histo[20]->GetXaxis()->SetTitle("Energy (MeV)");
histo[20]->GetYaxis()->SetTitle("Entries");

 histo[21] = new TH1D("Energy 11", "LaBr3:Ce Detector 11 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[21]) G4cout << "\n can't create histo 21" << G4endl;

histo[21]->GetXaxis()->SetTitle("Energy (MeV)");
histo[21]->GetYaxis()->SetTitle("Entries");

 histo[22] = new TH1D("Energy Resolution 11", "LaBr3:Ce Detector 11 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[22]) G4cout << "\n can't create histo 22" << G4endl;
 
histo[22]->GetXaxis()->SetTitle("Energy (MeV)");
histo[22]->GetYaxis()->SetTitle("Entries");

 histo[23] = new TH1D("Energy 12", "LaBr3:Ce Detector 12 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[23]) G4cout << "\n can't create histo 23" << G4endl;

histo[23]->GetXaxis()->SetTitle("Energy (MeV)");
histo[23]->GetYaxis()->SetTitle("Entries");

 histo[24] = new TH1D("Energy Resolution 12", "LaBr3:Ce Detector 12 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[24]) G4cout << "\n can't create histo 24" << G4endl;
 
histo[24]->GetXaxis()->SetTitle("Energy (MeV)");
histo[24]->GetYaxis()->SetTitle("Entries");

 histo[25] = new TH1D("Energy 13", "LaBr3:Ce Detector 13 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[25]) G4cout << "\n can't create histo 25" << G4endl;

histo[25]->GetXaxis()->SetTitle("Energy (MeV)");
histo[25]->GetYaxis()->SetTitle("Entries");

 histo[26] = new TH1D("Energy Resolution 13", "LaBr3:Ce Detector 13 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[26]) G4cout << "\n can't create histo 26" << G4endl;
 
histo[26]->GetXaxis()->SetTitle("Energy (MeV)");
histo[26]->GetYaxis()->SetTitle("Entries");

 histo[27] = new TH1D("Energy 14", "LaBr3:Ce Detector 14 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[27]) G4cout << "\n can't create histo 27" << G4endl;

histo[27]->GetXaxis()->SetTitle("Energy (MeV)");
histo[27]->GetYaxis()->SetTitle("Entries");

 histo[28] = new TH1D("Energy Resolution 14", "LaBr3:Ce Detector 14 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[28]) G4cout << "\n can't create histo 28" << G4endl;
 
histo[28]->GetXaxis()->SetTitle("Energy (MeV)");
histo[28]->GetYaxis()->SetTitle("Entries");

 histo[29] = new TH1D("Energy 15", "LaBr3:Ce Detector 15 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[29]) G4cout << "\n can't create histo 29" << G4endl;

histo[29]->GetXaxis()->SetTitle("Energy (MeV)");
histo[29]->GetYaxis()->SetTitle("Entries");

 histo[30] = new TH1D("Energy Resolution 15", "LaBr3:Ce Detector 15 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[30]) G4cout << "\n can't create histo 30" << G4endl;
 
histo[30]->GetXaxis()->SetTitle("Energy (MeV)");
histo[30]->GetYaxis()->SetTitle("Entries");

 histo[31] = new TH1D("Energy 16", "LaBr3:Ce Detector 16 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[31]) G4cout << "\n can't create histo 31" << G4endl;

histo[31]->GetXaxis()->SetTitle("Energy (MeV)");
histo[31]->GetYaxis()->SetTitle("Entries");

 histo[32] = new TH1D("Energy Resolution 16", "LaBr3:Ce Detector 16 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[32]) G4cout << "\n can't create histo 32" << G4endl;
 
histo[32]->GetXaxis()->SetTitle("Energy (MeV)");
histo[32]->GetYaxis()->SetTitle("Entries");

 histo[33] = new TH1D("Energy 17", "LaBr3:Ce Detector 17 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[33]) G4cout << "\n can't create histo 33" << G4endl;

histo[33]->GetXaxis()->SetTitle("Energy (MeV)");
histo[33]->GetYaxis()->SetTitle("Entries");

 histo[34] = new TH1D("Energy Resolution 17", "LaBr3:Ce Detector 17 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[34]) G4cout << "\n can't create histo 34" << G4endl;
 
histo[34]->GetXaxis()->SetTitle("Energy (MeV)");
histo[34]->GetYaxis()->SetTitle("Entries");

 histo[35] = new TH1D("Energy 18", "LaBr3:Ce Detector 18 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[35]) G4cout << "\n can't create histo 35" << G4endl;

histo[35]->GetXaxis()->SetTitle("Energy (MeV)");
histo[35]->GetYaxis()->SetTitle("Entries");

 histo[36] = new TH1D("Energy Resolution 18", "LaBr3:Ce Detector 18 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[36]) G4cout << "\n can't create histo 36" << G4endl;
 
histo[36]->GetXaxis()->SetTitle("Energy (MeV)");
histo[36]->GetYaxis()->SetTitle("Entries");

 histo[37] = new TH1D("Energy 19", "LaBr3:Ce Detector 19 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[37]) G4cout << "\n can't create histo 37" << G4endl;

histo[37]->GetXaxis()->SetTitle("Energy (MeV)");
histo[37]->GetYaxis()->SetTitle("Entries");

 histo[38] = new TH1D("Energy Resolution 19", "LaBr3:Ce Detector 19 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[38]) G4cout << "\n can't create histo 38" << G4endl;
 
histo[38]->GetXaxis()->SetTitle("Energy (MeV)");
histo[38]->GetYaxis()->SetTitle("Entries");

 histo[39] = new TH1D("Energy 20", "LaBr3:Ce Detector 20 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[39]) G4cout << "\n can't create histo 39" << G4endl;

histo[39]->GetXaxis()->SetTitle("Energy (MeV)");
histo[39]->GetYaxis()->SetTitle("Entries");

 histo[40] = new TH1D("Energy Resolution 20", "LaBr3:Ce Detector 20 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[40]) G4cout << "\n can't create histo 40" << G4endl;
 
histo[40]->GetXaxis()->SetTitle("Energy (MeV)");
histo[40]->GetYaxis()->SetTitle("Entries");

 histo[41] = new TH1D("Energy 21", "LaBr3:Ce Detector 21 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[41]) G4cout << "\n can't create histo 41" << G4endl;

histo[41]->GetXaxis()->SetTitle("Energy (MeV)");
histo[41]->GetYaxis()->SetTitle("Entries");

 histo[42] = new TH1D("Energy Resolution 21", "LaBr3:Ce Detector 21 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[42]) G4cout << "\n can't create histo 42" << G4endl;
 
histo[42]->GetXaxis()->SetTitle("Energy (MeV)");
histo[42]->GetYaxis()->SetTitle("Entries");

 histo[43] = new TH1D("Energy 22", "LaBr3:Ce Detector 22 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[43]) G4cout << "\n can't create histo 43" << G4endl;

histo[43]->GetXaxis()->SetTitle("Energy (MeV)");
histo[43]->GetYaxis()->SetTitle("Entries");

 histo[44] = new TH1D("Energy Resolution 22", "LaBr3:Ce Detector 22 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[44]) G4cout << "\n can't create histo 44" << G4endl;
 
histo[44]->GetXaxis()->SetTitle("Energy (MeV)");
histo[44]->GetYaxis()->SetTitle("Entries");

 histo[45] = new TH1D("Energy 23", "LaBr3:Ce Detector 23 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[45]) G4cout << "\n can't create histo 45" << G4endl;

histo[45]->GetXaxis()->SetTitle("Energy (MeV)");
histo[45]->GetYaxis()->SetTitle("Entries");

 histo[46] = new TH1D("Energy Resolution 23", "LaBr3:Ce Detector 23 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[46]) G4cout << "\n can't create histo 46" << G4endl;
 
histo[46]->GetXaxis()->SetTitle("Energy (MeV)");
histo[46]->GetYaxis()->SetTitle("Entries");

 histo[47] = new TH1D("Energy 24", "LaBr3:Ce Detector 24 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[47]) G4cout << "\n can't create histo 47" << G4endl;

histo[47]->GetXaxis()->SetTitle("Energy (MeV)");
histo[47]->GetYaxis()->SetTitle("Entries");

 histo[48] = new TH1D("Energy Resolution 24", "LaBr3:Ce Detector 24 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[48]) G4cout << "\n can't create histo 48" << G4endl;
 
histo[48]->GetXaxis()->SetTitle("Energy (MeV)");
histo[48]->GetYaxis()->SetTitle("Entries");

 histo[49] = new TH1D("Energy 25", "LaBr3:Ce Detector 25 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[49]) G4cout << "\n can't create histo 49" << G4endl;

histo[49]->GetXaxis()->SetTitle("Energy (MeV)");
histo[49]->GetYaxis()->SetTitle("Entries");

 histo[50] = new TH1D("Energy Resolution 25", "LaBr3:Ce Detector 25 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[50]) G4cout << "\n can't create histo 50" << G4endl;
 
histo[50]->GetXaxis()->SetTitle("Energy (MeV)");
histo[50]->GetYaxis()->SetTitle("Entries");

 histo[51] = new TH1D("Energy 26", "LaBr3:Ce Detector 26 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[51]) G4cout << "\n can't create histo 51" << G4endl;

histo[51]->GetXaxis()->SetTitle("Energy (MeV)");
histo[51]->GetYaxis()->SetTitle("Entries");

 histo[52] = new TH1D("Energy Resolution 26", "LaBr3:Ce Detector 26 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[52]) G4cout << "\n can't create histo 52" << G4endl;
 
histo[52]->GetXaxis()->SetTitle("Energy (MeV)");
histo[52]->GetYaxis()->SetTitle("Entries");

 histo[53] = new TH1D("Energy 27", "LaBr3:Ce Detector 27 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[53]) G4cout << "\n can't create histo 53" << G4endl;

histo[53]->GetXaxis()->SetTitle("Energy (MeV)");
histo[53]->GetYaxis()->SetTitle("Entries");

 histo[54] = new TH1D("Energy Resolution 27", "LaBr3:Ce Detector 27 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[54]) G4cout << "\n can't create histo 54" << G4endl;
 
histo[54]->GetXaxis()->SetTitle("Energy (MeV)");
histo[54]->GetYaxis()->SetTitle("Entries");

 histo[55] = new TH1D("Energy 28", "LaBr3:Ce Detector 28 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[55]) G4cout << "\n can't create histo 55" << G4endl;

histo[55]->GetXaxis()->SetTitle("Energy (MeV)");
histo[55]->GetYaxis()->SetTitle("Entries");

 histo[56] = new TH1D("Energy Resolution 28", "LaBr3:Ce Detector 28 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[56]) G4cout << "\n can't create histo 56" << G4endl;
 
histo[56]->GetXaxis()->SetTitle("Energy (MeV)");
histo[56]->GetYaxis()->SetTitle("Entries");

 histo[57] = new TH1D("Energy 29", "LaBr3:Ce Detector 29 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[57]) G4cout << "\n can't create histo 57" << G4endl;

histo[57]->GetXaxis()->SetTitle("Energy (MeV)");
histo[57]->GetYaxis()->SetTitle("Entries");

 histo[58] = new TH1D("Energy Resolution 29", "LaBr3:Ce Detector 29 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[58]) G4cout << "\n can't create histo 58" << G4endl;
 
histo[58]->GetXaxis()->SetTitle("Energy (MeV)");
histo[58]->GetYaxis()->SetTitle("Entries");

 histo[59] = new TH1D("Energy 30", "LaBr3:Ce Detector 30 Spectrum", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[59]) G4cout << "\n can't create histo 59" << G4endl;

histo[59]->GetXaxis()->SetTitle("Energy (MeV)");
histo[59]->GetYaxis()->SetTitle("Entries");

 histo[60] = new TH1D("Energy Resolution 30", "LaBr3:Ce Detector 30 Spectrum with Energy Resolution", 950, 0.5*CLHEP::MeV, 10*CLHEP::MeV);
 if (!histo[60]) G4cout << "\n can't create histo 60" << G4endl;
 
histo[60]->GetXaxis()->SetTitle("Energy (MeV)");
histo[60]->GetYaxis()->SetTitle("Entries");
 G4cout << "\n---> Histogram file is opened in " << fileName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{ 
  if (rootFile) {
    rootFile->Write();       // Writing the histograms to the file
    rootFile->Close();        // and closing the tree (and the file).
    G4cout << "\n----> Histogram Tree is saved \n" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
           << " does not exist. (xbin=" << xbin << " weight=" << weight << ")"
           << G4endl;
    return;
  }
 if  (histo[ih]) { histo[ih]->Fill(xbin, weight); }
}

void HistoManager::Fill2Histo(G4int ih, G4double ybin, G4double zbin, G4double weight)
{
  if (ih >= Max2Histo) {
    G4cout << "---> warning from HistoManager::Fill2Histo() : histo " << ih
           << " does not exist. (ybin=" << ybin << "zbin=" << zbin <<" weight=" << weight << ")"
           << G4endl;
    return;
  }
 if  (histo2[ih]) { histo2[ih]->Fill(ybin, zbin, weight); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double fac)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::Normalize() : histo " << ih
           << " does not exist. (fac=" << fac << ")" << G4endl;
    return;
  }
  if (histo[ih]) histo[ih]->Scale(fac);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::PrintStatistic()
{
  if(histo[1]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[1]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[1]->GetRMS(),  "Energy") << G4endl;
		}
            
  if(histo[3]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[3]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[3]->GetRMS(),  "Energy") << G4endl;
        }
        
  if(histo[5]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[5]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[5]->GetRMS(),  "Energy") << G4endl;
        }
        
  if(histo[7]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[7]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[7]->GetRMS(),  "Energy") << G4endl;
        }
        
  if(histo[9]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[9]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[9]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[11]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[11]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[11]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[13]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[13]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[13]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[15]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[15]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[15]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[17]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[17]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[17]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[19]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[19]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[19]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[21]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[21]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[21]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[23]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[23]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[23]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[25]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[25]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[25]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[27]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[27]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[27]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[29]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[29]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[29]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[31]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[31]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[31]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[33]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[33]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[33]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[35]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[35]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[35]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[37]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[37]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[37]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[39]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[39]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[39]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[41]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[41]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[41]->GetRMS(),  "Energy") << G4endl;
        }
        
   if(histo[43]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[43]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[43]->GetRMS(),  "Energy") << G4endl;
        }
        
    if(histo[45]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[45]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[45]->GetRMS(),  "Energy") << G4endl;
        }
        
    if(histo[47]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[47]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[47]->GetRMS(),  "Energy") << G4endl;
        }
        
    if(histo[49]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[49]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[49]->GetRMS(),  "Energy") << G4endl;
        }
        
    if(histo[51]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[51]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[51]->GetRMS(),  "Energy") << G4endl;
        }
        
    if(histo[53]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[53]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[53]->GetRMS(),  "Energy") << G4endl;
        }
        
    if(histo[55]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[55]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[55]->GetRMS(),  "Energy") << G4endl;
        }
        
    if(histo[57]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[57]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[57]->GetRMS(),  "Energy") << G4endl;
        }
        
    if(histo[59]) {
    G4cout << "\n ----> print histograms statistic \n" << G4endl;
    
    G4cout 
    << " ELaBr3 : mean = " << G4BestUnit(histo[59]->GetMean(), "Energy") 
            << " rms = " << G4BestUnit(histo[59]->GetRMS(),  "Energy") << G4endl;
        }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
