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
// 
//
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(HistoManager* histo)
:fHistoManager(histo)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << "\n--------------------Start of Run " << aRun->GetRunID() 
         << "--------------------------\n" << G4endl;
    
  // Initialize Cumulative Quantities
  //
  fSumELaBr3 = fSum2ELaBr3 = 0.;
  fSumELaBr32 = fSum2ELaBr32 = 0.;
  fSumELaBr33 = fSum2ELaBr33 = 0.;
  fSumELaBr34 = fSum2ELaBr34 = 0.;
  fSumELaBr35 = fSum2ELaBr35 = 0.;
  fSumELaBr36 = fSum2ELaBr36 = 0.;
  fSumELaBr37 = fSum2ELaBr37 = 0.;
  fSumELaBr38 = fSum2ELaBr38 = 0.;
  fSumELaBr39 = fSum2ELaBr39 = 0.;
  fSumELaBr310 = fSum2ELaBr310 = 0.;
  fSumELaBr311 = fSum2ELaBr311 = 0.;
  fSumELaBr312 = fSum2ELaBr312 = 0.;
  fSumELaBr313 = fSum2ELaBr313 = 0.;
  fSumELaBr314 = fSum2ELaBr314 = 0.;
  fSumELaBr315 = fSum2ELaBr315 = 0.;
  fSumELaBr316 = fSum2ELaBr316 = 0.;
  fSumELaBr317 = fSum2ELaBr317 = 0.;
  fSumELaBr318 = fSum2ELaBr318 = 0.;
  fSumELaBr319 = fSum2ELaBr319 = 0.;
  fSumELaBr320 = fSum2ELaBr320 = 0.;
  fSumELaBr321 = fSum2ELaBr321 = 0.;
  fSumELaBr322 = fSum2ELaBr322 = 0.;
  fSumELaBr323 = fSum2ELaBr323 = 0.;
  fSumELaBr324 = fSum2ELaBr324 = 0.;
  fSumELaBr325 = fSum2ELaBr325 = 0.;
  fSumELaBr326 = fSum2ELaBr326 = 0.;
  fSumELaBr327 = fSum2ELaBr327 = 0.;
  fSumELaBr328 = fSum2ELaBr328 = 0.;
  fSumELaBr329 = fSum2ELaBr329 = 0.;
  fSumELaBr330 = fSum2ELaBr330 = 0.;
  
  // Histograms
  //
  fHistoManager->book(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillPerEvent(G4double EDet)
{
  // Accumulate Statistics
  //
  fSumELaBr3 += EDet;
  fSum2ELaBr3 += EDet+EDet;
}

void RunAction::fillPerEvent2(G4double EDet2)
{
  // Accumulate Statistics
  //
  fSumELaBr32 += EDet2;
  fSum2ELaBr32 += EDet2+EDet2;
}

void RunAction::fillPerEvent3(G4double EDet3)
{
  // Accumulate Statistics
  //
  fSumELaBr33 += EDet3;
  fSum2ELaBr33 += EDet3+EDet3;
}

void RunAction::fillPerEvent4(G4double EDet4)
{
  // Accumulate Statistics
  //
  fSumELaBr34 += EDet4;
  fSum2ELaBr34 += EDet4+EDet4;
}

void RunAction::fillPerEvent5(G4double EDet5)
{
  // Accumulate Statistics
  //
  fSumELaBr35 += EDet5;
  fSum2ELaBr35 += EDet5+EDet5;
}

void RunAction::fillPerEvent6(G4double EDet6)
{
  // Accumulate Statistics
  //
  fSumELaBr36 += EDet6;
  fSum2ELaBr36 += EDet6+EDet6;
}

void RunAction::fillPerEvent7(G4double EDet7)
{
  // Accumulate Statistics
  //
  fSumELaBr37 += EDet7;
  fSum2ELaBr37 += EDet7+EDet7;
}

void RunAction::fillPerEvent8(G4double EDet8)
{
  // Accumulate Statistics
  //
  fSumELaBr38 += EDet8;
  fSum2ELaBr38 += EDet8+EDet8;
}

void RunAction::fillPerEvent9(G4double EDet9)
{
  // Accumulate Statistics
  //
  fSumELaBr39 += EDet9;
  fSum2ELaBr39 += EDet9+EDet9;
}

void RunAction::fillPerEvent10(G4double EDet10)
{
  // Accumulate Statistics
  //
  fSumELaBr310 += EDet10;
  fSum2ELaBr310 += EDet10+EDet10;
}

void RunAction::fillPerEvent11(G4double EDet11)
{
  // Accumulate Statistics
  //
  fSumELaBr311 += EDet11;
  fSum2ELaBr311 += EDet11+EDet11;
}

void RunAction::fillPerEvent12(G4double EDet12)
{
  // Accumulate Statistics
  //
  fSumELaBr312 += EDet12;
  fSum2ELaBr312 += EDet12+EDet12;
}

void RunAction::fillPerEvent13(G4double EDet13)
{
  // Accumulate Statistics
  //
  fSumELaBr313 += EDet13;
  fSum2ELaBr313 += EDet13+EDet13;
}

void RunAction::fillPerEvent14(G4double EDet14)
{
  // Accumulate Statistics
  //
  fSumELaBr314 += EDet14;
  fSum2ELaBr314 += EDet14+EDet14;
}

void RunAction::fillPerEvent15(G4double EDet15)
{
  // Accumulate Statistics
  //
  fSumELaBr315 += EDet15;
  fSum2ELaBr315 += EDet15+EDet15;
}

void RunAction::fillPerEvent16(G4double EDet16)
{
  // Accumulate Statistics
  //
  fSumELaBr316 += EDet16;
  fSum2ELaBr316 += EDet16+EDet16;
}

void RunAction::fillPerEvent17(G4double EDet17)
{
  // Accumulate Statistics
  //
  fSumELaBr317 += EDet17;
  fSum2ELaBr317 += EDet17+EDet17;
}

void RunAction::fillPerEvent18(G4double EDet18)
{
  // Accumulate Statistics
  //
  fSumELaBr318 += EDet18;
  fSum2ELaBr318 += EDet18+EDet18;
}

void RunAction::fillPerEvent19(G4double EDet19)
{
  // Accumulate Statistics
  //
  fSumELaBr319 += EDet19;
  fSum2ELaBr319 += EDet19+EDet19;
}

void RunAction::fillPerEvent20(G4double EDet20)
{
  // Accumulate Statistics
  //
  fSumELaBr320 += EDet20;
  fSum2ELaBr320 += EDet20+EDet20;
}

void RunAction::fillPerEvent21(G4double EDet21)
{
  // Accumulate Statistics
  //
  fSumELaBr321 += EDet21;
  fSum2ELaBr321 += EDet21+EDet21;
}

void RunAction::fillPerEvent22(G4double EDet22)
{
  // Accumulate Statistics
  //
  fSumELaBr322 += EDet22;
  fSum2ELaBr322 += EDet22+EDet22;
}

void RunAction::fillPerEvent23(G4double EDet23)
{
  // Accumulate Statistics
  //
  fSumELaBr323 += EDet23;
  fSum2ELaBr323 += EDet23+EDet23;
}

void RunAction::fillPerEvent24(G4double EDet24)
{
  // Accumulate Statistics
  //
  fSumELaBr324 += EDet24;
  fSum2ELaBr324 += EDet24+EDet24;
}

void RunAction::fillPerEvent25(G4double EDet25)
{
  // Accumulate Statistics
  //
  fSumELaBr325 += EDet25;
  fSum2ELaBr325 += EDet25+EDet25;
}

void RunAction::fillPerEvent26(G4double EDet26)
{
  // Accumulate Statistics
  //
  fSumELaBr326 += EDet26;
  fSum2ELaBr326 += EDet26+EDet26;
}

void RunAction::fillPerEvent27(G4double EDet27)
{
  // Accumulate Statistics
  //
  fSumELaBr327 += EDet27;
  fSum2ELaBr327 += EDet27+EDet27;
}

void RunAction::fillPerEvent28(G4double EDet28)
{
  // Accumulate Statistics
  //
  fSumELaBr328 += EDet28;
  fSum2ELaBr328 += EDet28+EDet28;
}

void RunAction::fillPerEvent29(G4double EDet29)
{
  // Accumulate Statistics
  //
  fSumELaBr329 += EDet29;
  fSum2ELaBr329 += EDet29+EDet29;
}

void RunAction::fillPerEvent30(G4double EDet30)
{
  // Accumulate Statistics
  //
  fSumELaBr330 += EDet30;
  fSum2ELaBr330 += EDet30+EDet30;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  // Compute Statistics: Mean and RMS
  //
  fSumELaBr3 /= NbOfEvents; fSum2ELaBr3 /= NbOfEvents;
  G4double rmsELaBr3 = fSum2ELaBr3 - fSumELaBr3*fSumELaBr3;
  if (rmsELaBr3 >0.) rmsELaBr3 = std::sqrt(rmsELaBr3); else rmsELaBr3 = 0.;
  
  fSumELaBr32 /= NbOfEvents; fSum2ELaBr32 /= NbOfEvents;
  G4double rmsELaBr32 = fSum2ELaBr32 - fSumELaBr32*fSumELaBr32;
  if (rmsELaBr32 >0.) rmsELaBr32 = std::sqrt(rmsELaBr32); else rmsELaBr32 = 0.;
  
  fSumELaBr33 /= NbOfEvents; fSum2ELaBr33 /= NbOfEvents;
  G4double rmsELaBr33 = fSum2ELaBr33 - fSumELaBr33*fSumELaBr33;
  if (rmsELaBr33 >0.) rmsELaBr33 = std::sqrt(rmsELaBr33); else rmsELaBr33 = 0.;
  
  fSumELaBr34 /= NbOfEvents; fSum2ELaBr34 /= NbOfEvents;
  G4double rmsELaBr34 = fSum2ELaBr34 - fSumELaBr34*fSumELaBr34;
  if (rmsELaBr34 >0.) rmsELaBr34 = std::sqrt(rmsELaBr34); else rmsELaBr34 = 0.;
  
  fSumELaBr35 /= NbOfEvents; fSum2ELaBr35 /= NbOfEvents;
  G4double rmsELaBr35 = fSum2ELaBr35 - fSumELaBr35*fSumELaBr35;
  if (rmsELaBr35 >0.) rmsELaBr35 = std::sqrt(rmsELaBr35); else rmsELaBr35 = 0.;
  
  fSumELaBr36 /= NbOfEvents; fSum2ELaBr36 /= NbOfEvents;
  G4double rmsELaBr36 = fSum2ELaBr36 - fSumELaBr36*fSumELaBr36;
  if (rmsELaBr36 >0.) rmsELaBr36 = std::sqrt(rmsELaBr36); else rmsELaBr36 = 0.;
  
  fSumELaBr37 /= NbOfEvents; fSum2ELaBr37 /= NbOfEvents;
  G4double rmsELaBr37 = fSum2ELaBr37 - fSumELaBr37*fSumELaBr37;
  if (rmsELaBr37 >0.) rmsELaBr37 = std::sqrt(rmsELaBr37); else rmsELaBr37 = 0.;
  
  fSumELaBr38 /= NbOfEvents; fSum2ELaBr38 /= NbOfEvents;
  G4double rmsELaBr38 = fSum2ELaBr38 - fSumELaBr38*fSumELaBr38;
  if (rmsELaBr38 >0.) rmsELaBr38 = std::sqrt(rmsELaBr38); else rmsELaBr38 = 0.;
  
  fSumELaBr39 /= NbOfEvents; fSum2ELaBr39 /= NbOfEvents;
  G4double rmsELaBr39 = fSum2ELaBr39 - fSumELaBr39*fSumELaBr39;
  if (rmsELaBr39 >0.) rmsELaBr39 = std::sqrt(rmsELaBr39); else rmsELaBr39 = 0.;
  
  fSumELaBr310 /= NbOfEvents; fSum2ELaBr310 /= NbOfEvents;
  G4double rmsELaBr310 = fSum2ELaBr310 - fSumELaBr310*fSumELaBr310;
  if (rmsELaBr310 >0.) rmsELaBr310 = std::sqrt(rmsELaBr310); else rmsELaBr310 = 0.;
  
  fSumELaBr311 /= NbOfEvents; fSum2ELaBr311 /= NbOfEvents;
  G4double rmsELaBr311 = fSum2ELaBr311 - fSumELaBr311*fSumELaBr311;
  if (rmsELaBr311 >0.) rmsELaBr311 = std::sqrt(rmsELaBr311); else rmsELaBr311 = 0.;
  
  fSumELaBr312 /= NbOfEvents; fSum2ELaBr312 /= NbOfEvents;
  G4double rmsELaBr312 = fSum2ELaBr312 - fSumELaBr312*fSumELaBr312;
  if (rmsELaBr312 >0.) rmsELaBr312 = std::sqrt(rmsELaBr312); else rmsELaBr312 = 0.;
  
  fSumELaBr313 /= NbOfEvents; fSum2ELaBr313 /= NbOfEvents;
  G4double rmsELaBr313 = fSum2ELaBr313 - fSumELaBr313*fSumELaBr313;
  if (rmsELaBr313 >0.) rmsELaBr313 = std::sqrt(rmsELaBr313); else rmsELaBr313 = 0.;
  
  fSumELaBr314 /= NbOfEvents; fSum2ELaBr314 /= NbOfEvents;
  G4double rmsELaBr314 = fSum2ELaBr314 - fSumELaBr314*fSumELaBr314;
  if (rmsELaBr314 >0.) rmsELaBr314 = std::sqrt(rmsELaBr314); else rmsELaBr314 = 0.;
  
  fSumELaBr315 /= NbOfEvents; fSum2ELaBr315 /= NbOfEvents;
  G4double rmsELaBr315 = fSum2ELaBr315 - fSumELaBr315*fSumELaBr315;
  if (rmsELaBr315 >0.) rmsELaBr315 = std::sqrt(rmsELaBr315); else rmsELaBr315 = 0.;
  
  fSumELaBr316 /= NbOfEvents; fSum2ELaBr316 /= NbOfEvents;
  G4double rmsELaBr316 = fSum2ELaBr316 - fSumELaBr316*fSumELaBr316;
  if (rmsELaBr316 >0.) rmsELaBr316 = std::sqrt(rmsELaBr316); else rmsELaBr316 = 0.;
  
  fSumELaBr317 /= NbOfEvents; fSum2ELaBr317 /= NbOfEvents;
  G4double rmsELaBr317 = fSum2ELaBr317 - fSumELaBr317*fSumELaBr317;
  if (rmsELaBr317 >0.) rmsELaBr317 = std::sqrt(rmsELaBr317); else rmsELaBr317 = 0.;
  
  fSumELaBr318 /= NbOfEvents; fSum2ELaBr318 /= NbOfEvents;
  G4double rmsELaBr318 = fSum2ELaBr318 - fSumELaBr318*fSumELaBr318;
  if (rmsELaBr318 >0.) rmsELaBr318 = std::sqrt(rmsELaBr318); else rmsELaBr318 = 0.;
  
  fSumELaBr319 /= NbOfEvents; fSum2ELaBr319 /= NbOfEvents;
  G4double rmsELaBr319 = fSum2ELaBr319 - fSumELaBr319*fSumELaBr319;
  if (rmsELaBr319 >0.) rmsELaBr319 = std::sqrt(rmsELaBr319); else rmsELaBr319 = 0.;
  
  fSumELaBr320 /= NbOfEvents; fSum2ELaBr320 /= NbOfEvents;
  G4double rmsELaBr320 = fSum2ELaBr320 - fSumELaBr320*fSumELaBr320;
  if (rmsELaBr320 >0.) rmsELaBr320 = std::sqrt(rmsELaBr320); else rmsELaBr320 = 0.;
  
  fSumELaBr321 /= NbOfEvents; fSum2ELaBr321 /= NbOfEvents;
  G4double rmsELaBr321 = fSum2ELaBr321 - fSumELaBr321*fSumELaBr321;
  if (rmsELaBr321 >0.) rmsELaBr321 = std::sqrt(rmsELaBr321); else rmsELaBr321 = 0.;
  
  fSumELaBr322 /= NbOfEvents; fSum2ELaBr322 /= NbOfEvents;
  G4double rmsELaBr322 = fSum2ELaBr322 - fSumELaBr322*fSumELaBr322;
  if (rmsELaBr322 >0.) rmsELaBr322 = std::sqrt(rmsELaBr322); else rmsELaBr322 = 0.;
  
  fSumELaBr323 /= NbOfEvents; fSum2ELaBr323 /= NbOfEvents;
  G4double rmsELaBr323 = fSum2ELaBr323 - fSumELaBr323*fSumELaBr323;
  if (rmsELaBr323 >0.) rmsELaBr323 = std::sqrt(rmsELaBr323); else rmsELaBr323 = 0.;
  
  fSumELaBr324 /= NbOfEvents; fSum2ELaBr324 /= NbOfEvents;
  G4double rmsELaBr324 = fSum2ELaBr324 - fSumELaBr324*fSumELaBr324;
  if (rmsELaBr324 >0.) rmsELaBr324 = std::sqrt(rmsELaBr324); else rmsELaBr324 = 0.;
  
  fSumELaBr325 /= NbOfEvents; fSum2ELaBr325 /= NbOfEvents;
  G4double rmsELaBr325 = fSum2ELaBr325 - fSumELaBr325*fSumELaBr325;
  if (rmsELaBr325 >0.) rmsELaBr325 = std::sqrt(rmsELaBr325); else rmsELaBr325 = 0.;
  
  fSumELaBr326 /= NbOfEvents; fSum2ELaBr326 /= NbOfEvents;
  G4double rmsELaBr326 = fSum2ELaBr326 - fSumELaBr326*fSumELaBr326;
  if (rmsELaBr326 >0.) rmsELaBr326 = std::sqrt(rmsELaBr326); else rmsELaBr326 = 0.;
  
  fSumELaBr327 /= NbOfEvents; fSum2ELaBr327 /= NbOfEvents;
  G4double rmsELaBr327 = fSum2ELaBr327 - fSumELaBr327*fSumELaBr327;
  if (rmsELaBr327 >0.) rmsELaBr327 = std::sqrt(rmsELaBr327); else rmsELaBr327 = 0.;
  
  fSumELaBr328 /= NbOfEvents; fSum2ELaBr328 /= NbOfEvents;
  G4double rmsELaBr328 = fSum2ELaBr328 - fSumELaBr328*fSumELaBr328;
  if (rmsELaBr328 >0.) rmsELaBr328 = std::sqrt(rmsELaBr328); else rmsELaBr328 = 0.;
  
  fSumELaBr329 /= NbOfEvents; fSum2ELaBr329 /= NbOfEvents;
  G4double rmsELaBr329 = fSum2ELaBr329 - fSumELaBr329*fSumELaBr329;
  if (rmsELaBr329 >0.) rmsELaBr329 = std::sqrt(rmsELaBr329); else rmsELaBr329 = 0.;
  
  fSumELaBr330 /= NbOfEvents; fSum2ELaBr330 /= NbOfEvents;
  G4double rmsELaBr330 = fSum2ELaBr330 - fSumELaBr330*fSumELaBr330;
  if (rmsELaBr330 >0.) rmsELaBr330 = std::sqrt(rmsELaBr330); else rmsELaBr330 = 0.;
  
  // Print
  //
  G4cout
     << "\n--------------------End of Run------------------------------\n"
     << "\n mean Energy in LaBr3:Ce     : " << G4BestUnit(fSumELaBr3,"Energy")
     << " +- "                          << G4BestUnit(rmsELaBr3,"Energy")
     << G4endl;
          
  // Save Histograms
  //
  fHistoManager->PrintStatistic();
  fHistoManager->save();   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
