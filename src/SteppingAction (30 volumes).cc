//
// 
//  License and Disclaimer                                           
//                                                                   
//  The  Geant4 software  is  copyright of the Copyright Holders  of 
//  the Geant4 Collaboration.  It is provided  under  the terms  and 
//  conditions of the Geant4 Software License,  included in the file 
//  LICENSE and available at  http://cern.ch/geant4/license .  These 
//  include a list of copyright holders.                             
//                                                                   
//  Neither the authors of this software system, nor their employing 
//  institutes,nor the agencies providing financial support for this 
//  work  make  any representation or  warranty, express or implied, 
//  regarding  this  software system or assume any liability for its 
//  use.  Please see the license in the file  LICENSE  and URL above 
//  for the full disclaimer and the limitation of liability.         
//                                                                   
//  This  code  implementation is the result of  the  scientific and 
//  technical work of the GEANT4 collaboration.                      
//  By using,  copying,  modifying or  distributing the software (or 
//  any work based  on the software)  you  agree  to acknowledge its 
//  use  in  resulting  scientific  publications,  and indicate your 
//  acceptance of all terms of the Geant4 Software license.          
// 
//
//
// $Id$
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"

#include <CLHEP/Units/SystemOfUnits.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det,
                                         EventAction* evt)
:fDetector(det), fEventAction(evt)                                         
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // get volume of the current step
  G4VPhysicalVolume* volume = 
  aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  // collect energy and track length step by step
  
  // following lines are to add resolution effects
  time_t seed = time( NULL );
  fRandomEngine = new CLHEP::HepJamesRandom( static_cast< long >( seed ) );
  fRandomGauss = new CLHEP::RandGauss( fRandomEngine );
  //resolution function taken from Saint Gobain document "BrilLanCeTM Scintillators Performance Summary"
  // G4double fERes = fEOrig + fEOrig  fRandomGauss->fire(0.0, (76.3/2.35)/std::sqrt(fEOrig*1000));
  // LaBr3 = BGO (9.7)/3.34 = 0.097*(sqrt(662)) and fEOrig is expressed in units of GeV
  
  if (volume == fDetector->GetScint()){
	G4double fEOrig = aStep->GetTotalEnergyDeposit();
    G4double fERes = fEOrig + fEOrig * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig*1000000));
    fEventAction->AddScint(fEOrig);
    fEventAction->AddResScint(fERes);
   }
  // Single GRIFFIN LaBr3 Detector
  if (volume == fDetector->GetScint1()){
	G4double fEOrig = aStep->GetTotalEnergyDeposit();
    G4double fERes = fEOrig + fEOrig * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig*1000000));
    fEventAction->AddScint(fEOrig);
    fEventAction->AddResScint(fERes);
  }
  
   if (volume == fDetector->GetScint2()){
    G4double fEOrig2 = aStep->GetTotalEnergyDeposit();
    G4double fERes2 = fEOrig2 + fEOrig2 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig2*1000000));
    fEventAction->AddScint2(fEOrig2);
    fEventAction->AddResScint2(fERes2);
  }
  if (volume == fDetector->GetScint3()){
	G4double fEOrig3 = aStep->GetTotalEnergyDeposit();
	G4double fERes3 = fEOrig3 + fEOrig3 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig3*1000000));
    fEventAction->AddScint3(fEOrig3);
    fEventAction->AddResScint3(fERes3);
  }
  if (volume == fDetector->GetScint4()){
	G4double fEOrig4 = aStep->GetTotalEnergyDeposit();
	G4double fERes4 = fEOrig4 + fEOrig4 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig4*1000000));
    fEventAction->AddScint4(fEOrig4);
    fEventAction->AddResScint4(fERes4);
  }
  if (volume == fDetector->GetScint5()){
	G4double fEOrig5 = aStep->GetTotalEnergyDeposit();
	G4double fERes5 = fEOrig5 + fEOrig5 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig5*1000000));
    fEventAction->AddScint5(fEOrig5);
    fEventAction->AddResScint5(fERes5);
  }
  if (volume == fDetector->GetScint6()){
	G4double fEOrig6 = aStep->GetTotalEnergyDeposit();
	G4double fERes6 = fEOrig6 + fEOrig6 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig6*1000000));
    fEventAction->AddScint6(fEOrig6);
    fEventAction->AddResScint6(fERes6);
  }
  if (volume == fDetector->GetScint7()){
	G4double fEOrig7 = aStep->GetTotalEnergyDeposit();
	G4double fERes7 = fEOrig7 + fEOrig7 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig7*1000000));
    fEventAction->AddScint7(fEOrig7);
    fEventAction->AddResScint7(fERes7);
  }
  if (volume == fDetector->GetScint8()){
	G4double fEOrig8 = aStep->GetTotalEnergyDeposit();
	G4double fERes8 = fEOrig8 + fEOrig8 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig8*1000000));
    fEventAction->AddScint8(fEOrig8);
    fEventAction->AddResScint8(fERes8);
  }
  if (volume == fDetector->GetScint9()){
	G4double fEOrig9 = aStep->GetTotalEnergyDeposit();
	G4double fERes9 = fEOrig9 + fEOrig9 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig9*1000000));
    fEventAction->AddScint9(fEOrig9);
    fEventAction->AddResScint9(fERes9);
  }
  if (volume == fDetector->GetScint10()){
	G4double fEOrig10 = aStep->GetTotalEnergyDeposit();
	G4double fERes10 = fEOrig10 + fEOrig10 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig10*1000000));
    fEventAction->AddScint10(fEOrig10);
    fEventAction->AddResScint10(fERes10);
  }
  if (volume == fDetector->GetScint11()){
	G4double fEOrig11 = aStep->GetTotalEnergyDeposit();
	G4double fERes11 = fEOrig11 + fEOrig11 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig11*1000000));
    fEventAction->AddScint11(fEOrig11);
    fEventAction->AddResScint11(fERes11);
  }
  if (volume == fDetector->GetScint12()){
	G4double fEOrig12 = aStep->GetTotalEnergyDeposit();
	G4double fERes12 = fEOrig12 + fEOrig12 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig12*1000000));
    fEventAction->AddScint12(fEOrig12);
    fEventAction->AddResScint12(fERes12);
  }
  if (volume == fDetector->GetScint13()){
	G4double fEOrig13 = aStep->GetTotalEnergyDeposit();
	G4double fERes13 = fEOrig13 + fEOrig13 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig13*1000000));
    fEventAction->AddScint13(fEOrig13);
    fEventAction->AddResScint13(fERes13);
  }
  if (volume == fDetector->GetScint14()){
	G4double fEOrig14 = aStep->GetTotalEnergyDeposit();
	G4double fERes14 = fEOrig14 + fEOrig14 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig14*1000000));
    fEventAction->AddScint14(fEOrig14);
    fEventAction->AddResScint14(fERes14);
  }
  if (volume == fDetector->GetScint15()){
	G4double fEOrig15 = aStep->GetTotalEnergyDeposit();
	G4double fERes15 = fEOrig15 + fEOrig15 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig15*1000000));
    fEventAction->AddScint15(fEOrig15);
    fEventAction->AddResScint15(fERes15);
  }
  if (volume == fDetector->GetScint16()){
	G4double fEOrig16 = aStep->GetTotalEnergyDeposit();
	G4double fERes16 = fEOrig16 + fEOrig16 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig16*1000000));
    fEventAction->AddScint16(fEOrig16);
    fEventAction->AddResScint16(fERes16);
  }
  if (volume == fDetector->GetScint17()){
	G4double fEOrig17 = aStep->GetTotalEnergyDeposit();
	G4double fERes17 = fEOrig17 + fEOrig17 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig17*1000000));
    fEventAction->AddScint17(fEOrig17);
    fEventAction->AddResScint17(fERes17);
  }
  if (volume == fDetector->GetScint18()){
	G4double fEOrig18 = aStep->GetTotalEnergyDeposit();
	G4double fERes18 = fEOrig18 + fEOrig18 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig18*1000000));
    fEventAction->AddScint18(fEOrig18);
    fEventAction->AddResScint18(fERes18);
  }
  if (volume == fDetector->GetScint19()){
	G4double fEOrig19 = aStep->GetTotalEnergyDeposit();
	G4double fERes19 = fEOrig19 + fEOrig19 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig19*1000000));
    fEventAction->AddScint19(fEOrig19);
    fEventAction->AddResScint19(fERes19);
  }
  if (volume == fDetector->GetScint20()){
	G4double fEOrig20 = aStep->GetTotalEnergyDeposit();
	G4double fERes20 = fEOrig20 + fEOrig20 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig20*1000000));
    fEventAction->AddScint20(fEOrig20);
    fEventAction->AddResScint20(fERes20);
  }
  if (volume == fDetector->GetScint21()){
	G4double fEOrig21 = aStep->GetTotalEnergyDeposit();
	G4double fERes21 = fEOrig21 + fEOrig21 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig21*1000000));
    fEventAction->AddScint21(fEOrig21);
    fEventAction->AddResScint21(fERes21);
  }
  if (volume == fDetector->GetScint22()){
	G4double fEOrig22 = aStep->GetTotalEnergyDeposit();
	G4double fERes22 = fEOrig22 + fEOrig22 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig22*1000000));
    fEventAction->AddScint22(fEOrig22);
    fEventAction->AddResScint22(fERes22);
  }
  if (volume == fDetector->GetScint23()){
	G4double fEOrig23 = aStep->GetTotalEnergyDeposit();
	G4double fERes23 = fEOrig23 + fEOrig23 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig23*1000000));
    fEventAction->AddScint23(fEOrig23);
    fEventAction->AddResScint23(fERes23);
  }
  if (volume == fDetector->GetScint24()){
	G4double fEOrig24 = aStep->GetTotalEnergyDeposit();
	G4double fERes24 = fEOrig24 + fEOrig24 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig24*1000000));
    fEventAction->AddScint24(fEOrig24);
    fEventAction->AddResScint24(fERes24);
  }
  if (volume == fDetector->GetScint25()){
	G4double fEOrig25 = aStep->GetTotalEnergyDeposit();
	G4double fERes25 = fEOrig25 + fEOrig25 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig25*1000000));
    fEventAction->AddScint25(fEOrig25);
    fEventAction->AddResScint25(fERes25);
  }
  if (volume == fDetector->GetScint26()){
	G4double fEOrig26 = aStep->GetTotalEnergyDeposit();
	G4double fERes26 = fEOrig26 + fEOrig26 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig26*1000000));
    fEventAction->AddScint26(fEOrig26);
    fEventAction->AddResScint26(fERes26);
  }
  if (volume == fDetector->GetScint27()){
	G4double fEOrig27 = aStep->GetTotalEnergyDeposit();
	G4double fERes27 = fEOrig27 + fEOrig27 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig27*1000000));
    fEventAction->AddScint27(fEOrig27);
    fEventAction->AddResScint27(fERes27);
  }
  if (volume == fDetector->GetScint28()){
	G4double fEOrig28 = aStep->GetTotalEnergyDeposit();
	G4double fERes28 = fEOrig28 + fEOrig28 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig28*1000000));
    fEventAction->AddScint28(fEOrig28);
    fEventAction->AddResScint28(fERes28);
  }
  if (volume == fDetector->GetScint29()){
	G4double fEOrig29 = aStep->GetTotalEnergyDeposit();
	G4double fERes29 = fEOrig29 + fEOrig29 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig29*1000000));
    fEventAction->AddScint29(fEOrig29);
    fEventAction->AddResScint29(fERes29);
  }
  if (volume == fDetector->GetScint30()){
	G4double fEOrig30 = aStep->GetTotalEnergyDeposit();
	G4double fERes30 = fEOrig30 + fEOrig30 * fRandomGauss->fire(0.0, 0.097*(sqrt(662))/std::sqrt(fEOrig30*1000000));
    fEventAction->AddScint30(fEOrig30);
    fEventAction->AddResScint30(fERes30);
  }
 
  //if (fEOrig > 0.){
    //G4cout << "fEOrig: " << fEOrig << G4endl;
    //G4cout << "Res Func: " << (76.3/std::sqrt(fEOrig*1000000))/100 << G4endl;
  //}

//  G4double fERes = fEOrig  fRandomGauss->fire(1.0, 0.01); // for checking

  // G4double stepLength = 0.;
  // if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
  //   stepLength = aStep->GetStepLength();
     
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
