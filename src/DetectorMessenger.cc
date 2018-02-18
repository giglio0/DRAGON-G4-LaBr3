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
/// \file electromagnetic/TestEm5/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:fDetector(Det)
{ 
  fTestemDir = new G4UIdirectory("/testem/");
  fTestemDir->SetGuidance("UI commands specific to this example.");
  
  fDetDir = new G4UIdirectory("/testem/det/");
  fDetDir->SetGuidance("detector construction commands");
  
  fSetNoDetsCmd = new G4UIcmdWithAnInteger("/testem/det/setNoDets",this);
  fSetNoDetsCmd->SetGuidance("Set the number of detectors for this geometry.");
  fSetNoDetsCmd->SetGuidance("A value > 1 will create an array from an input file");
  fSetNoDetsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fCrystalDiamCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setCrystalDiam",this);
  fCrystalDiamCmd->SetGuidance("Set the diameter of the Crystal crystal");
  fCrystalDiamCmd->SetParameterName("CrystalDiam",false);
  fCrystalDiamCmd->SetRange("CrystalDiam >0.");
  fCrystalDiamCmd->SetUnitCategory("Length");
  fCrystalDiamCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fCrystalLengthCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setCrystalLength",this);
  fCrystalLengthCmd->SetGuidance("Set the length of the Crystal crystal");
  fCrystalLengthCmd->SetParameterName("CrystalLength",false);
  fCrystalLengthCmd->SetRange("CrystalLength >0.");
  fCrystalLengthCmd->SetUnitCategory("Length");
  fCrystalLengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fGapThicknessCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setGapThickness",this);
  fGapThicknessCmd->SetGuidance("Set the thickness of the gap between Al and the Crystal crystal");
  fGapThicknessCmd->SetParameterName("GapThickness",false);
  fGapThicknessCmd->SetRange("GapThickness >0.");
  fGapThicknessCmd->SetUnitCategory("Length");
  fGapThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fGapMaterialCmd = new G4UIcmdWithAString("/testem/det/setGapMaterial",this);
  fGapMaterialCmd->SetGuidance("Select Material of the Gap.");
  fGapMaterialCmd->SetParameterName("GapMaterial",false);
  fGapMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAlCaseThicknessCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setAlCaseThickness",this);
  fAlCaseThicknessCmd->SetGuidance("Set the thickness of the Al casing");
  fAlCaseThicknessCmd->SetParameterName("CaseThickness",false);
  fAlCaseThicknessCmd->SetRange("CaseThickness >0.");
  fAlCaseThicknessCmd->SetUnitCategory("Length");
  fAlCaseThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
 
  fAlFaceThicknessCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setAlFaceThickness",this);
  fAlFaceThicknessCmd->SetGuidance("Set the thickness of the Al casing");
  fAlFaceThicknessCmd->SetParameterName("FaceThickness",false);
  fAlFaceThicknessCmd->SetRange("FaceThickness >0.");
  fAlFaceThicknessCmd->SetUnitCategory("Length");
  fAlFaceThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fPbCaseThicknessCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setPbCaseThickness",this);
  fPbCaseThicknessCmd->SetGuidance("Set the thickness of the Pb shielding surround");
  fPbCaseThicknessCmd->SetParameterName("PbThickness",false);
  fPbCaseThicknessCmd->SetRange("PbThickness >0.");
  fPbCaseThicknessCmd->SetUnitCategory("Length");
  fPbCaseThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fPMTDiameterCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setPMTDiameter",this);
  fPMTDiameterCmd->SetGuidance("Set the diamter of the PMT");
  fPMTDiameterCmd->SetParameterName("PMTDiameter",false);
  fPMTDiameterCmd->SetRange("PMTDiameter >0.");
  fPMTDiameterCmd->SetUnitCategory("Length");
  fPMTDiameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fPMTLengthCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setPMTLength",this);
  fPMTLengthCmd->SetGuidance("Set the Length of the PMT");
  fPMTLengthCmd->SetParameterName("PMTLength",false);
  fPMTLengthCmd->SetRange("PMTLength >0.");
  fPMTLengthCmd->SetUnitCategory("Length");
  fPMTLengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fUpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  fUpdateCmd->SetGuidance("Update calorimeter geometry.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
  fUpdateCmd->AvailableForStates(G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fDetDir;  
  delete fTestemDir;
  delete fCrystalDiamCmd;
  delete fCrystalLengthCmd;
  delete fGapThicknessCmd;
  delete fGapMaterialCmd;
  delete fAlCaseThicknessCmd;
  delete fAlFaceThicknessCmd;
  delete fPbCaseThicknessCmd;
  delete fPMTDiameterCmd;
  delete fPMTLengthCmd;
  delete fSetNoDetsCmd;
  delete fUpdateCmd;
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if ( command == fCrystalDiamCmd )
   {fDetector->SetCrystalDiameter(fCrystalDiamCmd->GetNewDoubleValue(newValue));}

  if ( command == fCrystalLengthCmd )
   {fDetector->SetCrystalLength(fCrystalLengthCmd->GetNewDoubleValue(newValue));}

  if ( command == fGapThicknessCmd )
   {fDetector->SetGapThickness(fGapThicknessCmd->GetNewDoubleValue(newValue));}
      
  if ( command == fGapMaterialCmd)
   {fDetector->SetGapMaterial(newValue);}

  if ( command == fAlCaseThicknessCmd )
   {fDetector->SetAlCaseThickness(fAlCaseThicknessCmd->GetNewDoubleValue(newValue));}

  if ( command == fAlFaceThicknessCmd )
   {fDetector->SetAlFaceThickness(fAlFaceThicknessCmd->GetNewDoubleValue(newValue));}

  if ( command == fPbCaseThicknessCmd)
   {fDetector->SetPbCaseThickness(fPbCaseThicknessCmd->GetNewDoubleValue(newValue));}

  if ( command == fPMTDiameterCmd)
   {fDetector->SetPMTDiameter(fPMTDiameterCmd->GetNewDoubleValue(newValue));}

  if ( command == fPMTLengthCmd)
   {fDetector->SetPMTLength(fPMTLengthCmd->GetNewDoubleValue(newValue));}

  if ( command == fSetNoDetsCmd)
   {fDetector->SetNoOfDetectors(fSetNoDetsCmd->GetNewIntValue(newValue));}

  if  ( command == fUpdateCmd )
   {fDetector->UpdateGeometry(); }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
