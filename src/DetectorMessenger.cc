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

  fWorldMaterialCmd = new G4UIcmdWithAString("/testem/det/setWorldMaterial",this);
  fWorldMaterialCmd->SetGuidance("Select Material of the World.");
  fWorldMaterialCmd->SetParameterName("WorldMaterial",true);
  fWorldMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fDetectorMaterialCmd = new G4UIcmdWithAString("/testem/det/setDetectorMaterial",this);
  fDetectorMaterialCmd->SetGuidance("Select Material of the Detector.");
  fDetectorMaterialCmd->SetParameterName("DetectorMaterial",true);
  fDetectorMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fDetectorGeometryCmd = new G4UIcmdWithAnInteger("/testem/det/setDetectorGeometry",this);
  fDetectorGeometryCmd->SetGuidance("Select Geometry of the Detector.");
  fDetectorGeometryCmd->SetParameterName("DetectorGeometry", true);
  fDetectorGeometryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fDetectorDiameterCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setDetectorDiameter",this);
  fDetectorDiameterCmd->SetGuidance("Set the diameter of the crystal");
  fDetectorDiameterCmd->SetParameterName("DetectorDiameter",true);
  fDetectorDiameterCmd->SetRange("DetectorDiameter >0.");
  fDetectorDiameterCmd->SetUnitCategory("Length");
  fDetectorDiameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fDetectorLengthCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setDetectorLength",this);
  fDetectorLengthCmd->SetGuidance("Set the length of the crystal");
  fDetectorLengthCmd->SetParameterName("DetectorLength",true);
  fDetectorLengthCmd->SetRange("DetectorLength >0.");
  fDetectorLengthCmd->SetUnitCategory("Length");
  fDetectorLengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fGapThicknessCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setGapThickness",this);
  fGapThicknessCmd->SetGuidance("Set the thickness of the gap between Al and the crystal");
  fGapThicknessCmd->SetParameterName("GapThickness",true);
  fGapThicknessCmd->SetRange("GapThickness >0.");
  fGapThicknessCmd->SetUnitCategory("Length");
  fGapThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fGapMaterialCmd = new G4UIcmdWithAString("/testem/det/setGapMaterial",this);
  fGapMaterialCmd->SetGuidance("Select Material of the Gap.");
  fGapMaterialCmd->SetParameterName("GapMaterial",true);
  fGapMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAlCaseThicknessCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setAlCaseThickness",this);
  fAlCaseThicknessCmd->SetGuidance("Set the thickness of the Al casing");
  fAlCaseThicknessCmd->SetParameterName("CaseThickness",true);
  fAlCaseThicknessCmd->SetRange("CaseThickness >0.");
  fAlCaseThicknessCmd->SetUnitCategory("Length");
  fAlCaseThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
 
  fPbCaseThicknessCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setPbCaseThickness",this);
  fPbCaseThicknessCmd->SetGuidance("Set the thickness of the Pb shielding surround");
  fPbCaseThicknessCmd->SetParameterName("PbThickness",true);
  fPbCaseThicknessCmd->SetRange("PbThickness >0.");
  fPbCaseThicknessCmd->SetUnitCategory("Length");
  fPbCaseThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fPMTDiameterCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setPMTDiameter",this);
  fPMTDiameterCmd->SetGuidance("Set the diamter of the PMT");
  fPMTDiameterCmd->SetParameterName("PMTDiameter",true);
  fPMTDiameterCmd->SetRange("PMTDiameter >0.");
  fPMTDiameterCmd->SetUnitCategory("Length");
  fPMTDiameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fPMTLengthCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setPMTLength",this);
  fPMTLengthCmd->SetGuidance("Set the Length of the PMT");
  fPMTLengthCmd->SetParameterName("PMTLength",true);
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
  delete fWorldMaterialCmd;
  delete fDetectorMaterialCmd;
  delete fDetectorGeometryCmd;
  delete fDetectorDiameterCmd;
  delete fDetectorLengthCmd;
  delete fGapThicknessCmd;
  delete fGapMaterialCmd;
  delete fAlCaseThicknessCmd;
  delete fPbCaseThicknessCmd;
  delete fPMTDiameterCmd;
  delete fPMTLengthCmd;
  delete fUpdateCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if ( command == fWorldMaterialCmd )
   {fDetector->SetWorldMaterial(newValue);}
  
  if ( command == fDetectorMaterialCmd )
   {fDetector->SetDetectorMaterial(newValue);}

  if ( command == fDetectorDiameterCmd )
   {fDetector->SetDetectorDiameter(fDetectorDiameterCmd->GetNewDoubleValue(newValue));}

  if ( command == fDetectorGeometryCmd )
   {fDetector->SetDetectorGeometry(fDetectorGeometryCmd->GetNewIntValue(newValue));}

  if ( command == fDetectorLengthCmd )
   {fDetector->SetDetectorLength(fDetectorLengthCmd->GetNewDoubleValue(newValue));}

  if ( command == fGapThicknessCmd )
   {fDetector->SetGapThickness(fGapThicknessCmd->GetNewDoubleValue(newValue));}
      
  if ( command == fGapMaterialCmd)
   {fDetector->SetGapMaterial(newValue);}

  if ( command == fAlCaseThicknessCmd )
   {fDetector->SetAlCaseThickness(fAlCaseThicknessCmd->GetNewDoubleValue(newValue));}

  if ( command == fPbCaseThicknessCmd)
   {fDetector->SetPbCaseThickness(fPbCaseThicknessCmd->GetNewDoubleValue(newValue));}

  if ( command == fPMTDiameterCmd)
   {fDetector->SetPMTDiameter(fPMTDiameterCmd->GetNewDoubleValue(newValue));}

  if ( command == fPMTLengthCmd)
   {fDetector->SetPMTLength(fPMTLengthCmd->GetNewDoubleValue(newValue));}

  if  ( command == fUpdateCmd )
   {fDetector->UpdateGeometry(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
