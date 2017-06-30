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
  fDetDir->SetGuidance("Detector construction commands.");

  fWorldMaterialCmd = new G4UIcmdWithAString("/testem/det/setWorldMaterial",this);
  fWorldMaterialCmd->SetGuidance("Select the material of the world.");
  fWorldMaterialCmd->SetParameterName("WorldMaterial",true);
  fWorldMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fDetectorMaterialCmd = new G4UIcmdWithAString("/testem/det/setDetectorMaterial",this);
  fDetectorMaterialCmd->SetGuidance("Select the material of the detector.");
  fDetectorMaterialCmd->SetParameterName("DetectorMaterial",true);
  fDetectorMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fDetectorGeometryCmd = new G4UIcmdWithAnInteger("/testem/det/setDetectorGeometry",this);
  fDetectorGeometryCmd->SetGuidance("Select geometry of the detector.");
  fDetectorGeometryCmd->SetParameterName("DetectorGeometry", true);
  fDetectorGeometryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fDetectorDiameterCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setDetectorDiameter",this);
  fDetectorDiameterCmd->SetGuidance("Set the diameter of the crystal.");
  fDetectorDiameterCmd->SetParameterName("DetectorDiameter",true);
  fDetectorDiameterCmd->SetRange("DetectorDiameter >0.");
  fDetectorDiameterCmd->SetUnitCategory("Length");
  fDetectorDiameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fDetectorLengthCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setDetectorLength",this);
  fDetectorLengthCmd->SetGuidance("Set the length of the crystal.");
  fDetectorLengthCmd->SetParameterName("DetectorLength",true);
  fDetectorLengthCmd->SetRange("DetectorLength >0.");
  fDetectorLengthCmd->SetUnitCategory("Length");
  fDetectorLengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  TemperatureCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setTemperature",this);
  TemperatureCmd->SetGuidance("Set the temperature of the detector.");
  TemperatureCmd->SetParameterName("Temperature",true);
  TemperatureCmd->SetRange("Temperature >=0.");
  TemperatureCmd->SetUnitCategory("Temperature");
  TemperatureCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  PressureCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setPressure",this);
  PressureCmd->SetGuidance("Set the pressure of the detector.");
  PressureCmd->SetParameterName("Pressure",true);
  PressureCmd->SetRange("Pressure >=0.");
  PressureCmd->SetUnitCategory("Pressure");
  PressureCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fGapThicknessCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setGapThickness",this);
  fGapThicknessCmd->SetGuidance("Set the thickness of the gap between the aluminum case and the crystal.");
  fGapThicknessCmd->SetParameterName("GapThickness",true);
  fGapThicknessCmd->SetRange("GapThickness >0.");
  fGapThicknessCmd->SetUnitCategory("Length");
  fGapThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  fCrystalMaterialCmd = new G4UIcmdWithAString("/testem/det/setCrystalMaterial",this);
  fCrystalMaterialCmd->SetGuidance("Select the material of the crystal.");
  fCrystalMaterialCmd->SetParameterName("CrystalMaterial",true);
  fCrystalMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fGapMaterialCmd = new G4UIcmdWithAString("/testem/det/setGapMaterial",this);
  fGapMaterialCmd->SetGuidance("Select the material of the gap.");
  fGapMaterialCmd->SetParameterName("GapMaterial",true);
  fGapMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fFaceGapMaterialCmd = new G4UIcmdWithAString("/testem/det/setFaceGapMaterial",this);
  fFaceGapMaterialCmd->SetGuidance("Select the material of the face of the gap.");
  fFaceGapMaterialCmd->SetParameterName("FaceGapMaterial",true);
  fFaceGapMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fAlCaseMaterialCmd = new G4UIcmdWithAString("/testem/det/setAlCaseMaterial",this);
  fAlCaseMaterialCmd->SetGuidance("Select the material of the aluminum case.");
  fAlCaseMaterialCmd->SetParameterName("AlCaseMaterial",true);
  fAlCaseMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fFaceAlCaseMaterialCmd = new G4UIcmdWithAString("/testem/det/setFaceAlCaseMaterial",this);
  fFaceAlCaseMaterialCmd->SetGuidance("Select the material of the face of the aluminum case.");
  fFaceAlCaseMaterialCmd->SetParameterName("FaceAlCaseMaterial",true);
  fFaceAlCaseMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fPbCaseMaterialCmd = new G4UIcmdWithAString("/testem/det/setPbCaseMaterial",this);
  fPbCaseMaterialCmd->SetGuidance("Select the material of the lead case.");
  fPbCaseMaterialCmd->SetParameterName("PbCaseMaterial",true);
  fPbCaseMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fPbCollarMaterialCmd = new G4UIcmdWithAString("/testem/det/setPbCollarMaterial",this);
  fPbCollarMaterialCmd->SetGuidance("Select the material of the lead collar.");
  fPbCollarMaterialCmd->SetParameterName("PbCollarMaterial",true);
  fPbCollarMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fPMTMaterialCmd = new G4UIcmdWithAString("/testem/det/setPMTMaterial",this);
  fPMTMaterialCmd->SetGuidance("Select the material of the photomultiplier tube.");
  fPMTMaterialCmd->SetParameterName("PMTMaterial",true);
  fPMTMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fPMTWinMaterialCmd = new G4UIcmdWithAString("/testem/det/setPMTWinMaterial",this);
  fPMTWinMaterialCmd->SetGuidance("Select the material of the photomultiplier tube window.");
  fPMTWinMaterialCmd->SetParameterName("PMTWinMaterial",true);
  fPMTWinMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAlCaseThicknessCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setAlCaseThickness",this);
  fAlCaseThicknessCmd->SetGuidance("Set the thickness of the aluminum casing.");
  fAlCaseThicknessCmd->SetParameterName("CaseThickness",true);
  fAlCaseThicknessCmd->SetRange("CaseThickness >0.");
  fAlCaseThicknessCmd->SetUnitCategory("Length");
  fAlCaseThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
 
  fPbCaseThicknessCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setPbCaseThickness",this);
  fPbCaseThicknessCmd->SetGuidance("Set the thickness of the surrounding lead shield.");
  fPbCaseThicknessCmd->SetParameterName("PbThickness",true);
  fPbCaseThicknessCmd->SetRange("PbThickness >0.");
  fPbCaseThicknessCmd->SetUnitCategory("Length");
  fPbCaseThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fPMTDiameterCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setPMTDiameter",this);
  fPMTDiameterCmd->SetGuidance("Set the diameter of the photomultiplier tube.");
  fPMTDiameterCmd->SetParameterName("PMTDiameter",true);
  fPMTDiameterCmd->SetRange("PMTDiameter >0.");
  fPMTDiameterCmd->SetUnitCategory("Length");
  fPMTDiameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  fPMTLengthCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setPMTLength",this);
  fPMTLengthCmd->SetGuidance("Set the length of the photomultiplier tube.");
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
  delete TemperatureCmd;
  delete PressureCmd;
  delete fGapThicknessCmd;
  delete fCrystalMaterialCmd;
  delete fGapMaterialCmd;
  delete fFaceGapMaterialCmd;
  delete fAlCaseMaterialCmd;
  delete fFaceAlCaseMaterialCmd;
  delete fPbCaseMaterialCmd;
  delete fPbCollarMaterialCmd;
  delete fPMTMaterialCmd;
  delete fPMTWinMaterialCmd;
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
   
  if ( command == TemperatureCmd )
   {fDetector->SetTemperature(TemperatureCmd->GetNewDoubleValue(newValue));}
   
  if ( command == PressureCmd )
   {fDetector->SetPressure(PressureCmd->GetNewDoubleValue(newValue));}
   
  if ( command == fGapThicknessCmd )
   {fDetector->SetGapThickness(fGapThicknessCmd->GetNewDoubleValue(newValue));}
      
  if ( command == fDetectorMaterialCmd )
   {fDetector->SetDetectorMaterial(newValue);}
   
  if ( command == fCrystalMaterialCmd )
   {fDetector->SetCrystalMaterial(newValue);}
   
  if ( command == fGapMaterialCmd)
   {fDetector->SetGapMaterial(newValue);}
   
  if ( command == fFaceGapMaterialCmd)
   {fDetector->SetFaceGapMaterial(newValue);}
   
  if ( command == fAlCaseMaterialCmd )
   {fDetector->SetAlCaseMaterial(newValue);}
   
  if ( command == fFaceAlCaseMaterialCmd)
   {fDetector->SetFaceAlCaseMaterial(newValue);}
   
  if ( command == fPbCaseMaterialCmd )
   {fDetector->SetPbCaseMaterial(newValue);}
   
  if ( command == fPbCollarMaterialCmd)
   {fDetector->SetPbCollarMaterial(newValue);}
   
  if ( command == fPMTMaterialCmd )
   {fDetector->SetPMTMaterial(newValue);}
   
  if ( command == fPMTWinMaterialCmd )
   {fDetector->SetPMTWinMaterial(newValue);}

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
