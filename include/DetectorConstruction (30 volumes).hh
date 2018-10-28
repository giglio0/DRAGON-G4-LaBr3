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
/// \file electromagnetic/TestEm5/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Trd.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4AssemblyVolume.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

class G4Box;
class G4SubtractionSolid;
class G4Tubs;
class G4CutTubs;
class G4Trd;
class G4Polyhedra;
class G4VPhysicalVolume;
class G4Material;
class G4String;
class G4MaterialPropertiesTable;
class G4MaterialCutsCouple;
class G4UniformMagField;
class DetectorMessenger;
class G4AssemblyVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
   ~DetectorConstruction();

  public:
  
     void SetWorldMaterial           (G4String);
     void SetDetectorMaterial        (G4String);
     void SetGasTargetMaterial       (G4String);
     void SetInnerTrapezoidMaterial  (G4String);
     void SetDetectorGeometry           (G4int);
     void SetTemperature             (G4double);
     void SetPressure                (G4double);
     void SetTotalDetectorDiameter   (G4double);
     void SetDetectorDiameter        (G4double);
     void SetTotalDetectorLength     (G4double);
     void SetDetectorLength	         (G4double);
     void SetGapThickness	         (G4double);
     void SetCrystalMaterial         (G4String);
     void SetGapMaterial             (G4String);
     void SetFaceGapMaterial         (G4String);
     void SetAlCaseMaterial          (G4String);   
     void SetFaceAlCaseMaterial      (G4String);
     void SetfAlFaceThickness        (G4String);
     void SetPbCaseMaterial          (G4String);
     void SetPbCollarMaterial        (G4String);
     void SetPMTMaterial             (G4String);
     void SetPMTWinMaterial          (G4String);
     void SetPhotoCathodeMaterial    (G4String);
     void SetPMTIntMaterial          (G4String);
     void SetAlCaseThickness         (G4double);
     void SetPbCaseThickness         (G4double);
     void SetPMTDiameter	         (G4double);
     void SetPMTLength		         (G4double);
     void SetMaterialPropertiesTable (G4String);
     
     virtual G4VPhysicalVolume* Construct();
     //virtual void ConstructSDandField();
     
     void UpdateGeometry();

  public:
  // Variable for BGO, Variable1 for LaBr3
     void PrintCalorParameters();
     
     G4Material* GetWorldMaterial()     {return fWorldMaterial;};
     G4double    GetWorldSizeX()        {return fWorldSizeX;};
     G4double    GetWorldSizeYZ()       {return fWorldSizeYZ;};

     G4Material* GetDetectorMaterial()           {return fDetectorMaterial;};
     G4Material* GetGasTargetMaterial()          {return fGasTargetMaterial;};
     G4Material* GetInnerTrapezoidMaterial()     {return fInnerTrapezoidMaterial;};
     G4double    GetTemperature()                {return Temperature;};
     G4double    GetPressure()                   {return Pressure;};
     G4int       GetDetectorGeometry()           {return fDetectorGeometry;};
     G4double    GetTotalDetectorDiameter()      {return fTotalDetectorDiameter;};
     G4double	 GetDetectorDiameter()		     {return fDetectorDiameter;};
     G4double    GetTotalDetectorLength()        {return fTotalDetectorLength;}; 
     G4double	 GetDetectorLength()		     {return fDetectorLength;};
     G4double	 GetGapThickness()		         {return fGapThickness;}; 
     G4double    GetGapFaceThickness()           {return fGapFaceThickness;};
     G4Material* GetCrystalMaterial()            {return fCrystalMaterial;};
     G4Material* GetGapMaterial()    	         {return fGapMaterial;};
     G4Material* GetFaceGapMaterial()            {return fFaceGapMaterial;};
     G4Material* GetAlCaseMaterial()             {return fAlCaseMaterial;};
     G4Material* GetFaceAlCaseMaterial()         {return fFaceAlCaseMaterial;};
     G4double    GetAlFaceThickness()            {return fAlFaceThickness;};
     G4Material* GetPbCaseMaterial()             {return fPbCaseMaterial;};
     G4Material* GetPbCollarMaterial()           {return fAlCaseMaterial;};
     G4Material* GetPMTMaterial()                {return fPMTMaterial;};
     G4Material* GetPMTWinMaterial()             {return fPMTWinMaterial;};
     G4Material* GetPhotoCathodeMaterial()       {return fPhotoCathodeMaterial;};
     G4Material* GetPMTIntMaterial()             {return fPMTIntMaterial;};
     G4double	 GetAlCaseThickness()	         {return fAlCaseThickness;}; 
     G4double	 GetPbCaseThickness()	         {return fPbCaseThickness;}; 
     G4double	 GetPMTDiameter()		         {return fPMTDiameter;}; 
     G4double	 GetPMTLength()			         {return fPMTLength;};
     G4VPhysicalVolume*    GetCrystal()           {return fPhysiCrystal;};
     G4VPhysicalVolume*    GetCrystal1()          {return fPhysiCrystal1;};
	 G4VPhysicalVolume*    GetCrystal2()          {return fPhysiCrystal2;};  
	 G4VPhysicalVolume*    GetCrystal3()          {return fPhysiCrystal3;};
	 G4VPhysicalVolume*    GetCrystal4()          {return fPhysiCrystal4;};  
	 G4VPhysicalVolume*    GetCrystal5()          {return fPhysiCrystal5;};
	 G4VPhysicalVolume*    GetCrystal6()          {return fPhysiCrystal6;};  
	 G4VPhysicalVolume*    GetCrystal7()          {return fPhysiCrystal7;};
	 G4VPhysicalVolume*    GetCrystal8()          {return fPhysiCrystal8;};  
	 G4VPhysicalVolume*    GetCrystal9()          {return fPhysiCrystal9;};
	 G4VPhysicalVolume*    GetCrystal10()          {return fPhysiCrystal10;};  
	 G4VPhysicalVolume*    GetCrystal11()          {return fPhysiCrystal11;};
	 G4VPhysicalVolume*    GetCrystal12()          {return fPhysiCrystal12;};  
	 G4VPhysicalVolume*    GetCrystal13()          {return fPhysiCrystal13;};
	 G4VPhysicalVolume*    GetCrystal14()          {return fPhysiCrystal14;};  
	 G4VPhysicalVolume*    GetCrystal15()          {return fPhysiCrystal15;};
	 G4VPhysicalVolume*    GetCrystal16()          {return fPhysiCrystal16;};  
	 G4VPhysicalVolume*    GetCrystal17()          {return fPhysiCrystal17;};
	 G4VPhysicalVolume*    GetCrystal18()          {return fPhysiCrystal18;};  
	 G4VPhysicalVolume*    GetCrystal19()          {return fPhysiCrystal19;};
	 G4VPhysicalVolume*    GetCrystal20()          {return fPhysiCrystal20;};  
	 G4VPhysicalVolume*    GetCrystal21()          {return fPhysiCrystal21;};
	 G4VPhysicalVolume*    GetCrystal22()          {return fPhysiCrystal22;};  
	 G4VPhysicalVolume*    GetCrystal23()          {return fPhysiCrystal23;};
	 G4VPhysicalVolume*    GetCrystal24()          {return fPhysiCrystal24;};  
	 G4VPhysicalVolume*    GetCrystal25()          {return fPhysiCrystal25;};
	 G4VPhysicalVolume*    GetCrystal26()          {return fPhysiCrystal26;};  
	 G4VPhysicalVolume*    GetCrystal27()          {return fPhysiCrystal27;};
	 G4VPhysicalVolume*    GetCrystal28()          {return fPhysiCrystal28;};  
	 G4VPhysicalVolume*    GetCrystal29()          {return fPhysiCrystal29;};
	 G4VPhysicalVolume*    GetCrystal30()          {return fPhysiCrystal30;};  

     const G4VPhysicalVolume* GetphysiWorld() {return fPhysiWorld;};
     const G4VPhysicalVolume* GetScint()  {return fPhysiCrystal;};
     const G4VPhysicalVolume* GetScint1() {return fPhysiCrystal1;};
     const G4VPhysicalVolume* GetScint2() {return fPhysiCrystal2;};
     const G4VPhysicalVolume* GetScint3() {return fPhysiCrystal3;};
     const G4VPhysicalVolume* GetScint4() {return fPhysiCrystal4;};
     const G4VPhysicalVolume* GetScint5() {return fPhysiCrystal5;};
     const G4VPhysicalVolume* GetScint6() {return fPhysiCrystal6;};
     const G4VPhysicalVolume* GetScint7() {return fPhysiCrystal7;};
     const G4VPhysicalVolume* GetScint8() {return fPhysiCrystal8;};
     const G4VPhysicalVolume* GetScint9() {return fPhysiCrystal9;};
     const G4VPhysicalVolume* GetScint10() {return fPhysiCrystal10;};
     const G4VPhysicalVolume* GetScint11() {return fPhysiCrystal11;};
     const G4VPhysicalVolume* GetScint12() {return fPhysiCrystal12;};
     const G4VPhysicalVolume* GetScint13() {return fPhysiCrystal13;};
     const G4VPhysicalVolume* GetScint14() {return fPhysiCrystal14;};
     const G4VPhysicalVolume* GetScint15() {return fPhysiCrystal15;};
     const G4VPhysicalVolume* GetScint16() {return fPhysiCrystal16;};
     const G4VPhysicalVolume* GetScint17() {return fPhysiCrystal17;};
     const G4VPhysicalVolume* GetScint18() {return fPhysiCrystal18;};
     const G4VPhysicalVolume* GetScint19() {return fPhysiCrystal19;};
     const G4VPhysicalVolume* GetScint20() {return fPhysiCrystal20;};
     const G4VPhysicalVolume* GetScint21() {return fPhysiCrystal21;};
     const G4VPhysicalVolume* GetScint22() {return fPhysiCrystal22;};
     const G4VPhysicalVolume* GetScint23() {return fPhysiCrystal23;};
     const G4VPhysicalVolume* GetScint24() {return fPhysiCrystal24;};
     const G4VPhysicalVolume* GetScint25() {return fPhysiCrystal25;};
     const G4VPhysicalVolume* GetScint26() {return fPhysiCrystal26;};
     const G4VPhysicalVolume* GetScint27() {return fPhysiCrystal27;};
     const G4VPhysicalVolume* GetScint28() {return fPhysiCrystal28;};
     const G4VPhysicalVolume* GetScint29() {return fPhysiCrystal29;};
     const G4VPhysicalVolume* GetScint30() {return fPhysiCrystal30;};

  private:
     G4VisAttributes* fCyanVisAtt;
     G4VisAttributes* fYellowVisAtt;
     G4VisAttributes* fMagnetaVisAtt;
     G4VisAttributes* fOrangeVisAtt;
     G4VisAttributes* fBlueVisAtt;
     G4VisAttributes* fGreenVisAtt;
     G4VisAttributes* fRedVisAtt;
     G4VisAttributes* fGreyVisAtt;
     G4VisAttributes* fWhiteVisAtt;
     G4VisAttributes* fBlackVisAtt;
     G4VisAttributes* fWireFrameVisAtt; 
     G4VisAttributes* fAuxEdgeVisAtt;

     // World
     G4Box*             fSolidWorld;
     G4LogicalVolume*   fLogicWorld;
     G4VPhysicalVolume* fPhysiWorld;
     G4Material*        fWorldMaterial;
     G4double           fWorldSizeX;
     G4double           fWorldSizeYZ;

     // Gas Target
     // Outer Rectangle Target Box
     G4Box*             fSolidOutTBox;
     G4LogicalVolume*   fLogicOutTBox;
     G4VPhysicalVolume* fPhysiOutTBox;
     G4double 			fOutTBoxX;
     G4double			fOutTBoxY;
     G4double			fOutTBoxZ;
     G4double 			fposEastArray;
     G4Box*             fSolidInTBox;
     G4LogicalVolume*   fLogicInTBox;
     G4VPhysicalVolume* fPhysiInTBox;
     G4Tubs*			fSolidEntHoleBox;
     G4LogicalVolume*	fLogicEntHoleBox;
     G4VPhysicalVolume*	fPhysiEntHoleBox;
     G4Tubs*			fSolidExitHoleBox;
     G4LogicalVolume*	fLogicExitHoleBox;
     G4VPhysicalVolume*	fPhysiExitHoleBox;
     G4Material*        fGasTargetMaterial;
          
     // Inner Trapezoidal Target
     G4Trd*		fSolidOutCell;
     G4LogicalVolume*	fLogicOutCell;
     G4VPhysicalVolume*	fPhysiOutCell;
     G4Trd*		fSolidInCell;
     G4LogicalVolume*	fLogicInCell;
     G4VPhysicalVolume*	fPhysiInCell;
     // Entrance Collimators
     G4Tubs*		fSolidEntHole;
     G4LogicalVolume*	fLogicEntHole;
     G4VPhysicalVolume*	fPhysiEntHole;
     G4Tubs*		fSolidEntDisk;
     G4LogicalVolume*	fLogicEntDisk;
     G4VPhysicalVolume*	fPhysiEntDisk;
     G4CutTubs*		fSolidEntCol;
     G4LogicalVolume*	fLogicEntCol;
     G4VPhysicalVolume*	fPhysiEntCol;
     // Exit Collimators
     G4Tubs*		fSolidExitHole;
     G4LogicalVolume*	fLogicExitHole;
     G4VPhysicalVolume*	fPhysiExitHole;
     G4Tubs*		fSolidExitDisk;
     G4LogicalVolume*	fLogicExitDisk;
     G4VPhysicalVolume*	fPhysiExitDisk;
     G4CutTubs*		fSolidExitCol;
     G4LogicalVolume*	fLogicExitCol;
     G4VPhysicalVolume*	fPhysiExitCol;
     G4Material*        fInnerTrapezoidMaterial;

     G4double               Temperature;
     G4double               Pressure;
     
     //LaBr3 Components
     //LaBr3 Detector Mother Holder
     G4Tubs* 		        fSolidDetector;
     G4Tubs* 		        fSolidDetector1;
     G4Polyhedra*           fSolidHDetector;
     G4Polyhedra*           fSolidHDetector1;
     G4Polyhedra*           fSolidHDetector2;
     G4Polyhedra*           fSolidHDetector3;
     G4Polyhedra*           fSolidHDetector4;
     G4Polyhedra*           fSolidHDetector5;
     G4Polyhedra*           fSolidHDetector6;
     G4Polyhedra*           fSolidHDetector7;
     G4Polyhedra*           fSolidHDetector8;
     G4Polyhedra*           fSolidHDetector9;
     G4Polyhedra*           fSolidHDetector10;
     G4Polyhedra*           fSolidHDetector11;
     G4Polyhedra*           fSolidHDetector12;
     G4Polyhedra*           fSolidHDetector13;
     G4Polyhedra*           fSolidHDetector14;
     G4Polyhedra*           fSolidHDetector15;
     G4Polyhedra*           fSolidHDetector16;
     G4Polyhedra*           fSolidHDetector17;
     G4Polyhedra*           fSolidHDetector18;
     G4Polyhedra*           fSolidHDetector19;
     G4Polyhedra*           fSolidHDetector20;
     G4Polyhedra*           fSolidHDetector21;
     G4Polyhedra*           fSolidHDetector22;
     G4Polyhedra*           fSolidHDetector23;
     G4Polyhedra*           fSolidHDetector24;
     G4Polyhedra*           fSolidHDetector25;
     G4Polyhedra*           fSolidHDetector26;
     G4Polyhedra*           fSolidHDetector27;
     G4Polyhedra*           fSolidHDetector28;
     G4Polyhedra*           fSolidHDetector29;
     G4Polyhedra*           fSolidHDetector30;
     G4LogicalVolume*	    fLogicDetector;
     G4LogicalVolume*	    fLogicDetector1;
     G4LogicalVolume*       fLogicDetector2;
     G4LogicalVolume*	    fLogicDetector3;
     G4LogicalVolume*	    fLogicDetector4;
     G4LogicalVolume*       fLogicDetector5;
     G4LogicalVolume*	    fLogicDetector6;
     G4LogicalVolume*	    fLogicDetector7;
     G4LogicalVolume*       fLogicDetector8;
     G4LogicalVolume*	    fLogicDetector9;
     G4LogicalVolume*       fLogicDetector10;
     G4LogicalVolume*	    fLogicDetector11;
     G4LogicalVolume*	    fLogicDetector12;
     G4LogicalVolume*       fLogicDetector13;
     G4LogicalVolume*	    fLogicDetector14;
     G4LogicalVolume*	    fLogicDetector15;
     G4LogicalVolume*       fLogicDetector16;
     G4LogicalVolume*	    fLogicDetector17;
     G4LogicalVolume*       fLogicDetector18;
     G4LogicalVolume*	    fLogicDetector19;
     G4LogicalVolume*	    fLogicDetector20;
     G4LogicalVolume*       fLogicDetector21;
     G4LogicalVolume*	    fLogicDetector22;
     G4LogicalVolume*	    fLogicDetector23;
     G4LogicalVolume*       fLogicDetector24;
     G4LogicalVolume*	    fLogicDetector25;
     G4LogicalVolume*       fLogicDetector26;
     G4LogicalVolume*	    fLogicDetector27;
     G4LogicalVolume*	    fLogicDetector28;
     G4LogicalVolume*       fLogicDetector29;
     G4LogicalVolume*	    fLogicDetector30;
     G4VPhysicalVolume*     fPhysiDetector;
     G4VPhysicalVolume*     fPhysiDetector1;
     G4VPhysicalVolume*     fPhysiDetector2;
     G4VPhysicalVolume*     fPhysiDetector3;
     G4VPhysicalVolume*     fPhysiDetector4;
     G4VPhysicalVolume*     fPhysiDetector5;
     G4VPhysicalVolume*     fPhysiDetector6;
     G4VPhysicalVolume*     fPhysiDetector7;
     G4VPhysicalVolume*     fPhysiDetector8;
     G4VPhysicalVolume*     fPhysiDetector9;
     G4VPhysicalVolume*     fPhysiDetector10;
     G4VPhysicalVolume*     fPhysiDetector11;
     G4VPhysicalVolume*     fPhysiDetector12;
     G4VPhysicalVolume*     fPhysiDetector13;
     G4VPhysicalVolume*     fPhysiDetector14;
     G4VPhysicalVolume*     fPhysiDetector15;
     G4VPhysicalVolume*     fPhysiDetector16;
     G4VPhysicalVolume*     fPhysiDetector17;
     G4VPhysicalVolume*     fPhysiDetector18;
     G4VPhysicalVolume*     fPhysiDetector19;
     G4VPhysicalVolume*     fPhysiDetector20;
     G4VPhysicalVolume*     fPhysiDetector21;
     G4VPhysicalVolume*     fPhysiDetector22;
     G4VPhysicalVolume*     fPhysiDetector23;
     G4VPhysicalVolume*     fPhysiDetector24;
     G4VPhysicalVolume*     fPhysiDetector25;
     G4VPhysicalVolume*     fPhysiDetector26;
     G4VPhysicalVolume*     fPhysiDetector27;
     G4VPhysicalVolume*     fPhysiDetector28;
     G4VPhysicalVolume*     fPhysiDetector29;
     G4VPhysicalVolume*     fPhysiDetector30;
     G4double               fTotalDetectorDiameter;
     G4double 		        fDetectorDiameter;
     G4double               fTotalDetectorLength;
     G4double 		        fDetectorLength;
     G4double               fTotalDetectorDiameter1;
     G4double 		        fDetectorDiameter1;
     G4double               fTotalDetectorLength1;
     G4double 		        fDetectorLength1;
     G4double		        fZposDetector;
     G4double 				fAirGap;
     G4int                  fDetectorGeometry;
     G4Material*            fDetectorMaterial;
     
     // Crystal
     G4Tubs* 		        fSolidCrystal;
     G4Tubs* 		        fSolidCrystal1;
     G4Polyhedra*           fSolidHCrystal;
     G4Polyhedra*           fSolidHCrystal1;
     G4Polyhedra*           fSolidHCrystal2;
     G4Polyhedra*           fSolidHCrystal3;
     G4Polyhedra*           fSolidHCrystal4;
     G4Polyhedra*           fSolidHCrystal5;
     G4Polyhedra*           fSolidHCrystal6;
     G4Polyhedra*           fSolidHCrystal7;
     G4Polyhedra*           fSolidHCrystal8;
     G4Polyhedra*           fSolidHCrystal9;
     G4Polyhedra*           fSolidHCrystal10;
     G4Polyhedra*           fSolidHCrystal11;
     G4Polyhedra*           fSolidHCrystal12;
     G4Polyhedra*           fSolidHCrystal13;
     G4Polyhedra*           fSolidHCrystal14;
     G4Polyhedra*           fSolidHCrystal15;
     G4Polyhedra*           fSolidHCrystal16;
     G4Polyhedra*           fSolidHCrystal17;
     G4Polyhedra*           fSolidHCrystal18;
     G4Polyhedra*           fSolidHCrystal19;
     G4Polyhedra*           fSolidHCrystal20;
     G4Polyhedra*           fSolidHCrystal21;
     G4Polyhedra*           fSolidHCrystal22;
     G4Polyhedra*           fSolidHCrystal23;
     G4Polyhedra*           fSolidHCrystal24;
     G4Polyhedra*           fSolidHCrystal25;
     G4Polyhedra*           fSolidHCrystal26;
     G4Polyhedra*           fSolidHCrystal27;
     G4Polyhedra*           fSolidHCrystal28;
     G4Polyhedra*           fSolidHCrystal29;
     G4Polyhedra*           fSolidHCrystal30;
     G4LogicalVolume*	    fLogicCrystal;
     G4LogicalVolume*       fLogicCrystal1;
     G4LogicalVolume*	    fLogicCrystal2;
     G4LogicalVolume*	    fLogicCrystal3;
     G4LogicalVolume*	    fLogicCrystal4;
     G4LogicalVolume*	    fLogicCrystal5;
     G4LogicalVolume*	    fLogicCrystal6;
     G4LogicalVolume*	    fLogicCrystal7;
     G4LogicalVolume*	    fLogicCrystal8;
     G4LogicalVolume*	    fLogicCrystal9;
     G4LogicalVolume*	    fLogicCrystal10;
     G4LogicalVolume*	    fLogicCrystal11;
     G4LogicalVolume*	    fLogicCrystal12;
     G4LogicalVolume*	    fLogicCrystal13;
     G4LogicalVolume*	    fLogicCrystal14;
     G4LogicalVolume*	    fLogicCrystal15;
     G4LogicalVolume*	    fLogicCrystal16;
     G4LogicalVolume*	    fLogicCrystal17;
     G4LogicalVolume*	    fLogicCrystal18;
     G4LogicalVolume*	    fLogicCrystal19;
     G4LogicalVolume*	    fLogicCrystal20;
     G4LogicalVolume*	    fLogicCrystal21;
     G4LogicalVolume*	    fLogicCrystal22;
     G4LogicalVolume*	    fLogicCrystal23;
     G4LogicalVolume*	    fLogicCrystal24;
     G4LogicalVolume*	    fLogicCrystal25;
     G4LogicalVolume*	    fLogicCrystal26;
     G4LogicalVolume*	    fLogicCrystal27;
     G4LogicalVolume*	    fLogicCrystal28;
     G4LogicalVolume*	    fLogicCrystal29;
     G4LogicalVolume*	    fLogicCrystal30;
     G4VPhysicalVolume*     fPhysiCrystal;
     G4VPhysicalVolume*     fPhysiCrystal1;
     G4VPhysicalVolume*     fPhysiCrystal2;
     G4VPhysicalVolume*     fPhysiCrystal3;
     G4VPhysicalVolume*     fPhysiCrystal4;
     G4VPhysicalVolume*     fPhysiCrystal5;
     G4VPhysicalVolume*     fPhysiCrystal6;
     G4VPhysicalVolume*     fPhysiCrystal7;
     G4VPhysicalVolume*     fPhysiCrystal8;
     G4VPhysicalVolume*     fPhysiCrystal9;
     G4VPhysicalVolume*     fPhysiCrystal10;
     G4VPhysicalVolume*     fPhysiCrystal11;
     G4VPhysicalVolume*     fPhysiCrystal12;
     G4VPhysicalVolume*     fPhysiCrystal13;
     G4VPhysicalVolume*     fPhysiCrystal14;
     G4VPhysicalVolume*     fPhysiCrystal15;
     G4VPhysicalVolume*     fPhysiCrystal16;
     G4VPhysicalVolume*     fPhysiCrystal17;
     G4VPhysicalVolume*     fPhysiCrystal18;
     G4VPhysicalVolume*     fPhysiCrystal19;
     G4VPhysicalVolume*     fPhysiCrystal20;
     G4VPhysicalVolume*     fPhysiCrystal21;
     G4VPhysicalVolume*     fPhysiCrystal22;
     G4VPhysicalVolume*     fPhysiCrystal23;
     G4VPhysicalVolume*     fPhysiCrystal24;
     G4VPhysicalVolume*     fPhysiCrystal25;
     G4VPhysicalVolume*     fPhysiCrystal26;
     G4VPhysicalVolume*     fPhysiCrystal27;
     G4VPhysicalVolume*     fPhysiCrystal28;
     G4VPhysicalVolume*     fPhysiCrystal29;
     G4VPhysicalVolume*     fPhysiCrystal30;
     
     G4double 		        fDiameter;
     G4double 		        fLength;
     G4double		        fZpos;
     G4double		        fZpos1;
     G4Material*            fCrystalMaterial;
     G4Material*            fCrystalMaterial1;
     
     //Gap between scintillator and aluminum casing (reflector?)
     //around crystal
     G4Tubs* 		        fSolidGap;
     G4Tubs* 		        fSolidGap1;
     G4Polyhedra*           fSolidHGap;
     G4Polyhedra*           fSolidHGap1;
     G4Polyhedra*           fSolidHGap2;
     G4Polyhedra*           fSolidHGap3;
     G4Polyhedra*           fSolidHGap4;
     G4Polyhedra*           fSolidHGap5;
     G4Polyhedra*           fSolidHGap6;
     G4Polyhedra*           fSolidHGap7;
     G4Polyhedra*           fSolidHGap8;
     G4Polyhedra*           fSolidHGap9;
     G4Polyhedra*           fSolidHGap10;
     G4Polyhedra*           fSolidHGap11;
     G4Polyhedra*           fSolidHGap12;
     G4Polyhedra*           fSolidHGap13;
     G4Polyhedra*           fSolidHGap14;
     G4Polyhedra*           fSolidHGap15;
     G4Polyhedra*           fSolidHGap16;
     G4Polyhedra*           fSolidHGap17;
     G4Polyhedra*           fSolidHGap18;
     G4Polyhedra*           fSolidHGap19;
     G4Polyhedra*           fSolidHGap20;
     G4Polyhedra*           fSolidHGap21;
     G4Polyhedra*           fSolidHGap22;
     G4Polyhedra*           fSolidHGap23;
     G4Polyhedra*           fSolidHGap24;
     G4Polyhedra*           fSolidHGap25;
     G4Polyhedra*           fSolidHGap26;
     G4Polyhedra*           fSolidHGap27;
     G4Polyhedra*           fSolidHGap28;
     G4Polyhedra*           fSolidHGap29;
     G4Polyhedra*           fSolidHGap30;
     G4LogicalVolume*	    fLogicGap;
     G4LogicalVolume*	    fLogicGap1;
     G4LogicalVolume*       fLogicGap2;
     G4LogicalVolume*       fLogicGap3;
     G4LogicalVolume*	    fLogicGap4;
     G4LogicalVolume*       fLogicGap5;
     G4LogicalVolume*       fLogicGap6;
     G4LogicalVolume*	    fLogicGap7;
     G4LogicalVolume*       fLogicGap8;
     G4LogicalVolume*       fLogicGap9;
     G4LogicalVolume*	    fLogicGap10;
     G4LogicalVolume*       fLogicGap11;
     G4LogicalVolume*       fLogicGap12;
     G4LogicalVolume*	    fLogicGap13;
     G4LogicalVolume*       fLogicGap14;
     G4LogicalVolume*       fLogicGap15;
     G4LogicalVolume*	    fLogicGap16;
     G4LogicalVolume*       fLogicGap17;
     G4LogicalVolume*       fLogicGap18;
     G4LogicalVolume*	    fLogicGap19;
     G4LogicalVolume*       fLogicGap20;
     G4LogicalVolume*       fLogicGap21;
     G4LogicalVolume*	    fLogicGap22;
     G4LogicalVolume*       fLogicGap23;
     G4LogicalVolume*       fLogicGap24;
     G4LogicalVolume*       fLogicGap25;
     G4LogicalVolume*	    fLogicGap26;
     G4LogicalVolume*       fLogicGap27;
     G4LogicalVolume*       fLogicGap28;
     G4LogicalVolume*	    fLogicGap29;
     G4LogicalVolume*       fLogicGap30;
     G4VPhysicalVolume*     fPhysiGap;
     G4VPhysicalVolume*     fPhysiGap1;
     G4VPhysicalVolume*     fPhysiGap2;
     G4VPhysicalVolume*     fPhysiGap3;
     G4VPhysicalVolume*     fPhysiGap4;
     G4VPhysicalVolume*     fPhysiGap5;
     G4VPhysicalVolume*     fPhysiGap6;
     G4VPhysicalVolume*     fPhysiGap7;
     G4VPhysicalVolume*     fPhysiGap8;
     G4VPhysicalVolume*     fPhysiGap9;
     G4VPhysicalVolume*     fPhysiGap10;
     G4VPhysicalVolume*     fPhysiGap11;
     G4VPhysicalVolume*     fPhysiGap12;
     G4VPhysicalVolume*     fPhysiGap13;
     G4VPhysicalVolume*     fPhysiGap14;
     G4VPhysicalVolume*     fPhysiGap15;
     G4VPhysicalVolume*     fPhysiGap16;
     G4VPhysicalVolume*     fPhysiGap17;
     G4VPhysicalVolume*     fPhysiGap18;
     G4VPhysicalVolume*     fPhysiGap19;
     G4VPhysicalVolume*     fPhysiGap20;
     G4VPhysicalVolume*     fPhysiGap21;
     G4VPhysicalVolume*     fPhysiGap22;
     G4VPhysicalVolume*     fPhysiGap23;
     G4VPhysicalVolume*     fPhysiGap24;
     G4VPhysicalVolume*     fPhysiGap25;
     G4VPhysicalVolume*     fPhysiGap26;
     G4VPhysicalVolume*     fPhysiGap27;
     G4VPhysicalVolume*     fPhysiGap28;
     G4VPhysicalVolume*     fPhysiGap29;
     G4VPhysicalVolume*     fPhysiGap30;
     G4double 		        fGapLength;
     G4double		        fZposGap;
     G4double 		        fGapLength1;
     G4double		        fZposGap1;
     G4Material* 	        fGapMaterial;
     //over face of crystal
     G4Tubs* 		        fSolidFaceGap;
     G4Tubs* 		        fSolidFaceGap1;
     G4Polyhedra*           fSolidHFaceGap;
     G4Polyhedra*           fSolidHFaceGap1;
     G4Polyhedra*           fSolidHFaceGap2;
     G4Polyhedra*           fSolidHFaceGap3;
     G4Polyhedra*           fSolidHFaceGap4;
     G4Polyhedra*           fSolidHFaceGap5;
     G4Polyhedra*           fSolidHFaceGap6;
     G4Polyhedra*           fSolidHFaceGap7;
     G4Polyhedra*           fSolidHFaceGap8;
     G4Polyhedra*           fSolidHFaceGap9;
     G4Polyhedra*           fSolidHFaceGap10;
     G4Polyhedra*           fSolidHFaceGap11;
     G4Polyhedra*           fSolidHFaceGap12;
     G4Polyhedra*           fSolidHFaceGap13;
     G4Polyhedra*           fSolidHFaceGap14;
     G4Polyhedra*           fSolidHFaceGap15;
     G4Polyhedra*           fSolidHFaceGap16;
     G4Polyhedra*           fSolidHFaceGap17;
     G4Polyhedra*           fSolidHFaceGap18;
     G4Polyhedra*           fSolidHFaceGap19;
     G4Polyhedra*           fSolidHFaceGap20;
     G4Polyhedra*           fSolidHFaceGap21;
     G4Polyhedra*           fSolidHFaceGap22;
     G4Polyhedra*           fSolidHFaceGap23;
     G4Polyhedra*           fSolidHFaceGap24;
     G4Polyhedra*           fSolidHFaceGap25;
     G4Polyhedra*           fSolidHFaceGap26;
     G4Polyhedra*           fSolidHFaceGap27;
     G4Polyhedra*           fSolidHFaceGap28;
     G4Polyhedra*           fSolidHFaceGap29;
     G4Polyhedra*           fSolidHFaceGap30;
     G4LogicalVolume*	    fLogicFaceGap;
     G4LogicalVolume*	    fLogicFaceGap1;
     G4LogicalVolume*       fLogicFaceGap2;
     G4LogicalVolume*       fLogicFaceGap3;
     G4LogicalVolume*	    fLogicFaceGap4;
     G4LogicalVolume*       fLogicFaceGap5;
     G4LogicalVolume*       fLogicFaceGap6;
     G4LogicalVolume*	    fLogicFaceGap7;
     G4LogicalVolume*       fLogicFaceGap8;
     G4LogicalVolume*       fLogicFaceGap9;
     G4LogicalVolume*	    fLogicFaceGap10;
     G4LogicalVolume*       fLogicFaceGap11;
     G4LogicalVolume*       fLogicFaceGap12;
     G4LogicalVolume*       fLogicFaceGap13;
     G4LogicalVolume*       fLogicFaceGap14;
     G4LogicalVolume*	    fLogicFaceGap15;
     G4LogicalVolume*       fLogicFaceGap16;
     G4LogicalVolume*       fLogicFaceGap17;
     G4LogicalVolume*	    fLogicFaceGap18;
     G4LogicalVolume*       fLogicFaceGap19;
     G4LogicalVolume*       fLogicFaceGap20;
     G4LogicalVolume*       fLogicFaceGap21;
     G4LogicalVolume*       fLogicFaceGap22;
     G4LogicalVolume*	    fLogicFaceGap23;
     G4LogicalVolume*       fLogicFaceGap24;
     G4LogicalVolume*       fLogicFaceGap25;
     G4LogicalVolume*	    fLogicFaceGap26;
     G4LogicalVolume*       fLogicFaceGap27;
     G4LogicalVolume*       fLogicFaceGap28;
     G4LogicalVolume*       fLogicFaceGap29;
     G4LogicalVolume*       fLogicFaceGap30;
     G4VPhysicalVolume*     fPhysiFaceGap;
     G4VPhysicalVolume*     fPhysiFaceGap1;
     G4VPhysicalVolume*     fPhysiFaceGap2;
     G4VPhysicalVolume*     fPhysiFaceGap3;
     G4VPhysicalVolume*     fPhysiFaceGap4;
     G4VPhysicalVolume*     fPhysiFaceGap5;
     G4VPhysicalVolume*     fPhysiFaceGap6;
     G4VPhysicalVolume*     fPhysiFaceGap7;
     G4VPhysicalVolume*     fPhysiFaceGap8;
     G4VPhysicalVolume*     fPhysiFaceGap9;
     G4VPhysicalVolume*     fPhysiFaceGap10;
     G4VPhysicalVolume*     fPhysiFaceGap11;
     G4VPhysicalVolume*     fPhysiFaceGap12;
     G4VPhysicalVolume*     fPhysiFaceGap13;
     G4VPhysicalVolume*     fPhysiFaceGap14;
     G4VPhysicalVolume*     fPhysiFaceGap15;
     G4VPhysicalVolume*     fPhysiFaceGap16;
     G4VPhysicalVolume*     fPhysiFaceGap17;
     G4VPhysicalVolume*     fPhysiFaceGap18;
     G4VPhysicalVolume*     fPhysiFaceGap19;
     G4VPhysicalVolume*     fPhysiFaceGap20;
     G4VPhysicalVolume*     fPhysiFaceGap21;
     G4VPhysicalVolume*     fPhysiFaceGap22;
     G4VPhysicalVolume*     fPhysiFaceGap23;
     G4VPhysicalVolume*     fPhysiFaceGap24;
     G4VPhysicalVolume*     fPhysiFaceGap25;
     G4VPhysicalVolume*     fPhysiFaceGap26;
     G4VPhysicalVolume*     fPhysiFaceGap27;
     G4VPhysicalVolume*     fPhysiFaceGap28;
     G4VPhysicalVolume*     fPhysiFaceGap29;
     G4VPhysicalVolume*     fPhysiFaceGap30;
     G4double		        fZposFaceGap;
     G4double		        fZposFaceGap1;
     G4Material* 	        fFaceGapMaterial;
     //common to both components of gap
     G4double 		        fGapDiameter;
     G4double 		        fGapThickness;
     G4double 				fGapFaceThickness;
     G4double 				fGapFaceThickness1;
     G4double 		        fGapDiameter1;
     G4double 		        fGapThickness1;
     
     //Aluminum casing surrounding scintillator crystal
     //around crystal
     G4Tubs* 		        fSolidAlCase;
     G4Tubs* 		        fSolidAlCase1;
     G4Polyhedra*           fSolidHAlCase;
     G4Polyhedra*           fSolidHAlCase1;
     G4Polyhedra*           fSolidHAlCase2;
     G4Polyhedra*           fSolidHAlCase3;
     G4Polyhedra*           fSolidHAlCase4;
     G4Polyhedra*           fSolidHAlCase5;
     G4Polyhedra*           fSolidHAlCase6;
     G4Polyhedra*           fSolidHAlCase7;
     G4Polyhedra*           fSolidHAlCase8;
     G4Polyhedra*           fSolidHAlCase9;
     G4Polyhedra*           fSolidHAlCase10;
     G4Polyhedra*           fSolidHAlCase11;
     G4Polyhedra*           fSolidHAlCase12;
     G4Polyhedra*           fSolidHAlCase13;
     G4Polyhedra*           fSolidHAlCase14;
     G4Polyhedra*           fSolidHAlCase15;
     G4Polyhedra*           fSolidHAlCase16;
     G4Polyhedra*           fSolidHAlCase17;
     G4Polyhedra*           fSolidHAlCase18;
     G4Polyhedra*           fSolidHAlCase19;
     G4Polyhedra*           fSolidHAlCase20;
     G4Polyhedra*           fSolidHAlCase21;
     G4Polyhedra*           fSolidHAlCase22;
     G4Polyhedra*           fSolidHAlCase23;
     G4Polyhedra*           fSolidHAlCase24;
     G4Polyhedra*           fSolidHAlCase25;
     G4Polyhedra*           fSolidHAlCase26;
     G4Polyhedra*           fSolidHAlCase27;
     G4Polyhedra*           fSolidHAlCase28;
     G4Polyhedra*           fSolidHAlCase29;
     G4Polyhedra*           fSolidHAlCase30;
     G4LogicalVolume*	    fLogicAlCase;
     G4LogicalVolume*	    fLogicAlCase1;
     G4LogicalVolume*       fLogicAlCase2;
     G4LogicalVolume*       fLogicAlCase3;
     G4LogicalVolume*	    fLogicAlCase4;
     G4LogicalVolume*       fLogicAlCase5;
     G4LogicalVolume*       fLogicAlCase6;
     G4LogicalVolume*	    fLogicAlCase7;
     G4LogicalVolume*       fLogicAlCase8;
     G4LogicalVolume*       fLogicAlCase9;
     G4LogicalVolume*	    fLogicAlCase10;
     G4LogicalVolume*       fLogicAlCase11;
     G4LogicalVolume*       fLogicAlCase12;
     G4LogicalVolume*	    fLogicAlCase13;
     G4LogicalVolume*       fLogicAlCase14;
     G4LogicalVolume*       fLogicAlCase15;
     G4LogicalVolume*	    fLogicAlCase16;
     G4LogicalVolume*       fLogicAlCase17;
     G4LogicalVolume*       fLogicAlCase18;
     G4LogicalVolume*	    fLogicAlCase19;
     G4LogicalVolume*       fLogicAlCase20;
     G4LogicalVolume*       fLogicAlCase21;
     G4LogicalVolume*	    fLogicAlCase22;
     G4LogicalVolume*       fLogicAlCase23;
     G4LogicalVolume*       fLogicAlCase24;
     G4LogicalVolume*	    fLogicAlCase25;
     G4LogicalVolume*       fLogicAlCase26;
     G4LogicalVolume*       fLogicAlCase27;
     G4LogicalVolume*	    fLogicAlCase28;
     G4LogicalVolume*       fLogicAlCase29;
     G4LogicalVolume*       fLogicAlCase30;
     G4VPhysicalVolume*     fPhysiAlCase;
     G4VPhysicalVolume*     fPhysiAlCase1;
     G4VPhysicalVolume*     fPhysiAlCase2;
     G4VPhysicalVolume*     fPhysiAlCase3;
     G4VPhysicalVolume*     fPhysiAlCase4;
     G4VPhysicalVolume*     fPhysiAlCase5;
     G4VPhysicalVolume*     fPhysiAlCase6;
     G4VPhysicalVolume*     fPhysiAlCase7;
     G4VPhysicalVolume*     fPhysiAlCase8;
     G4VPhysicalVolume*     fPhysiAlCase9;
     G4VPhysicalVolume*     fPhysiAlCase10;
     G4VPhysicalVolume*     fPhysiAlCase11;
     G4VPhysicalVolume*     fPhysiAlCase12;
     G4VPhysicalVolume*     fPhysiAlCase13;
     G4VPhysicalVolume*     fPhysiAlCase14;
     G4VPhysicalVolume*     fPhysiAlCase15;
     G4VPhysicalVolume*     fPhysiAlCase16;
     G4VPhysicalVolume*     fPhysiAlCase17;
     G4VPhysicalVolume*     fPhysiAlCase18;
     G4VPhysicalVolume*     fPhysiAlCase19;
     G4VPhysicalVolume*     fPhysiAlCase20;
     G4VPhysicalVolume*     fPhysiAlCase21;
     G4VPhysicalVolume*     fPhysiAlCase22;
     G4VPhysicalVolume*     fPhysiAlCase23;
     G4VPhysicalVolume*     fPhysiAlCase24;
     G4VPhysicalVolume*     fPhysiAlCase25;
     G4VPhysicalVolume*     fPhysiAlCase26;
     G4VPhysicalVolume*     fPhysiAlCase27;
     G4VPhysicalVolume*     fPhysiAlCase28;
     G4VPhysicalVolume*     fPhysiAlCase29;
     G4VPhysicalVolume*     fPhysiAlCase30;
     G4double		        fAlCaseLength;
     G4double		        fZposAlCase;
     G4double		        fAlCaseLength1;
     G4double		        fZposAlCase1;
     G4Material* 	        fAlCaseMaterial;
     //over face of crystal
     G4Tubs* 		        fSolidFaceAlCase;
     G4Tubs* 		        fSolidFaceAlCase1;
     G4Polyhedra*           fSolidHFaceAlCase;
     G4Polyhedra*           fSolidHFaceAlCase1;
     G4Polyhedra*           fSolidHFaceAlCase2;
     G4Polyhedra*           fSolidHFaceAlCase3;
     G4Polyhedra*           fSolidHFaceAlCase4;
     G4Polyhedra*           fSolidHFaceAlCase5;
     G4Polyhedra*           fSolidHFaceAlCase6;
     G4Polyhedra*           fSolidHFaceAlCase7;
     G4Polyhedra*           fSolidHFaceAlCase8;
     G4Polyhedra*           fSolidHFaceAlCase9;
     G4Polyhedra*           fSolidHFaceAlCase10;
     G4Polyhedra*           fSolidHFaceAlCase11;
     G4Polyhedra*           fSolidHFaceAlCase12;
     G4Polyhedra*           fSolidHFaceAlCase13;
     G4Polyhedra*           fSolidHFaceAlCase14;
     G4Polyhedra*           fSolidHFaceAlCase15;
     G4Polyhedra*           fSolidHFaceAlCase16;
     G4Polyhedra*           fSolidHFaceAlCase17;
     G4Polyhedra*           fSolidHFaceAlCase18;
     G4Polyhedra*           fSolidHFaceAlCase19;
     G4Polyhedra*           fSolidHFaceAlCase20;
     G4Polyhedra*           fSolidHFaceAlCase21;
     G4Polyhedra*           fSolidHFaceAlCase22;
     G4Polyhedra*           fSolidHFaceAlCase23;
     G4Polyhedra*           fSolidHFaceAlCase24;
     G4Polyhedra*           fSolidHFaceAlCase25;
     G4Polyhedra*           fSolidHFaceAlCase26;
     G4Polyhedra*           fSolidHFaceAlCase27;
     G4Polyhedra*           fSolidHFaceAlCase28;
     G4Polyhedra*           fSolidHFaceAlCase29;
     G4Polyhedra*           fSolidHFaceAlCase30;
     G4LogicalVolume*	    fLogicFaceAlCase;
     G4LogicalVolume*	    fLogicFaceAlCase1;
     G4LogicalVolume*       fLogicFaceAlCase2;
     G4LogicalVolume*	    fLogicFaceAlCase3;
     G4LogicalVolume*       fLogicFaceAlCase4;
     G4LogicalVolume*	    fLogicFaceAlCase5;
     G4LogicalVolume*       fLogicFaceAlCase6;
     G4LogicalVolume*	    fLogicFaceAlCase7;
     G4LogicalVolume*       fLogicFaceAlCase8;
     G4LogicalVolume*	    fLogicFaceAlCase9;
     G4LogicalVolume*       fLogicFaceAlCase10;
     G4LogicalVolume*	    fLogicFaceAlCase11;
     G4LogicalVolume*       fLogicFaceAlCase12;
     G4LogicalVolume*	    fLogicFaceAlCase13;
     G4LogicalVolume*       fLogicFaceAlCase14;
     G4LogicalVolume*	    fLogicFaceAlCase15;
     G4LogicalVolume*       fLogicFaceAlCase16;
     G4LogicalVolume*	    fLogicFaceAlCase17;
     G4LogicalVolume*       fLogicFaceAlCase18;
     G4LogicalVolume*	    fLogicFaceAlCase19;
     G4LogicalVolume*       fLogicFaceAlCase20;
     G4LogicalVolume*	    fLogicFaceAlCase21;
     G4LogicalVolume*       fLogicFaceAlCase22;
     G4LogicalVolume*	    fLogicFaceAlCase23;
     G4LogicalVolume*       fLogicFaceAlCase24;
     G4LogicalVolume*	    fLogicFaceAlCase25;
     G4LogicalVolume*       fLogicFaceAlCase26;
     G4LogicalVolume*	    fLogicFaceAlCase27;
     G4LogicalVolume*       fLogicFaceAlCase28;
     G4LogicalVolume*	    fLogicFaceAlCase29;
     G4LogicalVolume*       fLogicFaceAlCase30;
     G4VPhysicalVolume*     fPhysiFaceAlCase;
     G4VPhysicalVolume*     fPhysiFaceAlCase1;
     G4VPhysicalVolume*     fPhysiFaceAlCase2;
     G4VPhysicalVolume*     fPhysiFaceAlCase3;
     G4VPhysicalVolume*     fPhysiFaceAlCase4;
     G4VPhysicalVolume*     fPhysiFaceAlCase5;
     G4VPhysicalVolume*     fPhysiFaceAlCase6;
     G4VPhysicalVolume*     fPhysiFaceAlCase7;
     G4VPhysicalVolume*     fPhysiFaceAlCase8;
     G4VPhysicalVolume*     fPhysiFaceAlCase9;
     G4VPhysicalVolume*     fPhysiFaceAlCase10;
     G4VPhysicalVolume*     fPhysiFaceAlCase11;
     G4VPhysicalVolume*     fPhysiFaceAlCase12;
     G4VPhysicalVolume*     fPhysiFaceAlCase13;
     G4VPhysicalVolume*     fPhysiFaceAlCase14;
     G4VPhysicalVolume*     fPhysiFaceAlCase15;
     G4VPhysicalVolume*     fPhysiFaceAlCase16;
     G4VPhysicalVolume*     fPhysiFaceAlCase17;
     G4VPhysicalVolume*     fPhysiFaceAlCase18;
     G4VPhysicalVolume*     fPhysiFaceAlCase19;
     G4VPhysicalVolume*     fPhysiFaceAlCase20;
     G4VPhysicalVolume*     fPhysiFaceAlCase21;
     G4VPhysicalVolume*     fPhysiFaceAlCase22;
     G4VPhysicalVolume*     fPhysiFaceAlCase23;
     G4VPhysicalVolume*     fPhysiFaceAlCase24;
     G4VPhysicalVolume*     fPhysiFaceAlCase25;
     G4VPhysicalVolume*     fPhysiFaceAlCase26;
     G4VPhysicalVolume*     fPhysiFaceAlCase27;
     G4VPhysicalVolume*     fPhysiFaceAlCase28;
     G4VPhysicalVolume*     fPhysiFaceAlCase29;
     G4VPhysicalVolume*     fPhysiFaceAlCase30;
     G4double		        fZposFaceAlCase;
     G4double		        fZposFaceAlCase1;
     G4Material* 	        fFaceAlCaseMaterial;
     //common to both components of aluminum case
     G4double 		        fAlCaseDiameter;
     G4double 		        fAlCaseThickness;
     G4double 				fAlFaceThickness;
     G4double 		        fAlCaseDiameter1;
     G4double 		        fAlCaseThickness1;
     G4double 				fAlFaceThickness1;

     //Lead casing surrounding aluminum casing
     //around casing surround
     G4Tubs* 		        fSolidPbCase;
     G4Polyhedra*           fSolidHPbCase;
     G4LogicalVolume*	    fLogicPbCase;
     G4VPhysicalVolume*     fPhysiPbCase;
     G4Material* 	        fPbCaseMaterial;
      //at back in front of PMT
     G4Tubs* 		        fSolidPbCollar;
     G4Polyhedra*           fSolidHPbCollar;
     G4LogicalVolume*	    fLogicPbCollar;
     G4VPhysicalVolume*     fPhysiPbCollar;
     G4double		        fZposPbCollar;
     G4Material* 	        fPbCollarMaterial;
     //common to both components of lead surround
     G4double 		        fPbCaseDiameter;
     G4double 		        fPbCaseThickness;

     // PMT
     // PMT Body
     G4Tubs* 		        fSolidPMT;
     G4Tubs* 		        fSolidPMT1;
     G4Tubs*                fSolidPMT2;
     G4Tubs*                fSolidPMT3;
     G4Tubs*                fSolidPMT4;
     G4Tubs*                fSolidPMT5;
     G4Tubs*                fSolidPMT6;
     G4Tubs*                fSolidPMT7;
     G4Tubs*                fSolidPMT8;
     G4Tubs*                fSolidPMT9;
     G4Tubs*                fSolidPMT10;
     G4Tubs*                fSolidPMT11;
     G4Tubs*                fSolidPMT12;
     G4Tubs*                fSolidPMT13;
     G4Tubs*                fSolidPMT14;
     G4Tubs*                fSolidPMT15;
     G4Tubs*                fSolidPMT16;
     G4Tubs*                fSolidPMT17;
     G4Tubs*                fSolidPMT18;
     G4Tubs*                fSolidPMT19;
     G4Tubs*                fSolidPMT20;
     G4Tubs*                fSolidPMT21;
     G4Tubs*                fSolidPMT22;
     G4Tubs*                fSolidPMT23;
     G4Tubs*                fSolidPMT24;
     G4Tubs*                fSolidPMT25;
     G4Tubs*                fSolidPMT26;
     G4Tubs*                fSolidPMT27;
     G4Tubs*                fSolidPMT28;
     G4Tubs*                fSolidPMT29;
     G4Tubs*                fSolidPMT30;
     G4Tubs*                fSolidPMTI;
     G4Tubs*				fSolidPMTK;
     G4Polyhedra*           fSolidHPMT;
     G4Polyhedra*           fSolidHPMT1;
     
     G4LogicalVolume*	    fLogicPMT;
     G4LogicalVolume*	    fLogicPMT1;
     G4LogicalVolume*	    fLogicPMT2;
     G4LogicalVolume*	    fLogicPMT3;
     G4LogicalVolume*	    fLogicPMT4;
     G4LogicalVolume*	    fLogicPMT5;
     G4LogicalVolume*	    fLogicPMT6;
     G4LogicalVolume*	    fLogicPMT7;
     G4LogicalVolume*	    fLogicPMT8;
     G4LogicalVolume*	    fLogicPMT9;
     G4LogicalVolume*	    fLogicPMT10;
     G4LogicalVolume*	    fLogicPMT11;
     G4LogicalVolume*	    fLogicPMT12;
     G4LogicalVolume*	    fLogicPMT13;
     G4LogicalVolume*	    fLogicPMT14;
     G4LogicalVolume*	    fLogicPMT15;
     G4LogicalVolume*	    fLogicPMT16;
     G4LogicalVolume*	    fLogicPMT17;
     G4LogicalVolume*	    fLogicPMT18;
     G4LogicalVolume*	    fLogicPMT19;
     G4LogicalVolume*	    fLogicPMT20;
     G4LogicalVolume*	    fLogicPMT21;
     G4LogicalVolume*	    fLogicPMT22;
     G4LogicalVolume*	    fLogicPMT23;
     G4LogicalVolume*	    fLogicPMT24;
     G4LogicalVolume*	    fLogicPMT25;
     G4LogicalVolume*	    fLogicPMT26;
     G4LogicalVolume*	    fLogicPMT27;
     G4LogicalVolume*	    fLogicPMT28;
     G4LogicalVolume*	    fLogicPMT29;
     G4LogicalVolume*	    fLogicPMT30;
     G4VPhysicalVolume*     fPhysiPMT;
     G4VPhysicalVolume*     fPhysiPMT1;
     G4VPhysicalVolume*     fPhysiPMT2;
     G4VPhysicalVolume*     fPhysiPMT3;
     G4VPhysicalVolume*     fPhysiPMT4;
     G4VPhysicalVolume*     fPhysiPMT5;
     G4VPhysicalVolume*     fPhysiPMT6;
     G4VPhysicalVolume*     fPhysiPMT7;
     G4VPhysicalVolume*     fPhysiPMT8;
     G4VPhysicalVolume*     fPhysiPMT9;
     G4VPhysicalVolume*     fPhysiPMT10;
     G4VPhysicalVolume*     fPhysiPMT11;
     G4VPhysicalVolume*     fPhysiPMT12;
     G4VPhysicalVolume*     fPhysiPMT13;
     G4VPhysicalVolume*     fPhysiPMT14;
     G4VPhysicalVolume*     fPhysiPMT15;
     G4VPhysicalVolume*     fPhysiPMT16;
     G4VPhysicalVolume*     fPhysiPMT17;
     G4VPhysicalVolume*     fPhysiPMT18;
     G4VPhysicalVolume*     fPhysiPMT19;
     G4VPhysicalVolume*     fPhysiPMT20;
     G4VPhysicalVolume*     fPhysiPMT21;
     G4VPhysicalVolume*     fPhysiPMT22;
     G4VPhysicalVolume*     fPhysiPMT23;
     G4VPhysicalVolume*     fPhysiPMT24;
     G4VPhysicalVolume*     fPhysiPMT25;
     G4VPhysicalVolume*     fPhysiPMT26;
     G4VPhysicalVolume*     fPhysiPMT27;
     G4VPhysicalVolume*     fPhysiPMT28;
     G4VPhysicalVolume*     fPhysiPMT29;
     G4VPhysicalVolume*     fPhysiPMT30;
     G4LogicalVolume*	    fLogicPMTI;
     G4VPhysicalVolume*     fPhysiPMTI;
     G4Tubs* 				fSolidPMTJ;
     G4Tubs*                fSolidPMTJ2;
     G4Tubs* 				fSolidPMTJ3;
     G4Tubs*                fSolidPMTJ4;
     G4Tubs* 				fSolidPMTJ5;
     G4Tubs*                fSolidPMTJ6;
     G4Tubs* 				fSolidPMTJ7;
     G4Tubs*                fSolidPMTJ8;
     G4Tubs* 				fSolidPMTJ9;
     G4Tubs*                fSolidPMTJ10;
     G4Tubs* 				fSolidPMTJ11;
     G4Tubs*                fSolidPMTJ12;
     G4Tubs* 				fSolidPMTJ13;
     G4Tubs*                fSolidPMTJ14;
     G4Tubs* 				fSolidPMTJ15;
     G4Tubs*                fSolidPMTJ16;
     G4Tubs*                fSolidPMTJ17;
     G4Tubs* 				fSolidPMTJ18;
     G4Tubs*                fSolidPMTJ19;
     G4Tubs* 				fSolidPMTJ20;
     G4Tubs*                fSolidPMTJ21;
     G4Tubs* 				fSolidPMTJ22;
     G4Tubs*                fSolidPMTJ23;
     G4Tubs* 				fSolidPMTJ24;
     G4Tubs*                fSolidPMTJ25;
     G4Tubs* 				fSolidPMTJ26;
     G4Tubs*                fSolidPMTJ27;
     G4Tubs*                fSolidPMTJ28;
     G4Tubs* 				fSolidPMTJ29;
     G4Tubs*                fSolidPMTJ30;
     G4LogicalVolume*	    fLogicPMTJ;
     G4LogicalVolume*	    fLogicPMTJ2;
     G4LogicalVolume*	    fLogicPMTJ3;
     G4LogicalVolume*	    fLogicPMTJ4;
     G4LogicalVolume*	    fLogicPMTJ5;
     G4LogicalVolume*	    fLogicPMTJ6;
     G4LogicalVolume*	    fLogicPMTJ7;
     G4LogicalVolume*	    fLogicPMTJ8;
     G4LogicalVolume*	    fLogicPMTJ9;
     G4LogicalVolume*	    fLogicPMTJ10;
     G4LogicalVolume*	    fLogicPMTJ11;
     G4LogicalVolume*	    fLogicPMTJ12;
     G4LogicalVolume*	    fLogicPMTJ13;
     G4LogicalVolume*	    fLogicPMTJ14;
     G4LogicalVolume*	    fLogicPMTJ15;
     G4LogicalVolume*	    fLogicPMTJ16;
     G4LogicalVolume*	    fLogicPMTJ17;
     G4LogicalVolume*	    fLogicPMTJ18;
     G4LogicalVolume*	    fLogicPMTJ19;
     G4LogicalVolume*	    fLogicPMTJ20;
     G4LogicalVolume*	    fLogicPMTJ21;
     G4LogicalVolume*	    fLogicPMTJ22;
     G4LogicalVolume*	    fLogicPMTJ23;
     G4LogicalVolume*	    fLogicPMTJ24;
     G4LogicalVolume*	    fLogicPMTJ25;
     G4LogicalVolume*	    fLogicPMTJ26;
     G4LogicalVolume*	    fLogicPMTJ27;
     G4LogicalVolume*	    fLogicPMTJ28;
     G4LogicalVolume*	    fLogicPMTJ29;
     G4LogicalVolume*	    fLogicPMTJ30;
     G4VPhysicalVolume*     fPhysiPMTJ;
     G4VPhysicalVolume*     fPhysiPMTJ2;
     G4VPhysicalVolume*     fPhysiPMTJ3;
     G4VPhysicalVolume*     fPhysiPMTJ4;
     G4VPhysicalVolume*     fPhysiPMTJ5;
     G4VPhysicalVolume*     fPhysiPMTJ6;
     G4VPhysicalVolume*     fPhysiPMTJ7;
     G4VPhysicalVolume*     fPhysiPMTJ8;
     G4VPhysicalVolume*     fPhysiPMTJ9;
     G4VPhysicalVolume*     fPhysiPMTJ10;
     G4VPhysicalVolume*     fPhysiPMTJ11;
     G4VPhysicalVolume*     fPhysiPMTJ12;
     G4VPhysicalVolume*     fPhysiPMTJ13;
     G4VPhysicalVolume*     fPhysiPMTJ14;
     G4VPhysicalVolume*     fPhysiPMTJ15;
     G4VPhysicalVolume*     fPhysiPMTJ16;
     G4VPhysicalVolume*     fPhysiPMTJ17;
     G4VPhysicalVolume*     fPhysiPMTJ18;
     G4VPhysicalVolume*     fPhysiPMTJ19;
     G4VPhysicalVolume*     fPhysiPMTJ20;
     G4VPhysicalVolume*     fPhysiPMTJ21;
     G4VPhysicalVolume*     fPhysiPMTJ22;
     G4VPhysicalVolume*     fPhysiPMTJ23;
     G4VPhysicalVolume*     fPhysiPMTJ24;
     G4VPhysicalVolume*     fPhysiPMTJ25;
     G4VPhysicalVolume*     fPhysiPMTJ26;
     G4VPhysicalVolume*     fPhysiPMTJ27;
     G4VPhysicalVolume*     fPhysiPMTJ28;
     G4VPhysicalVolume*     fPhysiPMTJ29;
     G4VPhysicalVolume*     fPhysiPMTJ30;
     G4LogicalVolume*	    fLogicPMTK;
     G4VPhysicalVolume*     fPhysiPMTK;
     G4Tubs*				fSolidPMTL;
     G4Tubs*                fSolidPMTL2;
     G4Tubs*				fSolidPMTL3;
     G4Tubs*                fSolidPMTL4;
     G4Tubs*				fSolidPMTL5;
     G4Tubs*                fSolidPMTL6;
     G4Tubs*				fSolidPMTL7;
     G4Tubs*                fSolidPMTL8;
     G4Tubs*				fSolidPMTL9;
     G4Tubs*                fSolidPMTL10;
     G4Tubs*				fSolidPMTL11;
     G4Tubs*                fSolidPMTL12;
     G4Tubs*				fSolidPMTL13;
     G4Tubs*                fSolidPMTL14;
     G4Tubs*				fSolidPMTL15;
     G4Tubs*                fSolidPMTL16;
     G4Tubs*				fSolidPMTL17;
     G4Tubs*                fSolidPMTL18;
     G4Tubs*				fSolidPMTL19;
     G4Tubs*                fSolidPMTL20;
     G4Tubs*				fSolidPMTL21;
     G4Tubs*                fSolidPMTL22;
     G4Tubs*				fSolidPMTL23;
     G4Tubs*                fSolidPMTL24;
     G4Tubs*				fSolidPMTL25;
     G4Tubs*                fSolidPMTL26;
     G4Tubs*				fSolidPMTL27;
     G4Tubs*                fSolidPMTL28;
     G4Tubs*				fSolidPMTL29;
     G4Tubs*                fSolidPMTL30;
     G4LogicalVolume*	    fLogicPMTL;
     G4LogicalVolume*	    fLogicPMTL2;
     G4LogicalVolume*	    fLogicPMTL3;
     G4LogicalVolume*	    fLogicPMTL4;
     G4LogicalVolume*	    fLogicPMTL5;
     G4LogicalVolume*	    fLogicPMTL6;
     G4LogicalVolume*	    fLogicPMTL7;
     G4LogicalVolume*	    fLogicPMTL8;
     G4LogicalVolume*	    fLogicPMTL9;
     G4LogicalVolume*	    fLogicPMTL10;
     G4LogicalVolume*	    fLogicPMTL11;
     G4LogicalVolume*	    fLogicPMTL12;
     G4LogicalVolume*	    fLogicPMTL13;
     G4LogicalVolume*	    fLogicPMTL14;
     G4LogicalVolume*	    fLogicPMTL15;
     G4LogicalVolume*	    fLogicPMTL16;
     G4LogicalVolume*	    fLogicPMTL17;
     G4LogicalVolume*	    fLogicPMTL18;
     G4LogicalVolume*	    fLogicPMTL19;
     G4LogicalVolume*	    fLogicPMTL20;
     G4LogicalVolume*	    fLogicPMTL21;
     G4LogicalVolume*	    fLogicPMTL22;
     G4LogicalVolume*	    fLogicPMTL23;
     G4LogicalVolume*	    fLogicPMTL24;
     G4LogicalVolume*	    fLogicPMTL25;
     G4LogicalVolume*	    fLogicPMTL26;
     G4LogicalVolume*	    fLogicPMTL27;
     G4LogicalVolume*	    fLogicPMTL28;
     G4LogicalVolume*	    fLogicPMTL29;
     G4LogicalVolume*	    fLogicPMTL30;
     G4VPhysicalVolume*     fPhysiPMTL;
     G4VPhysicalVolume*     fPhysiPMTL2;
     G4VPhysicalVolume*     fPhysiPMTL3;
     G4VPhysicalVolume*     fPhysiPMTL4;
     G4VPhysicalVolume*     fPhysiPMTL5;
     G4VPhysicalVolume*     fPhysiPMTL6;
     G4VPhysicalVolume*     fPhysiPMTL7;
     G4VPhysicalVolume*     fPhysiPMTL8;
     G4VPhysicalVolume*     fPhysiPMTL9;
     G4VPhysicalVolume*     fPhysiPMTL10;
     G4VPhysicalVolume*     fPhysiPMTL11;
     G4VPhysicalVolume*     fPhysiPMTL12;
     G4VPhysicalVolume*     fPhysiPMTL13;
     G4VPhysicalVolume*     fPhysiPMTL14;
     G4VPhysicalVolume*     fPhysiPMTL15;
     G4VPhysicalVolume*     fPhysiPMTL16;
     G4VPhysicalVolume*     fPhysiPMTL17;
     G4VPhysicalVolume*     fPhysiPMTL18;
     G4VPhysicalVolume*     fPhysiPMTL19;
     G4VPhysicalVolume*     fPhysiPMTL20;
     G4VPhysicalVolume*     fPhysiPMTL21;
     G4VPhysicalVolume*     fPhysiPMTL22;
     G4VPhysicalVolume*     fPhysiPMTL23;
     G4VPhysicalVolume*     fPhysiPMTL24;
     G4VPhysicalVolume*     fPhysiPMTL25;
     G4VPhysicalVolume*     fPhysiPMTL26;
     G4VPhysicalVolume*     fPhysiPMTL27;
     G4VPhysicalVolume*     fPhysiPMTL28;
     G4VPhysicalVolume*     fPhysiPMTL29;
     G4VPhysicalVolume*     fPhysiPMTL30;
     
     G4double 		        fPMTDiameter;
     G4double 		        fPMTLength;
     G4double		        fZposPMT;
     G4double 		        fPMTDiameter1;
     G4double 		        fPMTLength1;
     G4double		        fZposPMT1;
     G4Material* 	        fPMTMaterial;
     
     // PMT Optical Window
     G4Tubs* 		        fSolidPMTWin;
     G4Polyhedra*           fSolidHPMTWin;
     G4Tubs* 		        fSolidPMTWin1;
     G4Polyhedra*           fSolidHPMTWin1;
     G4LogicalVolume*	    fLogicPMTWin;
     G4VPhysicalVolume*     fPhysiPMTWin;
     G4LogicalVolume*	    fLogicPMTWin1;
     G4VPhysicalVolume*     fPhysiPMTWin1;
     G4double 		        fPMTWinThickness;
     G4double		        fZposPMTWin;
     G4double 		        fPMTWinThickness1;
     G4double		        fZposPMTWin1;
     G4Material* 	        fPMTWinMaterial;
     G4Material*            fPhotoCathodeMaterial;
     G4Material*            fPMTIntMaterial;
     
     DetectorMessenger* fDetectorMessenger;
      
  private:
    
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter(); 
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif
