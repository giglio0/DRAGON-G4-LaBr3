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
     G4VPhysicalVolume*    GetCrystal()          {return fPhysiCrystal;};
	 G4VPhysicalVolume*    GetCrystal1()          {return fPhysiCrystal1;}; 
	 
     const G4VPhysicalVolume* GetphysiWorld() {return fPhysiWorld;};
     const G4VPhysicalVolume* GetScint() {return fPhysiCrystal;};
     const G4VPhysicalVolume* GetScint1() {return fPhysiCrystal1;};

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
     G4LogicalVolume*	    fLogicDetector;
     G4LogicalVolume*	    fLogicDetector1;
     G4VPhysicalVolume*     fPhysiDetector;
     G4VPhysicalVolume*     fPhysiDetector1;
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
     G4Polyhedra*           fSolidHCrystal;
     G4Tubs* 		        fSolidCrystal1;
     G4Polyhedra*           fSolidHCrystal1;
     G4LogicalVolume*	    fLogicCrystal;
     G4VPhysicalVolume*     fPhysiCrystal;
     G4LogicalVolume*	    fLogicCrystal1;
     G4VPhysicalVolume*     fPhysiCrystal1;
     G4double 		        fDiameter;
     G4double 		        fLength;
     G4double		        fZpos;
     G4double		        fZpos1;
     G4Material*            fCrystalMaterial;
     G4Material*            fCrystalMaterial1;
     
     //Gap between scintillator and aluminum casing (reflector?)
     //around crystal
     G4Tubs* 		        fSolidGap;
     G4Polyhedra*           fSolidHGap;
     G4Tubs* 		        fSolidGap1;
     G4Polyhedra*           fSolidHGap1;
     G4LogicalVolume*	    fLogicGap;
     G4VPhysicalVolume*     fPhysiGap;
     G4LogicalVolume*	    fLogicGap1;
     G4VPhysicalVolume*     fPhysiGap1;
     G4double 		        fGapLength;
     G4double		        fZposGap;
     G4double 		        fGapLength1;
     G4double		        fZposGap1;
     G4Material* 	        fGapMaterial;
     //over face of crystal
     G4Tubs* 		        fSolidFaceGap;
     G4Polyhedra*           fSolidHFaceGap;
     G4Tubs* 		        fSolidFaceGap1;
     G4Polyhedra*           fSolidHFaceGap1;
     G4LogicalVolume*	    fLogicFaceGap;
     G4VPhysicalVolume*     fPhysiFaceGap;
     G4LogicalVolume*	    fLogicFaceGap1;
     G4VPhysicalVolume*     fPhysiFaceGap1;
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
     G4Polyhedra*           fSolidHAlCase;
     G4Tubs* 		        fSolidAlCase1;
     G4Polyhedra*           fSolidHAlCase1;
     G4LogicalVolume*	    fLogicAlCase;
     G4VPhysicalVolume*     fPhysiAlCase;
     G4LogicalVolume*	    fLogicAlCase1;
     G4VPhysicalVolume*     fPhysiAlCase1;
     G4double		        fAlCaseLength;
     G4double		        fZposAlCase;
     G4double		        fAlCaseLength1;
     G4double		        fZposAlCase1;
     G4Material* 	        fAlCaseMaterial;
     //over face of crystal
     G4Tubs* 		        fSolidFaceAlCase;
     G4Polyhedra*           fSolidHFaceAlCase;
     G4Tubs* 		        fSolidFaceAlCase1;
     G4Polyhedra*           fSolidHFaceAlCase1;
     G4LogicalVolume*	    fLogicFaceAlCase;
     G4VPhysicalVolume*     fPhysiFaceAlCase;
     G4LogicalVolume*	    fLogicFaceAlCase1;
     G4VPhysicalVolume*     fPhysiFaceAlCase1;
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
     G4Tubs*                fSolidPMTI;
     G4Tubs* 				fSolidPMTJ;
     G4Tubs*				fSolidPMTK;
     G4Tubs*				fSolidPMTL;
     G4Polyhedra*           fSolidHPMT;
     G4Tubs* 		        fSolidPMT1;
     G4Polyhedra*           fSolidHPMT1;
     
     G4LogicalVolume*	    fLogicPMT;
     G4VPhysicalVolume*     fPhysiPMT;
     G4LogicalVolume*	    fLogicPMT1;
     G4VPhysicalVolume*     fPhysiPMT1;
     G4LogicalVolume*	    fLogicPMTI;
     G4VPhysicalVolume*     fPhysiPMTI;
     G4LogicalVolume*	    fLogicPMTJ;
     G4VPhysicalVolume*     fPhysiPMTJ;
     G4LogicalVolume*	    fLogicPMTK;
     G4VPhysicalVolume*     fPhysiPMTK;
     G4LogicalVolume*	    fLogicPMTL;
     G4VPhysicalVolume*     fPhysiPMTL;
     
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
