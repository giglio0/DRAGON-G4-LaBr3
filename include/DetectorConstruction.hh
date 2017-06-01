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
#include "G4Tubs.hh"
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
class G4Tubs;
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
     void SetDetectorGeometry           (G4int);
     void SetTotalDetectorDiameter   (G4double);
     void SetDetectorDiameter        (G4double);
     void SetTotalDetectorLength     (G4double);
     void SetDetectorLength	     (G4double);
     void SetGapThickness	     (G4double);
     void SetGapMaterial             (G4String);
     void SetAlCaseThickness         (G4double);
     void SetPbCaseThickness         (G4double);
     void SetPMTDiameter	     (G4double);
     void SetPMTLength		     (G4double);
     void SetMaterialPropertiesTable (G4String);
     
     virtual G4VPhysicalVolume* Construct();
     virtual void ConstructSDandField();

     void UpdateGeometry();

  public:


     void PrintCalorParameters();

     G4Material* GetWorldMaterial()     {return fWorldMaterial;};
     G4double    GetWorldSizeX()        {return fWorldSizeX;};
     G4double    GetWorldSizeYZ()       {return fWorldSizeYZ;};

     G4Material* GetDetectorMaterial()           {return fDetectorMaterial;};
     G4Material* GetResolutionDetectorMaterial() {return fResolutionDetectorMaterial;};
     G4int       GetDetectorGeometry()           {return fDetectorGeometry;};
     G4double    GetTotalDetectorDiameter()      {return fTotalDetectorDiameter;};
     G4double	 GetDetectorDiameter()		 {return fDetectorDiameter;};
     G4double    GetTotalDetectorLength()        {return fTotalDetectorLength;}; 
     G4double	 GetDetectorLength()		 {return fDetectorLength;};
     G4double	 GetGapThickness()		 {return fGapThickness;}; 
     G4Material* GetGapMaterial()    	         {return fGapMaterial;};
     G4double	 GetAlCaseThickness()	         {return fAlCaseThickness;}; 
     G4double	 GetPbCaseThickness()	         {return fPbCaseThickness;}; 
     G4double	 GetPMTDiameter()		 {return fPMTDiameter;}; 
     G4double	 GetPMTLength()			 {return fPMTLength;};
     G4VPhysicalVolume*    GetCrystal()          {return fPhysiCrystal;};
		 
     const G4VPhysicalVolume* GetphysiWorld() {return fPhysiWorld;};
     const G4VPhysicalVolume* GetScint() {return fPhysiCrystal;};

  private:
     G4VisAttributes* fBlueVisAtt;
     G4VisAttributes* fGreenVisAtt;
     G4VisAttributes* fRedVisAtt;
     G4VisAttributes* fGreyVisAtt;
     G4VisAttributes* fWireFrameVisAtt; 
     G4VisAttributes* fAuxEdgeVisAtt;

     // World
     G4Box*             fSolidWorld;
     G4LogicalVolume*   fLogicWorld;
     G4VPhysicalVolume* fPhysiWorld;
     G4Material*        fWorldMaterial;
     G4double           fWorldSizeX;
     G4double           fWorldSizeYZ;

     //LaBr3 Components
     //LaBr3 Detector Mother Holder
     G4Tubs* 		    fSolidDetector;
     G4Polyhedra*           fSolidDetector1;
     G4LogicalVolume*	    fLogicDetector;
     G4VPhysicalVolume*     fPhysiDetector;
     G4double               fTotalDetectorDiameter;
     G4double 		    fDetectorDiameter;
     G4double               fTotalDetectorLength;
     G4double 		    fDetectorLength;
     G4double		    fZposDetector;
     G4int                  fDetectorGeometry;

     // Crystal
     G4Tubs* 		    fSolidCrystal;
     G4Polyhedra*           fSolidCrystal1;
     G4LogicalVolume*	    fLogicCrystal;
     G4VPhysicalVolume*     fPhysiCrystal;
     G4double 		    fDiameter;
     G4double 		    fLength;
     G4double		    fZpos;
     G4double               fEnergyCrystal;
     G4Material* 	    fDetectorMaterial;
     G4Material*            fResolutionDetectorMaterial;

     //Gap between scintillator and aluminum casing (reflector?)
     //around crystal
     G4Tubs* 		    fSolidGap;
     G4Polyhedra*           fSolidGap1;
     G4LogicalVolume*	    fLogicGap;
     G4VPhysicalVolume*     fPhysiGap;
     G4double 		    fGapLength;
     G4double		    fZposGap;
     //over face of crystal
     G4Tubs* 		    fSolidFaceGap;
     G4Polyhedra*           fSolidFaceGap1;
     G4LogicalVolume*	    fLogicFaceGap;
     G4VPhysicalVolume*     fPhysiFaceGap;
     G4double		    fZposFaceGap;
     //common to both components of gap
     G4double 		    fGapDiameter;
     G4double 		    fGapThickness;
     G4Material* 	    fGapMaterial;

     //Aluminum casing surrounding scintillator crystal
     //around crystal
     G4Tubs* 		    fSolidAlCase;
     G4Polyhedra*           fSolidAlCase1;
     G4LogicalVolume*	    fLogicAlCase;
     G4VPhysicalVolume*     fPhysiAlCase;
     G4double		    fAlCaseLength;
     G4double		    fZposAlCase;
     //over face of crystal
     G4Tubs* 		    fSolidFaceAlCase;
     G4Polyhedra*           fSolidFaceAlCase1;
     G4LogicalVolume*	    fLogicFaceAlCase;
     G4VPhysicalVolume*     fPhysiFaceAlCase;
     G4double		    fZposFaceAlCase;
     //common to both components of aluminum case
     G4double 		    fAlCaseDiameter;
     G4double 		    fAlCaseThickness;
     G4Material* 	    fAlCaseMaterial;

     //Lead casing surrounding aluminum casing
     //around casing surround
     G4Tubs* 		    fSolidPbCase;
     G4Polyhedra*           fSolidPbCase1;
     G4LogicalVolume*	    fLogicPbCase;
     G4VPhysicalVolume*     fPhysiPbCase;
      //at back in front of PMT
     G4Tubs* 		    fSolidPbCollar;
     G4Polyhedra*           fSolidPbCollar1;
     G4LogicalVolume*	    fLogicPbCollar;
     G4VPhysicalVolume*     fPhysiPbCollar;
     G4double		    fZposPbCollar;
     //common to both components of lead surround
     G4double 		    fPbCaseDiameter;
     G4double 		    fPbCaseThickness;
     G4Material* 	    fPbCaseMaterial;

     // PMT
     // PMT Body
     G4Tubs* 		    fSolidPMT;
     G4Polyhedra*           fSolidPMT1;
     G4LogicalVolume*	    fLogicPMT;
     G4VPhysicalVolume*     fPhysiPMT;
     G4double 		    fPMTDiameter;
     G4double 		    fPMTLength;
     G4double		    fZposPMT;
     G4Material* 	    fPMTMaterial;
     
     // PMT Optical Window
     G4Tubs* 		    fSolidPMTWin;
     G4Polyhedra*           fSolidPMTWin1;
     G4LogicalVolume*	    fLogicPMTWin;
     G4VPhysicalVolume*     fPhysiPMTWin;
     G4double 		    fPMTWinThickness;
     G4double		    fZposPMTWin;
     G4Material* 	    fPMTWinMaterial;
     
     DetectorMessenger* fDetectorMessenger;
      
  private:
    
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

