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

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

class G4Box;
class G4Tubs;
class G4CutTubs;
class G4Trd;
class G4VPhysicalVolume;
class G4Material;
class G4MaterialCutsCouple;
class G4UniformMagField;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
   ~DetectorConstruction();

  public:

     void SetCrystalDiameter 	(G4double);
     void SetCrystalLength	(G4double);
     void SetGapThickness	(G4double);
     void SetGapMaterial    (G4String);
     void SetAlCaseThickness(G4double);
     void SetAlFaceThickness(G4double);
     void SetPbCaseThickness(G4double);
     void SetPMTDiameter	(G4double);
     void SetPMTLength		(G4double);
     void SetNoOfDetectors(G4int);
     
     virtual G4VPhysicalVolume* Construct();
//     virtual void ConstructSDandField();

     void PrintGeometryParameters();
     void UpdateGeometry();

  public:

     G4Material* GetWorldMaterial()     {return fWorldMaterial;};
     G4double    GetWorldSizeX()        {return fWorldSizeX;};
     G4double    GetWorldSizeYZ()       {return fWorldSizeYZ;};

     G4double	 GetDetectorDiameter()		{return fDetectorDiameter;}; 
     G4double	 GetDetectorLength()		{return fDetectorLength;}; 
     G4double	 GetGapThickness()		{return fGapThickness;}; 
     G4Material* GetGapMaterial()    	{return fGapMaterial;};
     G4double	 GetAlCaseThickness()	{return fAlCaseThickness;}; 
     G4double	 GetAlFaceThickness()	{return fAlFaceThickness;}; 
     G4double	 GetPbCaseThickness()	{return fPbCaseThickness;}; 
     G4double	 GetPMTDiameter()		{return fPMTDiameter;}; 
     G4double	 GetPMTLength()			{return fPMTLength;};
     G4int		 GetSetNoOfDetectors()  {return fNoOfDets;};
		 
     const G4VPhysicalVolume* GetphysiWorld() {return fPhysiWorld;};
     const G4VPhysicalVolume* GetScint() {return fPhysiDetectorCrystal;};

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


     //Detector Components
     
     G4int			fNoOfDets;
     //Detector Detector Mother Holder
     G4Tubs* 		fSolidDetector;
     G4LogicalVolume*	fLogicDetector;
     G4VPhysicalVolume* fPhysiDetector;
     G4double 		fDetectorDiameter;
     G4double 		fDetectorLength;
     G4double		fZposDetector;

     //Detector Crystal
     G4Tubs* 		fSolidDetectorCrystal;
     G4LogicalVolume*	fLogicDetectorCrystal;
     G4VPhysicalVolume* fPhysiDetectorCrystal;
     G4double 		fCrystalDiameter;
     G4double 		fCrystalLength;
     G4double		fZposCrystal;
     G4Material* 	fCrystalMaterial;

     //Gap between scintillator and aluminum casing (reflector?)
     //around crystal
     G4Tubs* 		fSolidGap;
     G4LogicalVolume*	fLogicGap;
     G4VPhysicalVolume* fPhysiGap;
     G4double 		fGapLength;
     G4double		fZposGap;
     //over face of crystal
     G4Tubs* 		fSolidFaceGap;
     G4LogicalVolume*	fLogicFaceGap;
     G4VPhysicalVolume* fPhysiFaceGap;
     G4double		fZposFaceGap;
     //common to both components of gap
     G4double 		fGapDiameter;
     G4double 		fGapThickness;
     G4Material* 	fGapMaterial;

     //Aluminum casing surrounding scintillator crystal
     //around crystal
     G4Tubs* 		fSolidAlCase;
     G4LogicalVolume*	fLogicAlCase;
     G4VPhysicalVolume* fPhysiAlCase;
     G4double		fAlCaseLength;
     G4double		fZposAlCase;
     G4double 		fAlCaseThickness;     
     //over face of crystal
     G4Tubs* 		fSolidFaceAlCase;
     G4LogicalVolume*	fLogicFaceAlCase;
     G4VPhysicalVolume* fPhysiFaceAlCase;
     G4double		fZposFaceAlCase;
     G4double		fAlFaceThickness;
     //common to both components of aluminum case
     G4double 		fAlCaseDiameter;
     G4Material* 	fAlCaseMaterial;

     //Lead casing surrounding aluminum casing
     //around casing surround
     G4Tubs* 		fSolidPbCase;
     G4LogicalVolume*	fLogicPbCase;
     G4VPhysicalVolume* fPhysiPbCase;
      //at back in front of PMT
     G4Tubs* 		fSolidPbCollar;
     G4LogicalVolume*	fLogicPbCollar;
     G4VPhysicalVolume* fPhysiPbCollar;
     G4double		fZposPbCollar;
     //common to both components of lead surround
     G4double 		fPbCaseDiameter;
     G4double 		fPbCaseThickness;
     G4Material* 	fPbCaseMaterial;

     //Detector PMT
     // PMT Body
     G4Tubs* 		fSolidPMT;
     G4LogicalVolume*	fLogicPMT;
     G4VPhysicalVolume* fPhysiPMT;
     G4double 		fPMTDiameter;
     G4double 		fPMTLength;
     G4double		fZposPMT;
     G4Material* 	fPMTMaterial;
     // PMT Optical Window
     G4Tubs* 		fSolidPMTWin;
     G4LogicalVolume*	fLogicPMTWin;
     G4VPhysicalVolume* fPhysiPMTWin;
     G4double 		fPMTWinDiameter;
     G4double 		fPMTWinThickness;
     G4double		fZposPMTWin;
     G4Material* 	fPMTWinMaterial;
     
     DetectorMessenger* fDetectorMessenger;
      
  private:
    
     void DefineMaterials();
     void ComputeGeometryParameters();
     G4VPhysicalVolume* ConstructGeometry();
     void ReadInGeometryFile();
     void ConstructTarget();
     void ConstructDetectorComponents();    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

