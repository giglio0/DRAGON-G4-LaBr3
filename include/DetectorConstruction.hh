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

#include "G4Trap.hh"

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

class G4Trap;

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

     void SetPbCaseMaterial          (G4String);

     void SetPbCollarMaterial        (G4String);

     void SetPMTMaterial             (G4String);

     void SetPMTWinMaterial          (G4String);

     void SetAlCaseThickness         (G4double);

     void SetPbCaseThickness         (G4double);

     void SetPMTDiameter	         (G4double);

     void SetPMTLength		         (G4double);

     void SetMaterialPropertiesTable (G4String);

     

     virtual G4VPhysicalVolume* Construct();

     //virtual void ConstructSDandField();

     

     void UpdateGeometry();



  public:



     void PrintCalorParameters();

     

     G4Material* GetWorldMaterial()     {return fWorldMaterial;};

     G4double    GetWorldSizeX()        {return fWorldSizeX;};

     G4double    GetWorldSizeYZ()       {return fWorldSizeYZ;};



     G4Material* GetDetectorMaterial()           {return fDetectorMaterial;};

     G4Material* GetGasTargetMaterial()          {return fGasTargetMaterial;};

     G4double    GetTemperature()                {return Temperature;};

     G4double    GetPressure()                   {return Pressure;};

     G4int       GetDetectorGeometry()           {return fDetectorGeometry;};

     G4double    GetTotalDetectorDiameter()      {return fTotalDetectorDiameter;};

     G4double	 GetDetectorDiameter()		     {return fDetectorDiameter;};

     G4double    GetTotalDetectorLength()        {return fTotalDetectorLength;}; 

     G4double	 GetDetectorLength()		     {return fDetectorLength;};

     G4double	 GetGapThickness()		         {return fGapThickness;}; 

     G4Material* GetCrystalMaterial()            {return fCrystalMaterial;};

     G4Material* GetGapMaterial()    	         {return fGapMaterial;};

     G4Material* GetFaceGapMaterial()            {return fFaceGapMaterial;};

     G4Material* GetAlCaseMaterial()             {return fAlCaseMaterial;};

     G4Material* GetFaceAlCaseMaterial()         {return fFaceAlCaseMaterial;};

     G4Material* GetPbCaseMaterial()             {return fPbCaseMaterial;};

     G4Material* GetPbCollarMaterial()           {return fAlCaseMaterial;};

     G4Material* GetPMTMaterial()                {return fPMTMaterial;};

     G4Material* GetPMTWinMaterial()             {return fPMTWinMaterial;};

     G4double	 GetAlCaseThickness()	         {return fAlCaseThickness;}; 

     G4double	 GetPbCaseThickness()	         {return fPbCaseThickness;}; 

     G4double	 GetPMTDiameter()		         {return fPMTDiameter;}; 

     G4double	 GetPMTLength()			         {return fPMTLength;};

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



     G4SubtractionSolid*    fSolidGasTarget;

     G4LogicalVolume*       fLogicGasTarget;

     G4VPhysicalVolume*     fPhysiGasTarget;

     G4SubtractionSolid*    fSolidGasTargetSideCuts;

     G4LogicalVolume*       fLogicGasTargetSideCuts;

     G4VPhysicalVolume*     fPhysiGasTargetSideCuts;

     G4Tubs*                fSolidInnerTube;

     G4LogicalVolume*       fLogicInnerTube;

     G4VPhysicalVolume*     fPhysiInnerTube;

     G4SubtractionSolid*    fSolidGasTargetSides;

     G4LogicalVolume*       fLogicGasTargetSides;

     G4VPhysicalVolume*     fPhysiGasTargetSides;

     

     // Collimator 

     G4SubtractionSolid*    fSolidLeftRing;

     G4LogicalVolume*       fLogicLeftRing;

     G4VPhysicalVolume*     fPhysiLeftRing;

     G4SubtractionSolid*    fSolidRightRing;

     G4LogicalVolume*       fLogicRightRing;

     G4VPhysicalVolume*     fPhysiRightRing;

     G4SubtractionSolid*    fSolidLeftConnectingRing;

     G4LogicalVolume*       fLogicLeftConnectingRing;

     G4VPhysicalVolume*     fPhysiLeftConnectingRing;

     G4SubtractionSolid*    fSolidRightConnectingRing;

     G4LogicalVolume*       fLogicRightConnectingRing;

     G4VPhysicalVolume*     fPhysiRightConnectingRing;

     G4SubtractionSolid*    fSolidLeftExternalRing;

     G4LogicalVolume*       fLogicLeftExternalRing;

     G4VPhysicalVolume*     fPhysiLeftExternalRing;

     G4SubtractionSolid*    fSolidRightExternalRing;

     G4LogicalVolume*       fLogicRightExternalRing;

     G4VPhysicalVolume*     fPhysiRightExternalRing;

     G4SubtractionSolid*    fSolidLeftOuterRing;

     G4LogicalVolume*       fLogicLeftOuterRing;

     G4VPhysicalVolume*     fPhysiLeftOuterRing;

     G4SubtractionSolid*    fSolidRightOuterRing;

     G4LogicalVolume*       fLogicRightOuterRing;

     G4VPhysicalVolume*     fPhysiRightOuterRing;

    

     G4CutTubs*             fSolidOblique;

     G4LogicalVolume*       fLogicOblique;

     G4VPhysicalVolume*     fPhysiOblique;

     G4SubtractionSolid*    fSolidTrapezoidCut;

     G4LogicalVolume*       fLogicTrapezoidCut;

     G4VPhysicalVolume*     fPhysiTrapezoidCut;

     G4SubtractionSolid*    fSolidTube;

     G4LogicalVolume*       fLogicTube;

     G4VPhysicalVolume*     fPhysiTube;

     G4Tubs*                fSolidDisk;

     G4LogicalVolume*       fLogicDisk;

     G4VPhysicalVolume*     fPhysiDisk;

     G4Material*            fGasTargetMaterial;

     G4Material*            fTubeMaterial;

     G4double               Temperature;

     G4double               Pressure;

     

     //LaBr3 Components

     //LaBr3 Detector Mother Holder

     G4Tubs* 		        fSolidDetector;

     G4Polyhedra*           fSolidHDetector;

     G4LogicalVolume*	    fLogicDetector;

     G4VPhysicalVolume*     fPhysiDetector;

     G4double               fTotalDetectorDiameter;

     G4double 		        fDetectorDiameter;

     G4double               fTotalDetectorLength;

     G4double 		        fDetectorLength;

     G4double		        fZposDetector;

     G4int                  fDetectorGeometry;

     G4Material*            fDetectorMaterial;

     

     // Crystal

     G4Tubs* 		        fSolidCrystal;

     G4Polyhedra*           fSolidHCrystal;

     G4LogicalVolume*	    fLogicCrystal;

     G4VPhysicalVolume*     fPhysiCrystal;

     G4double 		        fDiameter;

     G4double 		        fLength;

     G4double		        fZpos;

     G4Material*            fCrystalMaterial;

     

     //Gap between scintillator and aluminum casing (reflector?)

     //around crystal

     G4Tubs* 		        fSolidGap;

     G4Polyhedra*           fSolidHGap;

     G4LogicalVolume*	    fLogicGap;

     G4VPhysicalVolume*     fPhysiGap;

     G4double 		        fGapLength;

     G4double		        fZposGap;

     G4Material* 	        fGapMaterial;

     //over face of crystal

     G4Tubs* 		        fSolidFaceGap;

     G4Polyhedra*           fSolidHFaceGap;

     G4LogicalVolume*	    fLogicFaceGap;

     G4VPhysicalVolume*     fPhysiFaceGap;

     G4double		        fZposFaceGap;

     G4Material* 	        fFaceGapMaterial;

     //common to both components of gap

     G4double 		        fGapDiameter;

     G4double 		        fGapThickness;

     

      //Aluminum casing surrounding scintillator crystal

     //around crystal

     G4Tubs* 		        fSolidAlCase;

     G4Polyhedra*           fSolidHAlCase;

     G4LogicalVolume*	    fLogicAlCase;

     G4VPhysicalVolume*     fPhysiAlCase;

     G4double		        fAlCaseLength;

     G4double		        fZposAlCase;

     G4Material* 	        fAlCaseMaterial;

     //over face of crystal

     G4Tubs* 		        fSolidFaceAlCase;

     G4Polyhedra*           fSolidHFaceAlCase;

     G4LogicalVolume*	    fLogicFaceAlCase;

     G4VPhysicalVolume*     fPhysiFaceAlCase;

     G4double		        fZposFaceAlCase;

     G4Material* 	        fFaceAlCaseMaterial;

     //common to both components of aluminum case

     G4double 		        fAlCaseDiameter;

     G4double 		        fAlCaseThickness;



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

     G4Polyhedra*           fSolidHPMT;

     G4LogicalVolume*	    fLogicPMT;

     G4VPhysicalVolume*     fPhysiPMT;

     G4double 		        fPMTDiameter;

     G4double 		        fPMTLength;

     G4double		        fZposPMT;

     G4Material* 	        fPMTMaterial;

     

     // PMT Optical Window

     G4Tubs* 		        fSolidPMTWin;

     G4Polyhedra*           fSolidHPMTWin;

     G4LogicalVolume*	    fLogicPMTWin;

     G4VPhysicalVolume*     fPhysiPMTWin;

     G4double 		        fPMTWinThickness;

     G4double		        fZposPMTWin;

     G4Material* 	        fPMTWinMaterial;

     

     DetectorMessenger* fDetectorMessenger;

      

  private:

    

     void DefineMaterials();

     void ComputeCalorParameters();

     G4VPhysicalVolume* ConstructCalorimeter(); 

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
