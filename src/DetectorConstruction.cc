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
// Detector_v4
// src/DetectorConstruction.cc
// \brief Implementation of the DetectorConstruction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
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

#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: fSolidWorld(0),fLogicWorld(0),fPhysiWorld(0),
fSolidOutTBox(0), fLogicOutTBox(0),fPhysiOutTBox(0),
fSolidInTBox(0),fLogicInTBox(0),fPhysiInTBox(0),
fSolidEntHoleBox(0),fLogicEntHoleBox(0),fPhysiEntHoleBox(0),
fSolidExitHoleBox(0),fLogicExitHoleBox(0),fPhysiExitHoleBox(0),
fSolidOutCell(0),fLogicOutCell(0),fPhysiOutCell(0),
fSolidInCell(0),fLogicInCell(0),fPhysiInCell(0),
fSolidEntHole(0),fLogicEntHole(0),fPhysiEntHole(0),
fSolidEntDisk(0),fLogicEntDisk(0),fPhysiEntDisk(0),
fSolidEntCol(0),fLogicEntCol(0),fPhysiEntCol(0),
fSolidExitHole(0),fLogicExitHole(0),fPhysiExitHole(0),
fSolidExitDisk(0),fLogicExitDisk(0),fPhysiExitDisk(0),
fSolidExitCol(0),fLogicExitCol(0),fPhysiExitCol(0),
fSolidDetector(0),fLogicDetector(0),fPhysiDetector(0),
fSolidDetectorCrystal(0),fLogicDetectorCrystal(0),fPhysiDetectorCrystal(0),
fSolidGap(0),fLogicGap(0),fPhysiGap(0),
fSolidFaceGap(0),fLogicFaceGap(0),fPhysiFaceGap(0),
fSolidAlCase(0),fLogicAlCase(0),fPhysiAlCase(0),
fSolidFaceAlCase(0),fLogicFaceAlCase(0),fPhysiFaceAlCase(0),
fSolidPbCase(0),fLogicPbCase(0),fPhysiPbCase(0),
fSolidPbCollar(0),fLogicPbCollar(0),fPhysiPbCollar(0),
fSolidPMT(0),fLogicPMT(0),fPhysiPMT(0),
fSolidPMTWin(0),fLogicPMTWin(0),fPhysiPMTWin(0)
     
{
  //Setting the default for the world for the LaBr Sim
  fWorldSizeX = 100*cm;
  fWorldSizeYZ = 100*cm;

  //Set Detector Defaults
  fPMTDiameter = 5.1*cm;
  fPMTLength = 19.4*cm;
  fPMTWinThickness = 0.508*cm;
  fCrystalLength = 7.6*cm;
  fCrystalDiameter = 5.58*cm;
  fGapThickness = 0.01*cm;
  fAlCaseThickness = 0.05*cm;
  fAlFaceThickness = 0.025*cm;
  fPbCaseThickness = 0.5*cm;
  fOutTBoxX = 8.573*cm;
  fOutTBoxY = 2.381*cm;
  fOutTBoxZ = 12.859*cm;
  fNoOfDets = 10;
  ComputeGeometryParameters();
  
  // materials  
  DefineMaterials();
  // Default Materials
  fWorldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Galactic");
  fCrystalMaterial	= G4NistManager::Instance()->FindOrBuildMaterial("LaBr");
  fAlCaseMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  fPbCaseMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
  fPMTWinMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pyrex_Glass");
  fPMTMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  //Set Default Gap Material i.e. no reflector
  fGapMaterial 	= G4NistManager::Instance()->FindOrBuildMaterial("Galactic");

 //Set Visualization Attributes

  // Blue
  fBlueVisAtt = new G4VisAttributes(G4Colour(0.,0.,0.75));
  // Green
  fGreenVisAtt = new G4VisAttributes(G4Colour(0.,0.75,0.));
  // Red
  fRedVisAtt = new G4VisAttributes(G4Colour(0.75,0.,0.));
  // Grey
  fGreyVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  // WireFrame
  fWireFrameVisAtt = new G4VisAttributes(G4Colour(1.,1.,1.));
  fWireFrameVisAtt->SetForceWireframe(true);
  // Force Aux Edge Visible (Blue colour)
  fAuxEdgeVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.0));
  fAuxEdgeVisAtt->SetForceAuxEdgeVisible(true);

  // create commands for interactive definition of the calorimeter  
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  //This function illustrates the possible ways to define materials
 
  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;  

  G4int ncomponents, natoms;
  G4double fractionmass;
  G4double temperature, pressure;

  // define Elements
  G4Element* H  = new G4Element("Hydrogen",symbol="H",  z= 1, a=   1.01*g/mole);
  G4Element* C  = new G4Element("Carbon",  symbol="C",  z= 6, a=  12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N",  z= 7, a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",  symbol="O",  z= 8, a=  16.00*g/mole);
  G4Element* Na = new G4Element("Sodium",  symbol="Na", z=11, a=  22.99*g/mole);
  G4Element* I  = new G4Element("Iodine",  symbol="I" , z=53, a= 126.90*g/mole);
  G4Element* La = new G4Element("Lanthanum",  symbol="La" , z=57, a= 138.91*g/mole);
  G4Element* Br = new G4Element("Bromine",   symbol="Br", z=35, a= 79.90*g/mole);

  // define simple materials
  new G4Material("H2Liq"    , z= 1, a= 1.01*g/mole, density= 70.8*mg/cm3);
  new G4Material("Beryllium", z= 4, a= 9.01*g/mole, density= 1.848*g/cm3);
  new G4Material("Aluminium", z=13, a=26.98*g/mole, density= 2.700*g/cm3);
  new G4Material("Silicon"  , z=14, a=28.09*g/mole, density= 2.330*g/cm3);
  new G4Material("Iron",     z=26, a= 55.85*g/mole, density= 7.870*g/cm3);
  new G4Material("Copper",   z=29, a= 63.55*g/mole, density= 8.960*g/cm3);
  new G4Material("Germanium",z=32, a= 72.61*g/mole, density= 5.323*g/cm3);
  new G4Material("Silver",   z=47, a=107.87*g/mole, density= 10.50*g/cm3);
  new G4Material("Tungsten", z=74, a=183.85*g/mole, density= 19.30*g/cm3);
  new G4Material("Gold",     z=79, a=196.97*g/mole, density= 19.32*g/cm3);
  new G4Material("Lead",     z=82, a=207.19*g/mole, density= 11.35*g/cm3);

  // define a material from elements.   case 1: chemical molecule
  G4Material* CH = new G4Material("Plastic", density= 1.04*g/cm3, ncomponents=2);
  CH->AddElement(C, natoms=1);
  CH->AddElement(H, natoms=1);

  G4Material* NaI = new G4Material("NaI", density= 3.67*g/cm3, ncomponents=2);
  NaI->AddElement(Na, natoms=1);
  NaI->AddElement(I , natoms=1);
  NaI->GetIonisation()->SetMeanExcitationEnergy(452*eV);

  // define a material from elements.   case 2: mixture by fractional mass
  G4Material* LaBr = new G4Material("LaBr", density= 5.06*g/cm3, ncomponents=2);
  LaBr->AddElement(La, fractionmass = 0.366875);
  LaBr->AddElement(Br ,fractionmass = 0.633124);
  LaBr->GetIonisation()->SetMeanExcitationEnergy(454.5*eV);

  G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  // example of vacuum
  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  new G4Material("Galactic", z=1, a=1.01*g/mole,density,
                 kStateGas,temperature,pressure);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ComputeGeometryParameters()
{
  // Compute derived parameters of the calorimeter
  fGapLength = fCrystalLength + fPMTWinThickness + fGapThickness;
  fGapDiameter = fCrystalDiameter + 2*fGapThickness;
  fAlCaseLength = fGapLength + fAlFaceThickness;
  fAlCaseDiameter = fGapDiameter + 2*fAlCaseThickness;
  fPbCaseDiameter = fAlCaseDiameter + 2*fPbCaseThickness;
  fDetectorLength = fPMTWinThickness + fCrystalLength + fPMTLength + fGapThickness + fAlFaceThickness;
  //if (fPbCaseDiameter>fPMTDiameter){
	  //fPMTDiameter = fPbCaseDiameter;
	  //G4cout << "\n" << "\n" 
	  //<< "WARNING: PMT Diameter has been adjusted to compensate for large crystal/lead/aluminum size"  
	  //<< G4endl;
	  //fDetectorDiameter = fPMTDiameter;
  //}
  //else
  fDetectorDiameter = fAlCaseDiameter;

  fZposPbCollar = 0.5*fDetectorLength - fAlFaceThickness - fGapThickness - fCrystalLength - fPMTWinThickness + 0.5*fPbCaseThickness;
  fZposFaceAlCase = 0.5*fDetectorLength - 0.5*fAlFaceThickness;
  fZposFaceGap = 0.5*fDetectorLength - fAlFaceThickness - 0.5*fGapThickness ;
  fZposCrystal = 0.5*fDetectorLength - fGapThickness - fAlFaceThickness - 0.5*fCrystalLength;

  fZposAlCase = 0.5*fDetectorLength - 0.5*fAlCaseLength;
  fZposGap = 0.5*fDetectorLength - fAlFaceThickness - 0.5*fGapLength ;
  fZposPMTWin = 0.5*fDetectorLength - fAlFaceThickness - fGapThickness - fCrystalLength - 0.5*fPMTWinThickness;
  fZposPMT = 0.5*fDetectorLength - fAlFaceThickness - fGapThickness - fCrystalLength - fPMTWinThickness - 0.5*fPMTLength;

  fposEastArray = -0.5*fDetectorLength-fOutTBoxY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* DetectorConstruction::ConstructGeometry()
{ 
  // *********************************************************************************
  // Cleanup old geometry
  // *********************************************************************************

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // *********************************************************************************  
  // Compute the geometry parameters 
  // *********************************************************************************

  ComputeGeometryParameters();
        
  // *********************************************************************************
  // World Volume
  // *********************************************************************************
  
  fSolidWorld = new G4Box("World",                                //its name
                   fWorldSizeX/2,fWorldSizeYZ/2,fWorldSizeYZ/2);   //its size
                         
  fLogicWorld = new G4LogicalVolume(fSolidWorld,                //its solid
                                   fWorldMaterial,        //its material
                                   "World");                //its name
                                   
  fPhysiWorld = new G4PVPlacement(0,                        //no rotation
                                   G4ThreeVector(),        //at (0,0,0)
                                 fLogicWorld,                //its logical volume
                                 "World",                //its name
                                 0,                        //its mother  volume
                                 false,                        //no boolean operation
                                 0);                        //copy number
  fLogicWorld->SetVisAttributes(fWireFrameVisAtt);

  // *********************************************************************************
  // DRAGON Gas Target Components
  // *********************************************************************************

  ConstructTarget();

  // *********************************************************************************
  // Detector Components
  // *********************************************************************************

   G4RotationMatrix rotm  = G4RotationMatrix(0,0,0);
   G4ThreeVector position = G4ThreeVector(0.,0.,0.);
   G4Transform3D transform = G4Transform3D(rotm,position); 
   
  // Decide from user input whether to create a geometry with a single detector or 
  // read in an input file for parameters of a multi detector array
  if(fNoOfDets > 1) {
	ConstructDetectorComponents();

	G4int fRotX, fRotY, fRotZ;
	G4double fPosX, fPosY, fPosZ;
	std::ifstream fGeomDataFile ("geom_data.txt");
	if(fGeomDataFile.is_open()){
		G4int fDetectorCopy = 0;
		while(fGeomDataFile >> fRotX >> fRotY >> fRotZ 
						>> fPosX >> fPosY >> fPosZ){
			rotm = G4RotationMatrix(0,0,0);
			rotm.rotateX(fRotX*deg);
			rotm.rotateY(fRotY*deg);
			rotm.rotateZ(fRotZ*deg);
			position = G4ThreeVector(fPosX*cm,fPosY*cm,fPosZ*cm);
			transform = G4Transform3D(rotm,position);
			fPhysiDetector = new G4PVPlacement(transform,
								fLogicDetector, 
								"Detector", 
								fLogicWorld, 
								false, 
								fDetectorCopy,
								false);
			fDetectorCopy++;
		}
		fGeomDataFile.close();
	}
	else
		G4cout << "There was a problem opening the file";	
  }
  else {
  	ConstructDetectorComponents();
    rotm  = G4RotationMatrix(0,0,0);    
	rotm.rotateY(90*deg);
	position = G4ThreeVector(fposEastArray,0.,0.);
	transform = G4Transform3D(rotm,position);
	fPhysiDetector = new G4PVPlacement(transform,
						fLogicDetector, 
						"Detector", 
						fLogicWorld, 
						false, 
						0,
						false);
  }
  // *********************************************************************************
  // Print out all parameters for the geometry components
  // *********************************************************************************

  PrintGeometryParameters();         
  
  //always return the physical World
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructTarget()
{
  // DRAGON Outer Gas Target Box
  // Outer Wall 
  G4RotationMatrix rotm  = G4RotationMatrix(0,0,0);     
  rotm.rotateX(90*deg); 
  rotm.rotateY(-90*deg);
  rotm.rotateZ(0*deg);
  // Position target box so z-axis at beam height
  G4ThreeVector position = G4ThreeVector(0.,-9.684*cm,0.);
  G4Transform3D transform = G4Transform3D(rotm,position);  
  fSolidOutTBox = new G4Box("OutTBox", fOutTBoxX, fOutTBoxY, fOutTBoxZ);
  fLogicOutTBox = new G4LogicalVolume(fSolidOutTBox,fAlCaseMaterial,"OutTBox");
  fPhysiOutTBox = new G4PVPlacement(transform,
        			fLogicOutTBox, 
       				"OutTBox", 
       				fLogicWorld, 
       				false, 
    			    	0,
				    false);			    

  // DRAGON Outer Gas Target Box
  // Inner Wall
  rotm  = G4RotationMatrix(0,0,0);     
  position = G4ThreeVector(0.,0.,0.);
  transform = G4Transform3D(rotm,position);	
  fSolidInTBox = new G4Box("InTBox", 8.256*cm, 2.064*cm, 12.542*cm);
  fLogicInTBox = new G4LogicalVolume(fSolidInTBox,fWorldMaterial,"InTBox");
  fPhysiInTBox = new G4PVPlacement(transform,
        			fLogicInTBox, 
       				"InTBox", 
       				fLogicOutTBox, 
       				false, 
    			    	0,
				    false);
  				    
  // Entrance Hole on Target Box				    
  rotm  = G4RotationMatrix(0,0,0);
  rotm.rotateY(90*deg);    
  position = G4ThreeVector(-8.414*cm,0.,-9.684*cm);
  transform = G4Transform3D(rotm,position);  
  fSolidEntHoleBox = new G4Tubs("EntHoleBox", 0., 0.95*cm, 0.16*cm, 0.*deg, 360.*deg);
  fLogicEntHoleBox = new G4LogicalVolume(fSolidEntHoleBox, fWorldMaterial, "EntHoleBox");
  fPhysiEntHoleBox = new G4PVPlacement(transform,
        			fLogicEntHoleBox, 
       				"EntHoleBox", 
       				fLogicOutTBox, 
       				false, 
    			    	0);
   			    	
  // Exit Hole on Target Box    
  position = G4ThreeVector(8.414*cm,0.,-9.684*cm);
  transform = G4Transform3D(rotm,position);	 
  fSolidExitHoleBox = new G4Tubs("ExitHoleBox", 0., 0.95*cm, 0.16*cm, 0.*deg, 360.*deg);
  fLogicExitHoleBox = new G4LogicalVolume(fSolidExitHoleBox, fWorldMaterial, "ExitHoleBox");
  fPhysiExitHoleBox = new G4PVPlacement(transform,
        			fLogicExitHoleBox, 
       				"ExitHoleBox", 
       				fLogicOutTBox, 
       				false, 
    			    	0);   

  // DRAGON Inner Gas Cell
  // Outer Wall
  rotm  = G4RotationMatrix(0,0,0);
  // position inner gas cell at beam height relative to inner gas target box volume     
  position = G4ThreeVector(0.,0.,-7.655*cm);
  transform = G4Transform3D(rotm,position);	  
  fSolidOutCell = new G4Trd("OutCell",6.759*cm,1.901*cm,1.905*cm,1.905*cm,4.208*cm);
  fLogicOutCell = new G4LogicalVolume(fSolidOutCell,fAlCaseMaterial, "OutCell");
  fPhysiOutCell = new G4PVPlacement(transform,
        			fLogicOutCell, 
       				"OutCell", 
       				fLogicInTBox, 
       				false, 
    			    	0,
				    false);

  // DRAGON Inner Gas Cell
  // Inner Wall			    
  fSolidInCell = new G4Trd("InCell",6.208*cm,1.717*cm,1.588*cm,1.588*cm,3.890*cm);
  fLogicInCell = new G4LogicalVolume(fSolidInCell,fWorldMaterial, "InCell");
  fPhysiInCell = new G4PVPlacement(0,			//no rotation
            		      	G4ThreeVector(),        //no transform, at (0,0,0)
        			fLogicInCell, 
       				"InCell", 
       				fLogicOutCell, 
       				false, 
    			    	0,
				false);

  // Entrance Collimator
  // Entrance Collimator Hole
  rotm  = G4RotationMatrix(0,0,0);     
  rotm.rotateY(-60*deg);
//  position = G4ThreeVector(5.306*cm,0.,-2.008*cm);
  position = G4ThreeVector(-5.351*cm,0.,-2.087*cm);
  transform = G4Transform3D(rotm,position); 
  fSolidEntHole = new G4Tubs("EntHole", 0., 0.80*cm, 0.16*cm, 0.*deg, 360.*deg);
  fLogicEntHole = new G4LogicalVolume(fSolidEntHole, fWorldMaterial, "EntHole");
  fPhysiEntHole = new G4PVPlacement(transform,
        			fLogicEntHole, 
       				"EntHole", 
       				fLogicOutCell, 
       				false, 
    			    	0);   
    			    	
  // Entrance Collimator Thin Metal Disk without hole
  rotm  = G4RotationMatrix(0,0,0);     
  //rotm.rotateY(60*deg);
  position = G4ThreeVector(0.,0.,0.1196*cm);
  //position = G4ThreeVector(0.,0.,0.);
  transform = G4Transform3D(rotm,position); 
  fSolidEntDisk = new G4Tubs("EntDisk", 0., 0.79*cm, 0.039*cm, 0.*deg, 360.*deg);
  fLogicEntDisk = new G4LogicalVolume(fSolidEntDisk, fAlCaseMaterial, "EntDisk");
  fPhysiEntDisk = new G4PVPlacement(transform,
        			fLogicEntDisk, 
       				"EntDisk", 
       				fLogicEntHole, 
       				false, 
    			    	0); 
    			    	  
  // Entrance Collimator 6mm Hole
  rotm  = G4RotationMatrix(0,0,0);    
  rotm.rotateY(30*deg);
  rotm.rotateZ(180*deg);
  position = G4ThreeVector(0.,0.,0.);
  transform = G4Transform3D(rotm,position); 
  G4ThreeVector pLowNorm = G4ThreeVector(1.0,0.,-1.73);
  G4ThreeVector pHighNorm = G4ThreeVector(-1.0,0.,1.73);
  fSolidEntCol = new G4CutTubs("EntCol", 0.*cm, 0.3*cm, 0.046*cm, 0.*deg, 360.*deg,pLowNorm,pHighNorm);
  fLogicEntCol = new G4LogicalVolume(fSolidEntCol, fWorldMaterial, "EntCol");
  fPhysiEntCol = new G4PVPlacement(transform,
				fLogicEntCol, 
       				"EntCol", 
       				fLogicEntDisk, 
       				false, 
    			    	0); 

  // Exit Collimator
  // Exit Collimator Hole
  rotm  = G4RotationMatrix(0,0,0);     
  rotm.rotateY(60*deg);
//  position = G4ThreeVector(-5.306*cm,0.,-2.008*cm);
  position = G4ThreeVector(5.351*cm,0.,-2.087*cm);
  transform = G4Transform3D(rotm,position); 
  fSolidExitHole = new G4Tubs("ExitHole", 0., 0.80*cm, 0.16*cm, 0.*deg, 360.*deg);
  fLogicExitHole = new G4LogicalVolume(fSolidExitHole, fWorldMaterial, "ExitHole");
  fPhysiExitHole = new G4PVPlacement(transform,
        			fLogicExitHole, 
       				"ExitHole", 
       				fLogicOutCell, 
       				false, 
    			    	0);   
    			    	
  // Exit Collimator Thin Metal Disk without hole
  rotm  = G4RotationMatrix(0,0,0);     
  //rotm.rotateY(60*deg);
  position = G4ThreeVector(0.,0.,0.1196*cm);
  //position = G4ThreeVector(0.,0.,0.);
  transform = G4Transform3D(rotm,position); 
  fSolidExitDisk = new G4Tubs("ExitDisk", 0., 0.79*cm, 0.039*cm, 0.*deg, 360.*deg);
  fLogicExitDisk = new G4LogicalVolume(fSolidExitDisk, fAlCaseMaterial, "ExitDisk");
  fPhysiExitDisk = new G4PVPlacement(transform,
        			fLogicExitDisk, 
       				"ExitDisk", 
       				fLogicExitHole, 
       				false, 
    			    	0); 
    			    	  
  // Exit Collimator 6mm Hole
  rotm  = G4RotationMatrix(0,0,0);     
  rotm.rotateY(30*deg);
  position = G4ThreeVector(0.,0.,0.);
  transform = G4Transform3D(rotm,position); 
  pLowNorm = G4ThreeVector(1.0,0.,-1.73);
  pHighNorm = G4ThreeVector(-1.0,0.,1.73);
  fSolidExitCol = new G4CutTubs("ExitCol", 0.*cm, 0.3*cm, 0.046*cm, 0.*deg, 360.*deg,pLowNorm,pHighNorm);
  fLogicExitCol = new G4LogicalVolume(fSolidExitCol, fWorldMaterial, "ExitCol");
  fPhysiExitCol = new G4PVPlacement(transform,
				fLogicExitCol, 
       				"ExitCol", 
       				fLogicExitDisk, 
       				false, 
    			    	0); 
	
}

void DetectorConstruction::ConstructDetectorComponents()
{
  //Detector Holder Volume for Components 

  fSolidDetector = new G4Tubs("Detector", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicDetector = new G4LogicalVolume(fSolidDetector,fWorldMaterial, "Detector");
  fLogicDetector->SetVisAttributes(fWireFrameVisAtt);

  //Detector Crystal
  //
  fSolidDetectorCrystal = new G4Tubs("DetectorCrystal", 0., 0.5*fCrystalDiameter, 0.5*fCrystalLength, 0.*deg, 360.*deg);
  fLogicDetectorCrystal = new G4LogicalVolume(fSolidDetectorCrystal, fCrystalMaterial, "DetectorCrystal");
  fPhysiDetectorCrystal = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposCrystal), 
        			fLogicDetectorCrystal, 
       				"DetectorCrystal", 
       				fLogicDetector, 
       				false, 
    			    	0);
  fLogicDetectorCrystal->SetVisAttributes(fRedVisAtt);

  //Detector Gap (gap surrounding crystal and casing or a reflector as required)
  //
  fSolidGap = new G4Tubs("Gap", 0.5*fCrystalDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap = new G4LogicalVolume(fSolidGap, fGapMaterial, "Gap");
  fPhysiGap = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap, 
       				"Gap", 
       				fLogicDetector, 
       				false, 
    			    	0);

  //Detector FaceGap (gap between face of crystal and casing or a reflector as required)
  //
  fSolidFaceGap = new G4Tubs("FaceGap", 0., 0.5*fCrystalDiameter, 0.5*fGapThickness, 0.*deg, 360.*deg);
  fLogicFaceGap = new G4LogicalVolume(fSolidFaceGap, fGapMaterial, "FaceGap");
  fPhysiFaceGap = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap, 
       				"FaceGap", 
       				fLogicDetector, 
       				false, 
    			    	0);

  //fLogicGap->SetVisAttributes(fAuxEdgeVisAtt);
  //fLogicFaceGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap->SetVisAttributes(fBlueVisAtt);
  fLogicFaceGap->SetVisAttributes(fBlueVisAtt);

  //Detector Aluminum Casing (casing surrounding crystal)
  //
  fSolidAlCase = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase = new G4LogicalVolume(fSolidAlCase, fAlCaseMaterial, "AlCase");
  fPhysiAlCase = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase, 
       				"AlCase", 
       				fLogicDetector, 
       				false, 
    			    	0);

  //Detector Aluminum Face Casing (casing covering face of crystal)
  //
  fSolidFaceAlCase = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase = new G4LogicalVolume(fSolidFaceAlCase, fAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase, 
       				"FaceAlCase", 
       				fLogicDetector, 
       				false, 
    			    	0);
  //Lead casing surrounding aluminum casing (around casing surround)
  //
  //fSolidPbCase = new G4Tubs("PbCase", 0.5*fAlCaseDiameter, 0.5*fPbCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  //fLogicPbCase = new G4LogicalVolume(fSolidPbCase, fPbCaseMaterial, "PbCase");
  //fPhysiPbCase = new G4PVPlacement(0, 
        			//G4ThreeVector(0.,0.,fZposAlCase), 
        			//fLogicPbCase, 
       				//"PbCase", 
       				//fLogicDetector, 
       				//false, 
    			    	//0);

  //Lead Casing Ring (at back in front of PMT face)
  //
  //Check to see if collar is needed first. Large enough crystal and small PMT will not require collar.
  //if (fPbCaseDiameter<fPMTDiameter){  
	//fSolidPbCollar = new G4Tubs("PbCollar", 0.5*fPbCaseDiameter, 0.5*fPMTDiameter, 0.5*fPbCaseThickness, 0.*deg, 360.*deg);
	//fLogicPbCollar = new G4LogicalVolume(fSolidPbCollar, fPbCaseMaterial, "PbCollar");
	//fPhysiPbCollar = new G4PVPlacement(0, 
						//G4ThreeVector(0.,0.,fZposPbCollar), 
						//fLogicPbCollar, 
						//"PbCollar", 
						//fLogicDetector, 
						//false, 
    			    	//0);
	//fLogicPbCase->SetVisAttributes(fGreyVisAtt);
	//fLogicPbCollar->SetVisAttributes(fGreyVisAtt);
  //}
  //Detector PMT
  //
  fSolidPMT = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT = new G4LogicalVolume(fSolidPMT, fPMTMaterial,"PMT");
  fPhysiPMT = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT, 
       				"PMT", 
       				fLogicDetector, 
       				false, 
    			    	0);
  //PMT Optical Window
  //
  fSolidPMTWin = new G4Tubs("PMTWin", 0., 0.5*fCrystalDiameter, 0.5*fPMTWinThickness, 0.*deg, 360.*deg);
  fLogicPMTWin = new G4LogicalVolume(fSolidPMTWin, fPMTWinMaterial, "PMTWin");
  fPhysiPMTWin = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMTWin), 
        			fLogicPMTWin, 
       				"PMTWin", 
       				fLogicDetector, 
       				false, 
    			    	0);
  fLogicPMTWin->SetVisAttributes(fGreenVisAtt);
}

//void DetectorConstruction::ConstructSDandField()
//{
//}

void DetectorConstruction::PrintGeometryParameters()
{
  G4cout << "\n" << fWorldMaterial    << G4endl;
  G4cout << "\n" << "World Size X is: " << fWorldSizeX    << G4endl;
  G4cout << "\n" << "World Size YZ is: " << fWorldSizeYZ    << G4endl;  ;
  G4cout << "\n" << "Detector Diameter is: " << fDetectorDiameter    << G4endl;
  G4cout << "\n" << "Detector Length is: " << fDetectorLength    << G4endl;
  G4cout << "\n" << "Gap Thickness is: " << fGapThickness    << G4endl;
  G4cout << "\n" << "Gap  Material is: " << fGapMaterial    << G4endl;
  G4cout << "\n" << "Al Thickness is: " << fAlCaseThickness    << G4endl;
  G4cout << "\n" << "Pb Thickness is: " << fPbCaseThickness    << G4endl;
  G4cout << "\n" << "PMT Diameter is: " << fPMTDiameter    << G4endl;
  G4cout << "\n" << "PMT Length is: " << fPMTLength    << G4endl;
}
/*
void DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fWorldMaterial != pttoMaterial) {
    fWorldMaterial = pttoMaterial;                  
    if(fLogicWorld) fLogicWorld->SetMaterial(fWorldMaterial);
  }
}
*/    

void DetectorConstruction::SetCrystalDiameter(G4double val)
{
  fCrystalDiameter = val;
}

void DetectorConstruction::SetCrystalLength(G4double val)
{
  fCrystalLength = val;
}

void DetectorConstruction::SetGapThickness(G4double val)
{
  fGapThickness = val;
}

void DetectorConstruction::SetGapMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fGapMaterial != pttoMaterial) {
    fGapMaterial = pttoMaterial;                  
    if(fLogicGap) fLogicGap->SetMaterial(fGapMaterial);
  }
}

void DetectorConstruction::SetAlCaseThickness(G4double val)
{
  fAlCaseThickness = val;
}
void DetectorConstruction::SetAlFaceThickness(G4double val)
{
  fAlFaceThickness = val;
}
void DetectorConstruction::SetPbCaseThickness(G4double val)
{
  fPbCaseThickness = val;
}
void DetectorConstruction::SetPMTDiameter(G4double val)
{
  fPMTDiameter = val;
}
void DetectorConstruction::SetPMTLength(G4double val)
{
  fPMTLength = val;
}
void DetectorConstruction::SetNoOfDetectors(G4int val)
{
  fNoOfDets = val;
}
void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructGeometry());
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
//  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

