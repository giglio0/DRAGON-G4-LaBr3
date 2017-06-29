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
// * institutes,nor the agencies providing finanial support for this *
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
/// Si_Ion_Chamber_v4																																																											/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4TransportationManager.hh" 
#include "G4Navigator.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: fSolidWorld(0),fLogicWorld(0),fPhysiWorld(0),
fSolidDetector(0),fLogicDetector(0),fPhysiDetector(0),
fSolidCrystal(0),fLogicCrystal(0),fPhysiCrystal(0),
fSolidGap(0),fLogicGap(0),fPhysiGap(0),
fSolidFaceGap(0),fLogicFaceGap(0),fPhysiFaceGap(0),
fSolidAlCase(0),fLogicAlCase(0),fPhysiAlCase(0),
fSolidFaceAlCase(0),fLogicFaceAlCase(0),fPhysiFaceAlCase(0),
fSolidPbCollar(0),fLogicPbCollar(0),fPhysiPbCollar(0),
fSolidPMT(0),fLogicPMT(0),fPhysiPMT(0),
fSolidPMTWin(0),fLogicPMTWin(0),fPhysiPMTWin(0)
{
  //Setting the default for the world for the Sim
  fWorldSizeX = 100*cm;
  fWorldSizeYZ = 100*cm;

  //Set Detector Defaults
  fPMTDiameter = 7.5*cm;
  fPMTLength = 22.0*cm;
  fPMTWinThickness = 0.508*cm;
  fDetectorLength = 5.08*cm;
  fDetectorDiameter = 5.08*cm;
  fGapThickness = 0.127*cm;
  fAlCaseThickness = 0.2*cm;
  fPbCaseThickness = 0.5*cm;
  
  ComputeCalorParameters();
  
  // materials  
  DefineMaterials();
  // Default Materials
  fWorldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Galactic");
  fDetectorMaterial = G4NistManager::Instance()->FindOrBuildMaterial("LaBr3");
  fBoxMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Plastic");
  fCrystalMaterial = G4NistManager::Instance()->FindOrBuildMaterial("LaBr3");
  fAlCaseMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  fFaceAlCaseMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  fPbCaseMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
  fPbCollarMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
  fPMTMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pyrex_Glass");
  fPMTWinMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pyrex_Glass");
  //Set Default Gap Material i.e. no reflector
  fGapMaterial 	= G4NistManager::Instance()->FindOrBuildMaterial("Galactic");
  fFaceGapMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Galactic");

  //G4MaterialPropertyVector* AddProperty('Energy Resolution', [1*MeV, 2*MeV, 3*MeV, 4*MeV, 5*MeV, 6*MeV, 7*MeV, 8*MeV, 9*MeV, 10*MeV],
  //[1000*eV, 1414*eV, 1732*eV, 2000*eV, 2236*eV, 2449*eV, 2646*eV, 2828*eV, 3000*eV, 3162*eV], 10);
  //AddProperty('Energy Resolution', G4MaterialPropertyVector*);
  //G4MaterialPropertiesTable* fResolutionDetectorMaterial = new G4MaterialPropertiesTable();
  //fResolutionDetectorMaterial->AddConstProperty("RESOLUTIONSCALE", 1.0);
  //fDetectorMaterial->SetMaterialPropertiesTable(fResolutionDetectorMaterial);

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
  return ConstructCalorimeter();
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
  G4Element* Bi = new G4Element("Bismuth", symbol="Bi", z=83, a= 208.980*g/mole);
  G4Element* Ge = new G4Element("Germanium", symbol="Ge", z=32, a= 72.63*g/mole);

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
  G4Material* LaBr3 = new G4Material("LaBr3", density= 5.06*g/cm3, ncomponents=2);
  LaBr3->AddElement(La, fractionmass = 0.366875);
  LaBr3->AddElement(Br,fractionmass = 0.633125);
  LaBr3->GetIonisation()->SetMeanExcitationEnergy(454.5*eV);

  G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);
  
  G4Material* BGO = new G4Material("BGO", density = 7.13*g/cm3, ncomponents=3);
  BGO->AddElement(Bi, fractionmass=0.67099);
  BGO->AddElement(Ge, fractionmass=0.17490);
  BGO->AddElement(O, fractionmass=0.15411);
  BGO->GetIonisation()->SetMeanExcitationEnergy(534.1*eV);

  // example of vacuum
  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  new G4Material("Galactic", z=1, a=1.01*g/mole,density,
                 kStateGas,temperature,pressure);
}

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
  fGapLength = fDetectorLength + fPMTWinThickness + fGapThickness;
  fGapDiameter = fDetectorDiameter + 2*fGapThickness;
  fAlCaseLength = fGapLength + fAlCaseThickness;
  fAlCaseDiameter = fGapDiameter + 2*fAlCaseThickness;
  fPbCaseDiameter = fAlCaseDiameter + 2*fPbCaseThickness;
  fTotalDetectorLength = fPMTWinThickness + fDetectorLength + fPMTLength + fGapThickness + fAlCaseThickness;
  if (fPbCaseDiameter>fPMTDiameter){
	  fPMTDiameter = fPbCaseDiameter;
	  G4cout << "\n" << "\n" 
	  << "WARNING: PMT Diameter has been adjusted to compensate for large crystal/lead/aluminum size"  
	  << G4endl;
	  fTotalDetectorDiameter = fPMTDiameter;
  }
  else
	fTotalDetectorDiameter = fPMTDiameter;

  fZposPbCollar = 0.5*fTotalDetectorLength - fAlCaseThickness - fGapThickness - fDetectorLength - fPMTWinThickness +  
  0.5*fPbCaseThickness;
  fZposFaceAlCase = 0.5*fTotalDetectorLength - 0.5*fAlCaseThickness;
  fZposFaceGap = 0.5*fTotalDetectorLength - fAlCaseThickness - 0.5*fGapThickness ;
  fZpos = 0.5*fTotalDetectorLength - fGapThickness - fAlCaseThickness - 0.5*fDetectorLength;

  fZposAlCase = 0.5*fTotalDetectorLength - 0.5*fAlCaseLength;
  fZposGap = 0.5*fTotalDetectorLength - fAlCaseThickness - 0.5*fGapLength ;
  fZposPMTWin = 0.5*fTotalDetectorLength - fAlCaseThickness - fGapThickness - fDetectorLength - 0.5*fPMTWinThickness;
  fZposPMT = 0.5*fTotalDetectorLength - fAlCaseThickness - fGapThickness - fDetectorLength - fPMTWinThickness - 
  0.5*fPMTLength;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{ 
  // Cleanup old geometry
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
 
  // complete the Calor parameters definition 
  ComputeCalorParameters();

  // World
  //
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
  //Detector Holder Volume for Components
  //

if (fDetectorGeometry == 1){  
//  If rotation required  
//  G4RotationMatrix rotm  = G4RotationMatrix(0,90*deg,-90*deg);     
//  No Rotation Now
  G4RotationMatrix rotm  = G4RotationMatrix(0,0,0);     
  G4ThreeVector position = G4ThreeVector(0.,0.,0.);
  G4Transform3D transform = G4Transform3D(rotm,position);
  G4ThreeVector boxposition = G4ThreeVector(0.,0.,fTotalDetectorLength);
  G4Transform3D boxtransform = G4Transform3D(rotm,boxposition);
  
  G4Box *outerBox = new G4Box("Outer Box", 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength);
  G4Box *innerBox = new G4Box("Inner Box", 0.5*fTotalDetectorDiameter - 0.5*cm, 0.5*fTotalDetectorDiameter - 0.5*cm, 0.5*fTotalDetectorLength);
  fSolidBox = new G4SubtractionSolid("Box",outerBox,innerBox); 
  fLogicBox = new G4LogicalVolume(fSolidBox, fBoxMaterial, "Box");
  fPhysiBox = new G4PVPlacement(boxtransform,
        			fLogicBox, 
       				"Box", 
       				fLogicWorld, 
       				false, 
    			    	0,
				    false); 

  fSolidDetector = new G4Tubs("Detector", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fLogicDetector = new G4LogicalVolume(fSolidDetector,fWorldMaterial, "Detector");
  fPhysiDetector = new G4PVPlacement(transform,
        			fLogicDetector, 
       				"Detector", 
       				fLogicWorld, 
       				false, 
    			    	0,
				    false);

  fLogicDetector->SetVisAttributes(fWireFrameVisAtt);

  // PMT Optical Window
  //
  fSolidPMTWin = new G4Tubs("PMTWin", 0., 0.5*fDetectorDiameter, 0.5*fPMTWinThickness, 0.*deg, 360.*deg);
  fLogicPMTWin = new G4LogicalVolume(fSolidPMTWin, fPMTWinMaterial, "PMTWin");
  fPhysiPMTWin = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMTWin), 
        			fLogicPMTWin, 
       				"PMTWin", 
       				fLogicDetector, 
       				false, 
    			    	0);
  fLogicPMTWin->SetVisAttributes(fGreenVisAtt);

  // Crystal
  //
  fSolidCrystal = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal = new G4LogicalVolume(fSolidCrystal, fDetectorMaterial, "Crystal");
  fPhysiCrystal = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal, 
       				"Crystal", 
       				fLogicDetector, 
       				false, 
    			    	0);
  fLogicCrystal->SetVisAttributes(fRedVisAtt);

  //Gap (gap surrounding crystal and casing or a reflector as required)
  //
  fSolidGap = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap = new G4LogicalVolume(fSolidGap, fGapMaterial, "Gap");
  fPhysiGap = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap, 
       				"Gap", 
       				fLogicDetector, 
       				false, 
    			    	0);

  //FaceGap (gap between face of crystal and casing or a reflector as required)
  //
  fSolidFaceGap = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapThickness, 0.*deg, 360.*deg);
  fLogicFaceGap = new G4LogicalVolume(fSolidFaceGap, fGapMaterial, "FaceGap");
  fPhysiFaceGap = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap, 
       				"FaceGap", 
       				fLogicDetector, 
       				false, 
    			    	0);

  fLogicGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap->SetVisAttributes(fBlueVisAtt);
  fLogicFaceGap->SetVisAttributes(fBlueVisAtt);

  //Aluminum Casing (casing surrounding crystal)
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

  //Aluminum Face Casing (casing covering face of crystal)
  //
  fSolidFaceAlCase = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlCaseThickness, 0.*deg, 360.*deg);
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
  fSolidPbCase = new G4Tubs("PbCase", 0.5*fAlCaseDiameter, 0.5*fPbCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicPbCase = new G4LogicalVolume(fSolidPbCase, fPbCaseMaterial, "PbCase");
  fPhysiPbCase = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicPbCase, 
       				"PbCase", 
       				fLogicDetector, 
       				false, 
    			    	0);

  //Lead Casing Ring (at back in front of PMT face)
  //
  //Check to see if collar is needed first. Large enough crystal and small PMT will not require collar.
  if (fPbCaseDiameter<fPMTDiameter){  
	fSolidPbCollar = new G4Tubs("PbCollar", 0.5*fPbCaseDiameter, 0.5*fPMTDiameter, 0.5*fPbCaseThickness, 0.*deg,    
        360.*deg);
	fLogicPbCollar = new G4LogicalVolume(fSolidPbCollar, fPbCaseMaterial, "PbCollar");
	fPhysiPbCollar = new G4PVPlacement(0, 
						G4ThreeVector(0.,0.,fZposPbCollar), 
						fLogicPbCollar, 
						"PbCollar", 
						fLogicDetector, 
						false, 
    			    	0);
	fLogicPbCase->SetVisAttributes(fGreyVisAtt);
	fLogicPbCollar->SetVisAttributes(fGreyVisAtt);
  }
  //PMT
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
  PrintCalorParameters();
}

if (fDetectorGeometry == 2){
//  If rotation required  
//  G4RotationMatrix rotm  = G4RotationMatrix(0,90*deg,-90*deg);     
//  No Rotation Now
  G4RotationMatrix rotm  = G4RotationMatrix(0,0,0);     
  G4ThreeVector position = G4ThreeVector(0.,0.,0.);
  G4Transform3D transform = G4Transform3D(rotm,position);
  G4ThreeVector boxposition = G4ThreeVector(0.,0.,fTotalDetectorLength);
  G4Transform3D boxtransform = G4Transform3D(rotm,boxposition);
  const G4double zPlane[2] = {-0.5*fTotalDetectorLength, 0.5*fTotalDetectorLength};
  const G4double zPlane1[2] = {-0.5*fPMTWinThickness, 0.5*fPMTWinThickness};
  const G4double zPlane2[2] = {-0.5*fDetectorLength, 0.5*fDetectorLength};
  const G4double zPlane3[2] = {-0.5*fGapLength, 0.5*fGapLength};
  const G4double zPlane4[2] = {-0.5*fGapThickness, 0.5*fGapThickness};
  const G4double zPlane5[2] = {-0.5*fAlCaseLength, 0.5*fAlCaseLength};
  const G4double zPlane6[2] = {-0.5*fAlCaseThickness, 0.5*fAlCaseThickness};
  const G4double zPlane7[2] = {-0.5*fAlCaseLength, 0.5*fAlCaseLength};
  const G4double zPlane8[2] = {-0.5*fPbCaseThickness, 0.5*fPbCaseThickness};
  const G4double zPlane9[2] = {-0.5*fPMTLength, 0.5*fPMTLength};
  const G4double rInner[2] = {0,0};
  const G4double rInner1[2] = {0.5*fDetectorDiameter, 0.5*fDetectorDiameter};
  const G4double rInner2[2] = {0.5*fGapDiameter, 0.5*fGapDiameter};
  const G4double rInner3[2] = {0.5*fAlCaseDiameter, 0.5*fAlCaseDiameter};
  const G4double rInner4[2] = {0.5*fPbCaseDiameter,0.5*fPbCaseDiameter};
  const G4double rOuter[2] = {0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorDiameter};
  const G4double rOuter1[2] = {0.5*fGapDiameter, 0.5*fGapDiameter};
  const G4double rOuter2[2] = {0.5*fAlCaseDiameter, 0.5*fAlCaseDiameter};
  const G4double rOuter3[2] = {0.5*fPbCaseDiameter, 0.5*fPbCaseDiameter};
  const G4double rOuter4[2] = {0.5*fPMTDiameter, 0.5*fPMTDiameter};
  
  G4Box *outerBox = new G4Box("Outer Box", 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength);
  G4Box *innerBox = new G4Box("Inner Box", 0.5*fTotalDetectorDiameter - 0.5*cm, 0.5*fTotalDetectorDiameter - 0.5*cm, 0.5*fTotalDetectorLength);
  fSolidBox = new G4SubtractionSolid("Box",outerBox,innerBox); 
  fLogicBox = new G4LogicalVolume(fSolidBox, fBoxMaterial, "Box");
  fPhysiBox = new G4PVPlacement(boxtransform,
        			fLogicBox, 
       				"Box", 
       				fLogicWorld, 
       				false, 
    			    	0,
				    false); 
  
  fSolidDetector1 = new G4Polyhedra{"Detector", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fLogicDetector = new G4LogicalVolume(fSolidDetector1,fWorldMaterial, "Detector");
  fPhysiDetector = new G4PVPlacement(transform,
        			fLogicDetector, 
       				"Detector", 
       				fLogicWorld, 
       				false, 
    			    	0,
				    false);

  fLogicDetector->SetVisAttributes(fWireFrameVisAtt);

  // PMT Optical Window
  //
  fSolidPMTWin1 = new G4Polyhedra("PMTWin", 0.*deg, 360.*deg, 6, 2, zPlane1, rInner, rInner1);
  fLogicPMTWin = new G4LogicalVolume(fSolidPMTWin1, fPMTWinMaterial, "PMTWin");
  fPhysiPMTWin = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMTWin), 
        			fLogicPMTWin, 
       				"PMTWin", 
       				fLogicDetector, 
       				false, 
    			    	0);
  fLogicPMTWin->SetVisAttributes(fGreenVisAtt);

  // Crystal
  //
  fSolidCrystal1 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal = new G4LogicalVolume(fSolidCrystal1, fDetectorMaterial, "Crystal");
  fPhysiCrystal = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal, 
       				"Crystal", 
       				fLogicDetector, 
       				false, 
    			    	0);
  fLogicCrystal->SetVisAttributes(fRedVisAtt);

  //Gap (gap surrounding crystal and casing or a reflector as required)
  //
  fSolidGap1 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap = new G4LogicalVolume(fSolidGap1, fGapMaterial, "Gap");
  fPhysiGap = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap, 
       				"Gap", 
       				fLogicDetector, 
       				false, 
    			    	0);

  //FaceGap (gap between face of crystal and casing or a reflector as required)
  //
  fSolidFaceGap1 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane4, rInner, rInner1);
  fLogicFaceGap = new G4LogicalVolume(fSolidFaceGap1, fGapMaterial, "FaceGap");
  fPhysiFaceGap = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap, 
       				"FaceGap", 
       				fLogicDetector, 
       				false, 
    			    	0);

  fLogicGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap->SetVisAttributes(fBlueVisAtt);
  fLogicFaceGap->SetVisAttributes(fBlueVisAtt);

  //Aluminum Casing (casing surrounding crystal)
  //
  fSolidAlCase1 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase = new G4LogicalVolume(fSolidAlCase1, fAlCaseMaterial, "AlCase");
  fPhysiAlCase = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase, 
       				"AlCase", 
       				fLogicDetector, 
       				false, 
    			    	0);

  //Aluminum Face Casing (casing covering face of crystal)
  //
  fSolidFaceAlCase1 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase = new G4LogicalVolume(fSolidFaceAlCase1, fAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase, 
       				"FaceAlCase", 
       				fLogicDetector, 
       				false, 
    			    	0);
  //Lead casing surrounding aluminum casing (around casing surround)
  //
  fSolidPbCase1 = new G4Polyhedra("PbCase", 0.*deg, 360.*deg, 6, 2, zPlane7, rInner3, rOuter3);
  fLogicPbCase = new G4LogicalVolume(fSolidPbCase1, fPbCaseMaterial, "PbCase");
  fPhysiPbCase = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicPbCase, 
       				"PbCase", 
       				fLogicDetector, 
       				false, 
    			    	0);

  //Lead Casing Ring (at back in front of PMT face)
  //
  //Check to see if collar is needed first. Large enough crystal and small PMT will not require collar.
  if (fPbCaseDiameter<fPMTDiameter){  
	fSolidPbCollar1 = new G4Polyhedra("PbCollar", 0.*deg, 360.*deg, 6, 2, zPlane8, rInner4, rOuter4);
	fLogicPbCollar = new G4LogicalVolume(fSolidPbCollar1, fPbCaseMaterial, "PbCollar");
	fPhysiPbCollar = new G4PVPlacement(0, 
						G4ThreeVector(0.,0.,fZposPbCollar), 
						fLogicPbCollar, 
						"PbCollar", 
						fLogicDetector, 
						false, 
    			    	0);
	fLogicPbCase->SetVisAttributes(fGreyVisAtt);
	fLogicPbCollar->SetVisAttributes(fGreyVisAtt);
  }
  //PMT
  //
  fSolidPMT1 = new G4Polyhedra("PMT", 0.*deg, 360.*deg, 6, 2, zPlane9, rInner, rOuter4);
  fLogicPMT = new G4LogicalVolume(fSolidPMT1, fPMTMaterial,"PMT");
  fPhysiPMT = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT, 
       				"PMT", 
       				fLogicDetector, 
       				false, 
    			    	0);
  PrintCalorParameters();

  //always return the physical World
}  

if (fDetectorGeometry == 3){  
//  If rotation required  
//  G4RotationMatrix rotm  = G4RotationMatrix(0,90*deg,-90*deg);     
//  No Rotation Now
  
  fSolidDetector = new G4Tubs("Detector", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fLogicDetector = new G4LogicalVolume(fSolidDetector,fWorldMaterial, "Detector");
  G4RotationMatrix rotm = G4RotationMatrix(0,0,0); 
  G4RotationMatrix rotm180 = G4RotationMatrix(0,180.*deg,0);    
  G4ThreeVector P;
  G4Transform3D Tr;
  
  G4AssemblyVolume* assemblyDetector = new G4AssemblyVolume();
  //10 detectors with 10 straddling detectors
  //bottom 3
  P.setX(-1.5*fTotalDetectorDiameter); P.setY((-0.875*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(-0.5*fTotalDetectorDiameter); P.setY((-0.875*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(0.5*fTotalDetectorDiameter); P.setY((-0.875*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  //4
  P.setX(-2*fTotalDetectorDiameter); P.setY((-0.375*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(-fTotalDetectorDiameter); P.setY((-0.375*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(0.); P.setY((-0.375*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(fTotalDetectorDiameter); P.setY((-0.375*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  //3
  P.setX(-1.5*fTotalDetectorDiameter); P.setY((0.125*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(-0.5*fTotalDetectorDiameter); P.setY((0.125*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(0.5*fTotalDetectorDiameter); P.setY((0.125*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  //10 straddling detectors
  //top 4
  P.setX(0.5*sqrt(3)*fTotalDetectorDiameter); P.setY((0.625*sqrt(3) + 1/(16*sqrt(3)) + 0.5)*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(0.); P.setY((0.625*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(-fTotalDetectorDiameter); P.setY((0.625*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(0.5*(-2 - sqrt(3))*fTotalDetectorDiameter); P.setY((0.625*sqrt(3) + 1/(16*sqrt(3)) + 0.5)*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  //below the top 4
  P.setX((0.5*sqrt(3) + 0.5)*fTotalDetectorDiameter); P.setY((0.125*sqrt(3) + 1/(16*sqrt(3)) + 0.5)*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(0.5*(-3 - sqrt(3))*fTotalDetectorDiameter); P.setY((0.125*sqrt(3) + 1/(16*sqrt(3)) + 0.5)*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  //4
  P.setX((1 + 0.5*sqrt(3))*fTotalDetectorDiameter); P.setY((-0.375*sqrt(3) + 1/(16*sqrt(3)) - 0.5)*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((0.5 + 0.5*sqrt(3))*fTotalDetectorDiameter); P.setY((-0.875*sqrt(3) + 1/(16*sqrt(3)) - 0.5)*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-2 - 0.5*sqrt(3))*fTotalDetectorDiameter); P.setY((-0.375*sqrt(3) + 1/(16*sqrt(3)) - 0.5)*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-1.5 - 0.5*sqrt(3))*fTotalDetectorDiameter); P.setY((-0.875*sqrt(3) + 1/(16*sqrt(3)) - 0.5)*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  //end straddling detectors
  //10 detectors
  //bottom 3
  P.setX(-1.5*fTotalDetectorDiameter); P.setY((-0.875*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(-0.5*fTotalDetectorDiameter); P.setY((-0.875*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(0.5*fTotalDetectorDiameter); P.setY((-0.875*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  //4
  P.setX(-2*fTotalDetectorDiameter); P.setY((-0.375*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(-fTotalDetectorDiameter); P.setY((-0.375*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(0.); P.setY((-0.375*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(fTotalDetectorDiameter); P.setY((-0.375*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  //3
  P.setX(-1.5*fTotalDetectorDiameter); P.setY((0.125*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(-0.5*fTotalDetectorDiameter); P.setY((0.125*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(0.5*fTotalDetectorDiameter); P.setY((0.125*sqrt(3) + 1/(16*sqrt(3)))*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);

  assemblyDetector->MakeImprint(fLogicWorld, Tr);
  fLogicDetector->SetVisAttributes(fWireFrameVisAtt);

  // PMT Optical Window
  //
  fSolidPMTWin = new G4Tubs("PMTWin", 0., 0.5*fDetectorDiameter, 0.5*fPMTWinThickness, 0.*deg, 360.*deg);
  fLogicPMTWin = new G4LogicalVolume(fSolidPMTWin, fPMTWinMaterial, "PMTWin");
  fPhysiPMTWin = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMTWin), 
        			fLogicPMTWin,
       				"PMTWin", 
       				fLogicDetector, 
       				false, 
    			    	0);
  fLogicPMTWin->SetVisAttributes(fGreenVisAtt);

  // Crystal
  //
  fSolidCrystal = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal = new G4LogicalVolume(fSolidCrystal, fDetectorMaterial, "Crystal");
  fPhysiCrystal = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal, 
       				"Crystal", 
       				fLogicDetector, 
       				false, 
    			    	0);
  fLogicCrystal->SetVisAttributes(fRedVisAtt);

  //Gap (gap surrounding crystal and casing or a reflector as required)
  //
  fSolidGap = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap = new G4LogicalVolume(fSolidGap, fGapMaterial, "Gap");
  fPhysiGap = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap, 
       				"Gap", 
       				fLogicDetector, 
       				false, 
    			    	0);

  //FaceGap (gap between face of crystal and casing or a reflector as required)
  //
  fSolidFaceGap = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapThickness, 0.*deg, 360.*deg);
  fLogicFaceGap = new G4LogicalVolume(fSolidFaceGap, fGapMaterial, "FaceGap");
  fPhysiFaceGap = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap, 
       				"FaceGap", 
       				fLogicDetector, 
       				false, 
    			    	0);

  fLogicGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap->SetVisAttributes(fBlueVisAtt);
  fLogicFaceGap->SetVisAttributes(fBlueVisAtt);

  //Aluminum Casing (casing surrounding crystal)
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

  //Aluminum Face Casing (casing covering face of crystal)
  //
  fSolidFaceAlCase = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlCaseThickness, 0.*deg, 360.*deg);
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
  fSolidPbCase = new G4Tubs("PbCase", 0.5*fAlCaseDiameter, 0.5*fPbCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicPbCase = new G4LogicalVolume(fSolidPbCase, fPbCaseMaterial, "PbCase");
  fPhysiPbCase = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicPbCase, 
       				"PbCase", 
       				fLogicDetector, 
       				false, 
    			    	0);

  //Lead Casing Ring (at back in front of PMT face)
  //
  //Check to see if collar is needed first. Large enough crystal and small PMT will not require collar.
  if (fPbCaseDiameter<fPMTDiameter){  
	fSolidPbCollar = new G4Tubs("PbCollar", 0.5*fPbCaseDiameter, 0.5*fPMTDiameter, 0.5*fPbCaseThickness, 0.*deg,    
        360.*deg);
	fLogicPbCollar = new G4LogicalVolume(fSolidPbCollar, fPbCaseMaterial, "PbCollar");
	fPhysiPbCollar = new G4PVPlacement(0, 
						G4ThreeVector(0.,0.,fZposPbCollar), 
						fLogicPbCollar, 
						"PbCollar", 
						fLogicDetector, 
						false, 
    			    	0);
	fLogicPbCase->SetVisAttributes(fGreyVisAtt);
	fLogicPbCollar->SetVisAttributes(fGreyVisAtt);
  }
  //PMT
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
  PrintCalorParameters();
}

if (fDetectorGeometry == 4){
//  If rotation required  
//  G4RotationMatrix rotm  = G4RotationMatrix(0,90*deg,-90*deg);     
//  No Rotation Now
  const G4double zPlane[2] = {-0.5*fTotalDetectorLength, 0.5*fTotalDetectorLength};
  const G4double zPlane1[2] = {-0.5*fPMTWinThickness, 0.5*fPMTWinThickness};
  const G4double zPlane2[2] = {-0.5*fDetectorLength, 0.5*fDetectorLength};
  const G4double zPlane3[2] = {-0.5*fGapLength, 0.5*fGapLength};
  const G4double zPlane4[2] = {-0.5*fGapThickness, 0.5*fGapThickness};
  const G4double zPlane5[2] = {-0.5*fAlCaseLength, 0.5*fAlCaseLength};
  const G4double zPlane6[2] = {-0.5*fAlCaseThickness, 0.5*fAlCaseThickness};
  const G4double zPlane7[2] = {-0.5*fAlCaseLength, 0.5*fAlCaseLength};
  const G4double zPlane8[2] = {-0.5*fPbCaseThickness, 0.5*fPbCaseThickness};
  const G4double zPlane9[2] = {-0.5*fPMTLength, 0.5*fPMTLength};
  const G4double rInner[2] = {0,0};
  const G4double rInner1[2] = {0.5*fDetectorDiameter, 0.5*fDetectorDiameter};
  const G4double rInner2[2] = {0.5*fGapDiameter, 0.5*fGapDiameter};
  const G4double rInner3[2] = {0.5*fAlCaseDiameter, 0.5*fAlCaseDiameter};
  const G4double rInner4[2] = {0.5*fPbCaseDiameter,0.5*fPbCaseDiameter};
  const G4double rOuter[2] = {0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorDiameter};
  const G4double rOuter1[2] = {0.5*fGapDiameter, 0.5*fGapDiameter};
  const G4double rOuter2[2] = {0.5*fAlCaseDiameter, 0.5*fAlCaseDiameter};
  const G4double rOuter3[2] = {0.5*fPbCaseDiameter, 0.5*fPbCaseDiameter};
  const G4double rOuter4[2] = {0.5*fPMTDiameter, 0.5*fPMTDiameter};
  
  fSolidDetector1 = new G4Polyhedra{"Detector", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fLogicDetector = new G4LogicalVolume(fSolidDetector1,fWorldMaterial, "Detector");
  G4RotationMatrix rotm = G4RotationMatrix(0,0,0); 
  G4RotationMatrix rotm180 = G4RotationMatrix(0,180.*deg,0); 
  //G4RotationMatrix rotm90 = G4RotationMatrix (0,180.*deg,90.*deg);  
  G4ThreeVector P;
  G4Transform3D Tr;
  
  G4AssemblyVolume* assemblyDetector = new G4AssemblyVolume();
  //10 detectors with 10 straddling detectors
  //leftmost 3
  P.setX(-1.5*fTotalDetectorDiameter); P.setY(-1.5*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(-1.5*fTotalDetectorDiameter); P.setY(-0.5*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(-1.5*fTotalDetectorDiameter); P.setY(0.5*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  //4
  P.setX((-1.5 + 0.5*sqrt(3))*fTotalDetectorDiameter); P.setY(-2*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-1.5 + 0.5*sqrt(3))*fTotalDetectorDiameter); P.setY(-fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-1.5 + 0.5*sqrt(3))*fTotalDetectorDiameter); P.setY(0.); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-1.5 + 0.5*sqrt(3))*fTotalDetectorDiameter); P.setY(fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  //3
  P.setX((-1.5 + sqrt(3))*fTotalDetectorDiameter); P.setY(-1.5*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-1.5 + sqrt(3))*fTotalDetectorDiameter); P.setY(-0.5*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-1.5 + sqrt(3))*fTotalDetectorDiameter); P.setY(0.5*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter + 0.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //10 straddling detectors
  //rightmost 4
  P.setX((-1 + 1.5*sqrt(3))*fTotalDetectorDiameter); P.setY(-2*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-1.5 + 1.5*sqrt(3))*fTotalDetectorDiameter); P.setY(-fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-1.5 + 1.5*sqrt(3))*fTotalDetectorDiameter); P.setY(0.); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-1 + 1.5*sqrt(3))*fTotalDetectorDiameter); P.setY(fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  //next to the rightmost 4
  P.setX((-1 + sqrt(3))*fTotalDetectorDiameter); P.setY(1.5*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-1 + sqrt(3))*fTotalDetectorDiameter); P.setY(-2.5*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  //4
  P.setX((-2 + 0.5*sqrt(3))*fTotalDetectorDiameter); P.setY(2*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(-2*fTotalDetectorDiameter); P.setY(1.5*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-2 + 0.5*sqrt(3))*fTotalDetectorDiameter); P.setY(-3*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(-2*fTotalDetectorDiameter); P.setY(-2.5*fTotalDetectorDiameter); P.setZ(fTotalDetectorLength - 0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  //end straddling detectors
  //10 detectors
  //leftmost 3
  P.setX(-1.5*fTotalDetectorDiameter); P.setY(-1.5*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(-1.5*fTotalDetectorDiameter); P.setY(-0.5*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX(-1.5*fTotalDetectorDiameter); P.setY(0.5*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  //4
  P.setX((-1.5 + 0.5*sqrt(3))*fTotalDetectorDiameter); P.setY(-2*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-1.5 + 0.5*sqrt(3))*fTotalDetectorDiameter); P.setY(-fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-1.5 + 0.5*sqrt(3))*fTotalDetectorDiameter); P.setY(0.); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-1.5 + 0.5*sqrt(3))*fTotalDetectorDiameter); P.setY(fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  //3
  P.setX((-1.5 + sqrt(3))*fTotalDetectorDiameter); P.setY(-1.5*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-1.5 + sqrt(3))*fTotalDetectorDiameter); P.setY(-0.5*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  P.setX((-1.5 + sqrt(3))*fTotalDetectorDiameter); P.setY(0.5*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorDiameter - 4.5*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);

  assemblyDetector->MakeImprint(fLogicWorld, Tr);
  fLogicDetector->SetVisAttributes(fWireFrameVisAtt);

  // PMT Optical Window
  //
  fSolidPMTWin1 = new G4Polyhedra("PMTWin", 0.*deg, 360.*deg, 6, 2, zPlane1, rInner, rInner1);
  fLogicPMTWin = new G4LogicalVolume(fSolidPMTWin1, fPMTWinMaterial, "PMTWin");
  fPhysiPMTWin = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMTWin), 
        			fLogicPMTWin, 
       				"PMTWin", 
       				fLogicDetector, 
       				false, 
    			    	0);
  fLogicPMTWin->SetVisAttributes(fGreenVisAtt);

  // Crystal
  //
  fSolidCrystal1 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal = new G4LogicalVolume(fSolidCrystal1, fDetectorMaterial, "Crystal");
  fPhysiCrystal = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal, 
       				"Crystal", 
       				fLogicDetector, 
       				false, 
    			    	0);
  fLogicCrystal->SetVisAttributes(fRedVisAtt);

  //Gap (gap surrounding crystal and casing or a reflector as required)
  //
  fSolidGap1 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap = new G4LogicalVolume(fSolidGap1, fGapMaterial, "Gap");
  fPhysiGap = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap, 
       				"Gap", 
       				fLogicDetector, 
       				false, 
    			    	0);

  //FaceGap (gap between face of crystal and casing or a reflector as required)
  //
  fSolidFaceGap1 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane4, rInner, rInner1);
  fLogicFaceGap = new G4LogicalVolume(fSolidFaceGap1, fGapMaterial, "FaceGap");
  fPhysiFaceGap = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap, 
       				"FaceGap", 
       				fLogicDetector, 
       				false, 
    			    	0);

  fLogicGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap->SetVisAttributes(fBlueVisAtt);
  fLogicFaceGap->SetVisAttributes(fBlueVisAtt);

  //Aluminum Casing (casing surrounding crystal)
  //
  fSolidAlCase1 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase = new G4LogicalVolume(fSolidAlCase1, fAlCaseMaterial, "AlCase");
  fPhysiAlCase = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase, 
       				"AlCase", 
       				fLogicDetector, 
       				false, 
    			    	0);

  //Aluminum Face Casing (casing covering face of crystal)
  //
  fSolidFaceAlCase1 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase = new G4LogicalVolume(fSolidFaceAlCase1, fAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase, 
       				"FaceAlCase", 
       				fLogicDetector, 
       				false, 
    			    	0);
    			    	
  //Lead casing surrounding aluminum casing (around casing surround)
  //
  fSolidPbCase1 = new G4Polyhedra("PbCase", 0.*deg, 360.*deg, 6, 2, zPlane7, rInner3, rOuter3);
  fLogicPbCase = new G4LogicalVolume(fSolidPbCase1, fPbCaseMaterial, "PbCase");
  fPhysiPbCase = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicPbCase, 
       				"PbCase", 
       				fLogicDetector, 
       				false, 
    			    	0);

  //Lead Casing Ring (at back in front of PMT face)
  //
  //Check to see if collar is needed first. Large enough crystal and small PMT will not require collar.
  if (fPbCaseDiameter<fPMTDiameter){  
	fSolidPbCollar1 = new G4Polyhedra("PbCollar", 0.*deg, 360.*deg, 6, 2, zPlane8, rInner4, rOuter4);
	fLogicPbCollar = new G4LogicalVolume(fSolidPbCollar1, fPbCaseMaterial, "PbCollar");
	fPhysiPbCollar = new G4PVPlacement(0, 
						G4ThreeVector(0.,0.,fZposPbCollar), 
						fLogicPbCollar, 
						"PbCollar", 
						fLogicDetector, 
						false, 
    			    	0);
	fLogicPbCase->SetVisAttributes(fGreyVisAtt);
	fLogicPbCollar->SetVisAttributes(fGreyVisAtt);
  }
  //PMT
  //
  fSolidPMT1 = new G4Polyhedra("PMT", 0.*deg, 360.*deg, 6, 2, zPlane9, rInner, rOuter4);
  fLogicPMT = new G4LogicalVolume(fSolidPMT1, fPMTMaterial,"PMT");
  fPhysiPMT = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT, 
       				"PMT", 
       				fLogicDetector, 
       				false, 
    			    	0); 			    	
  PrintCalorParameters();

  //always return the physical World
	}  
return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n" << "The world is made of" << fWorldMaterial << G4endl;
  G4cout << "\n" << "World Size X: " << G4BestUnit(fWorldSizeX, "Length") << G4endl;
  G4cout << "\n" << "World Size YZ: " << G4BestUnit(fWorldSizeYZ, "Length") << G4endl;  ;
if (fDetectorGeometry == 1){
  G4cout << "\n" << "The geometry is cylindrical (1)." << G4endl;
}
if (fDetectorGeometry == 2){
  G4cout << "\n" << "The geometry is hexagonal (2)." << G4endl;
}
if (fDetectorGeometry == 3){
  G4cout << "\n" << "The geometry is cylindrical. The detectors are stacked in an array (3)." << G4endl;
}
if (fDetectorGeometry == 4){
  G4cout << "\n" << "The geometry is hexagonal. The detectors are stacked in an array (4)." << G4endl;
}
  G4cout << "\n" << "The detector is made of " << fDetectorMaterial << G4endl;
  G4cout << "\n" << "The box is made of" << fBoxMaterial << G4endl;
  G4cout << "\n" << "The crystal is made of" << fCrystalMaterial << G4endl;
  G4cout << "\n" << "The gap is made of " << fGapMaterial << G4endl;
  G4cout << "\n" << "The face of the gap is made of " << fFaceGapMaterial << G4endl;
  G4cout << "\n" << "The aluminum case is made of " << fAlCaseMaterial << G4endl;
  G4cout << "\n" << "The face of the aluminum case is made of " << fFaceAlCaseMaterial << G4endl;
  G4cout << "\n" << "The lead case is made of " << fPbCaseMaterial << G4endl;
  G4cout << "\n" << "The lead collar is made of " << fPbCollarMaterial << G4endl;
  G4cout << "\n" << "The photomultiplier tube is made of" << fPMTMaterial << G4endl;
  G4cout << "\n" << "The photomultiplier tube window is made of" << fBoxMaterial << G4endl;
  G4cout << "\n" << "Total Detector Diameter: " << G4BestUnit(fTotalDetectorDiameter, "Length") << G4endl;
  G4cout << "\n" << "Detector Diameter: " << G4BestUnit(fDetectorDiameter, "Length") << G4endl;
  G4cout << "\n" << "Total Detector Length: " << G4BestUnit(fTotalDetectorLength, "Length") << G4endl;
  G4cout << "\n" << "Detector Length: " << G4BestUnit(fDetectorLength, "Length") << G4endl;
  G4cout << "\n" << "Gap Thickness: " << 0.1*fGapThickness << " cm" << G4endl;
  G4cout << "\n" << "Al Thickness: " << 0.1*fAlCaseThickness << " cm" << G4endl;
  G4cout << "\n" << "Pb Thickness: " << 0.1*fPbCaseThickness << " cm" << G4endl;
  G4cout << "\n" << "PMT Diameter: " << G4BestUnit(fPMTDiameter, "Length") << G4endl;
  G4cout << "\n" << "PMT Length: " << G4BestUnit(fPMTLength, "Length") << G4endl;
}

void DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
   //search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fWorldMaterial != pttoMaterial) {
    fWorldMaterial = pttoMaterial;                  
    if(fLogicWorld) fLogicWorld->SetMaterial(fWorldMaterial);
  }
}
   
void DetectorConstruction::SetTotalDetectorDiameter(G4double val)
{
  fTotalDetectorDiameter = val;
}

void DetectorConstruction::SetDetectorDiameter(G4double val)
{
  fDetectorDiameter = val;
}

void DetectorConstruction::SetTotalDetectorLength(G4double val)
{
  fTotalDetectorLength = val;
}

void DetectorConstruction::SetDetectorLength(G4double val)
{
  fDetectorLength = val;
}

void DetectorConstruction::SetGapThickness(G4double val)
{
  fGapThickness = val;
}

void DetectorConstruction::SetDetectorMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fDetectorMaterial != pttoMaterial) {
    fDetectorMaterial = pttoMaterial;                  
    if(fLogicDetector) fLogicDetector->SetMaterial(fDetectorMaterial);
  }
}

void DetectorConstruction::SetBoxMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fBoxMaterial != pttoMaterial) {
    fBoxMaterial = pttoMaterial;                  
    if(fLogicBox) fLogicBox->SetMaterial(fBoxMaterial);
  }
}

void DetectorConstruction::SetCrystalMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fCrystalMaterial != pttoMaterial) {
    fCrystalMaterial = pttoMaterial;                  
    if(fLogicCrystal) fLogicCrystal->SetMaterial(fCrystalMaterial);
  }
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

void DetectorConstruction::SetFaceGapMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fFaceGapMaterial != pttoMaterial) {
    fFaceGapMaterial = pttoMaterial;                  
    if(fLogicFaceGap) fLogicFaceGap->SetMaterial(fFaceGapMaterial);
  }
}

void DetectorConstruction::SetAlCaseMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fAlCaseMaterial != pttoMaterial) {
    fAlCaseMaterial = pttoMaterial;                  
    if(fLogicAlCase) fLogicAlCase->SetMaterial(fAlCaseMaterial);
  }
}

void DetectorConstruction::SetFaceAlCaseMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fFaceAlCaseMaterial != pttoMaterial) {
    fFaceAlCaseMaterial = pttoMaterial;                  
    if(fLogicFaceAlCase) fLogicFaceAlCase->SetMaterial(fFaceAlCaseMaterial);
  }
}

void DetectorConstruction::SetPbCaseMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fPbCaseMaterial != pttoMaterial) {
    fPbCaseMaterial = pttoMaterial;                  
    if(fLogicPbCase) fLogicPbCase->SetMaterial(fPbCaseMaterial);
  }
}

void DetectorConstruction::SetPbCollarMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fPbCollarMaterial != pttoMaterial) {
    fPbCollarMaterial = pttoMaterial;                  
    if(fLogicPbCollar) fLogicPbCollar->SetMaterial(fPbCollarMaterial);
  }
}

void DetectorConstruction::SetPMTMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fPMTMaterial != pttoMaterial) {
    fPMTMaterial = pttoMaterial;                  
    if(fLogicPMT) fLogicPMT->SetMaterial(fPMTMaterial);
  }
}

void DetectorConstruction::SetPMTWinMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fPMTWinMaterial != pttoMaterial) {
    fPMTWinMaterial = pttoMaterial;                  
    if(fLogicPMTWin) fLogicPMTWin->SetMaterial(fPMTWinMaterial);
  }
}

void DetectorConstruction::SetAlCaseThickness(G4double val)
{
  fAlCaseThickness = val;
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
void DetectorConstruction::SetDetectorGeometry(G4int val)
{
  fDetectorGeometry = val;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
