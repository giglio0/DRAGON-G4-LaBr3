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

  //Set BGO and LaBr3 Detector Array Defaults according to the Geant3 dimensions
  fPMTDiameter = 5.9*cm;
  fPMTLength = 20.499*cm;
  fDetectorLength = 7.62*cm;
  fDetectorDiameter = 5.58*cm;
  fGapThickness = 0.0355*cm;
  fGapFaceThickness = 0.3175*cm;
  fAlCaseThickness = 0.0635*cm;
  fAlFaceThickness = 0.0635*cm;
  
  //Set Single (not array) LaBr3 Detector Defaults 
  fPMTDiameter1 = 5.1*cm;
  fPMTLength1 = 6.2*cm;
  fPMTWinThickness1 = 0.2*cm;
  fDetectorLength1 = 5.08*cm;
  fDetectorDiameter1 = 5.08*cm;
  fGapFaceThickness1 = 0.42*cm;
  fAlCaseThickness1 = 0.05*cm;
  fAlFaceThickness1 = 0.05*cm;
  
  fPbCaseThickness = 0.01*cm;
  
  fOutTBoxX = 8.573*cm;
  fOutTBoxY = 2.381*cm;
  fOutTBoxZ = 12.859*cm;
  fAirGap = 0.;
  Temperature = 298.15*kelvin;
  Pressure = 101325*pascal;
  
  ComputeCalorParameters();
  
  // materials  
  DefineMaterials();
  // Default Materials
  fWorldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Air");
  fDetectorMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Air");
  fCrystalMaterial = G4NistManager::Instance()->FindOrBuildMaterial("LaBr3Ce");
  fCrystalMaterial1 = G4NistManager::Instance()->FindOrBuildMaterial("LaBr3Ce");
  fAlCaseMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Aluminium");
  fFaceAlCaseMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Aluminium");
  fPbCaseMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Lead");
  fPbCollarMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Lead");
  fPMTMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Pyrex Glass");
  fPMTWinMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Optical");
  fPhotoCathodeMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Bialkali");
  fPMTIntMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Artificial Vacuum");
  fGasTargetMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Hydrogen Gas");
  fInnerTrapezoidMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Hydrogen Gas 2");
  //Set Default Gap Material i.e. no reflector
  fGapMaterial = G4NistManager::Instance()->FindOrBuildMaterial("MgO");
  fFaceGapMaterial = G4NistManager::Instance()->FindOrBuildMaterial("MgO");
  //Set Visualization Attributes

  fCyanVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
  fYellowVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0)) ;
  fMagnetaVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
  fOrangeVisAtt = new G4VisAttributes(G4Colour(0.5,0.75,0.75));
  // Blue
  fBlueVisAtt = new G4VisAttributes(G4Colour(0.,0.,0.75));
  // Green
  fGreenVisAtt = new G4VisAttributes(G4Colour(0.,0.75,0.));
  // Red
  fRedVisAtt = new G4VisAttributes(G4Colour(0.75,0.,0.));
  // Grey
  fGreyVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  fWhiteVisAtt = new G4VisAttributes(G4Colour(1.,1.,1., 0));
  fBlackVisAtt = new G4VisAttributes(G4Colour(0.,0.,0.));
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
  
  // define Elements
  G4Element* H  = new G4Element("Hydrogen",symbol="H", z=1, a= 1.01*g/mole);
  G4Element* B = new G4Element("Boron",symbol="B", z=5, a= 10.81*g/mole); 
  G4Element* C  = new G4Element("Carbon",  symbol="C", z=6, a= 12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N", z=7, a= 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",  symbol="O", z=8, a= 16.00*g/mole);
  G4Element* F = new G4Element("Flourine", symbol="F", z=9, a= 19.00*g/mole);
  G4Element* Na = new G4Element("Sodium",  symbol="Na", z=11, a= 22.99*g/mole);
  G4Element* Mg = new G4Element("Magnesium", symbol="Mg", z=12, a= 24.305*g/mole);
  G4Element* Al = new G4Element("Aluminium", symbol="Al", z=13, a=26.98*g/mole);
  G4Element* Si = new G4Element("Silicon", symbol="Si", z=14, a= 28.0855*g/mole);
  G4Element* Ar = new G4Element("Argon", symbol="Ar", z=18, a= 39.95*g/mole);
  G4Element* K = new G4Element("Potassium", symbol="K", z=19, a= 39.10*g/mole);
  G4Element* Ca = new G4Element("Calcium", symbol="Ca", z=20, a= 40.078*g/mole);
  G4Element* Ge = new G4Element("Germanium", symbol="Ge", z=32, a= 72.63*g/mole);
  G4Element* Br = new G4Element("Bromine",   symbol="Br", z=35, a= 79.90*g/mole);
  G4Element* Sr = new G4Element("Strontium",   symbol="Sr", z=38, a= 87.62*g/mole);
  G4Element* Y = new G4Element("Yttrium", symbol="Y", z=39, a= 88.906*g/mole);
  G4Element* Sb = new G4Element("Antimony", symbol="Sb", z=51, a = 121.76*g/mole);
  G4Element* I  = new G4Element("Iodine",  symbol="I", z=53, a= 126.90*g/mole);
  G4Element* Cs = new G4Element("Cesium", symbol="Cs", z=55, a= 132.91*g/mole);
  G4Element* La = new G4Element("Lanthanum",  symbol="La", z=57, a= 138.91*g/mole);
  G4Element* Ce = new G4Element("Cerium", symbol="Ce", z=58, a= 140.116*g/mole);
  G4Element* Lu = new G4Element("Lutetium", symbol="Lu", z=71, a= 174.967*g/mole);
  G4Element* Bi = new G4Element("Bismuth", symbol="Bi", z=83, a= 208.980*g/mole);
 
  // define simple materials
  new G4Material("H2Liq"    , z= 1, a= 1.01*g/mole, density= 70.8*mg/cm3);
  new G4Material("Beryllium", z= 4, a= 9.01*g/mole, density= 1.848*g/cm3);
  G4Material* Silicon = new G4Material("Silicon"  , z=14, a=28.09*g/mole, density= 2.330*g/cm3, kStateSolid, Temperature, Pressure);
  Silicon->GetIonisation()->SetMeanExcitationEnergy(173*eV);
  new G4Material("Iron",     z=26, a= 55.85*g/mole, density= 7.870*g/cm3);
  G4Material* Copper = new G4Material("Copper",   z=29, a= 63.55*g/mole, density= 8.960*g/cm3);
  Copper->GetIonisation()->SetMeanExcitationEnergy(322.*eV);
  new G4Material("Germanium",z=32, a= 72.61*g/mole, density= 5.323*g/cm3);
  G4Material* Silver = new G4Material("Silver",   z=47, a=107.87*g/mole, density= 10.50*g/cm3);
  Silver->GetIonisation()->SetMeanExcitationEnergy(470.*eV);
  new G4Material("Tungsten", z=74, a=183.85*g/mole, density= 19.30*g/cm3);
  new G4Material("Gold",     z=79, a=196.97*g/mole, density= 19.32*g/cm3);
  
  G4Material* Aluminium = new G4Material("Aluminium", z=13, a=26.98*g/mole, density= 2.699*g/cm3, kStateSolid, Temperature, Pressure);
  Aluminium->GetIonisation()->SetMeanExcitationEnergy(166*eV);
  G4Material* Lead = new G4Material("Lead", z=82, a=207.19*g/mole, density= 11.35*g/cm3, kStateSolid, Temperature, Pressure);
  Lead->GetIonisation()->SetMeanExcitationEnergy(823*eV);

  // define a material from elements.   case 1: chemical molecule
  // G4_SODIUM_IODIDE
  G4Material* NaI = new G4Material("NaI", density= 3.667*g/cm3, ncomponents=2, kStateSolid, Temperature, Pressure);
  NaI->AddElement(Na, natoms=1);
  NaI->AddElement(I , natoms=1);
  NaI->GetIonisation()->SetMeanExcitationEnergy(452*eV);

  // define a material from elements.   case 2: mixture by fractional mass
 
  G4Material* LaBr3 = new G4Material("LaBr3", density= 5.06*g/cm3, ncomponents=2, kStateSolid, Temperature, Pressure);
  LaBr3->AddElement(La, natoms=1);
  LaBr3->AddElement(Br, natoms=3);
  LaBr3->GetIonisation()->SetMeanExcitationEnergy(400.96625*eV);
  
  G4Material* SrBr2 = new G4Material("SrBr2", density= 4.216*g/cm3, ncomponents=2, kStateSolid, Temperature, Pressure);
  SrBr2->AddElement(Sr, natoms=1);
  SrBr2->AddElement(Br, natoms=2);
  SrBr2->GetIonisation()->SetMeanExcitationEnergy(400.96625*eV);
  
  G4Material* LaBr3Ce = new G4Material("LaBr3Ce", density= 5.08*g/cm3, ncomponents=2, kStateSolid, Temperature, Pressure);
  LaBr3Ce->AddMaterial(LaBr3, fractionmass=0.980894429);
  LaBr3Ce->AddElement(Ce, fractionmass=0.01910557);
  LaBr3Ce->GetIonisation()->SetMeanExcitationEnergy(403.298315*eV);
  
  G4Material* BGO = new G4Material("BGO", density = 7.13*g/cm3, ncomponents=3, kStateSolid, Temperature, Pressure);
  BGO->AddElement(Bi, natoms=1);
  BGO->AddElement(Ge, natoms=1);
  BGO->AddElement(O, natoms=1);
  BGO->GetIonisation()->SetMeanExcitationEnergy(534.1*eV);
 
  // G4_PYREX_GLASS
  G4Material* PMTWinMaterial = new G4Material("Pyrex Glass", density = 2.23*g/cm3, ncomponents=6, kStateSolid, Temperature, Pressure);
  PMTWinMaterial->AddElement(B, fractionmass=0.0400639);
  PMTWinMaterial->AddElement(O, fractionmass=0.539561);
  PMTWinMaterial->AddElement(Na, fractionmass=0.0281909);
  PMTWinMaterial->AddElement(Al, fractionmass=0.011644);
  PMTWinMaterial->AddElement(Si, fractionmass=0.377219);
  PMTWinMaterial->AddElement(K, fractionmass=0.00332099);
  PMTWinMaterial->GetIonisation()->SetMeanExcitationEnergy(134*eV);
  
  G4Material* Bialkali = new G4Material("Bialkali", density = 3.735419367*g/cm3, ncomponents=3, kStateSolid, Temperature, Pressure);
  Bialkali->AddElement(Cs, natoms=1);
  Bialkali->AddElement(K, natoms=1);
  Bialkali->AddElement(Sb, natoms=1);
  Bialkali->GetIonisation()->SetMeanExcitationEnergy(447.9233955*eV);
  
  G4Material* Optical = new G4Material("Optical", density = 0.965*g/cm3, ncomponents=4, kStateSolid, Temperature, Pressure);
  Optical->AddElement(C, natoms=2);
  Optical->AddElement(H, natoms=6);
  Optical->AddElement(O, natoms=1);
  Optical->AddElement(Si, natoms=1);
  Bialkali->GetIonisation()->SetMeanExcitationEnergy(113.8252413*eV);
  
  // G4_WATER
  G4Material* Water = new G4Material("Water", density= 1.*g/cm3, ncomponents=2, kStateLiquid, Temperature, Pressure);
  Water->AddElement(H, natoms=2);
  Water->AddElement(O, natoms=1);
  Water->GetIonisation()->SetMeanExcitationEnergy(78.*eV);
  
  // G4_WATER_VAPOR
  G4Material* WaterVapor = new G4Material("Water Vapor", density= 0.000756182*g/cm3, ncomponents=2, kStateGas, Temperature, Pressure);
  WaterVapor->AddElement(H, natoms=2);
  WaterVapor->AddElement(O, natoms=1);
  WaterVapor->GetIonisation()->SetMeanExcitationEnergy(71.6*eV);
  
  // G4_AIR
  G4Material* Air = new G4Material("Air", density= 0.00120479*g/cm3, ncomponents=4, kStateGas, 298.15*kelvin, 101325*pascal);
  Air->AddElement(C, fractionmass=0.000124);
  Air->AddElement(N, fractionmass=0.755268);
  Air->AddElement(O, fractionmass=0.231781);
  Air->AddElement(Ar, fractionmass=0.012827);
  Air->GetIonisation()->SetMeanExcitationEnergy(85.7*eV);
  
  // G4_MAGNESIUM_OXIDE
  G4Material* MgO = new G4Material("MgO", density= 3.58*g/cm3, ncomponents=2, kStateGas, Temperature, Pressure);
  MgO->AddElement(Mg, natoms=1);
  MgO->AddElement(O, natoms=1);
  MgO->GetIonisation()->SetMeanExcitationEnergy(143.8*eV);
  
  // G4_H
  G4Material* Hydrogen = new G4Material("Hydrogen Gas", density= 0.000083748*g/cm3, ncomponents=1, kStateGas, Temperature, (0.007/11400)*atmosphere);
  Hydrogen->AddElement(H, natoms=2);
  Hydrogen->GetIonisation()->SetMeanExcitationEnergy(19.2*eV);
  
  // G4_H
  G4Material* Hydrogen2 = new G4Material("Hydrogen Gas 2", density= 0.000083748*g/cm3, ncomponents=1, kStateGas, Temperature, (0.007/760)*atmosphere);
  Hydrogen2->AddElement(H, natoms=2);
  Hydrogen2->GetIonisation()->SetMeanExcitationEnergy(19.2*eV);
  
  // Artificial Vacuum
  new G4Material("Artificial Vacuum", z=1, a=1.01*g/mole, universe_mean_density, kStateGas, 298.15*kelvin,(3.e-18)*pascal);
  
  // example of vacuum
  //from PhysicalConstants.h
  new G4Material("Galactic", z=1, a=1.01*g/mole, universe_mean_density,
                 kStateGas,2.73*kelvin,(3.e-18)*pascal);
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
  // BGO and LaBr3 Array
  fGapLength = fDetectorLength + fGapFaceThickness;
  fGapDiameter = fDetectorDiameter + 2*fGapThickness;
  fAlCaseLength = fGapLength + fAlFaceThickness;
  fAlCaseDiameter = fGapDiameter + 2*fAlCaseThickness;
  fPbCaseDiameter = fAlCaseDiameter;
  fTotalDetectorLength = fDetectorLength + fPMTLength + fGapFaceThickness + fAlFaceThickness;
  
  // LaBr3
  fAlCaseLength1 = 5.75*cm;
  fAlCaseDiameter1 = 5.18*cm;
  fTotalDetectorLength1 = fPMTWinThickness1 + fDetectorLength1 + fPMTLength1 + fGapFaceThickness1 + fAlFaceThickness1;
  
  if (fPbCaseDiameter>fPMTDiameter){
	  fPMTDiameter = fPbCaseDiameter;
	  G4cout << "\n" << "\n" 
	  << "WARNING: PMT Diameter has been adjusted to compensate for large crystal/lead/aluminum size"  
	  << G4endl;
	  fTotalDetectorDiameter = fPMTDiameter;
	  fTotalDetectorDiameter1 = 5.18*cm;
  }
  else
	fTotalDetectorDiameter = fPMTDiameter;
	fTotalDetectorDiameter1 = 5.18*cm;
	
  // BGO Component Positions 
  fZposPbCollar = 0.5*fTotalDetectorLength - fAlFaceThickness - fGapThickness - fDetectorLength - fPMTWinThickness +  
  0.5*fPbCaseThickness;
  fZposFaceAlCase = 0.5*fTotalDetectorLength - 0.5*fAlFaceThickness;
  fZposFaceGap = 0.5*fTotalDetectorLength - fAlFaceThickness - 0.5*fGapFaceThickness ;
  fZpos = 0.5*fTotalDetectorLength - fGapFaceThickness - fAlFaceThickness - 0.5*fDetectorLength;
   
  fZposAlCase = 0.5*fTotalDetectorLength - 0.5*fAlCaseLength;
  fZposGap = 0.5*fTotalDetectorLength - fAlFaceThickness - 0.5*fGapLength ;
  fZposPMT = 0.5*fTotalDetectorLength - fAlFaceThickness - fGapFaceThickness - fDetectorLength - 0.5*fPMTLength;
  
  // LaBr3 Component Positions 
  fZposFaceAlCase1 = 0.5*fTotalDetectorLength1 - 0.5*fAlFaceThickness1;
  fZposFaceGap1 = 0.5*fTotalDetectorLength1 - fAlFaceThickness1 - 0.5*fGapFaceThickness1 ;
  fZpos1 = 0.5*fTotalDetectorLength1 - fAlFaceThickness1 - fGapFaceThickness1 - 0.5*fDetectorLength1;

  fZposAlCase1 = 0.5*fTotalDetectorLength1 - 0.5*fAlCaseLength1;
  fZposPMTWin1 = 0.5*fTotalDetectorLength1 - fAlFaceThickness1 - fGapFaceThickness1 - fDetectorLength1 - 0.5*fPMTWinThickness1;
  fZposPMT1 = 0.5*fTotalDetectorLength1 - fAlFaceThickness1 - fGapFaceThickness1 - fDetectorLength1 - fPMTWinThickness1 -
  0.5*fPMTLength1;
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
  fLogicWorld->SetVisAttributes(fWhiteVisAtt);
  //Detector Holder Volume for Components
  //
  //  If rotation required  
  //  G4RotationMatrix rotm  = G4RotationMatrix(0,90*deg,-90*deg);     
  //  No Rotation Now
  G4RotationMatrix rotm  = G4RotationMatrix(0,0,0.*deg);    
  G4ThreeVector position = G4ThreeVector(0.,0.,0.);
  G4Transform3D transform = G4Transform3D(rotm,position);
  G4ThreeVector P;
  G4Transform3D Tr;
  // The detectors on each side of the array face each other with this rotation matrix
  G4RotationMatrix rotm180 = G4RotationMatrix(0,180.*deg,0); 
  
if (fDetectorGeometry == 1 || fDetectorGeometry == 2 || fDetectorGeometry == 3 || fDetectorGeometry == 4){

if (fDetectorGeometry == 1 || fDetectorGeometry == 3){  
  // Single cylindrical detector 
  if (fDetectorGeometry == 1){
  fSolidDetector = new G4Tubs("Detector", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector1 = new G4Tubs("Detector", 0., 0.5*fTotalDetectorDiameter1, 0.5*fTotalDetectorLength1, 0.*deg, 360.*deg);
  fLogicDetector = new G4LogicalVolume(fSolidDetector, fDetectorMaterial, "BGO Detector");
  fLogicDetector1 = new G4LogicalVolume(fSolidDetector1, fDetectorMaterial, "LaBr3 Detector");
  //fPhysiDetector = new G4PVPlacement(transform,
        			//fLogicDetector, 
       				//"BGO Detector", 
       				//fLogicWorld, 
       				//false, 
    			    	//0,
				    //false);
  fPhysiDetector1 = new G4PVPlacement(transform,
        			fLogicDetector1, 
       				"LaBr3 Detector", 
       				fLogicWorld, 
       				false, 
    			    	0,
				    false);

  fLogicDetector->SetVisAttributes(fWireFrameVisAtt);
  fLogicDetector1->SetVisAttributes(fWireFrameVisAtt);
  }
  if (fDetectorGeometry == 3){  
  //  Cylindrical Detector Array	
 
  fSolidDetector = new G4Tubs("BGO Detector", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector1 = new G4Tubs("LaBr3 Detector", 0., 0.5*fTotalDetectorDiameter1, 0.5*fTotalDetectorLength1, 0.*deg, 360.*deg);
  fLogicDetector = new G4LogicalVolume(fSolidDetector,fDetectorMaterial, "BGO Detector");
  fLogicDetector1 = new G4LogicalVolume(fSolidDetector1,fDetectorMaterial, "LaBr3 Detector");
  G4AssemblyVolume* assemblyDetector = new G4AssemblyVolume();
  // 30 detectors are placed relative to the array volume 
  // Back - 24, 23, 22, 20, 18, 16, 14, 12, 26, 28, 30
  // 10.069 cm forward in Z for BGO detector positions, rotm180
  // 1
  P.setX(14.78*cm); P.setY(-4.96*cm); P.setZ(1.81*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(14.78*cm); P.setY(-4.96*cm); P.setZ(1.81*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 2
  P.setX(11.83*cm); P.setY(-10.08*cm); P.setZ(9.01*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(11.83*cm); P.setY(-10.08*cm); P.setZ(9.01*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 3
  P.setX(11.83*cm); P.setY(4.96*cm); P.setZ(2.309*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(11.83*cm); P.setY(4.96*cm); P.setZ(2.309*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 4
  P.setX(8.87*cm); P.setY(10.08*cm); P.setZ(9.01*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(8.87*cm); P.setY(10.08*cm); P.setZ(9.01*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 5
  P.setX(2.96*cm); P.setY(7.68*cm); P.setZ(9.01*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(2.96*cm); P.setY(7.68*cm); P.setZ(9.01*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 6
  P.setX(-2.96*cm); P.setY(7.68*cm); P.setZ(9.01*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(-2.96*cm); P.setY(7.68*cm); P.setZ(9.01*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 7
  P.setX(-8.87*cm); P.setY(10.08*cm); P.setZ(9.01*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(-8.87*cm); P.setY(10.08*cm); P.setZ(9.01*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 8
  P.setX(-11.83*cm); P.setY(-10.08*cm); P.setZ(9.01*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(-11.83*cm); P.setY(-10.08*cm); P.setZ(9.01*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 9
  P.setX(-11.83*cm); P.setY(4.96*cm); P.setZ(9.01*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(-11.83*cm); P.setY(4.96*cm); P.setZ(9.01*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 10
  P.setX(-14.78*cm); P.setY(-4.96*cm); P.setZ(9.01*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(-14.78*cm); P.setY(-4.96*cm); P.setZ(9.01*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 11
  P.setX(8.87*cm); P.setY(-2.56*cm); P.setZ(17.859*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(8.87*cm); P.setY(-2.56*cm); P.setZ(17.859*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 13
  P.setX(5.91*cm); P.setY(-7.68*cm); P.setZ(17.859*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(5.91*cm); P.setY(-7.68*cm); P.setZ(17.859*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 15
  P.setX(5.91*cm); P.setY(2.56*cm); P.setZ(17.859*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(5.91*cm); P.setY(2.56*cm); P.setZ(17.859*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 17
  P.setX(2.96*cm); P.setY(-2.56*cm); P.setZ(17.859*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(2.96*cm); P.setY(-2.56*cm); P.setZ(17.859*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 19
  P.setX(0.00*cm); P.setY(-7.68*cm); P.setZ(17.859*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(0.00*cm); P.setY(-7.68*cm); P.setZ(17.859*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 21
  P.setX(0.00*cm); P.setY(2.56*cm); P.setZ(17.859*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(0.00*cm); P.setY(2.56*cm); P.setZ(17.859*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 23
  P.setX(-2.96*cm); P.setY(-2.56*cm); P.setZ(17.859*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(-2.96*cm); P.setY(-2.56*cm); P.setZ(17.859*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 25
  P.setX(-5.91*cm); P.setY(-7.68*cm); P.setZ(17.859*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(-5.91*cm); P.setY(-7.68*cm); P.setZ(17.859*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 27
  P.setX(-5.91*cm); P.setY(2.56*cm); P.setZ(17.859*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(-5.91*cm); P.setY(2.56*cm); P.setZ(17.859*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 29
  P.setX(-8.87*cm); P.setY(-2.56*cm); P.setZ(17.859*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(-8.87*cm); P.setY(-2.56*cm); P.setZ(17.859*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 10.069 cm back in Z for BGO detector positions, rotm
  // 30
  P.setX(-8.87*cm); P.setY(-2.56*cm); P.setZ(-17.859*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(-8.87*cm); P.setY(-2.56*cm); P.setZ(-17.859*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 28
  P.setX(-5.91*cm); P.setY(2.56*cm); P.setZ(-17.859*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(-5.91*cm); P.setY(2.56*cm); P.setZ(-17.859*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 22
  P.setX(0.00*cm); P.setY(2.56*cm); P.setZ(-17.859*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(0.00*cm); P.setY(2.56*cm); P.setZ(-17.859*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 16
  P.setX(5.91*cm); P.setY(2.56*cm); P.setZ(-17.859*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(5.91*cm); P.setY(2.56*cm); P.setZ(-17.859*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 24
  P.setX(-2.96*cm); P.setY(-2.56*cm); P.setZ(-17.859*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(-2.96*cm); P.setY(-2.56*cm); P.setZ(-17.859*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 18
  P.setX(2.96*cm); P.setY(-2.56*cm); P.setZ(-17.859*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(2.96*cm); P.setY(-2.56*cm); P.setZ(-17.859*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 12
  P.setX(8.87*cm); P.setY(-2.56*cm); P.setZ(-17.859*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(8.87*cm); P.setY(-2.56*cm); P.setZ(-17.859*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 26
  P.setX(-5.91*cm); P.setY(-7.68*cm); P.setZ(-17.859*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(-5.91*cm); P.setY(-7.68*cm); P.setZ(-17.859*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 20
  P.setX(0.00*cm); P.setY(-7.68*cm); P.setZ(-17.859*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(0.00*cm); P.setY(-7.68*cm); P.setZ(-17.859*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 14
  P.setX(5.91*cm); P.setY(-7.68*cm); P.setZ(-17.859*cm);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(5.91*cm); P.setY(-7.68*cm); P.setZ(-17.859*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  
  // Assembly Placement in the World
  // Rotated 90° along the axes to be perpendicular to the gas target
  // Rotated 180° along the z-axis the place the array downside up instead of upside down
  // Coordinate transformation: (X, Y, Z) -> (Z, Y, X)
  G4ThreeVector WorldP(0,0,0);
  G4RotationMatrix worldrotm = G4RotationMatrix(90.*deg,90.*deg,270.*deg); 
  G4Transform3D WorldTr = G4Transform3D(worldrotm, WorldP);
  assemblyDetector->MakeImprint(fLogicWorld, WorldTr);
  fLogicDetector->SetVisAttributes(fWireFrameVisAtt);
  }
  // PMT Optical Window
  //
  //fSolidPMTWin = new G4Tubs("PMTWin", 0., 0.5*fDetectorDiameter, 0.5*fPMTWinThickness, 0.*deg, 360.*deg);
  //fLogicPMTWin = new G4LogicalVolume(fSolidPMTWin1, fPMTWinMaterial, "PMTWin");
  //fPhysiPMTWin = new G4PVPlacement(0, 
        			//G4ThreeVector(0.,0.,fZposPMTWin), 
        			//fLogicPMTWin, 
       				//"PMTWin", 
       				//fLogicDetector, 
       				//false, 
    			    	//0);
  //fLogicPMTWin->SetVisAttributes(fGreenVisAtt);
  
  fSolidPMTWin1 = new G4Tubs("PMTWin", 0., 0.5*fDetectorDiameter1, 0.5*fPMTWinThickness1, 0.*deg, 360.*deg);
  fLogicPMTWin1 = new G4LogicalVolume(fSolidPMTWin1, fPMTWinMaterial, "PMTWin");
  fPhysiPMTWin1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMTWin1), 
        			fLogicPMTWin1, 
       				"PMTWin", 
       				fLogicDetector1, 
       				false, 
    			    	0);
  fLogicPMTWin1->SetVisAttributes(fGreenVisAtt);

  //// Crystal
  //
  fSolidCrystal = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal = new G4LogicalVolume(fSolidCrystal, fCrystalMaterial, "Crystal");
  fPhysiCrystal = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal, 
       				"Crystal", 
       				fLogicDetector, 
       				false, 
    			    	0);
    			    	
  fSolidCrystal1 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter1, 0.5*fDetectorLength1, 0.*deg, 360.*deg);
  fLogicCrystal1 = new G4LogicalVolume(fSolidCrystal1, fCrystalMaterial1, "Crystal");
  fPhysiCrystal1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos1), 
        			fLogicCrystal1, 
       				"Crystal", 
       				fLogicDetector1, 
       				false, 
    			    	0);
    			    	
  fLogicCrystal->SetVisAttributes(fRedVisAtt);
  fLogicCrystal1->SetVisAttributes(fMagnetaVisAtt);

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
  fSolidFaceGap = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap = new G4LogicalVolume(fSolidFaceGap, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap, 
       				"FaceGap", 
       				fLogicDetector, 
       				false, 
    			    	0);
    
  fSolidFaceGap1 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter1, 0.5*fGapFaceThickness1, 0.*deg, 360.*deg);			    	
  fLogicFaceGap1 = new G4LogicalVolume(fSolidFaceGap1, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap1), 
        			fLogicFaceGap1, 
       				"FaceGap", 
       				fLogicDetector1, 
       				false, 
    			    	0);

  fLogicGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap->SetVisAttributes(fBlueVisAtt);
  fLogicFaceGap->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap1->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap1->SetVisAttributes(fCyanVisAtt);

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
  
  fSolidAlCase1 = new G4Tubs("AlCase", 0.5*fDetectorDiameter1, 0.5*fAlCaseDiameter1, 0.5*fAlCaseLength1, 0.*deg, 360.*deg); 			    	
  fLogicAlCase1 = new G4LogicalVolume(fSolidAlCase1, fAlCaseMaterial, "AlCase");
  fPhysiAlCase1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase1), 
        			fLogicAlCase1, 
       				"AlCase", 
       				fLogicDetector1, 
       				false, 
    			    	0);
  
  //Aluminum Face Casing (casing covering face of crystal)
  //
  fSolidFaceAlCase = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase = new G4LogicalVolume(fSolidFaceAlCase, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase, 
       				"FaceAlCase", 
       				fLogicDetector, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase1 = new G4Tubs("FaceAlCase", 0., 0.5*fDetectorDiameter1, 0.5*fAlFaceThickness1, 0.*deg, 360.*deg);  			    	
  fLogicFaceAlCase1 = new G4LogicalVolume(fSolidFaceAlCase1, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase1), 
        			fLogicFaceAlCase1, 
       				"FaceAlCase", 
       				fLogicDetector1, 
       				false, 
    			    	0);
    			    	
  fLogicAlCase->SetVisAttributes(fYellowVisAtt);
  fLogicFaceAlCase->SetVisAttributes(fYellowVisAtt);
  fLogicAlCase1->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase1->SetVisAttributes(fGreyVisAtt);
  
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
    			    	
  fSolidPMT1 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter1, 0.5*fPMTLength1, 0.*deg, 360.*deg); 			    	
  fLogicPMT1 = new G4LogicalVolume(fSolidPMT1, fPMTMaterial,"PMT");
  fPhysiPMT1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT1), 
        			fLogicPMT1, 
       				"PMT", 
       				fLogicDetector1, 
       				false, 
    			    	0);
    			    	
  fSolidPMTI = new G4Tubs("Vacuum", 0., 0.5*fPMTDiameter1 - 0.25*cm, 0.5*fPMTLength1 - 0.508*cm, 0.*deg, 360.*deg); 			    	
  fLogicPMTI = new G4LogicalVolume(fSolidPMTI, fPMTIntMaterial,"Vacuum");
  fPhysiPMTI = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTI, 
       				"Vacuum", 
       				fLogicPMT1, 
       				false, 
    			    	0);
    			    	
  fSolidPMTK = new G4Tubs("Photocathode", 0., 0.5*fPMTDiameter1 - 0.25*cm, 0.254*cm, 0.*deg, 360.*deg);			    	
  fLogicPMTK = new G4LogicalVolume(fSolidPMTK, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTK = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,2.846*cm), 
        			fLogicPMTK, 
       				"Photocathode", 
       				fLogicPMT1, 
       				false, 
    			    	0);
    			    	
  fLogicPMTK->SetVisAttributes(fYellowVisAtt);
  }
  
if (fDetectorGeometry == 2 || fDetectorGeometry == 4){
//  Single hexagonal prism detector	
//  If rotation required  
//  G4RotationMatrix rotm  = G4RotationMatrix(0,90*deg,-90*deg);     
//  No Rotation Now
// BGO Dimensions
  const G4double zPlane[2] = {-0.5*fTotalDetectorLength, 0.5*fTotalDetectorLength};
  const G4double zPlane1[2] = {-0.5*fPMTWinThickness, 0.5*fPMTWinThickness};
  const G4double zPlane2[2] = {-0.5*fDetectorLength, 0.5*fDetectorLength};
  const G4double zPlane3[2] = {-0.5*fGapLength, 0.5*fGapLength};
  const G4double zPlane4[2] = {-0.5*fGapThickness, 0.5*fGapThickness};
  const G4double zPlane5[2] = {-0.5*fAlCaseLength, 0.5*fAlCaseLength};
  const G4double zPlane6[2] = {-0.5*fAlFaceThickness, 0.5*fAlFaceThickness};
  const G4double zPlane7[2] = {-0.5*fAlCaseLength, 0.5*fAlCaseLength};
  const G4double zPlane8[2] = {-0.5*fPbCaseThickness, 0.5*fPbCaseThickness};
  const G4double zPlane9[2] = {-0.5*fPMTLength, 0.5*fPMTLength};
  const G4double zPlane10[2] = {-0.5*fGapFaceThickness, 0.5*fGapFaceThickness};
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
// LaBr3 Dimensions   
  const G4double zPlaneL[2] = {-0.5*fTotalDetectorLength1, 0.5*fTotalDetectorLength1};
  const G4double zPlaneL1[2] = {-0.5*fPMTWinThickness1, 0.5*fPMTWinThickness1};
  const G4double zPlaneL2[2] = {-0.5*fDetectorLength1, 0.5*fDetectorLength1};
  const G4double zPlaneL5[2] = {-0.5*fAlCaseLength1, 0.5*fAlCaseLength1};
  const G4double zPlaneL6[2] = {-0.5*fAlFaceThickness1, 0.5*fAlFaceThickness1};
  const G4double zPlaneL7[2] = {-0.5*fAlCaseLength1, 0.5*fAlCaseLength1};
  const G4double zPlaneL8[2] = {-0.5*fPbCaseThickness, 0.5*fPbCaseThickness};
  const G4double zPlaneL9[2] = {-0.5*fPMTLength1, 0.5*fPMTLength1};
  const G4double rInnerL[2] = {0,0};
  const G4double rInnerL1[2] = {0.5*fDetectorDiameter1, 0.5*fDetectorDiameter1};
  const G4double rInnerL3[2] = {0.5*fAlCaseDiameter1, 0.5*fAlCaseDiameter1};
  const G4double rInnerL4[2] = {0.5*fPbCaseDiameter,0.5*fPbCaseDiameter};
  const G4double rOuterL[2] = {0.5*fTotalDetectorDiameter1, 0.5*fTotalDetectorDiameter1};
  const G4double rOuterL2[2] = {0.5*fAlCaseDiameter1, 0.5*fAlCaseDiameter1};
  const G4double rOuterL3[2]= {0.5*fPbCaseDiameter, 0.5*fPbCaseDiameter};
  const G4double rOuterL4[2] = {0.5*fPMTDiameter1, 0.5*fPMTDiameter1};
  
  if (fDetectorGeometry == 2){
  fSolidHDetector = new G4Polyhedra{"BGO Detector", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector1 = new G4Polyhedra{"LaBr3 Detector", 0.*deg, 360.*deg, 6, 2, zPlaneL, rInnerL, rOuterL};
  fLogicDetector = new G4LogicalVolume(fSolidHDetector,fDetectorMaterial, "BGO Detector");
  fLogicDetector1 = new G4LogicalVolume(fSolidHDetector1,fDetectorMaterial, "LaBr3 Detector");
  fPhysiDetector = new G4PVPlacement(transform,
        			fLogicDetector, 
       				"BGO Detector", 
       				fLogicWorld, 
       				false, 
    			    	0,
				    false);
  //fPhysiDetector1 = new G4PVPlacement(transform,
        			//fLogicDetector1, 
       				//"LaBr3 Detector", 
       				//fLogicWorld, 
       				//false, 
    			    	//0,
				    //false);
  fLogicDetector->SetVisAttributes(fWireFrameVisAtt);
  //fLogicDetector1->SetVisAttributes(fWireFrameVisAtt);
  }
  if (fDetectorGeometry == 4){
  // Hexagonal Prism Detector Array 
   
  fSolidHDetector = new G4Polyhedra{"BGO Detector", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidDetector1 = new G4Tubs("Detector", 0., 0.5*fTotalDetectorDiameter1, 0.5*fTotalDetectorLength1, 0.*deg, 360.*deg);
  fLogicDetector = new G4LogicalVolume(fSolidHDetector,fDetectorMaterial, "BGO Detector");
  fLogicDetector1 = new G4LogicalVolume(fSolidDetector1,fDetectorMaterial, "LaBr3 Detector");
  G4AssemblyVolume* assemblyDetector = new G4AssemblyVolume();
  // Back - 24, 23, 22, 20, 18, 16, 14, 12, 26, 28, 30
  // 10.059 cm forward in Z for BGO detector positions, rotm180
  //1
  P.setX(-4.96*cm); P.setY(14.78*cm); P.setZ(1.8*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //P.setX(-4.96*cm); P.setY(14.78*cm); P.setZ(1.8*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  //2
  P.setX(-10.08*cm); P.setY(11.83*cm); P.setZ(9*cm);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(-10.08*cm); P.setY(11.83*cm); P.setZ(9*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////3
  //P.setX(4.96*cm); P.setY(11.83*cm); P.setZ(2.299*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(4.96*cm); P.setY(11.83*cm); P.setZ(2.299*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////4
  //P.setX(10.08*cm); P.setY(8.87*cm); P.setZ(9*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(10.08*cm); P.setY(8.87*cm); P.setZ(9*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////5
  //P.setX(7.68*cm); P.setY(2.96*cm); P.setZ(9*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(7.68*cm); P.setY(2.96*cm); P.setZ(9*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////6
  //P.setX(7.68*cm); P.setY(-2.96*cm); P.setZ(9*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(7.68*cm); P.setY(-2.96*cm); P.setZ(9*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////7
  //P.setX(10.08*cm); P.setY(-8.87*cm); P.setZ(9*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(10.08*cm); P.setY(-8.87*cm); P.setZ(9*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////8
  //P.setX(-10.08*cm); P.setY(-11.83*cm); P.setZ(9*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(-10.08*cm); P.setY(-11.83*cm); P.setZ(9*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////9
  //P.setX(4.96*cm); P.setY(-11.83*cm); P.setZ(9*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(4.96*cm); P.setY(-11.83*cm); P.setZ(9*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////10
  //P.setX(-4.96*cm); P.setY(-14.78*cm); P.setZ(9*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(-4.96*cm); P.setY(-14.78*cm); P.setZ(9*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////11
  //P.setX(-2.56*cm); P.setY(8.87*cm); P.setZ(17.849*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(-2.56*cm); P.setY(8.87*cm); P.setZ(17.849*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////13
  //P.setX(-7.68*cm); P.setY(5.91*cm); P.setZ(17.849*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(-7.68*cm); P.setY(5.91*cm); P.setZ(17.849*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////15
  //P.setX(2.56*cm); P.setY(5.91*cm); P.setZ(17.849*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(2.56*cm); P.setY(5.91*cm); P.setZ(17.849*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////17
  //P.setX(-2.56*cm); P.setY(2.96*cm); P.setZ(17.849*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(-2.56*cm); P.setY(2.96*cm); P.setZ(17.849*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////19
  //P.setX(-7.68*cm); P.setY(0.00*cm); P.setZ(17.849*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(-7.68*cm); P.setY(0.00*cm); P.setZ(17.849*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////21
  //P.setX(2.56*cm); P.setY(0.00*cm); P.setZ(17.849*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(2.56*cm); P.setY(0.00*cm); P.setZ(17.849*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////23
  //P.setX(-2.56*cm); P.setY(-2.96*cm); P.setZ(17.849*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(-2.56*cm); P.setY(-2.96*cm); P.setZ(17.849*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////25
  //P.setX(-7.68*cm); P.setY(-5.91*cm); P.setZ(17.849*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(-7.68*cm); P.setY(-5.91*cm); P.setZ(17.849*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////27
  //P.setX(2.56*cm); P.setY(-5.91*cm); P.setZ(17.849*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(2.56*cm); P.setY(-5.91*cm); P.setZ(17.849*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////29
  //P.setX(-2.56*cm); P.setY(-8.87*cm); P.setZ(17.849*cm);
  //Tr = G4Transform3D(rotm180,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(-2.56*cm); P.setY(-8.87*cm); P.setZ(17.849*cm);
  ////Tr = G4Transform3D(rotm180,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////10.059 cm back in Z for BGO detector positions, rotm
  ////30
  //P.setX(-2.56*cm); P.setY(-8.87*cm); P.setZ(-17.849*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(-2.56*cm); P.setY(-8.87*cm); P.setZ(-17.849*cm);
  ////Tr = G4Transform3D(rotm,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////28
  //P.setX(2.56*cm); P.setY(-5.91*cm); P.setZ(-17.849*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(2.56*cm); P.setY(-5.91*cm); P.setZ(-17.849*cm);
  ////Tr = G4Transform3D(rotm,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////22
  //P.setX(2.56*cm); P.setY(0.00*cm); P.setZ(-17.849*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(2.56*cm); P.setY(0.00*cm); P.setZ(-17.849*cm);
  ////Tr = G4Transform3D(rotm,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////16
  //P.setX(2.56*cm); P.setY(5.91*cm); P.setZ(-17.849*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(2.56*cm); P.setY(5.91*cm); P.setZ(-17.849*cm);
  ////Tr = G4Transform3D(rotm,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////24
  //P.setX(-2.56*cm); P.setY(-2.96*cm); P.setZ(-17.849*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(-2.56*cm); P.setY(-2.96*cm); P.setZ(-17.849*cm);
  ////Tr = G4Transform3D(rotm,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////18
  //P.setX(-2.56*cm); P.setY(2.96*cm); P.setZ(-17.849*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(-2.56*cm); P.setY(2.96*cm); P.setZ(-17.849*cm);
  ////Tr = G4Transform3D(rotm,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////12
  //P.setX(-2.56*cm); P.setY(8.87*cm); P.setZ(-17.849*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(-2.56*cm); P.setY(8.87*cm); P.setZ(-17.849*cm);
  ////Tr = G4Transform3D(rotm,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////26
  //P.setX(-7.68*cm); P.setY(-5.91*cm); P.setZ(-17.849*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(-7.68*cm); P.setY(-5.91*cm); P.setZ(-17.849*cm);
  ////Tr = G4Transform3D(rotm,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  ////20
  //P.setX(-7.68*cm); P.setY(0.00*cm); P.setZ(-17.849*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  //////P.setX(-7.68*cm); P.setY(0.00*cm); P.setZ(-17.849*cm);
  ////Tr = G4Transform3D(rotm,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  //// 14
  //P.setX(-7.68*cm); P.setY(5.91*cm); P.setZ(-17.849*cm);
  //Tr = G4Transform3D(rotm,P);
  //assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  
  ////P.setX(-7.68*cm); P.setY(5.91*cm); P.setZ(-17.849*cm);
  ////Tr = G4Transform3D(rotm,P);
  ////assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);

  // Assembly Placement in the World
  // Rotated 270° along the x-axis to be viewed from the same perspective as the cylindrical array
  // Rotated 90° along the z-axis to be perpendicular to the gas target
  // Rotated 180° along the z-axis to place the array downside up instead of upside down
  // Coordinate transformation: (X, Y, Z) -> (Y, Z, X)
  G4ThreeVector WorldP(0,0,0);
  G4RotationMatrix worldrotm = G4RotationMatrix(0.*deg,90.*deg,270.*deg);
  G4Transform3D WorldTr = G4Transform3D(worldrotm, WorldP);
  assemblyDetector->MakeImprint(fLogicWorld, WorldTr);
  fLogicDetector->SetVisAttributes(fWireFrameVisAtt);
  fLogicDetector1->SetVisAttributes(fWireFrameVisAtt);
  }
  // PMT Optical Window
  //
  //fSolidPMTWin = new G4Tubs("PMTWin", 0., 0.5*fDetectorDiameter, 0.5*fPMTWinThickness, 0.*deg, 360.*deg);
  //fLogicPMTWin = new G4LogicalVolume(fSolidPMTWin, fPMTWinMaterial, "PMTWin");
  //fPhysiPMTWin = new G4PVPlacement(0, 
        			//G4ThreeVector(0.,0.,fZposPMTWin), 
        			//fLogicPMTWin, 
       				//"PMTWin", 
       				//fLogicDetector, 
       				//false, 
    			    	//0);
    			    	
  //fLogicPMTWin->SetVisAttributes(fGreenVisAtt);
  
  fSolidPMTWin1 = new G4Tubs("PMTWin", 0., 0.5*fDetectorDiameter1, 0.5*fPMTWinThickness1, 0.*deg, 360.*deg);
  fLogicPMTWin1 = new G4LogicalVolume(fSolidPMTWin1, fPMTWinMaterial, "PMTWin");
  fPhysiPMTWin1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMTWin1), 
        			fLogicPMTWin1, 
       				"PMTWin", 
       				fLogicDetector1, 
       				false, 
    			    	0);
    			    	
  fLogicPMTWin1->SetVisAttributes(fGreenVisAtt);

  // Crystal
  //
  fSolidHCrystal = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal = new G4LogicalVolume(fSolidHCrystal, fCrystalMaterial, "Crystal");
  fPhysiCrystal = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal, 
       				"Crystal", 
       				fLogicDetector, 
       				false, 
    			    	0);
    			    	
  fSolidCrystal1 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter1, 0.5*fDetectorLength1, 0.*deg, 360.*deg);			    	
  fLogicCrystal1 = new G4LogicalVolume(fSolidCrystal1, fCrystalMaterial1, "Crystal");
  fPhysiCrystal1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos1), 
        			fLogicCrystal1, 
       				"Crystal", 
       				fLogicDetector1, 
       				false, 
    			    	0);
    			    	
  fLogicCrystal->SetVisAttributes(fRedVisAtt);
  fLogicCrystal1->SetVisAttributes(fMagnetaVisAtt);
  
  //Gap (gap surrounding crystal and casing or a reflector as required)
  
  fSolidHGap = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap = new G4LogicalVolume(fSolidHGap, fGapMaterial, "Gap");
  fPhysiGap = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap, 
       				"Gap", 
       				fLogicDetector, 
       				false, 
    			    	0);

  //FaceGap (gap between face of crystal and casing or a reflector as required)
  //
  fSolidHFaceGap = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);//zPlane4
  fLogicFaceGap = new G4LogicalVolume(fSolidHFaceGap, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap, 
       				"FaceGap", 
       				fLogicDetector, 
       				false, 
    			    	0);
  			    	
  fSolidFaceGap1 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter1, 0.5*fGapFaceThickness1, 0.*deg, 360.*deg);			    	
  fLogicFaceGap1 = new G4LogicalVolume(fSolidFaceGap1, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap1), 
        			fLogicFaceGap1,
       				"FaceGap", 
       				fLogicDetector1, 
       				false, 
    			    	0);

  fLogicGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap1->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap1->SetVisAttributes(fCyanVisAtt);
  
  //Aluminum Casing (casing surrounding crystal)
  //
  fSolidHAlCase = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase = new G4LogicalVolume(fSolidHAlCase, fAlCaseMaterial, "AlCase");
  fPhysiAlCase = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase, 
       				"AlCase", 
       				fLogicDetector, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase1 = new G4Tubs("AlCase", 0.5*fDetectorDiameter1, 0.5*fAlCaseDiameter1, 0.5*fAlCaseLength1, 0.*deg, 360.*deg); 			    				    	
  fLogicAlCase1 = new G4LogicalVolume(fSolidAlCase1, fAlCaseMaterial, "AlCase");
  fPhysiAlCase1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase1), 
        			fLogicAlCase1, 
       				"AlCase", 
       				fLogicDetector1, 
       				false, 
    			    	0);

  //Aluminum Face Casing (casing covering face of crystal)
  //
  fSolidHFaceAlCase = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase = new G4LogicalVolume(fSolidHFaceAlCase, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase, 
       				"FaceAlCase", 
       				fLogicDetector, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase1 = new G4Tubs("FaceAlCase", 0., 0.5*fDetectorDiameter1, 0.5*fAlFaceThickness1, 0.*deg, 360.*deg); 
  fLogicFaceAlCase1 = new G4LogicalVolume(fSolidFaceAlCase1, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase1), 
        			fLogicFaceAlCase1, 
       				"FaceAlCase", 
       				fLogicDetector1, 
       				false, 
    			    	0);
    			    	
  fLogicAlCase->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase1->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase1->SetVisAttributes(fGreyVisAtt);
  
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
  
  fSolidPMTJ = new G4Tubs("Vacuum", 0., 0.5*fPMTDiameter - 0.16*cm, 0.5*fPMTLength - 0.508*cm, 0.*deg, 360.*deg);		    	
  fLogicPMTJ = new G4LogicalVolume(fSolidPMTJ, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ, 
       				"Vacuum", 
       				fLogicPMT, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL = new G4Tubs("Photocathode", 0., 0.5*fPMTDiameter - 0.16*cm, 0.254*cm, 0.*deg, 360.*deg);		    	
  fLogicPMTL = new G4LogicalVolume(fSolidPMTL, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,9.9955*cm), 
        			fLogicPMTL, 
       				"Photocathode", 
       				fLogicPMT, 
       				false, 
    			    	0);
    			    	
    			    	
  fSolidPMT1 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter1, 0.5*fPMTLength1, 0.*deg, 360.*deg); 			    	
  fLogicPMT1 = new G4LogicalVolume(fSolidPMT1, fPMTMaterial,"PMT");
  fPhysiPMT1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT1), 
        			fLogicPMT1, 
       				"PMT", 
       				fLogicDetector1, 
       				false, 
    			    	0);
    			    	
  fSolidPMTI = new G4Tubs("Vacuum", 0., 0.5*fPMTDiameter1 - 0.25*cm, 0.5*fPMTLength1 - 0.508*cm, 0.*deg, 360.*deg);			    	
  fLogicPMTI = new G4LogicalVolume(fSolidPMTI, fPMTIntMaterial,"Vacuum");
  fPhysiPMTI = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTI, 
       				"Vacuum", 
       				fLogicPMT1, 
       				false, 
    			    	0);
  
  fSolidPMTK = new G4Tubs("Photocathode", 0., 0.5*fPMTDiameter1 - 0.25*cm, 0.254*cm, 0.*deg, 360.*deg);			    	
  fLogicPMTK = new G4LogicalVolume(fSolidPMTK, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTK = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,2.846*cm), 
        			fLogicPMTK, 
       				"Photocathode", 
       				fLogicPMT1, 
       				false, 
    			    	0);
    			    	
  fLogicPMTK->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL->SetVisAttributes(fYellowVisAtt);
  }
  if (fDetectorGeometry == 3 || fDetectorGeometry == 4){
 
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
  fLogicOutTBox->SetVisAttributes(fWireFrameVisAtt);
  // DRAGON Outer Gas Target Box
  // Inner Wall
  rotm  = G4RotationMatrix(0,0,0);     
  position = G4ThreeVector(0.,0.,0.);
  transform = G4Transform3D(rotm,position);	
  fSolidInTBox = new G4Box("InTBox", 8.256*cm, 2.064*cm, 12.542*cm);
  fLogicInTBox = new G4LogicalVolume(fSolidInTBox,fGasTargetMaterial,"InTBox");
  fPhysiInTBox = new G4PVPlacement(transform,
        			fLogicInTBox, 
       				"InTBox", 
       				fLogicOutTBox, 
       				false, 
    			    	0,
				    false);
  fLogicInTBox->SetVisAttributes(fWireFrameVisAtt);				    
  // Entrance Hole on Target Box				    
  rotm  = G4RotationMatrix(0,0,0);
  rotm.rotateY(90*deg);    
  position = G4ThreeVector(-8.4145*cm,0.,-9.684*cm);
  transform = G4Transform3D(rotm,position);  
  fSolidEntHoleBox = new G4Tubs("EntHoleBox", 0., 0.95*cm, 0.1585*cm, 0.*deg, 360.*deg);
  fLogicEntHoleBox = new G4LogicalVolume(fSolidEntHoleBox, fGasTargetMaterial, "EntHoleBox");
  fPhysiEntHoleBox = new G4PVPlacement(transform,
        			fLogicEntHoleBox, 
       				"EntHoleBox", 
       				fLogicOutTBox, 
       				false, 
    			    	0);
  fLogicEntHoleBox->SetVisAttributes(fWireFrameVisAtt);			    	
  // Exit Hole on Target Box    
  position = G4ThreeVector(8.4145*cm,0.,-9.684*cm);
  transform = G4Transform3D(rotm,position);	 
  fSolidExitHoleBox = new G4Tubs("ExitHoleBox", 0., 0.95*cm, 0.1585*cm, 0.*deg, 360.*deg);
  fLogicExitHoleBox = new G4LogicalVolume(fSolidExitHoleBox, fGasTargetMaterial, "ExitHoleBox");
  fPhysiExitHoleBox = new G4PVPlacement(transform,
        			fLogicExitHoleBox, 
       				"ExitHoleBox", 
       				fLogicOutTBox, 
       				false, 
    			    	0);   
  fLogicExitHoleBox->SetVisAttributes(fWireFrameVisAtt);
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
  fLogicOutCell->SetVisAttributes(fWireFrameVisAtt);

  // DRAGON Inner Gas Cell
  // Inner Wall			    
  fSolidInCell = new G4Trd("InCell",6.208*cm,1.717*cm,1.588*cm,1.588*cm,3.890*cm);
  fLogicInCell = new G4LogicalVolume(fSolidInCell,fInnerTrapezoidMaterial, "InCell");
  fPhysiInCell = new G4PVPlacement(0,			//no rotation
            		      	G4ThreeVector(),        //no transform, at (0,0,0)
        			fLogicInCell, 
       				"InCell", 
       				fLogicOutCell, 
       				false, 
    			    	0,
				    false);
  fLogicInCell->SetVisAttributes(fWireFrameVisAtt);
  // Entrance Collimator
  // Entrance Collimator Hole
  rotm  = G4RotationMatrix(0,0,0);     
  rotm.rotateY(-60*deg);
  // position = G4ThreeVector(5.306*cm,0.,-2.008*cm);
  position = G4ThreeVector(-5.351*cm,0.,-2.087*cm);
  transform = G4Transform3D(rotm,position); 
  fSolidEntHole = new G4Tubs("EntHole", 0.*cm, 0.80*cm, 0.1585*cm, 0.*deg, 360.*deg);
  fLogicEntHole = new G4LogicalVolume(fSolidEntHole, fInnerTrapezoidMaterial, "EntHole");
  fPhysiEntHole = new G4PVPlacement(transform,
        			fLogicEntHole, 
       				"EntHole", 
       				fLogicOutCell, 
       				false, 
    			    	0);   
  fLogicEntHole->SetVisAttributes(fWireFrameVisAtt);	
  
  // Entrance Collimator Thin Metal Disk without hole
  rotm  = G4RotationMatrix(0,0,0);     
  //rotm.rotateY(60*deg);
  //position = G4ThreeVector(0.,0.,0.1196*cm);
  position = G4ThreeVector(0.,0.,0.1195*cm);
  //position = G4ThreeVector(0.,0.,0.);
  transform = G4Transform3D(rotm,position); 
  fSolidEntDisk = new G4Tubs("EntDisk", 0.*cm, 0.79*cm, 0.039*cm, 0.*deg, 360.*deg);
  fLogicEntDisk = new G4LogicalVolume(fSolidEntDisk, fAlCaseMaterial, "EntDisk");
  fPhysiEntDisk = new G4PVPlacement(transform,
        			fLogicEntDisk, 
       				"EntDisk", 
       				fLogicEntHole, 
       				false, 
    			    	0); 
  fLogicEntDisk->SetVisAttributes(fWireFrameVisAtt); 
  
   // Entrance Collimator 6mm Hole
  rotm  = G4RotationMatrix(0,0,0);    
  rotm.rotateY(30*deg);
  rotm.rotateZ(180*deg);
  position = G4ThreeVector(0.,0.,0.);
  transform = G4Transform3D(rotm,position); 
  //G4ThreeVector pLowNorm = G4ThreeVector(1.0,0.,-1.73);
  //G4ThreeVector pHighNorm = G4ThreeVector(-1.0,0.,1.73);
  G4ThreeVector pLowNorm = G4ThreeVector(1.0,0.,-sqrt(3));
  G4ThreeVector pHighNorm = G4ThreeVector(-1.0,0.,sqrt(3));
  fSolidEntCol = new G4CutTubs("EntCol", 0.*cm, 0.3*cm, (0.039/cos(pi/6))*cm, 0.*deg, 360.*deg,pLowNorm,pHighNorm);
  fLogicEntCol = new G4LogicalVolume(fSolidEntCol, fInnerTrapezoidMaterial, "EntCol");
  fPhysiEntCol = new G4PVPlacement(transform,
				fLogicEntCol, 
       				"EntCol", 
       				fLogicEntDisk, 
       				false, 
    			    	0); 
  fLogicEntCol->SetVisAttributes(fWireFrameVisAtt); 			    	  
  
  // Exit Collimator
  // Exit Collimator Hole
  rotm  = G4RotationMatrix(0,0,0);     
  rotm.rotateY(60*deg);
//  position = G4ThreeVector(-5.306*cm,0.,-2.008*cm);
  position = G4ThreeVector(5.351*cm,0.,-2.087*cm);
  transform = G4Transform3D(rotm,position); 
  //fSolidExitHole = new G4Tubs("ExitHole", 0., 0.80*cm, 0.1585*cm, 0.*deg, 360.*deg);
  fSolidExitHole = new G4Tubs("ExitHole", 0.*cm, 0.80*cm, 0.1585*cm, 0.*deg, 360.*deg);
  fLogicExitHole = new G4LogicalVolume(fSolidExitHole, fInnerTrapezoidMaterial, "ExitHole");
  fPhysiExitHole = new G4PVPlacement(transform,
        			fLogicExitHole, 
       				"ExitHole", 
       				fLogicOutCell, 
       				false, 
    			    	0);   
  fLogicExitHole->SetVisAttributes(fWireFrameVisAtt); 	
  
  // Exit Collimator Thin Metal Disk without hole
  rotm  = G4RotationMatrix(0,0,0);     
  //rotm.rotateY(60*deg);
  //position = G4ThreeVector(0.,0.,0.1196*cm);
  position = G4ThreeVector(0.,0.,0.1195*cm);
  //position = G4ThreeVector(0.,0.,0.);
  transform = G4Transform3D(rotm,position); 
  fSolidExitDisk = new G4Tubs("ExitDisk", 0.*cm, 0.79*cm, 0.039*cm, 0.*deg, 360.*deg);
  fLogicExitDisk = new G4LogicalVolume(fSolidExitDisk, fAlCaseMaterial, "ExitDisk");
  fPhysiExitDisk = new G4PVPlacement(transform,
        			fLogicExitDisk, 
       				"ExitDisk", 
       				fLogicExitHole, 
       				false, 
    			    	0); 
  fLogicExitDisk->SetVisAttributes(fWireFrameVisAtt); 	
  
  // Exit Collimator 6mm/8mm Hole
  rotm  = G4RotationMatrix(0,0,0);     
  rotm.rotateY(30*deg);
  position = G4ThreeVector(0.,0.,0.);
  transform = G4Transform3D(rotm,position); 
  //pLowNorm = G4ThreeVector(1.0,0.,-1.73);
  //pHighNorm = G4ThreeVector(-1.0,0.,1.73);
  pLowNorm = G4ThreeVector(1.0,0.,-sqrt(3));
  pHighNorm = G4ThreeVector(-1.0,0.,sqrt(3));
  //fSolidExitCol = new G4CutTubs("ExitCol", 0.*cm, 0.3*cm, 0.046*cm, 0.*deg, 360.*deg,pLowNorm,pHighNorm);
  fSolidExitCol = new G4CutTubs("ExitCol", 0.*cm, 0.4*cm, (0.039/cos(pi/6))*cm, 0.*deg, 360.*deg,pLowNorm,pHighNorm);
  fLogicExitCol = new G4LogicalVolume(fSolidExitCol, fInnerTrapezoidMaterial, "ExitCol");
  fPhysiExitCol = new G4PVPlacement(transform,
				fLogicExitCol, 
       				"ExitCol", 
       				fLogicExitDisk, 
       				false, 
    			    	0); 
  fLogicExitCol->SetVisAttributes(fWireFrameVisAtt);		    	  
  
  }
  PrintCalorParameters();
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
  G4cout << "\n" << "The objects' temperature is: " << Temperature << " kelvin" << G4endl;
  G4cout << "\n" << "The objects' pressure is: " << G4BestUnit(Pressure, "Pressure") << G4endl;
  G4cout << "\n" << "The red crystal is made of" << fCrystalMaterial << G4endl;
  G4cout << "\n" << "The magneta crystal is made of" << fCrystalMaterial1 << G4endl;
  G4cout << "\n" << "The gap is made of " << fGapMaterial << G4endl;
  G4cout << "\n" << "The face of the gap is made of " << fFaceGapMaterial << G4endl;
  G4cout << "\n" << "The aluminum case is made of " << fAlCaseMaterial << G4endl;
  G4cout << "\n" << "The face of the aluminum case is made of " << fFaceAlCaseMaterial << G4endl;
  G4cout << "\n" << "The lead case, left lead ring, and connector are made of " << fPbCaseMaterial << G4endl;
  G4cout << "\n" << "The lead collar is made of " << fPbCollarMaterial << G4endl;
  G4cout << "\n" << "The photomultiplier tube is made of" << fPMTMaterial << G4endl;
  G4cout << "\n" << "The photomultiplier tube interior is made of" << fPMTIntMaterial << G4endl;
  G4cout << "\n" << "The photocathode is made of" << fPhotoCathodeMaterial << G4endl;
  G4cout << "\n" << "The photomultiplier tube window is made of" << fPMTWinMaterial << G4endl;
  G4cout << "\n" << "Total BGO Detector Diameter: " << G4BestUnit(fTotalDetectorDiameter, "Length") << G4endl;
  G4cout << "\n" << "BGO Scintillator Diameter: " << G4BestUnit(fDetectorDiameter, "Length") << G4endl;
  G4cout << "\n" << "Total BGO Detector Length: " << G4BestUnit(fTotalDetectorLength, "Length") << G4endl;
  G4cout << "\n" << "BGO Scintillator Length: " << G4BestUnit(fDetectorLength, "Length") << G4endl;
  G4cout << "\n" << "Gap Thickness: " << 0.1*fGapThickness << " cm" << G4endl;
  G4cout << "\n" << "Al Thickness: " << 0.1*fAlCaseThickness << " cm" << G4endl;
  G4cout << "\n" << "Pb Thickness: " << 0.1*fPbCaseThickness << " cm" << G4endl;
  G4cout << "\n" << "PMT Diameter: " << G4BestUnit(fPMTDiameter, "Length") << G4endl;
  G4cout << "\n" << "PMT Length: " << G4BestUnit(fPMTLength, "Length") << G4endl;
  G4cout << "\n" << "Total LaBr3 Detector Diameter: " << G4BestUnit(fTotalDetectorDiameter1, "Length") << G4endl;
  G4cout << "\n" << "LaBr3 Scintillator Diameter: " << G4BestUnit(fDetectorDiameter1, "Length") << G4endl;
  G4cout << "\n" << "Total LaBr3 Detector Length: " << G4BestUnit(fTotalDetectorLength1, "Length") << G4endl;
  G4cout << "\n" << "LaBr3 Scintillator Length: " << G4BestUnit(fDetectorLength1, "Length") << G4endl;
  G4cout << "\n" << "Gap Thickness: " << 0.1*fGapThickness1 << " cm" << G4endl;
  G4cout << "\n" << "Al Thickness: " << 0.1*fAlCaseThickness1 << " cm" << G4endl;
  G4cout << "\n" << "PMT Diameter: " << G4BestUnit(fPMTDiameter1, "Length") << G4endl;
  G4cout << "\n" << "PMT Length: " << G4BestUnit(fPMTLength1, "Length") << G4endl;
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

void DetectorConstruction::SetTemperature(G4double val)
{
  Temperature = val;
}

void DetectorConstruction::SetPressure(G4double val)
{
  Pressure = val;
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
