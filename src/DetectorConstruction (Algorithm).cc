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
fSolidPMT(0),fLogicPMT(0),fPhysiPMT(0),
fSolidPMTWin(0),fLogicPMTWin(0),fPhysiPMTWin(0)
{
  // Setting the default for the world for the Sim.
  fWorldSizeX = 100*cm;
  fWorldSizeYZ = 100*cm;

  // Set BGO and LaBr3:Ce Detector array defaults according to the Geant3 dimensions.
  // Adjust fPMTDiameter to change the detector size. 
  fPMTDiameter = 5.9*cm;
  fPMTLength = (20.499*fPMTDiameter)/5.9;
  fPMTWinThickness = (0.508*fPMTDiameter)/5.9;
  fDetectorLength = (5.08*fPMTDiameter)/5.9;
  fDetectorDiameter = (5.08*fPMTDiameter)/5.9;
  fGapThickness = (0.0355*fPMTDiameter)/5.9;
  fGapFaceThickness = (0.3175*fPMTDiameter)/5.9;
  fAlCaseThickness = (0.0635*fPMTDiameter)/5.9;
  fAlFaceThickness = (0.0635*fPMTDiameter)/5.9;
  fVacuumDiameter = fPMTDiameter - (0.32*fPMTDiameter)/5.9;
  fVacuumLength = fPMTLength - (1.016*fPMTDiameter)/5.9;
  fPhotoCathodeDiameter = fPMTDiameter - (0.32*fPMTDiameter)/5.9;
  fPhotoCathodeLength = (0.508*fPMTDiameter)/5.9;
  fAirGap = 0*cm;
  
  // Set Single (and the timing array) LaBr3:Ce detector defaults.
  // Adjust fAlCaseDiameter1 to change the detector size. 
  fAlCaseDiameter1 = 6.1*cm;
  fPMTDiameter1 = (5.1*fAlCaseDiameter1)/5.18;
  fPMTLength1 = (6.2*fAlCaseDiameter1)/5.18;
  fPMTWinThickness1 = (0.2*fAlCaseDiameter1)/5.18;
  fDetectorLength1 = (7.62*fAlCaseDiameter1)/5.18;
  fDetectorDiameter1 = (5*fAlCaseDiameter1)/5.18;
  fGapFaceThickness1 = (0.42*fAlCaseDiameter1)/5.18;
  fAlCaseThickness1 = (0.05*fAlCaseDiameter1)/5.18;
  fAlFaceThickness1 = (0.05*fAlCaseDiameter1)/5.18;
  fAlCaseLength1 = fDetectorLength1 + fPMTWinThickness1 + fGapFaceThickness1 + fAlFaceThickness1;
  
  fOutTBoxX = 8.573*cm;
  fOutTBoxY = 2.381*cm;
  fOutTBoxZ = 12.859*cm;
  Temperature = 298.15*kelvin;
  Pressure = 101325*pascal;
  
  ComputeCalorParameters();
  
  // Materials  
  DefineMaterials();
  // Default Materials
  fWorldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Air");
  fDetectorMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Air");
  fCrystalMaterial = G4NistManager::Instance()->FindOrBuildMaterial("BGO");
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
  // Set Default Gap Material i.e. no reflector
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
  // fWireFrameVisAtt = new G4VisAttributes(G4Colour(1.,1.,1.));
  fWireFrameVisAtt = new G4VisAttributes(G4Colour(0,0,0,1));
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
  // This function illustrates the possible ways to define materials.
 
  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;  

  G4int ncomponents, natoms;
  G4double fractionmass;
  
  // Define Elements
  // Atomic masses are taken from the National Library of Medicine and rounded to 2 decimal places.
  G4Element* H  = new G4Element("Hydrogen",symbol="H", z=1, a= 1.01*g/mole);
  G4Element* B = new G4Element("Boron",symbol="B", z=5, a= 10.81*g/mole); 
  G4Element* C  = new G4Element("Carbon",  symbol="C", z=6, a= 12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N", z=7, a= 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",  symbol="O", z=8, a= 16.00*g/mole);
  G4Element* F = new G4Element("Flourine", symbol="F", z=9, a= 19.00*g/mole);
  G4Element* Na = new G4Element("Sodium",  symbol="Na", z=11, a= 22.99*g/mole);
  G4Element* Mg = new G4Element("Magnesium", symbol="Mg", z=12, a= 24.31*g/mole);
  G4Element* Al = new G4Element("Aluminium", symbol="Al", z=13, a=26.98*g/mole);
  G4Element* Si = new G4Element("Silicon", symbol="Si", z=14, a= 28.09*g/mole);
  G4Element* Ar = new G4Element("Argon", symbol="Ar", z=18, a= 39.90*g/mole);
  G4Element* K = new G4Element("Potassium", symbol="K", z=19, a= 39.10*g/mole);
  G4Element* Ca = new G4Element("Calcium", symbol="Ca", z=20, a= 40.08*g/mole);
  G4Element* Ge = new G4Element("Germanium", symbol="Ge", z=32, a= 72.63*g/mole);
  G4Element* Br = new G4Element("Bromine",  symbol="Br", z=35, a= 79.90*g/mole);
  G4Element* Sr = new G4Element("Strontium",  symbol="Sr", z=38, a= 87.62*g/mole);
  G4Element* Y = new G4Element("Yttrium", symbol="Y", z=39, a= 88.91*g/mole);
  G4Element* Sb = new G4Element("Antimony", symbol="Sb", z=51, a= 121.76*g/mole);
  G4Element* I  = new G4Element("Iodine",  symbol="I", z=53, a= 126.90*g/mole);
  G4Element* Cs = new G4Element("Cesium", symbol="Cs", z=55, a= 132.91*g/mole);
  G4Element* La = new G4Element("Lanthanum", symbol="La", z=57, a= 138.91*g/mole);
  G4Element* Ce = new G4Element("Cerium", symbol="Ce", z=58, a= 140.12*g/mole);
  G4Element* Lu = new G4Element("Lutetium", symbol="Lu", z=71, a= 174.97*g/mole);
  G4Element* Bi = new G4Element("Bismuth", symbol="Bi", z=83, a= 208.98*g/mole);
 
  // Define Simple Materials
  new G4Material("H2Liq"    , z= 1, a= 1.01*g/mole, density= 70.8*mg/cm3);
  new G4Material("Beryllium", z= 4, a= 9.01*g/mole, density= 1.848*g/cm3);
  G4Material* Silicon = new G4Material("Silicon"  , z=14, a=28.09*g/mole, density= 2.330*g/cm3, kStateSolid, Temperature, Pressure);
  Silicon->GetIonisation()->SetMeanExcitationEnergy(173*eV);
  new G4Material("Iron",     z=26, a= 55.84*g/mole, density= 7.874*g/cm3);
  G4Material* Copper = new G4Material("Copper",   z=29, a= 63.55*g/mole, density= 8.960*g/cm3);
  Copper->GetIonisation()->SetMeanExcitationEnergy(322.*eV);
  new G4Material("Germanium",z=32, a= 72.61*g/mole, density= 5.323*g/cm3);
  G4Material* Silver = new G4Material("Silver",   z=47, a=107.87*g/mole, density= 10.50*g/cm3);
  Silver->GetIonisation()->SetMeanExcitationEnergy(470.*eV);
  new G4Material("Tungsten", z=74, a=183.84*g/mole, density= 19.30*g/cm3);
  new G4Material("Gold",     z=79, a=196.97*g/mole, density= 19.32*g/cm3);
  
  G4Material* Aluminium = new G4Material("Aluminium", z=13, a=26.98*g/mole, density= 2.699*g/cm3, kStateSolid, Temperature, Pressure);
  Aluminium->GetIonisation()->SetMeanExcitationEnergy(166*eV);
  G4Material* Lead = new G4Material("Lead", z=82, a=207*g/mole, density= 11.35*g/cm3, kStateSolid, Temperature, Pressure);
  Lead->GetIonisation()->SetMeanExcitationEnergy(823*eV);

  // Define a material from elements. Case 1: Chemical Molecule
  // G4_SODIUM_IODIDE
  G4Material* NaI = new G4Material("NaI", density= 3.667*g/cm3, ncomponents=2, kStateSolid, Temperature, Pressure);
  NaI->AddElement(Na, natoms=1);
  NaI->AddElement(I , natoms=1);
  NaI->GetIonisation()->SetMeanExcitationEnergy(452*eV);

  // Define a material from elements. Case 2: Mixture by fractional mass.
 
  G4Material* LaBr3 = new G4Material("LaBr3", density= 5.06*g/cm3, ncomponents=2, kStateSolid, Temperature, Pressure);
  LaBr3->AddElement(La, natoms=1);
  LaBr3->AddElement(Br, natoms=3);
  LaBr3->GetIonisation()->SetMeanExcitationEnergy(401*eV);
  
  G4Material* LaBr3Ce = new G4Material("LaBr3Ce", density= 5.08*g/cm3, ncomponents=2, kStateSolid, Temperature, Pressure);
  LaBr3Ce->AddMaterial(LaBr3, fractionmass=0.980894429);
  LaBr3Ce->AddElement(Ce, fractionmass=0.01910557);
  LaBr3Ce->GetIonisation()->SetMeanExcitationEnergy(403.3*eV);
  
  // G4Material* LaBr3Ce = new G4Material("LaBr3Ce", density= 5.08*g/cm3, ncomponents=2);
  // LaBr3Ce->AddMaterial(LaBr3, fractionmass=0.980894429);
  // LaBr3Ce->AddElement(Ce, fractionmass=0.01910557);
  // LaBr3Ce->GetIonisation()->SetMeanExcitationEnergy(403.298315*eV);
  
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
  
  // Only material where the density is estimated. 
  G4Material* Bialkali = new G4Material("Bialkali", density = 3.74*g/cm3, ncomponents=3, kStateSolid, Temperature, Pressure);
  Bialkali->AddElement(Cs, natoms=1);
  Bialkali->AddElement(K, natoms=1);
  Bialkali->AddElement(Sb, natoms=1);
  Bialkali->GetIonisation()->SetMeanExcitationEnergy(447.9*eV);
  
  // Known density of polydimethylsiloxane. 
  G4Material* Optical = new G4Material("Optical", density = 0.965*g/cm3, ncomponents=4, kStateSolid, Temperature, Pressure);
  Optical->AddElement(C, natoms=2);
  Optical->AddElement(H, natoms=6);
  Optical->AddElement(O, natoms=1);
  Optical->AddElement(Si, natoms=1);
  Bialkali->GetIonisation()->SetMeanExcitationEnergy(113.8*eV);
  
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
  
  // Standard temperature and pressure. 
  // G4Material* Air = new G4Material("Air", density= 0.00120479*g/cm3, ncomponents=4);
  // Air->AddElement(C, fractionmass=0.000124);
  // Air->AddElement(N, fractionmass=0.755268);
  // Air->AddElement(O, fractionmass=0.231781);
  // Air->AddElement(Ar, fractionmass=0.012827);
  // Air->GetIonisation()->SetMeanExcitationEnergy(85.7*eV);
  
  // G4_MAGNESIUM_OXIDE
  G4Material* MgO = new G4Material("MgO", density= 3.58*g/cm3, ncomponents=2, kStateGas, Temperature, Pressure);
  MgO->AddElement(Mg, natoms=1);
  MgO->AddElement(O, natoms=1);
  MgO->GetIonisation()->SetMeanExcitationEnergy(143.8*eV);
  
  // G4_H for 1 Torr.
  G4Material* Hydrogen = new G4Material("Hydrogen Gas", density= 0.000083748*g/cm3, ncomponents=1, kStateGas, Temperature, (1/760)*atmosphere);
  Hydrogen->AddElement(H, natoms=2);
  Hydrogen->GetIonisation()->SetMeanExcitationEnergy(19.2*eV);
  
  // G4_H for 2 x 10^-6 Torr.
  G4Material* Hydrogen2 = new G4Material("Hydrogen Gas 2", density= 0.000083748*g/cm3, ncomponents=1, kStateGas, Temperature, (0.000002/760)*atmosphere);
  Hydrogen2->AddElement(H, natoms=2);
  Hydrogen2->GetIonisation()->SetMeanExcitationEnergy(19.2*eV);
  
  // Artificial Vacuum
  new G4Material("Artificial Vacuum", z=1, a=1.01*g/mole, universe_mean_density, kStateGas, 298.15*kelvin, (3.e-18)*pascal);
  
  // Example of vacuum
  // from PhysicalConstants.h.
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
  // Compute derived parameters of the calorimeter.
  // BGO and LaBr3:Ce Array
  fGapLength = fDetectorLength + fGapFaceThickness;
  fGapDiameter = fDetectorDiameter + 2*fGapThickness;
  fAlCaseLength = fGapLength + fAlFaceThickness;
  fAlCaseDiameter = fGapDiameter + 2*fAlCaseThickness;
  fPbCaseDiameter = fAlCaseDiameter;
  fTotalDetectorLength = fDetectorLength + fPMTLength + fGapFaceThickness + fAlFaceThickness;
  
  // GRIFFIN LaBr3:Ce
  fVacuumDiameter1 = fPMTDiameter1 - (0.5*fAlCaseDiameter1)/5.18;
  fVacuumLength1 =  fPMTLength1 - (1.016*fAlCaseDiameter1)/5.18;
  fPhotoCathodeDiameter1 = fPMTDiameter1 - (0.5*fAlCaseDiameter1)/5.18;
  fPhotoCathodeLength1 = 0.508*fAlCaseDiameter1/5.18;
  fZposPhotoCathode1 = 0.5*fPMTLength1 - 0.5*fPhotoCathodeLength1;
  fTotalDetectorLength1 = fPMTWinThickness1 + fDetectorLength1 + fPMTLength1 + fGapFaceThickness1 + fAlFaceThickness1;
  
  // This has been commented out to allow fPMTDiameter and fAlCaseDiameter1 to be adjusted by DetectorMessenger.cc. 
  //if (fPbCaseDiameter>fPMTDiameter){
	  //fPMTDiameter = fPbCaseDiameter;
	  //G4cout << "\n" << "\n" 
	  //<< "WARNING: PMT Diameter has been adjusted to compensate for large crystal/lead/aluminum size."  
	  //<< G4endl;
	  //fTotalDetectorDiameter = fPMTDiameter;
	  //fTotalDetectorDiameter1 = fAlCaseDiameter1;
  //}
  //else
  fTotalDetectorDiameter = fPMTDiameter;
  fTotalDetectorDiameter1 = fAlCaseDiameter1;
	
  // BGO Component Positions 
  fZposFaceAlCase = 0.5*fTotalDetectorLength - 0.5*fAlFaceThickness;
  fZposFaceGap = 0.5*fTotalDetectorLength - fAlFaceThickness - 0.5*fGapFaceThickness ;
  fZpos = 0.5*fTotalDetectorLength - fGapFaceThickness - fAlFaceThickness - 0.5*fDetectorLength;
  fZposPMTWin = 0.5*fTotalDetectorLength - fAlFaceThickness - fGapFaceThickness - fDetectorLength - 0.5*fPMTWinThickness;
  fZposAlCase = 0.5*fTotalDetectorLength - 0.5*fAlCaseLength;
  fZposGap = 0.5*fTotalDetectorLength - fAlFaceThickness - 0.5*fGapLength ;
  fZposPMT = 0.5*fTotalDetectorLength - fAlFaceThickness - fGapFaceThickness - fDetectorLength - 0.5*fPMTLength;
  fZposPhotoCathode = 0.5*fPMTLength - 0.5*fPhotoCathodeLength;
  
  // LaBr3:Ce Component Positions 
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
  // Cleanup old geometry.
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
 
  // Complete the Calor parameters definition.
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
  // Detector Holder Volume for Components
  //
  //  If rotation required  
  //  G4RotationMatrix rotm  = G4RotationMatrix(0,90*deg,-90*deg);     
  //  No Rotation Now
  G4RotationMatrix rotm  = G4RotationMatrix(0,0,0.*deg);    
  G4ThreeVector position = G4ThreeVector(0.,0.,0.);
  G4Transform3D transform = G4Transform3D(rotm,position);
  G4ThreeVector P;
  G4Transform3D Tr;
  // The detectors on each side of the array face each other with this rotation matrix.
  G4RotationMatrix rotm180 = G4RotationMatrix(0,180.*deg,0); 
  
if (fDetectorGeometry == 1 || fDetectorGeometry == 2 || fDetectorGeometry == 3 || fDetectorGeometry == 4 || fDetectorGeometry == 5){

if (fDetectorGeometry == 1 || fDetectorGeometry == 3){  
  // Single Cylindrical Detector 
  if (fDetectorGeometry == 1){
	  
  fSolidDetector = new G4Tubs("LaBr3:Ce Detector 1", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector1 = new G4Tubs("GRIFFIN LaBr3:Ce Detector", 0., 0.5*fTotalDetectorDiameter1, 0.5*fTotalDetectorLength1, 0.*deg, 360.*deg);
  fSolidDetector2 = new G4Tubs("LaBr3:Ce Detector 2", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector3 = new G4Tubs("LaBr3:Ce Detector 3", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector4 = new G4Tubs("LaBr3:Ce Detector 4", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector5 = new G4Tubs("LaBr3:Ce Detector 5", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector6 = new G4Tubs("LaBr3:Ce Detector 6", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector7 = new G4Tubs("LaBr3:Ce Detector 7", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector8 = new G4Tubs("LaBr3:Ce Detector 8", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector9 = new G4Tubs("LaBr3:Ce Detector 9", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector10 = new G4Tubs("LaBr3:Ce Detector 10", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector11 = new G4Tubs("LaBr3:Ce Detector 11", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector12 = new G4Tubs("LaBr3:Ce Detector 12", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector13 = new G4Tubs("LaBr3:Ce Detector 13", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector14 = new G4Tubs("LaBr3:Ce Detector 14", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector15 = new G4Tubs("LaBr3:Ce Detector 15", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector16 = new G4Tubs("LaBr3:Ce Detector 16", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector17 = new G4Tubs("LaBr3:Ce Detector 17", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector18 = new G4Tubs("LaBr3:Ce Detector 18", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector19 = new G4Tubs("LaBr3:Ce Detector 19", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector20 = new G4Tubs("LaBr3:Ce Detector 20", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector21 = new G4Tubs("LaBr3:Ce Detector 21", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector22 = new G4Tubs("LaBr3:Ce Detector 22", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector23 = new G4Tubs("LaBr3:Ce Detector 23", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector24 = new G4Tubs("LaBr3:Ce Detector 24", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector25 = new G4Tubs("LaBr3:Ce Detector 25", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector26 = new G4Tubs("LaBr3:Ce Detector 26", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector27 = new G4Tubs("LaBr3:Ce Detector 27", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector28 = new G4Tubs("LaBr3:Ce Detector 28", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector29 = new G4Tubs("LaBr3:Ce Detector 29", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector30 = new G4Tubs("LaBr3:Ce Detector 30", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  
  fLogicDetector = new G4LogicalVolume(fSolidDetector, fDetectorMaterial, "LaBr3:Ce Detector 1");
  fLogicDetector1 = new G4LogicalVolume(fSolidDetector1, fDetectorMaterial, "GRIFFIN LaBr3:Ce Detector");
  fLogicDetector2  = new G4LogicalVolume(fSolidDetector2,fDetectorMaterial, "LaBr3:Ce Detector 2");
  fLogicDetector3 = new G4LogicalVolume(fSolidDetector3, fDetectorMaterial, "LaBr3:Ce Detector 3");
  fLogicDetector4  = new G4LogicalVolume(fSolidDetector4,fDetectorMaterial, "LaBr3:Ce Detector 4");
  fLogicDetector5 = new G4LogicalVolume(fSolidDetector5, fDetectorMaterial, "LaBr3:Ce Detector 5");
  fLogicDetector6  = new G4LogicalVolume(fSolidDetector6,fDetectorMaterial, "LaBr3:Ce Detector 6");
  fLogicDetector7 = new G4LogicalVolume(fSolidDetector7, fDetectorMaterial, "LaBr3:Ce Detector 7");
  fLogicDetector8  = new G4LogicalVolume(fSolidDetector8,fDetectorMaterial, "LaBr3:Ce Detector 8");
  fLogicDetector9 = new G4LogicalVolume(fSolidDetector9, fDetectorMaterial, "LaBr3:Ce Detector 9");
  fLogicDetector10  = new G4LogicalVolume(fSolidDetector10,fDetectorMaterial, "LaBr3:Ce Detector 10");
  fLogicDetector11 = new G4LogicalVolume(fSolidDetector11, fDetectorMaterial, "LaBr3:Ce Detector 11");
  fLogicDetector12  = new G4LogicalVolume(fSolidDetector12,fDetectorMaterial, "LaBr3:Ce Detector 12");
  fLogicDetector13 = new G4LogicalVolume(fSolidDetector13, fDetectorMaterial, "LaBr3:Ce Detector 13");
  fLogicDetector14  = new G4LogicalVolume(fSolidDetector14,fDetectorMaterial, "LaBr3:Ce Detector 14");
  fLogicDetector15 = new G4LogicalVolume(fSolidDetector15, fDetectorMaterial, "LaBr3:Ce Detector 15");
  fLogicDetector16  = new G4LogicalVolume(fSolidDetector16,fDetectorMaterial, "LaBr3:Ce Detector 16");
  fLogicDetector17 = new G4LogicalVolume(fSolidDetector17, fDetectorMaterial, "LaBr3:Ce Detector 17");
  fLogicDetector18  = new G4LogicalVolume(fSolidDetector18,fDetectorMaterial, "LaBr3:Ce Detector 18");
  fLogicDetector19 = new G4LogicalVolume(fSolidDetector19, fDetectorMaterial, "LaBr3:Ce Detector 19");
  fLogicDetector20  = new G4LogicalVolume(fSolidDetector20,fDetectorMaterial, "LaBr3:Ce Detector 20");
  fLogicDetector21 = new G4LogicalVolume(fSolidDetector21, fDetectorMaterial, "LaBr3:Ce Detector 21");
  fLogicDetector22  = new G4LogicalVolume(fSolidDetector22,fDetectorMaterial, "LaBr3:Ce Detector 22");
  fLogicDetector23 = new G4LogicalVolume(fSolidDetector23, fDetectorMaterial, "LaBr3:Ce Detector 23");
  fLogicDetector24  = new G4LogicalVolume(fSolidDetector24,fDetectorMaterial, "LaBr3:Ce Detector 24");
  fLogicDetector25 = new G4LogicalVolume(fSolidDetector25, fDetectorMaterial, "LaBr3:Ce Detector 25");
  fLogicDetector26  = new G4LogicalVolume(fSolidDetector26,fDetectorMaterial, "LaBr3:Ce Detector 26");
  fLogicDetector27 = new G4LogicalVolume(fSolidDetector27, fDetectorMaterial, "LaBr3:Ce Detector 27");
  fLogicDetector28  = new G4LogicalVolume(fSolidDetector28,fDetectorMaterial, "LaBr3:Ce Detector 28");
  fLogicDetector29 = new G4LogicalVolume(fSolidDetector29, fDetectorMaterial, "LaBr3:Ce Detector 29");
  fLogicDetector30  = new G4LogicalVolume(fSolidDetector30,fDetectorMaterial, "LaBr3:Ce Detector 30");
 
  //fPhysiDetector = new G4PVPlacement(transform,
        			//fLogicDetector, 
       				//"BGO Detector", 
       				//fLogicWorld, 
       				//false, 
    			    	//0,
				    //false);
  fPhysiDetector1 = new G4PVPlacement(transform,
        			fLogicDetector1, 
       				"LaBr3:Ce Detector", 
       				fLogicWorld, 
       				false, 
    			    	0,
				    false);

  fLogicDetector->SetVisAttributes(fWireFrameVisAtt);
  fLogicDetector1->SetVisAttributes(fWireFrameVisAtt);
  }
  if (fDetectorGeometry == 3){  
  //  Cylindrical Detector Array
  //  This array is treated as a single volume for gamma ray detection.	
 
  fSolidDetector = new G4Tubs("LaBr3:Ce Detector 1", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector1 = new G4Tubs("GRIFFIN LaBr3:Ce Detector", 0., 0.5*fTotalDetectorDiameter1, 0.5*fTotalDetectorLength1, 0.*deg, 360.*deg);
  fSolidDetector2 = new G4Tubs("LaBr3:Ce Detector 2", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector3 = new G4Tubs("LaBr3:Ce Detector 3", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector4 = new G4Tubs("LaBr3:Ce Detector 4", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector5 = new G4Tubs("LaBr3:Ce Detector 5", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector6 = new G4Tubs("LaBr3:Ce Detector 6", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector7 = new G4Tubs("LaBr3:Ce Detector 7", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector8 = new G4Tubs("LaBr3:Ce Detector 8", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector9 = new G4Tubs("LaBr3:Ce Detector 9", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector10 = new G4Tubs("LaBr3:Ce Detector 10", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector11 = new G4Tubs("LaBr3:Ce Detector 11", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector12 = new G4Tubs("LaBr3:Ce Detector 12", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector13 = new G4Tubs("LaBr3:Ce Detector 13", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector14 = new G4Tubs("LaBr3:Ce Detector 14", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector15 = new G4Tubs("LaBr3:Ce Detector 15", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector16 = new G4Tubs("LaBr3:Ce Detector 16", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector17 = new G4Tubs("LaBr3:Ce Detector 17", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector18 = new G4Tubs("LaBr3:Ce Detector 18", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector19 = new G4Tubs("LaBr3:Ce Detector 19", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector20 = new G4Tubs("LaBr3:Ce Detector 20", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector21 = new G4Tubs("LaBr3:Ce Detector 21", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector22 = new G4Tubs("LaBr3:Ce Detector 22", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector23 = new G4Tubs("LaBr3:Ce Detector 23", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector24 = new G4Tubs("LaBr3:Ce Detector 24", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector25 = new G4Tubs("LaBr3:Ce Detector 25", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector26 = new G4Tubs("LaBr3:Ce Detector 26", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector27 = new G4Tubs("LaBr3:Ce Detector 27", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector28 = new G4Tubs("LaBr3:Ce Detector 28", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector29 = new G4Tubs("LaBr3:Ce Detector 29", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  fSolidDetector30 = new G4Tubs("LaBr3:Ce Detector 30", 0., 0.5*fTotalDetectorDiameter, 0.5*fTotalDetectorLength, 0.*deg, 360.*deg);
  
  fLogicDetector = new G4LogicalVolume(fSolidDetector,fDetectorMaterial, "LaBr3:Ce Detector 1");
  fLogicDetector1 = new G4LogicalVolume(fSolidDetector1,fDetectorMaterial, "GRIFFIN LaBr3:Ce Detector");
  fLogicDetector2  = new G4LogicalVolume(fSolidDetector2,fDetectorMaterial, "LaBr3:Ce Detector 2");
  fLogicDetector3 = new G4LogicalVolume(fSolidDetector3, fDetectorMaterial, "LaBr3:Ce Detector 3");
  fLogicDetector4  = new G4LogicalVolume(fSolidDetector4,fDetectorMaterial, "LaBr3:Ce Detector 4");
  fLogicDetector5 = new G4LogicalVolume(fSolidDetector5, fDetectorMaterial, "LaBr3:Ce Detector 5");
  fLogicDetector6  = new G4LogicalVolume(fSolidDetector6,fDetectorMaterial, "LaBr3:Ce Detector 6");
  fLogicDetector7 = new G4LogicalVolume(fSolidDetector7, fDetectorMaterial, "LaBr3:Ce Detector 7");
  fLogicDetector8  = new G4LogicalVolume(fSolidDetector8,fDetectorMaterial, "LaBr3:Ce Detector 8");
  fLogicDetector9 = new G4LogicalVolume(fSolidDetector9, fDetectorMaterial, "LaBr3:Ce Detector 9");
  fLogicDetector10  = new G4LogicalVolume(fSolidDetector10,fDetectorMaterial, "LaBr3:Ce Detector 10");
  fLogicDetector11 = new G4LogicalVolume(fSolidDetector11, fDetectorMaterial, "LaBr3:Ce Detector 11");
  fLogicDetector12  = new G4LogicalVolume(fSolidDetector12,fDetectorMaterial, "LaBr3:Ce Detector 12");
  fLogicDetector13 = new G4LogicalVolume(fSolidDetector13, fDetectorMaterial, "LaBr3:Ce Detector 13");
  fLogicDetector14  = new G4LogicalVolume(fSolidDetector14,fDetectorMaterial, "LaBr3:Ce Detector 14");
  fLogicDetector15 = new G4LogicalVolume(fSolidDetector15, fDetectorMaterial, "LaBr3:Ce Detector 15");
  fLogicDetector16  = new G4LogicalVolume(fSolidDetector16,fDetectorMaterial, "LaBr3:Ce Detector 16");
  fLogicDetector17 = new G4LogicalVolume(fSolidDetector17, fDetectorMaterial, "LaBr3:Ce Detector 17");
  fLogicDetector18  = new G4LogicalVolume(fSolidDetector18,fDetectorMaterial, "LaBr3:Ce Detector 18");
  fLogicDetector19 = new G4LogicalVolume(fSolidDetector19, fDetectorMaterial, "LaBr3:Ce Detector 19");
  fLogicDetector20  = new G4LogicalVolume(fSolidDetector20,fDetectorMaterial, "LaBr3:Ce Detector 20");
  fLogicDetector21 = new G4LogicalVolume(fSolidDetector21, fDetectorMaterial, "LaBr3:Ce Detector 21");
  fLogicDetector22  = new G4LogicalVolume(fSolidDetector22,fDetectorMaterial, "LaBr3:Ce Detector 22");
  fLogicDetector23 = new G4LogicalVolume(fSolidDetector23, fDetectorMaterial, "LaBr3:Ce Detector 23");
  fLogicDetector24  = new G4LogicalVolume(fSolidDetector24,fDetectorMaterial, "LaBr3:Ce Detector 24");
  fLogicDetector25 = new G4LogicalVolume(fSolidDetector25, fDetectorMaterial, "LaBr3:Ce Detector 25");
  fLogicDetector26  = new G4LogicalVolume(fSolidDetector26,fDetectorMaterial, "LaBr3:Ce Detector 26");
  fLogicDetector27 = new G4LogicalVolume(fSolidDetector27, fDetectorMaterial, "LaBr3:Ce Detector 27");
  fLogicDetector28  = new G4LogicalVolume(fSolidDetector28,fDetectorMaterial, "LaBr3:Ce Detector 28");
  fLogicDetector29 = new G4LogicalVolume(fSolidDetector29, fDetectorMaterial, "LaBr3:Ce Detector 29");
  fLogicDetector30  = new G4LogicalVolume(fSolidDetector30,fDetectorMaterial, "LaBr3:Ce Detector 30");
  
  G4AssemblyVolume* assemblyDetector = new G4AssemblyVolume();
  // 30 detectors are placed relative to the array volume.
  // Left-Hand Array, rotm180
  // 13
  //P.setX(-5.91*cm); P.setY(-7.68*cm); P.setZ(17.849*cm);
  P.setX(-fTotalDetectorDiameter - fAirGap); P.setY((-5/(2*sqrt(3)))*fTotalDetectorDiameter - fAirGap); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector13,Tr);
  // 19
  //P.setX(0.00*cm); P.setY(-7.68*cm); P.setZ(17.849*cm);
  P.setX(0*cm); P.setY((-5/(2*sqrt(3)))*fTotalDetectorDiameter - fAirGap); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector19,Tr);
  // 25
  //P.setX(5.91*cm); P.setY(-7.68*cm); P.setZ(17.849*cm);
  P.setX(fTotalDetectorDiameter + fAirGap); P.setY((-5/(2*sqrt(3)))*fTotalDetectorDiameter - fAirGap); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector25,Tr);
  // 11
  //P.setX(-8.87*cm); P.setY(-2.56*cm); P.setZ(17.849*cm);
   P.setX(-1.5*fTotalDetectorDiameter - 1.5*fAirGap); P.setY((-1/sqrt(3))*fTotalDetectorDiameter); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector11,Tr);
  // 17
  //P.setX(-2.96*cm); P.setY(-2.56*cm); P.setZ(17.849*cm);
  P.setX(-0.5*fTotalDetectorDiameter - 0.5*fAirGap); P.setY((-1/sqrt(3))*fTotalDetectorDiameter); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector17,Tr);
  // 23
  //P.setX(2.96*cm); P.setY(-2.56*cm); P.setZ(17.849*cm);
  P.setX(0.5*fTotalDetectorDiameter + 0.5*fAirGap); P.setY((-1/sqrt(3))*fTotalDetectorDiameter); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector23,Tr);
  // 29
  //P.setX(8.87*cm); P.setY(-2.56*cm); P.setZ(17.849*cm);
   P.setX(1.5*fTotalDetectorDiameter + 1.5*fAirGap); P.setY((-1/sqrt(3))*fTotalDetectorDiameter); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector29,Tr);
  // 15
  //P.setX(-5.91*cm); P.setY(2.56*cm); P.setZ(17.849*cm);
  P.setX(-fTotalDetectorDiameter - fAirGap); P.setY((1/(2*sqrt(3)))*fTotalDetectorDiameter + fAirGap); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector15,Tr);
  // 21
  //P.setX(0.00*cm); P.setY(2.56*cm); P.setZ(17.849*cm);
  P.setX(0.); P.setY((1/(2*sqrt(3)))*fTotalDetectorDiameter + fAirGap); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector21,Tr);
  // 27
  //P.setX(5.91*cm); P.setY(2.56*cm); P.setZ(17.849*cm);
  P.setX(fTotalDetectorDiameter + fAirGap); P.setY((1/(2*sqrt(3)))*fTotalDetectorDiameter + fAirGap); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector27,Tr);
  // Straddling Detectors
  // 4
  //P.setX(-8.87*cm); P.setY(10.08*cm); P.setZ(9.*cm);
  if (fPMTDiameter < (6.35*sqrt(3)/(4 - sqrt(3)))*cm){
  P.setX((-0.5 - 0.5*sqrt(3))*fTotalDetectorDiameter - (0.5 + 0.5*sqrt(3))*fAirGap); P.setY(3.175*cm + fTotalDetectorDiameter + 2.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  else if (fPMTDiameter >= (6.35*sqrt(3)/(4 - sqrt(3)))*cm){
  P.setX((-0.5 - 0.5*sqrt(3))*fTotalDetectorDiameter - (0.5 + 0.5*sqrt(3))*fAirGap); P.setY(((4 + sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter + 2.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector4,Tr);
  // 5 
  //P.setX(-2.96*cm); P.setY(7.68*cm); P.setZ(9.*cm);
  if (fPMTDiameter < (6.35*sqrt(3)/(4 - sqrt(3)))*cm){
  P.setX(-0.5*fTotalDetectorDiameter - 0.5*fAirGap); P.setY(3.175*cm + 0.5*fTotalDetectorDiameter + 2*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  else if (fPMTDiameter >= (6.35*sqrt(3)/(4 - sqrt(3)))*cm){
  P.setX(-0.5*fTotalDetectorDiameter - 0.5*fAirGap); P.setY((2/(sqrt(3)))*fTotalDetectorDiameter + 2*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector5,Tr);
  // 6
  //P.setX(2.96*cm); P.setY(7.68*cm); P.setZ(9.*cm);
  if (fPMTDiameter < (6.35*sqrt(3)/(4 - sqrt(3)))*cm){ 
  P.setX(0.5*fTotalDetectorDiameter + 0.5*fAirGap); P.setY(3.175*cm + 0.5*fTotalDetectorDiameter + 2*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);  
  }
  else if (fPMTDiameter >= (6.35*sqrt(3)/(4 - sqrt(3)))*cm){
  P.setX(0.5*fTotalDetectorDiameter + 0.5*fAirGap); P.setY((2/(sqrt(3)))*fTotalDetectorDiameter + 2*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector6,Tr); 
  // 7
  //P.setX(8.87*cm); P.setY(10.08*cm); P.setZ(9.*cm);
  if (fPMTDiameter < (6.35*sqrt(3)/(4 - sqrt(3)))*cm){
  P.setX((0.5 + 0.5*sqrt(3))*fTotalDetectorDiameter + (0.5 + 0.5*sqrt(3))*fAirGap); P.setY(3.175*cm + fTotalDetectorDiameter + 2.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  else if (fPMTDiameter >= (6.35*sqrt(3)/(4 - sqrt(3)))*cm){
  P.setX((0.5 + 0.5*sqrt(3))*fTotalDetectorDiameter + (0.5 + 0.5*sqrt(3))*fAirGap); P.setY(((4 + sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter + 2.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector7,Tr);
  // 3
  //P.setX(-11.83*cm); P.setY(4.96*cm); P.setZ(2.299*cm); 
  if (fPMTDiameter < 2*fOutTBoxX/(1 + sqrt(3))){
  P.setX(-fOutTBoxX - 0.5*fTotalDetectorDiameter - (1 + 0.5*sqrt(3))*fAirGap); P.setY(((1 + sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter + 1.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY - 6.701*cm);
  }
  else if (fPMTDiameter >= 2*fOutTBoxX/(1 + sqrt(3))){
  P.setX((-1 - 0.5*sqrt(3))*fTotalDetectorDiameter - (1 + 0.5*sqrt(3))*fAirGap); P.setY(((1 + sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter + 1.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY - 6.701*cm);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector3,Tr);
  // 9
  //P.setX(11.83*cm); P.setY(4.96*cm); P.setZ(9.*cm);
  if (fPMTDiameter < 2*fOutTBoxX/(1 + sqrt(3))){
  P.setX(fOutTBoxX + 0.5*fTotalDetectorDiameter + (1 + 0.5*sqrt(3))*fAirGap); P.setY(((1 + sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter + 1.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  else if (fPMTDiameter >= 2*fOutTBoxX/(1 + sqrt(3))){
  P.setX((1 + 0.5*sqrt(3))*fTotalDetectorDiameter + (1 + 0.5*sqrt(3))*fAirGap); P.setY(((1 + sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter + 1.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector9,Tr);
  // 1
  //P.setX(-14.78*cm); P.setY(-4.96*cm); P.setZ(1.8*cm);
  if (fPMTDiameter < 2*fOutTBoxX/(1 + sqrt(3))){
  P.setX(-fOutTBoxX - fTotalDetectorDiameter - (1.5 + 0.5*sqrt(3))*fAirGap); P.setY(((-2 - sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter - 0.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY - 7.2*cm);  
  }
  else if (fPMTDiameter >= 2*fOutTBoxX/(1 + sqrt(3))){
  P.setX((-1.5 - 0.5*sqrt(3))*fTotalDetectorDiameter - (1.5 + 0.5*sqrt(3))*fAirGap); P.setY(((-2 - sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter - 0.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY - 7.2*cm);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  // 2
  //P.setX(-11.83*cm); P.setY(-10.08*cm); P.setZ(9.*cm);
  if (fPMTDiameter < 2*fOutTBoxX/(1 + sqrt(3))){
  P.setX(-fOutTBoxX - 0.5*fTotalDetectorDiameter - (1 + 0.5*sqrt(3))*fAirGap); P.setY(((-5 - sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter - 1.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  else if (fPMTDiameter >= 2*fOutTBoxX/(1 + sqrt(3))){
   P.setX((-1 - 0.5*sqrt(3))*fTotalDetectorDiameter - (1 + 0.5*sqrt(3))*fAirGap); P.setY(((-5 - sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter - 1.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector2,Tr);
  // 8
  //P.setX(11.83*cm); P.setY(-10.08*cm); P.setZ(9.*cm);
  if (fPMTDiameter < 2*fOutTBoxX/(1 + sqrt(3))){
  P.setX(fOutTBoxX + 0.5*fTotalDetectorDiameter + (1 + 0.5*sqrt(3))*fAirGap); P.setY(((-5 - sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter - 1.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  else if (fPMTDiameter >= 2*fOutTBoxX/(1 + sqrt(3))){
  P.setX((1 + 0.5*sqrt(3))*fTotalDetectorDiameter + (1 + 0.5*sqrt(3))*fAirGap); P.setY(((-5 - sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter - 1.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector8,Tr);
  // 10
  //P.setX(14.78*cm); P.setY(-4.96*cm); P.setZ(9.*cm);
  if (fPMTDiameter < 2*fOutTBoxX/(1 + sqrt(3))){
  P.setX(fOutTBoxX + fTotalDetectorDiameter + (1.5 + 0.5*sqrt(3))*fAirGap); P.setY(((-2 - sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter - 0.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  else if (fPMTDiameter >= 2*fOutTBoxX/(1 + sqrt(3))){
  P.setX((1.5 + 0.5*sqrt(3))*fTotalDetectorDiameter + (1.5 + 0.5*sqrt(3))*fAirGap); P.setY(((-2 - sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter - 0.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector10,Tr);
  // Right-Hand Array, rotm
  // 14
  //P.setX(-5.91*cm); P.setY(-7.68*cm); P.setZ(-17.849*cm);
  P.setX(-fTotalDetectorDiameter - fAirGap); P.setY((-5/(2*sqrt(3)))*fTotalDetectorDiameter - fAirGap); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector14,Tr);
  // 20
  //P.setX(0.00*cm); P.setY(-7.68*cm); P.setZ(-17.849*cm);
  P.setX(0.); P.setY((-5/(2*sqrt(3)))*fTotalDetectorDiameter - fAirGap); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector20,Tr);
  // 26
  //P.setX(5.91*cm); P.setY(-7.68*cm); P.setZ(-17.849*cm);
  P.setX(fTotalDetectorDiameter + fAirGap); P.setY((-5/(2*sqrt(3)))*fTotalDetectorDiameter - fAirGap); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector26,Tr);
  // 12
  //P.setX(-8.87*cm); P.setY(-2.56*cm); P.setZ(-17.849*cm);
  P.setX(-1.5*fTotalDetectorDiameter - 1.5*fAirGap); P.setY((-1/sqrt(3))*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector12,Tr);
  // 18
  //P.setX(-2.96*cm); P.setY(-2.56*cm); P.setZ(-17.849*cm);
  P.setX(-0.5*fTotalDetectorDiameter - 0.5*fAirGap); P.setY((-1/sqrt(3))*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector18,Tr);
  // 24
  //P.setX(2.96*cm); P.setY(-2.56*cm); P.setZ(-17.849*cm);
  P.setX(0.5*fTotalDetectorDiameter + 0.5*fAirGap); P.setY((-1/sqrt(3))*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector24,Tr);
  // 30
  //P.setX(8.87*cm); P.setY(-2.56*cm); P.setZ(-17.849*cm);
  P.setX(1.5*fTotalDetectorDiameter + 1.5*fAirGap); P.setY((-1/sqrt(3))*fTotalDetectorDiameter); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector30,Tr);
  // 16
  //P.setX(-5.91*cm); P.setY(2.56*cm); P.setZ(-17.849*cm);
  P.setX(-fTotalDetectorDiameter - fAirGap); P.setY((1/(2*sqrt(3)))*fTotalDetectorDiameter + fAirGap); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector16,Tr);
  // 22
  //P.setX(0.00*cm); P.setY(2.56*cm); P.setZ(-17.849*cm);
  P.setX(0.); P.setY((1/(2*sqrt(3)))*fTotalDetectorDiameter + fAirGap); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector22,Tr);
  // 28
  //P.setX(5.91*cm); P.setY(2.56*cm); P.setZ(-17.849*cm);
  P.setX(fTotalDetectorDiameter + fAirGap); P.setY((1/(2*sqrt(3)))*fTotalDetectorDiameter + fAirGap); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector28,Tr);
  // Assembly Placement in the World
  // Rotated 90° along the axes to be perpendicular to the gas target.
  // Rotated 180° along the z-axis the place the array downside up instead of upside down.
  // Coordinate transformation: (X, Y, Z) -> (Z, Y, -X)
  G4ThreeVector WorldP(0,0,0);
  G4RotationMatrix worldrotm = G4RotationMatrix(90.*deg,90.*deg,270.*deg); 
  G4Transform3D WorldTr = G4Transform3D(worldrotm, WorldP);
  assemblyDetector->MakeImprint(fLogicWorld, WorldTr);
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
  
  // Extra fPMTWin here which causes overlaps.
  //fSolidPMTWin = new G4Tubs("PMTWin", 0., 0.5*fDetectorDiameter, 0.5*fPMTWinThickness, 0.*deg, 360.*deg);
  //fLogicPMTWin = new G4LogicalVolume(fSolidPMTWin, fPMTWinMaterial, "PMTWin");
  //fPhysiPMTWin = new G4PVPlacement(0, 
        			//G4ThreeVector(0.,0.,fZposPMTWin), 
        			//fLogicPMTWin, 
       				//"PMTWin", 
       				//fLogicPMT, 
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
  fSolidCrystal = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal = new G4LogicalVolume(fSolidCrystal, fCrystalMaterial1, "Crystal");
  fPhysiCrystal = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal, 
       				"Crystal", 
       				fLogicDetector, 
       				false, 
    			    	1);
 			    	
  fSolidCrystal1 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter1, 0.5*fDetectorLength1, 0.*deg, 360.*deg);
  fLogicCrystal1 = new G4LogicalVolume(fSolidCrystal1, fCrystalMaterial1, "Crystal");
  fPhysiCrystal1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos1), 
        			fLogicCrystal1, 
       				"Crystal", 
       				fLogicDetector1, 
       				false, 
    			    	0);
    			    	
  fSolidCrystal2 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal2 = new G4LogicalVolume(fSolidCrystal2, fCrystalMaterial1, "Crystal");
  fPhysiCrystal2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal2, 
       				"Crystal", 
       				fLogicDetector2, 
       				false, 
    			    	2);
    			    	
  fSolidCrystal3 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal3 = new G4LogicalVolume(fSolidCrystal3, fCrystalMaterial1, "Crystal");
  fPhysiCrystal3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal3, 
       				"Crystal", 
       				fLogicDetector3, 
       				false, 
    			    	3);  
    			    	
  fSolidCrystal4 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal4 = new G4LogicalVolume(fSolidCrystal4, fCrystalMaterial1, "Crystal");
  fPhysiCrystal4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal4, 
       				"Crystal", 
       				fLogicDetector4, 
       				false, 
    			    	4); 
    			    	
  fSolidCrystal5 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal5 = new G4LogicalVolume(fSolidCrystal5, fCrystalMaterial1, "Crystal");
  fPhysiCrystal5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal5, 
       				"Crystal", 
       				fLogicDetector5, 
       				false, 
    			    	5); 
    			    	
  fSolidCrystal6 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal6 = new G4LogicalVolume(fSolidCrystal6, fCrystalMaterial1, "Crystal");
  fPhysiCrystal6 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal6, 
       				"Crystal", 
       				fLogicDetector6, 
       				false, 
    			    	6);			    	
    			    	
  fSolidCrystal7 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal7 = new G4LogicalVolume(fSolidCrystal7, fCrystalMaterial1, "Crystal");
  fPhysiCrystal7 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal7, 
       				"Crystal", 
       				fLogicDetector7, 
       				false, 
    			    	7);	
    			    	
  fSolidCrystal8 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal8 = new G4LogicalVolume(fSolidCrystal8, fCrystalMaterial1, "Crystal");
  fPhysiCrystal8 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal8, 
       				"Crystal", 
       				fLogicDetector8, 
       				false, 
    			    	8);	
    			    	
  fSolidCrystal9 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal9 = new G4LogicalVolume(fSolidCrystal9, fCrystalMaterial1, "Crystal");
  fPhysiCrystal9 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal9, 
       				"Crystal", 
       				fLogicDetector9, 
       				false, 
    			    	9);	
    			    	
  fSolidCrystal10 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal10 = new G4LogicalVolume(fSolidCrystal10, fCrystalMaterial1, "Crystal");
  fPhysiCrystal10 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal10, 
       				"Crystal", 
       				fLogicDetector10, 
       				false, 
    			    	10);
    			    	
  fSolidCrystal11 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal11 = new G4LogicalVolume(fSolidCrystal11, fCrystalMaterial1, "Crystal");
  fPhysiCrystal11 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal11, 
       				"Crystal", 
       				fLogicDetector11, 
       				false, 
    			    	11);
    			    	
  fSolidCrystal12 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal12 = new G4LogicalVolume(fSolidCrystal12, fCrystalMaterial1, "Crystal");
  fPhysiCrystal12 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal12, 
       				"Crystal", 
       				fLogicDetector12, 
       				false, 
    			    	12);	
    			    	
  fSolidCrystal13 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal13 = new G4LogicalVolume(fSolidCrystal13, fCrystalMaterial1, "Crystal");
  fPhysiCrystal13 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal13, 
       				"Crystal", 
       				fLogicDetector13, 
       				false, 
    			    	13);	
    			    		
  fSolidCrystal14 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal14 = new G4LogicalVolume(fSolidCrystal14, fCrystalMaterial1, "Crystal");
  fPhysiCrystal14 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal14, 
       				"Crystal", 
       				fLogicDetector14, 
       				false, 
    			    	14);
    			    	
  fSolidCrystal15 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal15 = new G4LogicalVolume(fSolidCrystal15, fCrystalMaterial1, "Crystal");
  fPhysiCrystal15 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal15, 
       				"Crystal", 
       				fLogicDetector15, 
       				false, 
    			    	15);
    			    	
  fSolidCrystal16 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal16 = new G4LogicalVolume(fSolidCrystal16, fCrystalMaterial1, "Crystal");
  fPhysiCrystal16 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal16, 
       				"Crystal", 
       				fLogicDetector16, 
       				false, 
    			    	16);
    			    	
  fSolidCrystal17 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal17 = new G4LogicalVolume(fSolidCrystal17, fCrystalMaterial1, "Crystal");
  fPhysiCrystal17 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal17, 
       				"Crystal", 
       				fLogicDetector17, 
       				false, 
    			    	17);
    			    	
    			    	  			    	
  fSolidCrystal18 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal18 = new G4LogicalVolume(fSolidCrystal18, fCrystalMaterial1, "Crystal");
  fPhysiCrystal18 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal18, 
       				"Crystal", 
       				fLogicDetector18, 
       				false, 
    			    	18);
    			    	
  fSolidCrystal19 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal19 = new G4LogicalVolume(fSolidCrystal19, fCrystalMaterial1, "Crystal");
  fPhysiCrystal19 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal19, 
       				"Crystal", 
       				fLogicDetector19, 
       				false, 
    			    	19);
    			    	
  fSolidCrystal20 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal20 = new G4LogicalVolume(fSolidCrystal20, fCrystalMaterial1, "Crystal");
  fPhysiCrystal20 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal20, 
       				"Crystal", 
       				fLogicDetector20, 
       				false, 
    			    	20);
    			    	
  fSolidCrystal21 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal21 = new G4LogicalVolume(fSolidCrystal21, fCrystalMaterial1, "Crystal");
  fPhysiCrystal21 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal21, 
       				"Crystal", 
       				fLogicDetector21, 
       				false, 
    			    	21);
    			    	
  fSolidCrystal22 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal22 = new G4LogicalVolume(fSolidCrystal22, fCrystalMaterial1, "Crystal");
  fPhysiCrystal22 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal22, 
       				"Crystal", 
       				fLogicDetector22, 
       				false, 
    			    	22);
    			    	
  fSolidCrystal23 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal23 = new G4LogicalVolume(fSolidCrystal23, fCrystalMaterial1, "Crystal");
  fPhysiCrystal23 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal23, 
       				"Crystal", 
       				fLogicDetector23, 
       				false, 
    			    	23);
    			    	
  fSolidCrystal24 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal24 = new G4LogicalVolume(fSolidCrystal24, fCrystalMaterial1, "Crystal");
  fPhysiCrystal24 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal24, 
       				"Crystal", 
       				fLogicDetector24, 
       				false, 
    			    	24);
    			    	
  fSolidCrystal25 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal25 = new G4LogicalVolume(fSolidCrystal25, fCrystalMaterial1, "Crystal");
  fPhysiCrystal25 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal25, 
       				"Crystal", 
       				fLogicDetector25, 
       				false, 
    			    	25);
    			    	
  fSolidCrystal26 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal26 = new G4LogicalVolume(fSolidCrystal26, fCrystalMaterial1, "Crystal");
  fPhysiCrystal26 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal26, 
       				"Crystal", 
       				fLogicDetector26, 
       				false, 
    			    	26);
    			    	
  fSolidCrystal27 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal27 = new G4LogicalVolume(fSolidCrystal27, fCrystalMaterial1, "Crystal");
  fPhysiCrystal27 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal27, 
       				"Crystal", 
       				fLogicDetector27, 
       				false, 
    			    	27);
    			    	
  fSolidCrystal28 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal28 = new G4LogicalVolume(fSolidCrystal28, fCrystalMaterial1, "Crystal");
  fPhysiCrystal28 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal28, 
       				"Crystal", 
       				fLogicDetector28, 
       				false, 
    			    	28);
    			    	
  fSolidCrystal29 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal29 = new G4LogicalVolume(fSolidCrystal29, fCrystalMaterial1, "Crystal");
  fPhysiCrystal29 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal29, 
       				"Crystal", 
       				fLogicDetector29, 
       				false, 
    			    	29);
    			    	
  fSolidCrystal30 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter, 0.5*fDetectorLength, 0.*deg, 360.*deg);
  fLogicCrystal30 = new G4LogicalVolume(fSolidCrystal30, fCrystalMaterial1, "Crystal");
  fPhysiCrystal30 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal30, 
       				"Crystal", 
       				fLogicDetector30, 
       				false, 
    			    	30);
    			    	
  fLogicCrystal->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal1->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal2->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal3->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal4->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal5->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal6->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal7->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal8->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal9->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal10->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal11->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal12->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal13->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal14->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal15->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal16->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal17->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal18->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal19->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal20->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal21->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal22->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal23->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal24->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal25->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal26->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal27->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal28->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal29->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal30->SetVisAttributes(fMagnetaVisAtt);

  // Gap (gap surrounding crystal and casing or a reflector as required)
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
    			    	
  fSolidGap2 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap2 = new G4LogicalVolume(fSolidGap2, fGapMaterial, "Gap");
  fPhysiGap2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap2, 
       				"Gap", 
       				fLogicDetector2, 
       				false, 
    			    	0);

  fSolidGap3 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap3 = new G4LogicalVolume(fSolidGap3, fGapMaterial, "Gap");
  fPhysiGap3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap3, 
       				"Gap", 
       				fLogicDetector3, 
       				false, 
    			    	0);
    			    	
  fSolidGap4 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap4 = new G4LogicalVolume(fSolidGap4, fGapMaterial, "Gap");
  fPhysiGap4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap4, 
       				"Gap", 
       				fLogicDetector4, 
       				false, 
    			    	0);
    			    	
  fSolidGap5 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap5 = new G4LogicalVolume(fSolidGap5, fGapMaterial, "Gap");
  fPhysiGap5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap5, 
       				"Gap", 
       				fLogicDetector5, 
       				false, 
    			    	0);
    			    	
  fSolidGap6 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap6 = new G4LogicalVolume(fSolidGap6, fGapMaterial, "Gap");
  fPhysiGap6 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap6, 
       				"Gap", 
       				fLogicDetector6, 
       				false, 
    			    	0);
    			    	
  fSolidGap7 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap7 = new G4LogicalVolume(fSolidGap7, fGapMaterial, "Gap");
  fPhysiGap7 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap7, 
       				"Gap", 
       				fLogicDetector7, 
       				false, 
    			    	0);
    			    	
  fSolidGap8 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap8 = new G4LogicalVolume(fSolidGap8, fGapMaterial, "Gap");
  fPhysiGap8 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap8, 
       				"Gap", 
       				fLogicDetector8, 
       				false, 
    			    	0);
    			    	
  fSolidGap9 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap9 = new G4LogicalVolume(fSolidGap9, fGapMaterial, "Gap");
  fPhysiGap9 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap9, 
       				"Gap", 
       				fLogicDetector9, 
       				false, 
    			    	0);
    			    	
  fSolidGap10 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap10 = new G4LogicalVolume(fSolidGap10, fGapMaterial, "Gap");
  fPhysiGap10 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap10, 
       				"Gap", 
       				fLogicDetector10, 
       				false, 
    			    	0);
    			    	
  fSolidGap11 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap11 = new G4LogicalVolume(fSolidGap11, fGapMaterial, "Gap");
  fPhysiGap11 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap11, 
       				"Gap", 
       				fLogicDetector11, 
       				false, 
    			    	0);
    			    	
  fSolidGap12 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap12 = new G4LogicalVolume(fSolidGap12, fGapMaterial, "Gap");
  fPhysiGap12 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap12, 
       				"Gap", 
       				fLogicDetector12, 
       				false, 
    			    	0);
    			    	
  fSolidGap13 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap13 = new G4LogicalVolume(fSolidGap13, fGapMaterial, "Gap");
  fPhysiGap13 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap13, 
       				"Gap", 
       				fLogicDetector13, 
       				false, 
    			    	0);
    			    	
  fSolidGap14 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap14 = new G4LogicalVolume(fSolidGap14, fGapMaterial, "Gap");
  fPhysiGap14 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap14, 
       				"Gap", 
       				fLogicDetector14, 
       				false, 
    			    	0);
    			    	
  fSolidGap15 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap15 = new G4LogicalVolume(fSolidGap15, fGapMaterial, "Gap");
  fPhysiGap15 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap15, 
       				"Gap", 
       				fLogicDetector15, 
       				false, 
    			    	0);
    			    	
  fSolidGap16 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap16 = new G4LogicalVolume(fSolidGap16, fGapMaterial, "Gap");
  fPhysiGap16 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap16, 
       				"Gap", 
       				fLogicDetector16, 
       				false, 
    			    	0);
    			    	
  fSolidGap17 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap17 = new G4LogicalVolume(fSolidGap17, fGapMaterial, "Gap");
  fPhysiGap17 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap17, 
       				"Gap", 
       				fLogicDetector17, 
       				false, 
    			    	0);
    			    	
  fSolidGap18 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap18 = new G4LogicalVolume(fSolidGap18, fGapMaterial, "Gap");
  fPhysiGap18 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap18, 
       				"Gap", 
       				fLogicDetector18, 
       				false, 
    			    	0);
    			    	
  fSolidGap19 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap19 = new G4LogicalVolume(fSolidGap19, fGapMaterial, "Gap");
  fPhysiGap19 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap19, 
       				"Gap", 
       				fLogicDetector19, 
       				false, 
    			    	0);
    			    	
  fSolidGap20 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap20 = new G4LogicalVolume(fSolidGap20, fGapMaterial, "Gap");
  fPhysiGap20 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap20, 
       				"Gap", 
       				fLogicDetector20, 
       				false, 
    			    	0);
    			    	
  fSolidGap21 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap21 = new G4LogicalVolume(fSolidGap21, fGapMaterial, "Gap");
  fPhysiGap21 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap21, 
       				"Gap", 
       				fLogicDetector21, 
       				false, 
    			    	0);
    			    	
  fSolidGap22 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap22 = new G4LogicalVolume(fSolidGap22, fGapMaterial, "Gap");
  fPhysiGap22 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap22, 
       				"Gap", 
       				fLogicDetector22, 
       				false, 
    			    	0);
    			    	
  fSolidGap23 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap23 = new G4LogicalVolume(fSolidGap23, fGapMaterial, "Gap");
  fPhysiGap23 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap23, 
       				"Gap", 
       				fLogicDetector23, 
       				false, 
    			    	0);
    			    	
  fSolidGap24 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap24 = new G4LogicalVolume(fSolidGap24, fGapMaterial, "Gap");
  fPhysiGap24 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap24, 
       				"Gap", 
       				fLogicDetector24, 
       				false, 
    			    	0);
    			    	
  fSolidGap25 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap25 = new G4LogicalVolume(fSolidGap25, fGapMaterial, "Gap");
  fPhysiGap25 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap25, 
       				"Gap", 
       				fLogicDetector25, 
       				false, 
    			    	0);
    			    	
  fSolidGap26 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap26 = new G4LogicalVolume(fSolidGap26, fGapMaterial, "Gap");
  fPhysiGap26 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap26, 
       				"Gap", 
       				fLogicDetector26, 
       				false, 
    			    	0);
    			    	
  fSolidGap27 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap27 = new G4LogicalVolume(fSolidGap27, fGapMaterial, "Gap");
  fPhysiGap27 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap27, 
       				"Gap", 
       				fLogicDetector27, 
       				false, 
    			    	0);
    			    	
  fSolidGap28 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap28 = new G4LogicalVolume(fSolidGap28, fGapMaterial, "Gap");
  fPhysiGap28 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap28, 
       				"Gap", 
       				fLogicDetector28, 
       				false, 
    			    	0);
    			    	
  fSolidGap29 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap29 = new G4LogicalVolume(fSolidGap29, fGapMaterial, "Gap");
  fPhysiGap29 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap29, 
       				"Gap", 
       				fLogicDetector29, 
       				false, 
    			    	0);
    			    	
  fSolidGap30 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg);
  fLogicGap30 = new G4LogicalVolume(fSolidGap30, fGapMaterial, "Gap");
  fPhysiGap30 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap30, 
       				"Gap", 
       				fLogicDetector30, 
       				false, 
    			    	0);

  // FaceGap (gap between face of crystal and casing or a reflector as required)
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
    			    	
  fSolidFaceGap2 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap2 = new G4LogicalVolume(fSolidFaceGap2, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap2, 
       				"FaceGap", 
       				fLogicDetector2, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap3 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap3 = new G4LogicalVolume(fSolidFaceGap3, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap3, 
       				"FaceGap", 
       				fLogicDetector3, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap4 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap4 = new G4LogicalVolume(fSolidFaceGap4, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap4, 
       				"FaceGap", 
       				fLogicDetector4, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap5 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap5 = new G4LogicalVolume(fSolidFaceGap5, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap5, 
       				"FaceGap", 
       				fLogicDetector5, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap6 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap6 = new G4LogicalVolume(fSolidFaceGap6, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap6 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap6, 
       				"FaceGap", 
       				fLogicDetector6, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap7 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap7 = new G4LogicalVolume(fSolidFaceGap7, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap7 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap7, 
       				"FaceGap", 
       				fLogicDetector7, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap8 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap8 = new G4LogicalVolume(fSolidFaceGap8, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap8 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap8, 
       				"FaceGap", 
       				fLogicDetector8, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap9 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap9 = new G4LogicalVolume(fSolidFaceGap9, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap9 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap9, 
       				"FaceGap", 
       				fLogicDetector9, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap10 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap10 = new G4LogicalVolume(fSolidFaceGap10, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap10 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap10, 
       				"FaceGap", 
       				fLogicDetector10, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap11 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap11 = new G4LogicalVolume(fSolidFaceGap11, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap11 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap11, 
       				"FaceGap", 
       				fLogicDetector11, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap12 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap12 = new G4LogicalVolume(fSolidFaceGap12, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap12 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap12, 
       				"FaceGap", 
       				fLogicDetector12, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap13 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap13 = new G4LogicalVolume(fSolidFaceGap13, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap13 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap13, 
       				"FaceGap", 
       				fLogicDetector13, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap14 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap14 = new G4LogicalVolume(fSolidFaceGap14, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap14 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap14, 
       				"FaceGap", 
       				fLogicDetector14, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap15 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap15 = new G4LogicalVolume(fSolidFaceGap15, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap15 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap15, 
       				"FaceGap", 
       				fLogicDetector15, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap16 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap16 = new G4LogicalVolume(fSolidFaceGap16, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap16 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap16, 
       				"FaceGap", 
       				fLogicDetector16, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap17 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap17 = new G4LogicalVolume(fSolidFaceGap17, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap17 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap17, 
       				"FaceGap", 
       				fLogicDetector17, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap18 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap18 = new G4LogicalVolume(fSolidFaceGap18, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap18 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap18, 
       				"FaceGap", 
       				fLogicDetector18, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap19 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap19 = new G4LogicalVolume(fSolidFaceGap19, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap19 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap19, 
       				"FaceGap", 
       				fLogicDetector19, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap20 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap20 = new G4LogicalVolume(fSolidFaceGap20, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap20 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap20, 
       				"FaceGap", 
       				fLogicDetector20, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap21 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap21 = new G4LogicalVolume(fSolidFaceGap21, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap21 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap21, 
       				"FaceGap", 
       				fLogicDetector21, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap22 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap22 = new G4LogicalVolume(fSolidFaceGap22, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap22 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap22, 
       				"FaceGap", 
       				fLogicDetector22, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap23 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap23 = new G4LogicalVolume(fSolidFaceGap23, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap23 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap23, 
       				"FaceGap", 
       				fLogicDetector23, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap24 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap24 = new G4LogicalVolume(fSolidFaceGap24, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap24 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap24, 
       				"FaceGap", 
       				fLogicDetector24, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap25 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap25 = new G4LogicalVolume(fSolidFaceGap25, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap25 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap25, 
       				"FaceGap", 
       				fLogicDetector25, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap26 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap26 = new G4LogicalVolume(fSolidFaceGap26, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap26 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap26, 
       				"FaceGap", 
       				fLogicDetector26, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap27 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap27 = new G4LogicalVolume(fSolidFaceGap27, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap27 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap27, 
       				"FaceGap", 
       				fLogicDetector27, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap28 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap28 = new G4LogicalVolume(fSolidFaceGap28, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap28 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap28, 
       				"FaceGap", 
       				fLogicDetector28, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap29 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap29 = new G4LogicalVolume(fSolidFaceGap29, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap29 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap29, 
       				"FaceGap", 
       				fLogicDetector29, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap30 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter, 0.5*fGapFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceGap30 = new G4LogicalVolume(fSolidFaceGap30, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap30 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap30, 
       				"FaceGap", 
       				fLogicDetector30, 
       				false, 
    			    	0);
    			    	
  fLogicGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap2->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap3->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap4->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap5->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap6->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap7->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap8->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap9->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap10->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap11->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap12->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap13->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap14->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap15->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap16->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap17->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap18->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap19->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap20->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap21->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap22->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap23->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap24->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap25->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap26->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap27->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap28->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap29->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap30->SetVisAttributes(fAuxEdgeVisAtt);
  
  fLogicGap->SetVisAttributes(fBlueVisAtt);
  fLogicGap2->SetVisAttributes(fBlueVisAtt);
  fLogicGap3->SetVisAttributes(fBlueVisAtt);
  fLogicGap4->SetVisAttributes(fBlueVisAtt);
  fLogicGap5->SetVisAttributes(fBlueVisAtt);
  fLogicGap6->SetVisAttributes(fBlueVisAtt);
  fLogicGap7->SetVisAttributes(fBlueVisAtt);
  fLogicGap8->SetVisAttributes(fBlueVisAtt);
  fLogicGap9->SetVisAttributes(fBlueVisAtt);
  fLogicGap10->SetVisAttributes(fBlueVisAtt);
  fLogicGap11->SetVisAttributes(fBlueVisAtt);
  fLogicGap12->SetVisAttributes(fBlueVisAtt);
  fLogicGap13->SetVisAttributes(fBlueVisAtt);
  fLogicGap14->SetVisAttributes(fBlueVisAtt);
  fLogicGap15->SetVisAttributes(fBlueVisAtt);
  fLogicGap16->SetVisAttributes(fBlueVisAtt);
  fLogicGap17->SetVisAttributes(fBlueVisAtt);
  fLogicGap18->SetVisAttributes(fBlueVisAtt);
  fLogicGap19->SetVisAttributes(fBlueVisAtt);
  fLogicGap20->SetVisAttributes(fBlueVisAtt);
  fLogicGap21->SetVisAttributes(fBlueVisAtt);
  fLogicGap22->SetVisAttributes(fBlueVisAtt);
  fLogicGap23->SetVisAttributes(fBlueVisAtt);
  fLogicGap24->SetVisAttributes(fBlueVisAtt);
  fLogicGap25->SetVisAttributes(fBlueVisAtt);
  fLogicGap26->SetVisAttributes(fBlueVisAtt);
  fLogicGap27->SetVisAttributes(fBlueVisAtt);
  fLogicGap28->SetVisAttributes(fBlueVisAtt);
  fLogicGap29->SetVisAttributes(fBlueVisAtt);
  fLogicGap30->SetVisAttributes(fBlueVisAtt);
  
  fLogicFaceGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap2->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap3->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap4->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap5->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap6->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap7->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap8->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap9->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap10->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap11->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap12->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap13->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap14->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap15->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap16->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap17->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap18->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap19->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap20->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap21->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap22->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap23->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap24->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap25->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap26->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap27->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap28->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap29->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap30->SetVisAttributes(fAuxEdgeVisAtt);
  
  fLogicFaceGap->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap2->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap3->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap4->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap5->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap6->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap7->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap8->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap9->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap10->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap11->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap12->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap13->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap14->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap15->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap16->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap17->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap18->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap19->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap20->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap21->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap22->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap23->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap24->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap25->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap26->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap27->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap28->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap29->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap30->SetVisAttributes(fCyanVisAtt);

  // Aluminum Casing (casing surrounding crystal)
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
    			    	
  fSolidAlCase2 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase2 = new G4LogicalVolume(fSolidAlCase2, fAlCaseMaterial, "AlCase");
  fPhysiAlCase2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase2, 
       				"AlCase", 
       				fLogicDetector2, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase3 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase3 = new G4LogicalVolume(fSolidAlCase3, fAlCaseMaterial, "AlCase");
  fPhysiAlCase3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase3, 
       				"AlCase", 
       				fLogicDetector3, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase4 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase4 = new G4LogicalVolume(fSolidAlCase4, fAlCaseMaterial, "AlCase");
  fPhysiAlCase4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase4, 
       				"AlCase", 
       				fLogicDetector4, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase5 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase5 = new G4LogicalVolume(fSolidAlCase5, fAlCaseMaterial, "AlCase");
  fPhysiAlCase5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase5, 
       				"AlCase", 
       				fLogicDetector5, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase6 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase6 = new G4LogicalVolume(fSolidAlCase6, fAlCaseMaterial, "AlCase");
  fPhysiAlCase6 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase6, 
       				"AlCase", 
       				fLogicDetector6, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase7 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase7 = new G4LogicalVolume(fSolidAlCase7, fAlCaseMaterial, "AlCase");
  fPhysiAlCase7 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase7, 
       				"AlCase", 
       				fLogicDetector7, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase8 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase8 = new G4LogicalVolume(fSolidAlCase8, fAlCaseMaterial, "AlCase");
  fPhysiAlCase8 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase8, 
       				"AlCase", 
       				fLogicDetector8, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase9 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase9 = new G4LogicalVolume(fSolidAlCase9, fAlCaseMaterial, "AlCase");
  fPhysiAlCase9 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase9, 
       				"AlCase", 
       				fLogicDetector9, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase10 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase10 = new G4LogicalVolume(fSolidAlCase10, fAlCaseMaterial, "AlCase");
  fPhysiAlCase10 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase10, 
       				"AlCase", 
       				fLogicDetector10, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase11 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase11 = new G4LogicalVolume(fSolidAlCase11, fAlCaseMaterial, "AlCase");
  fPhysiAlCase11 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase11, 
       				"AlCase", 
       				fLogicDetector11, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase12 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase12 = new G4LogicalVolume(fSolidAlCase12, fAlCaseMaterial, "AlCase");
  fPhysiAlCase12 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase12, 
       				"AlCase", 
       				fLogicDetector12, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase13 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase13 = new G4LogicalVolume(fSolidAlCase13, fAlCaseMaterial, "AlCase");
  fPhysiAlCase13 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase13, 
       				"AlCase", 
       				fLogicDetector13, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase14 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase14 = new G4LogicalVolume(fSolidAlCase14, fAlCaseMaterial, "AlCase");
  fPhysiAlCase14 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase14, 
       				"AlCase", 
       				fLogicDetector14, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase15 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase15 = new G4LogicalVolume(fSolidAlCase15, fAlCaseMaterial, "AlCase");
  fPhysiAlCase15 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase15, 
       				"AlCase", 
       				fLogicDetector15, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase16 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase16 = new G4LogicalVolume(fSolidAlCase16, fAlCaseMaterial, "AlCase");
  fPhysiAlCase16 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase16, 
       				"AlCase", 
       				fLogicDetector16, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase17 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase17 = new G4LogicalVolume(fSolidAlCase17, fAlCaseMaterial, "AlCase");
  fPhysiAlCase17 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase17, 
       				"AlCase", 
       				fLogicDetector17, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase18 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase18 = new G4LogicalVolume(fSolidAlCase18, fAlCaseMaterial, "AlCase");
  fPhysiAlCase18 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase18, 
       				"AlCase", 
       				fLogicDetector18, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase19 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase19 = new G4LogicalVolume(fSolidAlCase19, fAlCaseMaterial, "AlCase");
  fPhysiAlCase19 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase19, 
       				"AlCase", 
       				fLogicDetector19, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase20 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase20 = new G4LogicalVolume(fSolidAlCase20, fAlCaseMaterial, "AlCase");
  fPhysiAlCase20 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase20, 
       				"AlCase", 
       				fLogicDetector20, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase21 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase21 = new G4LogicalVolume(fSolidAlCase21, fAlCaseMaterial, "AlCase");
  fPhysiAlCase21 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase21, 
       				"AlCase", 
       				fLogicDetector21, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase22 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase22 = new G4LogicalVolume(fSolidAlCase22, fAlCaseMaterial, "AlCase");
  fPhysiAlCase22 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase22, 
       				"AlCase", 
       				fLogicDetector22, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase23 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase23 = new G4LogicalVolume(fSolidAlCase23, fAlCaseMaterial, "AlCase");
  fPhysiAlCase23 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase23, 
       				"AlCase", 
       				fLogicDetector23, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase24 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase24 = new G4LogicalVolume(fSolidAlCase24, fAlCaseMaterial, "AlCase");
  fPhysiAlCase24 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase24, 
       				"AlCase", 
       				fLogicDetector24, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase25 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase25 = new G4LogicalVolume(fSolidAlCase25, fAlCaseMaterial, "AlCase");
  fPhysiAlCase25 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase25, 
       				"AlCase", 
       				fLogicDetector25, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase26 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase26 = new G4LogicalVolume(fSolidAlCase26, fAlCaseMaterial, "AlCase");
  fPhysiAlCase26 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase26, 
       				"AlCase", 
       				fLogicDetector26, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase27 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase27 = new G4LogicalVolume(fSolidAlCase27, fAlCaseMaterial, "AlCase");
  fPhysiAlCase27 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase27, 
       				"AlCase", 
       				fLogicDetector27, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase28 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase28 = new G4LogicalVolume(fSolidAlCase28, fAlCaseMaterial, "AlCase");
  fPhysiAlCase28 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase28, 
       				"AlCase", 
       				fLogicDetector28, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase29 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase29 = new G4LogicalVolume(fSolidAlCase29, fAlCaseMaterial, "AlCase");
  fPhysiAlCase29 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase29, 
       				"AlCase", 
       				fLogicDetector29, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase30 = new G4Tubs("AlCase", 0.5*fGapDiameter, 0.5*fAlCaseDiameter, 0.5*fAlCaseLength, 0.*deg, 360.*deg);
  fLogicAlCase30 = new G4LogicalVolume(fSolidAlCase30, fAlCaseMaterial, "AlCase");
  fPhysiAlCase30 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase30, 
       				"AlCase", 
       				fLogicDetector30, 
       				false, 
    			    	0);
  
  // Aluminum Face Casing (casing covering face of crystal)
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
    			    	
  fSolidFaceAlCase2 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase2 = new G4LogicalVolume(fSolidFaceAlCase2, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase2, 
       				"FaceAlCase", 
       				fLogicDetector2, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase3 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase3 = new G4LogicalVolume(fSolidFaceAlCase3, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase3, 
       				"FaceAlCase", 
       				fLogicDetector3, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase4 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase4 = new G4LogicalVolume(fSolidFaceAlCase4, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase4, 
       				"FaceAlCase", 
       				fLogicDetector4, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase5 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase5 = new G4LogicalVolume(fSolidFaceAlCase5, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase5, 
       				"FaceAlCase", 
       				fLogicDetector5, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase6 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase6 = new G4LogicalVolume(fSolidFaceAlCase6, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase6 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase6, 
       				"FaceAlCase", 
       				fLogicDetector6, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase7 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase7 = new G4LogicalVolume(fSolidFaceAlCase7, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase7 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase7, 
       				"FaceAlCase", 
       				fLogicDetector7, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase8 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase8 = new G4LogicalVolume(fSolidFaceAlCase8, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase8 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase8, 
       				"FaceAlCase", 
       				fLogicDetector8, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase9 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase9 = new G4LogicalVolume(fSolidFaceAlCase9, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase9 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase9, 
       				"FaceAlCase", 
       				fLogicDetector9, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase10 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase10 = new G4LogicalVolume(fSolidFaceAlCase10, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase10 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase10, 
       				"FaceAlCase", 
       				fLogicDetector10, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase11 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase11 = new G4LogicalVolume(fSolidFaceAlCase11, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase11 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase11, 
       				"FaceAlCase", 
       				fLogicDetector11, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase12 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase12 = new G4LogicalVolume(fSolidFaceAlCase12, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase12 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase12, 
       				"FaceAlCase", 
       				fLogicDetector12, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase13 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase13 = new G4LogicalVolume(fSolidFaceAlCase13, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase13 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase13, 
       				"FaceAlCase", 
       				fLogicDetector13, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase14 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase14 = new G4LogicalVolume(fSolidFaceAlCase14, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase14 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase14, 
       				"FaceAlCase", 
       				fLogicDetector14, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase15 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase15 = new G4LogicalVolume(fSolidFaceAlCase15, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase15 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase15, 
       				"FaceAlCase", 
       				fLogicDetector15, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase16 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase16 = new G4LogicalVolume(fSolidFaceAlCase16, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase16 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase16, 
       				"FaceAlCase", 
       				fLogicDetector16, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase17 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase17 = new G4LogicalVolume(fSolidFaceAlCase17, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase17 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase17, 
       				"FaceAlCase", 
       				fLogicDetector17, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase18 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase18 = new G4LogicalVolume(fSolidFaceAlCase18, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase18 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase18, 
       				"FaceAlCase", 
       				fLogicDetector18, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase19 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase19 = new G4LogicalVolume(fSolidFaceAlCase19, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase19 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase19, 
       				"FaceAlCase", 
       				fLogicDetector19, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase20 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase20 = new G4LogicalVolume(fSolidFaceAlCase20, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase20 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase20, 
       				"FaceAlCase", 
       				fLogicDetector20, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase21 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase21 = new G4LogicalVolume(fSolidFaceAlCase21, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase21 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase21, 
       				"FaceAlCase", 
       				fLogicDetector21, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase22 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase22 = new G4LogicalVolume(fSolidFaceAlCase22, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase22 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase22, 
       				"FaceAlCase", 
       				fLogicDetector22, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase23 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase23 = new G4LogicalVolume(fSolidFaceAlCase23, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase23 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase23, 
       				"FaceAlCase", 
       				fLogicDetector23, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase24 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase24 = new G4LogicalVolume(fSolidFaceAlCase24, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase24 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase24, 
       				"FaceAlCase", 
       				fLogicDetector24, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase25 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase25 = new G4LogicalVolume(fSolidFaceAlCase25, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase25 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase25, 
       				"FaceAlCase", 
       				fLogicDetector25, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase26 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase26 = new G4LogicalVolume(fSolidFaceAlCase26, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase26 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase26, 
       				"FaceAlCase", 
       				fLogicDetector26, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase27 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase27 = new G4LogicalVolume(fSolidFaceAlCase27, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase27 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase27, 
       				"FaceAlCase", 
       				fLogicDetector27, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase28 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase28 = new G4LogicalVolume(fSolidFaceAlCase28, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase28 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase28, 
       				"FaceAlCase", 
       				fLogicDetector28, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase29 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase29 = new G4LogicalVolume(fSolidFaceAlCase29, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase29 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase29, 
       				"FaceAlCase", 
       				fLogicDetector29,
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase30 = new G4Tubs("FaceAlCase", 0., 0.5*fGapDiameter, 0.5*fAlFaceThickness, 0.*deg, 360.*deg);
  fLogicFaceAlCase30 = new G4LogicalVolume(fSolidFaceAlCase30, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase30 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase30, 
       				"FaceAlCase", 
       				fLogicDetector30, 
       				false, 
    			    	0);
    			    	
  fLogicAlCase->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase2->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase3->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase4->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase5->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase6->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase7->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase8->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase9->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase10->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase11->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase12->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase13->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase14->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase15->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase16->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase17->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase18->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase19->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase20->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase21->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase22->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase23->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase24->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase25->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase26->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase27->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase28->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase29->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase30->SetVisAttributes(fGreyVisAtt);
  
  fLogicFaceAlCase->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase2->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase3->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase4->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase5->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase6->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase7->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase8->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase9->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase10->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase11->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase12->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase13->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase14->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase15->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase16->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase17->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase18->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase19->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase20->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase21->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase22->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase23->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase24->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase25->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase26->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase27->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase28->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase29->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase30->SetVisAttributes(fGreyVisAtt);
  
  fLogicAlCase1->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase1->SetVisAttributes(fGreyVisAtt);
  
  // PMT
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
    			    	
  fSolidPMT2 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT2 = new G4LogicalVolume(fSolidPMT2, fPMTMaterial,"PMT");
  fPhysiPMT2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT2, 
       				"PMT", 
       				fLogicDetector2, 
       				false, 
    			    	0);
    			    	
  fSolidPMT3 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT3 = new G4LogicalVolume(fSolidPMT3, fPMTMaterial,"PMT");
  fPhysiPMT3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT3, 
       				"PMT", 
       				fLogicDetector3, 
       				false, 
    			    	0);
    			    	
  fSolidPMT4 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT4 = new G4LogicalVolume(fSolidPMT4, fPMTMaterial,"PMT");
  fPhysiPMT4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT4, 
       				"PMT", 
       				fLogicDetector4, 
       				false, 
    			    	0);
    			    	
  fSolidPMT5 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT5 = new G4LogicalVolume(fSolidPMT5, fPMTMaterial,"PMT");
  fPhysiPMT5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT5, 
       				"PMT", 
       				fLogicDetector5, 
       				false, 
    			    	0);
    			    	
  fSolidPMT6 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT6 = new G4LogicalVolume(fSolidPMT6, fPMTMaterial,"PMT");
  fPhysiPMT6 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT6, 
       				"PMT", 
       				fLogicDetector6, 
       				false, 
    			    	0);
    			    	
  fSolidPMT7 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT7 = new G4LogicalVolume(fSolidPMT7, fPMTMaterial,"PMT");
  fPhysiPMT7 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT7, 
       				"PMT", 
       				fLogicDetector7, 
       				false, 
    			    	0);
    			    	
  fSolidPMT8 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT8 = new G4LogicalVolume(fSolidPMT8, fPMTMaterial,"PMT");
  fPhysiPMT8 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT8, 
       				"PMT", 
       				fLogicDetector8, 
       				false, 
    			    	0);
    			    	
  fSolidPMT9 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT9 = new G4LogicalVolume(fSolidPMT9, fPMTMaterial,"PMT");
  fPhysiPMT9 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT9, 
       				"PMT", 
       				fLogicDetector9, 
       				false, 
    			    	0);
    			    	
  fSolidPMT10 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT10 = new G4LogicalVolume(fSolidPMT10, fPMTMaterial,"PMT");
  fPhysiPMT10 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT10, 
       				"PMT", 
       				fLogicDetector10, 
       				false, 
    			    	0);
    			    	
  fSolidPMT11 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT11 = new G4LogicalVolume(fSolidPMT11, fPMTMaterial,"PMT");
  fPhysiPMT11 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT11, 
       				"PMT", 
       				fLogicDetector11, 
       				false, 
    			    	0);
    			    	
  fSolidPMT12 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT12 = new G4LogicalVolume(fSolidPMT12, fPMTMaterial,"PMT");
  fPhysiPMT12 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT12, 
       				"PMT", 
       				fLogicDetector12, 
       				false, 
    			    	0);
    			    	
  fSolidPMT13 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT13 = new G4LogicalVolume(fSolidPMT13, fPMTMaterial,"PMT");
  fPhysiPMT13 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT13, 
       				"PMT", 
       				fLogicDetector13, 
       				false, 
    			    	0);
    			    	
  fSolidPMT14 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT14 = new G4LogicalVolume(fSolidPMT14, fPMTMaterial,"PMT");
  fPhysiPMT14 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT14, 
       				"PMT", 
       				fLogicDetector14, 
       				false, 
    			    	0);
    			    	
  fSolidPMT15 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT15 = new G4LogicalVolume(fSolidPMT15, fPMTMaterial,"PMT");
  fPhysiPMT15 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT15, 
       				"PMT", 
       				fLogicDetector15, 
       				false, 
    			    	0);
    			    	
  fSolidPMT16 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT16 = new G4LogicalVolume(fSolidPMT16, fPMTMaterial,"PMT");
  fPhysiPMT16 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT16, 
       				"PMT", 
       				fLogicDetector16, 
       				false, 
    			    	0);
    			    	
  fSolidPMT17 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT17 = new G4LogicalVolume(fSolidPMT17, fPMTMaterial,"PMT");
  fPhysiPMT17 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT17, 
       				"PMT", 
       				fLogicDetector17, 
       				false, 
    			    	0);
    			    	
  fSolidPMT18 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT18 = new G4LogicalVolume(fSolidPMT18, fPMTMaterial,"PMT");
  fPhysiPMT18 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT18, 
       				"PMT", 
       				fLogicDetector18, 
       				false, 
    			    	0);
    			    	
  fSolidPMT19 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT19 = new G4LogicalVolume(fSolidPMT19, fPMTMaterial,"PMT");
  fPhysiPMT19 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT19, 
       				"PMT", 
       				fLogicDetector19, 
       				false, 
    			    	0);
    			    	
  fSolidPMT20 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT20 = new G4LogicalVolume(fSolidPMT20, fPMTMaterial,"PMT");
  fPhysiPMT20 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT20, 
       				"PMT", 
       				fLogicDetector20, 
       				false, 
    			    	0);
    			    	
  fSolidPMT21 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT21 = new G4LogicalVolume(fSolidPMT21, fPMTMaterial,"PMT");
  fPhysiPMT21 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT21, 
       				"PMT", 
       				fLogicDetector21, 
       				false, 
    			    	0);
    			    	
  fSolidPMT22 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT22 = new G4LogicalVolume(fSolidPMT22, fPMTMaterial,"PMT");
  fPhysiPMT22 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT22, 
       				"PMT", 
       				fLogicDetector22, 
       				false, 
    			    	0);
    			    	
  fSolidPMT23 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT23 = new G4LogicalVolume(fSolidPMT23, fPMTMaterial,"PMT");
  fPhysiPMT23 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT23, 
       				"PMT", 
       				fLogicDetector23, 
       				false, 
    			    	0);
    			    	
  fSolidPMT24 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT24 = new G4LogicalVolume(fSolidPMT24, fPMTMaterial,"PMT");
  fPhysiPMT24 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT24, 
       				"PMT", 
       				fLogicDetector24, 
       				false, 
    			    	0);
    			    	
  fSolidPMT25 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT25 = new G4LogicalVolume(fSolidPMT25, fPMTMaterial,"PMT");
  fPhysiPMT25 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT25, 
       				"PMT", 
       				fLogicDetector25, 
       				false, 
    			    	0);
    			    	
  fSolidPMT26 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT26 = new G4LogicalVolume(fSolidPMT26, fPMTMaterial,"PMT");
  fPhysiPMT26 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT26, 
       				"PMT", 
       				fLogicDetector26, 
       				false, 
    			    	0);
    			    	
  fSolidPMT27 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT27 = new G4LogicalVolume(fSolidPMT27, fPMTMaterial,"PMT");
  fPhysiPMT27 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT27, 
       				"PMT", 
       				fLogicDetector27, 
       				false, 
    			    	0);
    			    	
  fSolidPMT28 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT28 = new G4LogicalVolume(fSolidPMT28, fPMTMaterial,"PMT");
  fPhysiPMT28 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT28, 
       				"PMT", 
       				fLogicDetector28, 
       				false, 
    			    	0);
    			    	
  fSolidPMT29 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT29 = new G4LogicalVolume(fSolidPMT29, fPMTMaterial,"PMT");
  fPhysiPMT29 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT29, 
       				"PMT", 
       				fLogicDetector29, 
       				false, 
    			    	0);
    			    	
  fSolidPMT30 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg);
  fLogicPMT30 = new G4LogicalVolume(fSolidPMT30, fPMTMaterial,"PMT");
  fPhysiPMT30 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT30, 
       				"PMT", 
       				fLogicDetector30, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI = new G4LogicalVolume(fSolidPMTAI, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI, 
       				"Vacuum", 
       				fLogicPMT, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI2 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI2 = new G4LogicalVolume(fSolidPMTAI2, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI2, 
       				"Vacuum", 
       				fLogicPMT2, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI3 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI3 = new G4LogicalVolume(fSolidPMTAI3, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI3, 
       				"Vacuum", 
       				fLogicPMT3, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI4 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI4 = new G4LogicalVolume(fSolidPMTAI4, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI4, 
       				"Vacuum", 
       				fLogicPMT4, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI5 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI5 = new G4LogicalVolume(fSolidPMTAI5, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI5, 
       				"Vacuum", 
       				fLogicPMT5, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI6 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI6 = new G4LogicalVolume(fSolidPMTAI6, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI6 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI6, 
       				"Vacuum", 
       				fLogicPMT6, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI7 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI7 = new G4LogicalVolume(fSolidPMTAI7, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI7 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI7, 
       				"Vacuum", 
       				fLogicPMT7, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI8 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI8 = new G4LogicalVolume(fSolidPMTAI8, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI8 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI8, 
       				"Vacuum", 
       				fLogicPMT8, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI9 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI9 = new G4LogicalVolume(fSolidPMTAI9, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI9 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI9, 
       				"Vacuum", 
       				fLogicPMT9, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI10 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI10 = new G4LogicalVolume(fSolidPMTAI10, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI10 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI10, 
       				"Vacuum", 
       				fLogicPMT10, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI11 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI11 = new G4LogicalVolume(fSolidPMTAI11, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI11 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI11, 
       				"Vacuum", 
       				fLogicPMT11, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI12 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI12 = new G4LogicalVolume(fSolidPMTAI12, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI12 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI12, 
       				"Vacuum", 
       				fLogicPMT12, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI13 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI13 = new G4LogicalVolume(fSolidPMTAI13, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI13 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI13, 
       				"Vacuum", 
       				fLogicPMT13, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI14 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI14 = new G4LogicalVolume(fSolidPMTAI14, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI14 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI14, 
       				"Vacuum", 
       				fLogicPMT14, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI15 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI15 = new G4LogicalVolume(fSolidPMTAI15, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI15 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI15, 
       				"Vacuum", 
       				fLogicPMT15, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI16 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI16 = new G4LogicalVolume(fSolidPMTAI16, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI16 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI16, 
       				"Vacuum", 
       				fLogicPMT16, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI17 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI17 = new G4LogicalVolume(fSolidPMTAI17, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI17 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI17, 
       				"Vacuum", 
       				fLogicPMT17, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI18 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI18 = new G4LogicalVolume(fSolidPMTAI18, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI18 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI18, 
       				"Vacuum", 
       				fLogicPMT18, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI19 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI19 = new G4LogicalVolume(fSolidPMTAI19, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI19 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI19, 
       				"Vacuum", 
       				fLogicPMT19, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI20 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI20 = new G4LogicalVolume(fSolidPMTAI20, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI20 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI20, 
       				"Vacuum", 
       				fLogicPMT20, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI21 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI21 = new G4LogicalVolume(fSolidPMTAI21, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI21 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI21, 
       				"Vacuum", 
       				fLogicPMT21, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI22 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI22 = new G4LogicalVolume(fSolidPMTAI22, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI22 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI22, 
       				"Vacuum", 
       				fLogicPMT22, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI23 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI23 = new G4LogicalVolume(fSolidPMTAI23, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI23 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI23, 
       				"Vacuum", 
       				fLogicPMT23, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI24 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI24 = new G4LogicalVolume(fSolidPMTAI24, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI24 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI24, 
       				"Vacuum", 
       				fLogicPMT24, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI25 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI25 = new G4LogicalVolume(fSolidPMTAI25, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI25 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI25, 
       				"Vacuum", 
       				fLogicPMT25, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI26 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI26 = new G4LogicalVolume(fSolidPMTAI26, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI26 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI26, 
       				"Vacuum", 
       				fLogicPMT26, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI27 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI27 = new G4LogicalVolume(fSolidPMTAI27, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI27 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI27, 
       				"Vacuum", 
       				fLogicPMT27, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI28 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI28 = new G4LogicalVolume(fSolidPMTAI28, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI28 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI28, 
       				"Vacuum", 
       				fLogicPMT28, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI29 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI29 = new G4LogicalVolume(fSolidPMTAI29, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI29 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI29, 
       				"Vacuum", 
       				fLogicPMT29, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAI30 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg); 			    	
  fLogicPMTAI30 = new G4LogicalVolume(fSolidPMTAI30, fPMTIntMaterial,"Vacuum");
  fPhysiPMTAI30 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTAI30, 
       				"Vacuum", 
       				fLogicPMT30, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK = new G4LogicalVolume(fSolidPMTAK, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK, 
       				"Photocathode", 
       				fLogicPMT, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK2 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK2 = new G4LogicalVolume(fSolidPMTAK2, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK2, 
       				"Photocathode", 
       				fLogicPMT2, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK3 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK3 = new G4LogicalVolume(fSolidPMTAK3, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK3, 
       				"Photocathode", 
       				fLogicPMT3, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK4 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK4 = new G4LogicalVolume(fSolidPMTAK4, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK4, 
       				"Photocathode", 
       				fLogicPMT4, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK5 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK5 = new G4LogicalVolume(fSolidPMTAK5, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK5, 
       				"Photocathode", 
       				fLogicPMT5, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK6 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK6 = new G4LogicalVolume(fSolidPMTAK6, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK6 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK6, 
       				"Photocathode", 
       				fLogicPMT6, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK7 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK7 = new G4LogicalVolume(fSolidPMTAK7, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK7 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK7, 
       				"Photocathode", 
       				fLogicPMT7, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK8 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK8 = new G4LogicalVolume(fSolidPMTAK8, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK8 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK8, 
       				"Photocathode", 
       				fLogicPMT8, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK9 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK9 = new G4LogicalVolume(fSolidPMTAK9, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK9 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK9, 
       				"Photocathode", 
       				fLogicPMT9, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK10 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK10 = new G4LogicalVolume(fSolidPMTAK10, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK10 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK10, 
       				"Photocathode", 
       				fLogicPMT10, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK11 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK11 = new G4LogicalVolume(fSolidPMTAK11, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK11 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK11, 
       				"Photocathode", 
       				fLogicPMT11, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK12 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK12 = new G4LogicalVolume(fSolidPMTAK12, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK12 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK12, 
       				"Photocathode", 
       				fLogicPMT12, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK13 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK13 = new G4LogicalVolume(fSolidPMTAK13, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK13 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK13, 
       				"Photocathode", 
       				fLogicPMT13, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK14 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK14 = new G4LogicalVolume(fSolidPMTAK14, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK14 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK14, 
       				"Photocathode", 
       				fLogicPMT14, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK15 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK15 = new G4LogicalVolume(fSolidPMTAK15, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK15 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK15, 
       				"Photocathode", 
       				fLogicPMT15, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK16 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK16 = new G4LogicalVolume(fSolidPMTAK16, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK16 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK16, 
       				"Photocathode", 
       				fLogicPMT16, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK17 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK17 = new G4LogicalVolume(fSolidPMTAK17, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK17 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK17, 
       				"Photocathode", 
       				fLogicPMT17, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK18 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK18 = new G4LogicalVolume(fSolidPMTAK18, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK18 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK18, 
       				"Photocathode", 
       				fLogicPMT18, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK19 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK19 = new G4LogicalVolume(fSolidPMTAK19, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK19 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK19, 
       				"Photocathode", 
       				fLogicPMT19, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK20 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK20 = new G4LogicalVolume(fSolidPMTAK20, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK20 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK20, 
       				"Photocathode", 
       				fLogicPMT20, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK21 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK21 = new G4LogicalVolume(fSolidPMTAK21, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK21 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK21, 
       				"Photocathode", 
       				fLogicPMT21, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK22 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK22 = new G4LogicalVolume(fSolidPMTAK22, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK22 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK22, 
       				"Photocathode", 
       				fLogicPMT22, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK23 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK23 = new G4LogicalVolume(fSolidPMTAK23, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK23 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK23, 
       				"Photocathode", 
       				fLogicPMT23, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK24 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK24 = new G4LogicalVolume(fSolidPMTAK24, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK24 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK24, 
       				"Photocathode", 
       				fLogicPMT24, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK25 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK25 = new G4LogicalVolume(fSolidPMTAK25, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK25 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK25, 
       				"Photocathode", 
       				fLogicPMT25, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK26 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK26 = new G4LogicalVolume(fSolidPMTAK26, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK26 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK26, 
       				"Photocathode", 
       				fLogicPMT26, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK27 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK27 = new G4LogicalVolume(fSolidPMTAK27, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK27 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK27, 
       				"Photocathode", 
       				fLogicPMT27, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK28 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK28 = new G4LogicalVolume(fSolidPMTAK28, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK28 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK28, 
       				"Photocathode", 
       				fLogicPMT28, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK29 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK29 = new G4LogicalVolume(fSolidPMTAK29, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK29 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK29, 
       				"Photocathode", 
       				fLogicPMT29, 
       				false, 
    			    	0);
    			    	
  fSolidPMTAK30 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);			    	
  fLogicPMTAK30 = new G4LogicalVolume(fSolidPMTAK30, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTAK30 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTAK30, 
       				"Photocathode", 
       				fLogicPMT30, 
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
    			    	
  fSolidPMTI = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter1, 0.5*fVacuumLength1, 0.*deg, 360.*deg); 			    	
  fLogicPMTI = new G4LogicalVolume(fSolidPMTI, fPMTIntMaterial,"Vacuum");
  fPhysiPMTI = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTI, 
       				"Vacuum", 
       				fLogicPMT1, 
       				false, 
    			    	0);
    			    	
  fSolidPMTK = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter1, 0.5*fPhotoCathodeLength1, 0.*deg, 360.*deg);			    	
  fLogicPMTK = new G4LogicalVolume(fSolidPMTK, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTK = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode1), 
        			fLogicPMTK, 
       				"Photocathode", 
       				fLogicPMT1, 
       				false, 
    			    	0);
    			    	
  fLogicPMTK->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK2->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK3->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK4->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK5->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK6->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK7->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK8->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK9->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK10->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK11->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK12->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK13->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK14->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK15->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK16->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK17->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK18->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK19->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK20->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK21->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK22->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK23->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK24->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK25->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK26->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK27->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK28->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK29->SetVisAttributes(fYellowVisAtt);
  fLogicPMTAK30->SetVisAttributes(fYellowVisAtt);
  }
  
if (fDetectorGeometry == 2 || fDetectorGeometry == 4){
  //  Single Hexagonal Prism Detector	
  //  If rotation required  
  //  G4RotationMatrix rotm  = G4RotationMatrix(0,90*deg,-90*deg);     
  //  No Rotation Now
  //  BGO Dimensions
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
  // GRIFFIN LaBr3:Ce Dimensions   
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
  fSolidDetector1 = new G4Tubs("Detector", 0., 0.5*fTotalDetectorDiameter1, 0.5*fTotalDetectorLength1, 0.*deg, 360.*deg);
  fSolidHDetector2 = new G4Polyhedra{"BGO Detector 2", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector3 = new G4Polyhedra{"BGO Detector 3", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector4 = new G4Polyhedra{"BGO Detector 4", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector5 = new G4Polyhedra{"BGO Detector 5", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector6 = new G4Polyhedra{"BGO Detector 6", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector7 = new G4Polyhedra{"BGO Detector 7", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector8 = new G4Polyhedra{"BGO Detector 8", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector9 = new G4Polyhedra{"BGO Detector 9", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector10 = new G4Polyhedra{"BGO Detector 10", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector11 = new G4Polyhedra{"BGO Detector 11", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector12 = new G4Polyhedra{"BGO Detector 12", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector13 = new G4Polyhedra{"BGO Detector 13", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector14 = new G4Polyhedra{"BGO Detector 14", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector15 = new G4Polyhedra{"BGO Detector 15", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector16 = new G4Polyhedra{"BGO Detector 16", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector17 = new G4Polyhedra{"BGO Detector 17", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector18 = new G4Polyhedra{"BGO Detector 18", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector19 = new G4Polyhedra{"BGO Detector 19", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector20 = new G4Polyhedra{"BGO Detector 20", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector21 = new G4Polyhedra{"BGO Detector 21", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector22 = new G4Polyhedra{"BGO Detector 22", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector23 = new G4Polyhedra{"BGO Detector 23", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector24 = new G4Polyhedra{"BGO Detector 24", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector25 = new G4Polyhedra{"BGO Detector 25", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector26 = new G4Polyhedra{"BGO Detector 26", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector27 = new G4Polyhedra{"BGO Detector 27", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector28 = new G4Polyhedra{"BGO Detector 28", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector29 = new G4Polyhedra{"BGO Detector 29", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector30 = new G4Polyhedra{"BGO Detector 30", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fLogicDetector = new G4LogicalVolume(fSolidHDetector,fDetectorMaterial, "BGO Detector 1");
  fLogicDetector1 = new G4LogicalVolume(fSolidDetector1,fDetectorMaterial, "LaBr3:Ce Detector");
  fLogicDetector2 = new G4LogicalVolume(fSolidHDetector2,fDetectorMaterial, "BGO Detector 2");
  fLogicDetector3 = new G4LogicalVolume(fSolidHDetector3,fDetectorMaterial, "BGO Detector 3");
  fLogicDetector4 = new G4LogicalVolume(fSolidHDetector4,fDetectorMaterial, "BGO Detector 4");
  fLogicDetector5 = new G4LogicalVolume(fSolidHDetector5,fDetectorMaterial, "BGO Detector 5");
  fLogicDetector6 = new G4LogicalVolume(fSolidHDetector6,fDetectorMaterial, "BGO Detector 6");
  fLogicDetector7 = new G4LogicalVolume(fSolidHDetector7,fDetectorMaterial, "BGO Detector 7");
  fLogicDetector8 = new G4LogicalVolume(fSolidHDetector8,fDetectorMaterial, "BGO Detector 8");
  fLogicDetector9 = new G4LogicalVolume(fSolidHDetector9,fDetectorMaterial, "BGO Detector 9");
  fLogicDetector10 = new G4LogicalVolume(fSolidHDetector10,fDetectorMaterial, "BGO Detector 10");
  fLogicDetector11 = new G4LogicalVolume(fSolidHDetector11,fDetectorMaterial, "BGO Detector 11");
  fLogicDetector12 = new G4LogicalVolume(fSolidHDetector12,fDetectorMaterial, "BGO Detector 12");
  fLogicDetector13 = new G4LogicalVolume(fSolidHDetector13,fDetectorMaterial, "BGO Detector 13");
  fLogicDetector14 = new G4LogicalVolume(fSolidHDetector14,fDetectorMaterial, "BGO Detector 14");
  fLogicDetector15 = new G4LogicalVolume(fSolidHDetector15,fDetectorMaterial, "BGO Detector 15");
  fLogicDetector16 = new G4LogicalVolume(fSolidHDetector16,fDetectorMaterial, "BGO Detector 16");
  fLogicDetector17 = new G4LogicalVolume(fSolidHDetector17,fDetectorMaterial, "BGO Detector 17");
  fLogicDetector18 = new G4LogicalVolume(fSolidHDetector18,fDetectorMaterial, "BGO Detector 18");
  fLogicDetector19 = new G4LogicalVolume(fSolidHDetector19,fDetectorMaterial, "BGO Detector 19");
  fLogicDetector20 = new G4LogicalVolume(fSolidHDetector20,fDetectorMaterial, "BGO Detector 20");
  fLogicDetector21 = new G4LogicalVolume(fSolidHDetector21,fDetectorMaterial, "BGO Detector 21");
  fLogicDetector22 = new G4LogicalVolume(fSolidHDetector22,fDetectorMaterial, "BGO Detector 22");
  fLogicDetector23 = new G4LogicalVolume(fSolidHDetector23,fDetectorMaterial, "BGO Detector 23");
  fLogicDetector24 = new G4LogicalVolume(fSolidHDetector24,fDetectorMaterial, "BGO Detector 24");
  fLogicDetector25 = new G4LogicalVolume(fSolidHDetector25,fDetectorMaterial, "BGO Detector 25");
  fLogicDetector26 = new G4LogicalVolume(fSolidHDetector26,fDetectorMaterial, "BGO Detector 26");
  fLogicDetector27 = new G4LogicalVolume(fSolidHDetector27,fDetectorMaterial, "BGO Detector 27");
  fLogicDetector28 = new G4LogicalVolume(fSolidHDetector28,fDetectorMaterial, "BGO Detector 28");
  fLogicDetector29 = new G4LogicalVolume(fSolidHDetector29,fDetectorMaterial, "BGO Detector 29");
  fLogicDetector30 = new G4LogicalVolume(fSolidHDetector30,fDetectorMaterial, "BGO Detector 30");
  fPhysiDetector = new G4PVPlacement(transform,
        			fLogicDetector, 
       				"BGO Detector", 
       				fLogicWorld, 
       				false, 
    			    	0,
				    false);
  //fPhysiDetector1 = new G4PVPlacement(transform,
        			//fLogicDetector1, 
       				//"LaBr3:Ce Detector", 
       				//fLogicWorld, 
       				//false, 
    			    	//0,
				    //false);
  fLogicDetector->SetVisAttributes(fWireFrameVisAtt);
  fLogicDetector1->SetVisAttributes(fWireFrameVisAtt);
  }
  if (fDetectorGeometry == 4){
  // Hexagonal Prism Detector Array 
  // Each detector is treated as an individual volume for gamma ray detection. 
  fSolidHDetector = new G4Polyhedra{"BGO Detector", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidDetector1 = new G4Tubs("Detector", 0., 0.5*fTotalDetectorDiameter1, 0.5*fTotalDetectorLength1, 0.*deg, 360.*deg);
  fSolidHDetector2 = new G4Polyhedra{"BGO Detector 2", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector3 = new G4Polyhedra{"BGO Detector 3", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector4 = new G4Polyhedra{"BGO Detector 4", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector5 = new G4Polyhedra{"BGO Detector 5", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector6 = new G4Polyhedra{"BGO Detector 6", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector7 = new G4Polyhedra{"BGO Detector 7", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector8 = new G4Polyhedra{"BGO Detector 8", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector9 = new G4Polyhedra{"BGO Detector 9", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector10 = new G4Polyhedra{"BGO Detector 10", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector11 = new G4Polyhedra{"BGO Detector 11", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector12 = new G4Polyhedra{"BGO Detector 12", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector13 = new G4Polyhedra{"BGO Detector 13", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector14 = new G4Polyhedra{"BGO Detector 14", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector15 = new G4Polyhedra{"BGO Detector 15", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector16 = new G4Polyhedra{"BGO Detector 16", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector17 = new G4Polyhedra{"BGO Detector 17", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector18 = new G4Polyhedra{"BGO Detector 18", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector19 = new G4Polyhedra{"BGO Detector 19", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector20 = new G4Polyhedra{"BGO Detector 20", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector21 = new G4Polyhedra{"BGO Detector 21", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector22 = new G4Polyhedra{"BGO Detector 22", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector23 = new G4Polyhedra{"BGO Detector 23", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector24 = new G4Polyhedra{"BGO Detector 24", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector25 = new G4Polyhedra{"BGO Detector 25", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector26 = new G4Polyhedra{"BGO Detector 26", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector27 = new G4Polyhedra{"BGO Detector 27", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector28 = new G4Polyhedra{"BGO Detector 28", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector29 = new G4Polyhedra{"BGO Detector 29", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fSolidHDetector30 = new G4Polyhedra{"BGO Detector 30", 0.*deg, 360.*deg, 6, 2, zPlane, rInner, rOuter};
  fLogicDetector = new G4LogicalVolume(fSolidHDetector,fDetectorMaterial, "BGO Detector 1");
  fLogicDetector1 = new G4LogicalVolume(fSolidDetector1,fDetectorMaterial, "LaBr3:Ce Detector");
  fLogicDetector2 = new G4LogicalVolume(fSolidHDetector2,fDetectorMaterial, "BGO Detector 2");
  fLogicDetector3 = new G4LogicalVolume(fSolidHDetector3,fDetectorMaterial, "BGO Detector 3");
  fLogicDetector4 = new G4LogicalVolume(fSolidHDetector4,fDetectorMaterial, "BGO Detector 4");
  fLogicDetector5 = new G4LogicalVolume(fSolidHDetector5,fDetectorMaterial, "BGO Detector 5");
  fLogicDetector6 = new G4LogicalVolume(fSolidHDetector6,fDetectorMaterial, "BGO Detector 6");
  fLogicDetector7 = new G4LogicalVolume(fSolidHDetector7,fDetectorMaterial, "BGO Detector 7");
  fLogicDetector8 = new G4LogicalVolume(fSolidHDetector8,fDetectorMaterial, "BGO Detector 8");
  fLogicDetector9 = new G4LogicalVolume(fSolidHDetector9,fDetectorMaterial, "BGO Detector 9");
  fLogicDetector10 = new G4LogicalVolume(fSolidHDetector10,fDetectorMaterial, "BGO Detector 10");
  fLogicDetector11 = new G4LogicalVolume(fSolidHDetector11,fDetectorMaterial, "BGO Detector 11");
  fLogicDetector12 = new G4LogicalVolume(fSolidHDetector12,fDetectorMaterial, "BGO Detector 12");
  fLogicDetector13 = new G4LogicalVolume(fSolidHDetector13,fDetectorMaterial, "BGO Detector 13");
  fLogicDetector14 = new G4LogicalVolume(fSolidHDetector14,fDetectorMaterial, "BGO Detector 14");
  fLogicDetector15 = new G4LogicalVolume(fSolidHDetector15,fDetectorMaterial, "BGO Detector 15");
  fLogicDetector16 = new G4LogicalVolume(fSolidHDetector16,fDetectorMaterial, "BGO Detector 16");
  fLogicDetector17 = new G4LogicalVolume(fSolidHDetector17,fDetectorMaterial, "BGO Detector 17");
  fLogicDetector18 = new G4LogicalVolume(fSolidHDetector18,fDetectorMaterial, "BGO Detector 18");
  fLogicDetector19 = new G4LogicalVolume(fSolidHDetector19,fDetectorMaterial, "BGO Detector 19");
  fLogicDetector20 = new G4LogicalVolume(fSolidHDetector20,fDetectorMaterial, "BGO Detector 20");
  fLogicDetector21 = new G4LogicalVolume(fSolidHDetector21,fDetectorMaterial, "BGO Detector 21");
  fLogicDetector22 = new G4LogicalVolume(fSolidHDetector22,fDetectorMaterial, "BGO Detector 22");
  fLogicDetector23 = new G4LogicalVolume(fSolidHDetector23,fDetectorMaterial, "BGO Detector 23");
  fLogicDetector24 = new G4LogicalVolume(fSolidHDetector24,fDetectorMaterial, "BGO Detector 24");
  fLogicDetector25 = new G4LogicalVolume(fSolidHDetector25,fDetectorMaterial, "BGO Detector 25");
  fLogicDetector26 = new G4LogicalVolume(fSolidHDetector26,fDetectorMaterial, "BGO Detector 26");
  fLogicDetector27 = new G4LogicalVolume(fSolidHDetector27,fDetectorMaterial, "BGO Detector 27");
  fLogicDetector28 = new G4LogicalVolume(fSolidHDetector28,fDetectorMaterial, "BGO Detector 28");
  fLogicDetector29 = new G4LogicalVolume(fSolidHDetector29,fDetectorMaterial, "BGO Detector 29");
  fLogicDetector30 = new G4LogicalVolume(fSolidHDetector30,fDetectorMaterial, "BGO Detector 30");
  
  G4AssemblyVolume* assemblyDetector = new G4AssemblyVolume();
  // 30 detectors are placed relative to the array volume.
  // Left-Hand Array, rotm180
  // 25
  //P.setX(-7.68*cm); P.setY(-5.91*cm); P.setZ(17.849*cm);
  P.setX((-5/(2*sqrt(3)))*fTotalDetectorDiameter - fAirGap); P.setY(-fTotalDetectorDiameter - fAirGap); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector25,Tr);
  // 19
  //P.setX(-7.68*cm); P.setY(0.00*cm); P.setZ(17.849*cm);
  P.setX((-5/(2*sqrt(3)))*fTotalDetectorDiameter - fAirGap); P.setY(0.); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector19,Tr);
  // 13
  //P.setX(-7.68*cm); P.setY(5.91*cm); P.setZ(17.849*cm);
  P.setX((-5/(2*sqrt(3)))*fTotalDetectorDiameter - fAirGap); P.setY(fTotalDetectorDiameter + fAirGap); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector13,Tr);
  // 29
  //P.setX(-2.56*cm); P.setY(-8.87*cm); P.setZ(17.849*cm);
  P.setX((-1/sqrt(3))*fTotalDetectorDiameter); P.setY(-1.5*fTotalDetectorDiameter - 1.5*fAirGap); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector29,Tr);
  // 23
  //P.setX(-2.56*cm); P.setY(-2.96*cm); P.setZ(17.849*cm);
  P.setX((-1/sqrt(3))*fTotalDetectorDiameter); P.setY(-0.5*fTotalDetectorDiameter - 0.5*fAirGap); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector23,Tr);
  // 17
  //P.setX(-2.56*cm); P.setY(2.96*cm); P.setZ(17.849*cm);
  P.setX((-1/sqrt(3))*fTotalDetectorDiameter); P.setY(0.5*fTotalDetectorDiameter + 0.5*fAirGap); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector17,Tr);
  // 11
  //P.setX(-2.56*cm); P.setY(8.87*cm); P.setZ(17.849*cm);
  P.setX((-1/sqrt(3))*fTotalDetectorDiameter); P.setY(1.5*fTotalDetectorDiameter + 1.5*fAirGap); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector11,Tr);
  // 27
  //P.setX(2.56*cm); P.setY(-5.91*cm); P.setZ(17.849*cm);
  P.setX((1/(2*sqrt(3)))*fTotalDetectorDiameter + fAirGap); P.setY(-fTotalDetectorDiameter - fAirGap); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector27,Tr);
  // 21
  //P.setX(2.56*cm); P.setY(0.00*cm); P.setZ(17.849*cm);
  P.setX((1/(2*sqrt(3)))*fTotalDetectorDiameter + fAirGap); P.setY(0.); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector21,Tr);
  // 15
  //P.setX(2.56*cm); P.setY(5.91*cm); P.setZ(17.849*cm);
  P.setX((1/(2*sqrt(3)))*fTotalDetectorDiameter + fAirGap); P.setY(fTotalDetectorDiameter + fAirGap); P.setZ(0.5*fTotalDetectorLength + fOutTBoxY);
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector15,Tr);
  // Straddling Detectors
  // 7
  //P.setX(10.08*cm); P.setY(-8.87*cm); P.setZ(9*cm);
  if (fPMTDiameter < 3.175*sqrt(3)*cm){
  P.setX(3.175*cm + fTotalDetectorDiameter*((2 + sqrt(3))/(2*sqrt(3))) + 2.5*fAirGap); P.setY(-1.5*fTotalDetectorDiameter - (0.5 + 0.5*sqrt(3))*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  else if (fPMTDiameter >= 3.175*sqrt(3)*cm){
  P.setX(((4 + sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter + 2.5*fAirGap); P.setY(-1.5*fTotalDetectorDiameter - (0.5 + 0.5*sqrt(3))*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector7,Tr);
  // 6
  //P.setX(7.68*cm); P.setY(-2.96*cm); P.setZ(9*cm);
  if (fPMTDiameter < 3.175*sqrt(3)*cm){
  P.setX(3.175*cm + fTotalDetectorDiameter/sqrt(3) + 2*fAirGap); P.setY(-0.5*fTotalDetectorDiameter - 0.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  else if (fPMTDiameter >= 3.175*sqrt(3)*cm){
  P.setX((2/sqrt(3))*fTotalDetectorDiameter + 2*fAirGap); P.setY(-0.5*fTotalDetectorDiameter - 0.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector6,Tr);
  // 5
  //P.setX(7.68*cm); P.setY(2.96*cm); P.setZ(9*cm);
  if (fPMTDiameter < 3.175*sqrt(3)*cm){ 
  P.setX(3.175*cm + fTotalDetectorDiameter/sqrt(3) + 2*fAirGap); P.setY(0.5*fTotalDetectorDiameter + 0.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);  
  }
  else if (fPMTDiameter >= 3.175*sqrt(3)*cm){
  P.setX((2/sqrt(3))*fTotalDetectorDiameter + 2*fAirGap); P.setY(0.5*fTotalDetectorDiameter + 0.5*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector5,Tr); 
  // 4
  //P.setX(10.08*cm); P.setY(8.87*cm); P.setZ(9*cm);
  if (fPMTDiameter < 3.175*sqrt(3)*cm){
  P.setX(3.175*cm + fTotalDetectorDiameter*((2 + sqrt(3))/(2*sqrt(3))) + 2.5*fAirGap); P.setY(1.5*fTotalDetectorDiameter + (0.5 + 0.5*sqrt(3))*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  else if (fPMTDiameter >= 3.175*sqrt(3)*cm){
  P.setX(((4 + sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter + 2.5*fAirGap); P.setY(1.5*fTotalDetectorDiameter + (0.5 + 0.5*sqrt(3))*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector4,Tr);
  // 3
  //P.setX(4.96*cm); P.setY(11.83*cm); P.setZ(2.299*cm); 
  if (fPMTDiameter < 2*fOutTBoxX/3){
  P.setX(((1 + sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter + 1.5*fAirGap); P.setY(fOutTBoxX + 0.5*fTotalDetectorDiameter + (1 + 0.5*sqrt(3))*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY - 6.701*cm);
  }
  else if (fPMTDiameter >= 2*fOutTBoxX/3){
  P.setX(((1 + sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter + 1.5*fAirGap); P.setY(2*fTotalDetectorDiameter + (1 + 0.5*sqrt(3))*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY - 6.701*cm);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector3,Tr);
  // 9
  //P.setX(4.96*cm); P.setY(-11.83*cm); P.setZ(9*cm);
  if (fPMTDiameter < 2*fOutTBoxX/3){
  P.setX(((1 + sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter + 1.5*fAirGap); P.setY(-fOutTBoxX - 0.5*fTotalDetectorDiameter - (1 + 0.5*sqrt(3))*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  else if (fPMTDiameter >= 2*fOutTBoxX/3){
  P.setX(((1 + sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter + 1.5*fAirGap); P.setY(-2*fTotalDetectorDiameter - (1 + 0.5*sqrt(3))*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector9,Tr);
  // 1
  //P.setX(-4.96*cm); P.setY(14.78*cm); P.setZ(1.8*cm);
  if (fPMTDiameter < 2*fOutTBoxX/3){
  P.setX(((-2 - sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter - 0.5*fAirGap); P.setY(fOutTBoxX + fTotalDetectorDiameter + (1.5 + 0.5*sqrt(3))*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY - 7.2*cm);  
  }
  else if (fPMTDiameter >= 2*fOutTBoxX/3){
  P.setX(((-2 - sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter - 0.5*fAirGap); P.setY(2.5*fTotalDetectorDiameter + (1.5 + 0.5*sqrt(3))*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY - 7.2*cm);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector,Tr);
  // 2
  //P.setX(-10.08*cm); P.setY(11.83*cm); P.setZ(9*cm);
  if (fPMTDiameter < 2*fOutTBoxX/3){
  P.setX(((-5 - sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter - 1.5*fAirGap); P.setY(fOutTBoxX + 0.5*fTotalDetectorDiameter + (1 + 0.5*sqrt(3))*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  else if (fPMTDiameter >= 2*fOutTBoxX/3){
  P.setX(((-5 - sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter - 1.5*fAirGap); P.setY(2*fTotalDetectorDiameter + (1 + 0.5*sqrt(3))*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector2,Tr);
  // 10
  //P.setX(-4.96*cm); P.setY(-14.78*cm); P.setZ(9*cm);
  if (fPMTDiameter < 2*fOutTBoxX/3){
  P.setX(((-2 - sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter - 0.5*fAirGap); P.setY(-fOutTBoxX - fTotalDetectorDiameter - (1.5 + 0.5*sqrt(3))*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  else if (fPMTDiameter >= 2*fOutTBoxX/3){
  P.setX(((-2 - sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter - 0.5*fAirGap); P.setY(-2.5*fTotalDetectorDiameter - (1.5 + 0.5*sqrt(3))*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector10,Tr);
  // 8
  //P.setX(-10.08*cm); P.setY(-11.83*cm); P.setZ(9*cm);
  if (fPMTDiameter < 2*fOutTBoxX/3){
  P.setX(((-5 - sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter - 1.5*fAirGap); P.setY(-fOutTBoxX - 0.5*fTotalDetectorDiameter - (1 + 0.5*sqrt(3))*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  else if (fPMTDiameter >= 2*fOutTBoxX/3){
  P.setX(((-5 - sqrt(3))/(2*sqrt(3)))*fTotalDetectorDiameter - 1.5*fAirGap); P.setY(-2*fTotalDetectorDiameter - (1 + 0.5*sqrt(3))*fAirGap); P.setZ(0.5*fTotalDetectorLength - fOutTBoxY);
  }
  Tr = G4Transform3D(rotm180,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector8,Tr);
  // Right-Hand Array, rotm 
  // 26
  //P.setX(-7.68*cm); P.setY(-5.91*cm); P.setZ(-17.849*cm);
  P.setX((-5/(2*sqrt(3)))*fTotalDetectorDiameter - fAirGap); P.setY(-fTotalDetectorDiameter - fAirGap); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector26,Tr);
  // 20
  //P.setX(-7.68*cm); P.setY(0.00*cm); P.setZ(-17.849*cm);
  P.setX((-5/(2*sqrt(3)))*fTotalDetectorDiameter - fAirGap); P.setY(0.); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector20,Tr);
  // 14
  //P.setX(-7.68*cm); P.setY(5.91*cm); P.setZ(-17.849*cm);
  P.setX((-5/(2*sqrt(3)))*fTotalDetectorDiameter - fAirGap); P.setY(fTotalDetectorDiameter + fAirGap); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector14,Tr);
  // 30
  //P.setX(-2.56*cm); P.setY(-8.87*cm); P.setZ(-17.849*cm);
  P.setX((-1/sqrt(3))*fTotalDetectorDiameter); P.setY(-1.5*fTotalDetectorDiameter - 1.5*fAirGap); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector30,Tr);
  // 24
  //P.setX(-2.56*cm); P.setY(-2.96*cm); P.setZ(-17.849*cm);
  P.setX((-1/sqrt(3))*fTotalDetectorDiameter); P.setY(-0.5*fTotalDetectorDiameter - 0.5*fAirGap); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector24,Tr);
  // 18
  //P.setX(-2.56*cm); P.setY(2.96*cm); P.setZ(-17.849*cm);
  P.setX((-1/sqrt(3))*fTotalDetectorDiameter); P.setY(0.5*fTotalDetectorDiameter + 0.5*fAirGap); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector18,Tr);
  // 12
  //P.setX(-2.56*cm); P.setY(8.87*cm); P.setZ(-17.849*cm);
  P.setX((-1/sqrt(3))*fTotalDetectorDiameter); P.setY(1.5*fTotalDetectorDiameter + 1.5*fAirGap); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector12,Tr);
  // 28
  //P.setX(2.56*cm); P.setY(-5.91*cm); P.setZ(-17.849*cm);
  P.setX((1/(2*sqrt(3)))*fTotalDetectorDiameter + fAirGap); P.setY(-fTotalDetectorDiameter - fAirGap); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector28,Tr);
  // 22
  //P.setX(2.56*cm); P.setY(0.00*cm); P.setZ(-17.849*cm);
  P.setX((1/(2*sqrt(3)))*fTotalDetectorDiameter + fAirGap); P.setY(0.); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector22,Tr);
  // 16
  //P.setX(2.56*cm); P.setY(5.91*cm); P.setZ(-17.849*cm);
  P.setX((1/(2*sqrt(3)))*fTotalDetectorDiameter + fAirGap); P.setY(fTotalDetectorDiameter + fAirGap); P.setZ(-0.5*fTotalDetectorLength - fOutTBoxY);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector16,Tr);
  // Assembly Placement in the World
  // Rotated 270° along the x-axis to be viewed from the same perspective as the cylindrical array
  // Rotated 90° along the z-axis to be perpendicular to the gas target.
  // Rotated 180° along the z-axis to place the array downside up instead of upside down.
  // Coordinate transformation: (X, Y, Z) -> (Y, -Z, -X)
  G4ThreeVector WorldP(0,0,0);
  G4RotationMatrix worldrotm = G4RotationMatrix(0.*deg,90.*deg,270.*deg);
  G4Transform3D WorldTr = G4Transform3D(worldrotm, WorldP);
  assemblyDetector->MakeImprint(fLogicWorld, WorldTr);
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
    			    	
  fSolidPMTWin1 = new G4Tubs("PMTWin", 0., 0.5*fDetectorDiameter1, 0.5*fPMTWinThickness1, 0.*deg, 360.*deg);
  fLogicPMTWin1 = new G4LogicalVolume(fSolidPMTWin1, fPMTWinMaterial, "PMTWin");
  fPhysiPMTWin1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMTWin1), 
        			fLogicPMTWin1, 
       				"PMTWin", 
       				fLogicDetector1, 
       				false, 
    			    	0);
    			    	
  //fLogicPMTWin->SetVisAttributes(fGreenVisAtt); 
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
    			    	
  fSolidHCrystal2 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal2 = new G4LogicalVolume(fSolidHCrystal2, fCrystalMaterial, "Crystal");
  fPhysiCrystal2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal2, 
       				"Crystal", 
       				fLogicDetector2, 
       				false, 
    			    	2);
    			    	
  fSolidHCrystal3 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal3 = new G4LogicalVolume(fSolidHCrystal3, fCrystalMaterial, "Crystal");
  fPhysiCrystal3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal3, 
       				"Crystal", 
       				fLogicDetector3, 
       				false, 
    			    	3);
    			    	
  fSolidHCrystal4 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal4 = new G4LogicalVolume(fSolidHCrystal4, fCrystalMaterial, "Crystal");
  fPhysiCrystal4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal4, 
       				"Crystal", 
       				fLogicDetector4, 
       				false, 
    			    	4);
    			    	
  fSolidHCrystal5 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal5 = new G4LogicalVolume(fSolidHCrystal5, fCrystalMaterial, "Crystal");
  fPhysiCrystal5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal5, 
       				"Crystal", 
       				fLogicDetector5, 
       				false, 
    			    	5);
    			    	
  fSolidHCrystal6 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal6 = new G4LogicalVolume(fSolidHCrystal6, fCrystalMaterial, "Crystal");
  fPhysiCrystal6 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal6, 
       				"Crystal", 
       				fLogicDetector6, 
       				false, 
    			    	6);
    			    	
  fSolidHCrystal7 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal7 = new G4LogicalVolume(fSolidHCrystal7, fCrystalMaterial, "Crystal");
  fPhysiCrystal7 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal7, 
       				"Crystal", 
       				fLogicDetector7, 
       				false, 
    			    	7);
    			    	
  fSolidHCrystal8 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal8 = new G4LogicalVolume(fSolidHCrystal8, fCrystalMaterial, "Crystal");
  fPhysiCrystal8 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal8, 
       				"Crystal", 
       				fLogicDetector8, 
       				false, 
    			    	8);
    			    	
  fSolidHCrystal9 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal9 = new G4LogicalVolume(fSolidHCrystal9, fCrystalMaterial, "Crystal");
  fPhysiCrystal9 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal9, 
       				"Crystal", 
       				fLogicDetector9, 
       				false, 
    			    	9);
    			    	
  fSolidHCrystal10 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal10 = new G4LogicalVolume(fSolidHCrystal10, fCrystalMaterial, "Crystal");
  fPhysiCrystal10 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal10, 
       				"Crystal", 
       				fLogicDetector10, 
       				false, 
    			    	10);
    			    	
  fSolidHCrystal11 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal11 = new G4LogicalVolume(fSolidHCrystal11, fCrystalMaterial, "Crystal");
  fPhysiCrystal11 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal11, 
       				"Crystal", 
       				fLogicDetector11, 
       				false, 
    			    	11);
    			    	
  fSolidHCrystal12 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal12 = new G4LogicalVolume(fSolidHCrystal12, fCrystalMaterial, "Crystal");
  fPhysiCrystal12 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal12, 
       				"Crystal", 
       				fLogicDetector12, 
       				false, 
    			    	12);
  
  fSolidHCrystal13 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal13 = new G4LogicalVolume(fSolidHCrystal13, fCrystalMaterial, "Crystal");
  fPhysiCrystal13 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal13, 
       				"Crystal", 
       				fLogicDetector13, 
       				false, 
    			    	13);
    			    	
  fSolidHCrystal14 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal14 = new G4LogicalVolume(fSolidHCrystal14, fCrystalMaterial, "Crystal");
  fPhysiCrystal14 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal14, 
       				"Crystal", 
       				fLogicDetector14, 
       				false, 
    			    	14); 	
    			    	
  fSolidHCrystal15 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal15 = new G4LogicalVolume(fSolidHCrystal15, fCrystalMaterial, "Crystal");
  fPhysiCrystal15 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal15, 
       				"Crystal", 
       				fLogicDetector15, 
       				false, 
    			    	15);	
    			    	
  fSolidHCrystal16 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal16 = new G4LogicalVolume(fSolidHCrystal16, fCrystalMaterial, "Crystal");
  fPhysiCrystal16 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal16, 
       				"Crystal", 
       				fLogicDetector16, 
       				false, 
    			    	16);	    
    			    	
  fSolidHCrystal17 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal17 = new G4LogicalVolume(fSolidHCrystal17, fCrystalMaterial, "Crystal");
  fPhysiCrystal17 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal17, 
       				"Crystal", 
       				fLogicDetector17, 
       				false, 
    			    	17);	
    			    	
  fSolidHCrystal18 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal18 = new G4LogicalVolume(fSolidHCrystal18, fCrystalMaterial, "Crystal");
  fPhysiCrystal18 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal18, 
       				"Crystal", 
       				fLogicDetector18, 
       				false, 
    			    	18);
    			    	
  fSolidHCrystal19 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal19 = new G4LogicalVolume(fSolidHCrystal19, fCrystalMaterial, "Crystal");
  fPhysiCrystal19 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal19, 
       				"Crystal", 
       				fLogicDetector19, 
       				false, 
    			    	19);
    			    	
  fSolidHCrystal20 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal20 = new G4LogicalVolume(fSolidHCrystal20, fCrystalMaterial, "Crystal");
  fPhysiCrystal20 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal20, 
       				"Crystal", 
       				fLogicDetector20, 
       				false, 
    			    	20);
    			    	
  fSolidHCrystal21 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal21 = new G4LogicalVolume(fSolidHCrystal21, fCrystalMaterial, "Crystal");
  fPhysiCrystal21 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal21, 
       				"Crystal", 
       				fLogicDetector21, 
       				false, 
    			    	21);
    			    	
  fSolidHCrystal22 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal22 = new G4LogicalVolume(fSolidHCrystal22, fCrystalMaterial, "Crystal");
  fPhysiCrystal22 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal22, 
       				"Crystal", 
       				fLogicDetector22, 
       				false, 
    			    	22);
    			    	
  fSolidHCrystal23 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal23 = new G4LogicalVolume(fSolidHCrystal23, fCrystalMaterial, "Crystal");
  fPhysiCrystal23 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal23, 
       				"Crystal", 
       				fLogicDetector23, 
       				false, 
    			    	23);
    			    	
  fSolidHCrystal24 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal24 = new G4LogicalVolume(fSolidHCrystal24, fCrystalMaterial, "Crystal");
  fPhysiCrystal24 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal24, 
       				"Crystal", 
       				fLogicDetector24, 
       				false, 
    			    	24);
    			    	
  fSolidHCrystal25 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal25 = new G4LogicalVolume(fSolidHCrystal25, fCrystalMaterial, "Crystal");
  fPhysiCrystal25 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal25, 
       				"Crystal", 
       				fLogicDetector25, 
       				false, 
    			    	25);
    			    	
  fSolidHCrystal26 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal26 = new G4LogicalVolume(fSolidHCrystal26, fCrystalMaterial, "Crystal");
  fPhysiCrystal26 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal26, 
       				"Crystal", 
       				fLogicDetector26, 
       				false, 
    			    	26);
    			    	
  fSolidHCrystal27 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal27 = new G4LogicalVolume(fSolidHCrystal27, fCrystalMaterial, "Crystal");
  fPhysiCrystal27 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal27, 
       				"Crystal", 
       				fLogicDetector27, 
       				false, 
    			    	27);
    			    	
  fSolidHCrystal28 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal28 = new G4LogicalVolume(fSolidHCrystal28, fCrystalMaterial, "Crystal");
  fPhysiCrystal28 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal28, 
       				"Crystal", 
       				fLogicDetector28, 
       				false, 
    			    	28);
    			    	
  fSolidHCrystal29 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal29 = new G4LogicalVolume(fSolidHCrystal29, fCrystalMaterial, "Crystal");
  fPhysiCrystal29 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal29, 
       				"Crystal", 
       				fLogicDetector29, 
       				false, 
    			    	29);
    			    	
  fSolidHCrystal30 = new G4Polyhedra("Crystal", 0.*deg, 360.*deg, 6, 2, zPlane2, rInner, rInner1);
  fLogicCrystal30 = new G4LogicalVolume(fSolidHCrystal30, fCrystalMaterial, "Crystal");
  fPhysiCrystal30 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos), 
        			fLogicCrystal30, 
       				"Crystal", 
       				fLogicDetector30, 
       				false, 
    			    	30);
    			    	
  fLogicCrystal->SetVisAttributes(fRedVisAtt);
  fLogicCrystal2->SetVisAttributes(fRedVisAtt);
  fLogicCrystal3->SetVisAttributes(fRedVisAtt);
  fLogicCrystal4->SetVisAttributes(fRedVisAtt);
  fLogicCrystal5->SetVisAttributes(fRedVisAtt);
  fLogicCrystal6->SetVisAttributes(fRedVisAtt);
  fLogicCrystal7->SetVisAttributes(fRedVisAtt);
  fLogicCrystal8->SetVisAttributes(fRedVisAtt);
  fLogicCrystal9->SetVisAttributes(fRedVisAtt);
  fLogicCrystal10->SetVisAttributes(fRedVisAtt);
  fLogicCrystal11->SetVisAttributes(fRedVisAtt);
  fLogicCrystal12->SetVisAttributes(fRedVisAtt);
  fLogicCrystal13->SetVisAttributes(fRedVisAtt);
  fLogicCrystal14->SetVisAttributes(fRedVisAtt);
  fLogicCrystal15->SetVisAttributes(fRedVisAtt);
  fLogicCrystal16->SetVisAttributes(fRedVisAtt);
  fLogicCrystal17->SetVisAttributes(fRedVisAtt);
  fLogicCrystal18->SetVisAttributes(fRedVisAtt);
  fLogicCrystal19->SetVisAttributes(fRedVisAtt);
  fLogicCrystal20->SetVisAttributes(fRedVisAtt);
  fLogicCrystal21->SetVisAttributes(fRedVisAtt);
  fLogicCrystal22->SetVisAttributes(fRedVisAtt);
  fLogicCrystal23->SetVisAttributes(fRedVisAtt);
  fLogicCrystal24->SetVisAttributes(fRedVisAtt);
  fLogicCrystal25->SetVisAttributes(fRedVisAtt);
  fLogicCrystal26->SetVisAttributes(fRedVisAtt);
  fLogicCrystal27->SetVisAttributes(fRedVisAtt);
  fLogicCrystal28->SetVisAttributes(fRedVisAtt);
  fLogicCrystal29->SetVisAttributes(fRedVisAtt);
  fLogicCrystal30->SetVisAttributes(fRedVisAtt);
  fLogicCrystal1->SetVisAttributes(fMagnetaVisAtt);
  
  // Gap (gap surrounding crystal and casing or a reflector as required)
  
  fSolidHGap = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap = new G4LogicalVolume(fSolidHGap, fGapMaterial, "Gap");
  fPhysiGap = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap, 
       				"Gap", 
       				fLogicDetector, 
       				false, 
    			    	0);
    			    	
  fSolidGap1 = new G4Tubs("Gap", 0.5*fDetectorDiameter, 0.5*fGapDiameter, 0.5*fGapLength, 0.*deg, 360.*deg); 			    	
  fLogicGap1 = new G4LogicalVolume(fSolidGap1, fGapMaterial, "Gap");
  fPhysiGap1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap1), 
        			fLogicGap1, 
       				"Gap", 
       				fLogicDetector1, 
       				false, 
    			    	0);
   			    	
  fSolidHGap2 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap2 = new G4LogicalVolume(fSolidHGap2, fGapMaterial, "Gap");
  fPhysiGap2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap2, 
       				"Gap", 
       				fLogicDetector2, 
       				false, 
    			    	0);
    			    	
  fSolidHGap3 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap3 = new G4LogicalVolume(fSolidHGap3, fGapMaterial, "Gap");
  fPhysiGap3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap3, 
       				"Gap", 
       				fLogicDetector3, 
       				false, 
    			    	0);
    			    	
  fSolidHGap4 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap4 = new G4LogicalVolume(fSolidHGap4, fGapMaterial, "Gap");
  fPhysiGap4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap4, 
       				"Gap", 
       				fLogicDetector4, 
       				false, 
    			    	0);
    			    	
  fSolidHGap5 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap5 = new G4LogicalVolume(fSolidHGap5, fGapMaterial, "Gap");
  fPhysiGap5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap5, 
       				"Gap", 
       				fLogicDetector5, 
       				false, 
    			    	0);
    			    	
  fSolidHGap6 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap6 = new G4LogicalVolume(fSolidHGap6, fGapMaterial, "Gap");
  fPhysiGap6 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap6,
       				"Gap", 
       				fLogicDetector6, 
       				false, 
    			    	0);
    			    	
  fSolidHGap7 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap7 = new G4LogicalVolume(fSolidHGap7, fGapMaterial, "Gap");
  fPhysiGap7 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap7, 
       				"Gap", 
       				fLogicDetector7, 
       				false, 
    			    	0);
    			    	
  fSolidHGap8 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap8 = new G4LogicalVolume(fSolidHGap8, fGapMaterial, "Gap");
  fPhysiGap8 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap8, 
       				"Gap", 
       				fLogicDetector8, 
       				false, 
    			    	0);
    			    	
  fSolidHGap9 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap9 = new G4LogicalVolume(fSolidHGap9, fGapMaterial, "Gap");
  fPhysiGap9 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap9, 
       				"Gap", 
       				fLogicDetector9, 
       				false, 
    			    	0);
    			    	
  fSolidHGap10 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap10 = new G4LogicalVolume(fSolidHGap10, fGapMaterial, "Gap");
  fPhysiGap10 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap10, 
       				"Gap", 
       				fLogicDetector10, 
       				false, 
    			    	0);
    			    	
  fSolidHGap11 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap11 = new G4LogicalVolume(fSolidHGap11, fGapMaterial, "Gap");
  fPhysiGap11 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap11, 
       				"Gap", 
       				fLogicDetector11, 
       				false, 
    			    	0);
    			    	
  fSolidHGap12 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap12 = new G4LogicalVolume(fSolidHGap12, fGapMaterial, "Gap");
  fPhysiGap12 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap12, 
       				"Gap", 
       				fLogicDetector12, 
       				false, 
    			    	0);
    			    	
  fSolidHGap13 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap13 = new G4LogicalVolume(fSolidHGap13, fGapMaterial, "Gap");
  fPhysiGap13 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap13, 
       				"Gap", 
       				fLogicDetector13, 
       				false, 
    			    	0);
    			    	
  fSolidHGap14 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap14 = new G4LogicalVolume(fSolidHGap14, fGapMaterial, "Gap");
  fPhysiGap14 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap14, 
       				"Gap", 
       				fLogicDetector14, 
       				false, 
    			    	0);
    			    	
  fSolidHGap15 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap15 = new G4LogicalVolume(fSolidHGap15, fGapMaterial, "Gap");
  fPhysiGap15 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap15, 
       				"Gap", 
       				fLogicDetector15, 
       				false, 
    			    	0);
    			    	
  fSolidHGap16 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap16 = new G4LogicalVolume(fSolidHGap16, fGapMaterial, "Gap");
  fPhysiGap16 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap16, 
       				"Gap", 
       				fLogicDetector16, 
       				false, 
    			    	0);
    			    	
  fSolidHGap17 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap17 = new G4LogicalVolume(fSolidHGap17, fGapMaterial, "Gap");
  fPhysiGap17 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap17, 
       				"Gap", 
       				fLogicDetector17, 
       				false, 
    			    	0);
    			    	
  fSolidHGap18 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap18 = new G4LogicalVolume(fSolidHGap18, fGapMaterial, "Gap");
  fPhysiGap18 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap18, 
       				"Gap", 
       				fLogicDetector18, 
       				false, 
    			    	0);
    			    	
  fSolidHGap19 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap19 = new G4LogicalVolume(fSolidHGap19, fGapMaterial, "Gap");
  fPhysiGap19 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap19, 
       				"Gap", 
       				fLogicDetector19, 
       				false, 
    			    	0);
    			    	
  fSolidHGap20 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap20 = new G4LogicalVolume(fSolidHGap20, fGapMaterial, "Gap");
  fPhysiGap20 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap20, 
       				"Gap", 
       				fLogicDetector20, 
       				false, 
    			    	0);
    			    	
  fSolidHGap21 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap21 = new G4LogicalVolume(fSolidHGap21, fGapMaterial, "Gap");
  fPhysiGap21 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap21, 
       				"Gap", 
       				fLogicDetector21, 
       				false, 
    			    	0);
    			    	
  fSolidHGap22 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap22 = new G4LogicalVolume(fSolidHGap22, fGapMaterial, "Gap");
  fPhysiGap22 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap22, 
       				"Gap", 
       				fLogicDetector22, 
       				false, 
    			    	0);
    			    	
  fSolidHGap23 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap23 = new G4LogicalVolume(fSolidHGap23, fGapMaterial, "Gap");
  fPhysiGap23 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap23, 
       				"Gap", 
       				fLogicDetector23, 
       				false, 
    			    	0);
    			    	
  fSolidHGap24 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap24 = new G4LogicalVolume(fSolidHGap24, fGapMaterial, "Gap");
  fPhysiGap24 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap24, 
       				"Gap", 
       				fLogicDetector24, 
       				false, 
    			    	0);
    			    	
  fSolidHGap25 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap25 = new G4LogicalVolume(fSolidHGap25, fGapMaterial, "Gap");
  fPhysiGap25 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap25, 
       				"Gap", 
       				fLogicDetector25, 
       				false, 
    			    	0);
    			    	
  fSolidHGap26 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap26 = new G4LogicalVolume(fSolidHGap26, fGapMaterial, "Gap");
  fPhysiGap26 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap26, 
       				"Gap", 
       				fLogicDetector26, 
       				false, 
    			    	0);
    			    	
  fSolidHGap27 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap27 = new G4LogicalVolume(fSolidHGap27, fGapMaterial, "Gap");
  fPhysiGap27 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap27, 
       				"Gap", 
       				fLogicDetector27, 
       				false, 
    			    	0);
    			    	
  fSolidHGap28 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap28 = new G4LogicalVolume(fSolidHGap28, fGapMaterial, "Gap");
  fPhysiGap28 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap28, 
       				"Gap", 
       				fLogicDetector28, 
       				false, 
    			    	0);
    			    	
  fSolidHGap29 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap29 = new G4LogicalVolume(fSolidHGap29, fGapMaterial, "Gap");
  fPhysiGap29 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap29, 
       				"Gap", 
       				fLogicDetector29, 
       				false, 
    			    	0);
    			    	
  fSolidHGap30 = new G4Polyhedra("Gap", 0.*deg, 360.*deg, 6, 2, zPlane3, rInner1, rOuter1);
  fLogicGap30 = new G4LogicalVolume(fSolidHGap30, fGapMaterial, "Gap");
  fPhysiGap30 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposGap), 
        			fLogicGap30, 
       				"Gap", 
       				fLogicDetector30, 
       				false, 
    			    	0);

  // FaceGap (gap between face of crystal and casing or a reflector as required)
  //
  fSolidHFaceGap = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
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
    			    	
  fSolidHFaceGap2 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap2 = new G4LogicalVolume(fSolidHFaceGap2, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap2, 
       				"FaceGap", 
       				fLogicDetector2, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap3 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap3 = new G4LogicalVolume(fSolidHFaceGap3, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap3, 
       				"FaceGap", 
       				fLogicDetector3, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap4 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap4 = new G4LogicalVolume(fSolidHFaceGap4, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap4, 
       				"FaceGap", 
       				fLogicDetector4, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap5 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap5 = new G4LogicalVolume(fSolidHFaceGap5, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap5, 
       				"FaceGap", 
       				fLogicDetector5, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap6 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap6 = new G4LogicalVolume(fSolidHFaceGap6, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap6 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap6, 
       				"FaceGap", 
       				fLogicDetector6, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap7 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap7 = new G4LogicalVolume(fSolidHFaceGap7, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap7 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap7, 
       				"FaceGap", 
       				fLogicDetector7, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap8 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap8 = new G4LogicalVolume(fSolidHFaceGap8, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap8 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap8, 
       				"FaceGap", 
       				fLogicDetector8, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap9 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap9 = new G4LogicalVolume(fSolidHFaceGap9, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap9 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap9, 
       				"FaceGap", 
       				fLogicDetector9, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap10 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap10 = new G4LogicalVolume(fSolidHFaceGap10, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap10 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap10, 
       				"FaceGap", 
       				fLogicDetector10, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap11 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap11 = new G4LogicalVolume(fSolidHFaceGap11, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap11 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap11, 
       				"FaceGap", 
       				fLogicDetector11, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap12 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap12 = new G4LogicalVolume(fSolidHFaceGap12, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap12 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap12, 
       				"FaceGap", 
       				fLogicDetector12, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap13 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap13 = new G4LogicalVolume(fSolidHFaceGap13, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap13 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap13, 
       				"FaceGap", 
       				fLogicDetector13, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap14 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap14 = new G4LogicalVolume(fSolidHFaceGap14, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap14 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap14, 
       				"FaceGap", 
       				fLogicDetector14, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap15 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap15 = new G4LogicalVolume(fSolidHFaceGap15, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap15 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap15, 
       				"FaceGap", 
       				fLogicDetector15, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap16 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap16 = new G4LogicalVolume(fSolidHFaceGap16, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap16 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap16, 
       				"FaceGap", 
       				fLogicDetector16, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap17 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap17 = new G4LogicalVolume(fSolidHFaceGap17, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap17 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap17, 
       				"FaceGap", 
       				fLogicDetector17, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap18 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap18 = new G4LogicalVolume(fSolidHFaceGap18, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap18 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap18, 
       				"FaceGap", 
       				fLogicDetector18, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap19 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap19 = new G4LogicalVolume(fSolidHFaceGap19, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap19 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap19, 
       				"FaceGap", 
       				fLogicDetector19, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap20 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap20 = new G4LogicalVolume(fSolidHFaceGap20, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap20 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap20, 
       				"FaceGap", 
       				fLogicDetector20, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap21 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap21 = new G4LogicalVolume(fSolidHFaceGap21, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap21 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap21, 
       				"FaceGap", 
       				fLogicDetector21, 
       				false, 
    			    	0);
    			    
  fSolidHFaceGap22 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap22 = new G4LogicalVolume(fSolidHFaceGap22, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap22 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap22, 
       				"FaceGap", 
       				fLogicDetector22, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap23 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap23 = new G4LogicalVolume(fSolidHFaceGap23, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap23 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap23, 
       				"FaceGap", 
       				fLogicDetector23, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap24 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap24 = new G4LogicalVolume(fSolidHFaceGap24, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap24 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap24, 
       				"FaceGap", 
       				fLogicDetector24, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap25 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap25 = new G4LogicalVolume(fSolidHFaceGap25, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap25 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap25, 
       				"FaceGap", 
       				fLogicDetector25, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap26 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap26 = new G4LogicalVolume(fSolidHFaceGap26, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap26 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap26, 
       				"FaceGap", 
       				fLogicDetector26, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap27 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap27 = new G4LogicalVolume(fSolidHFaceGap27, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap27 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap27, 
       				"FaceGap", 
       				fLogicDetector27, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap28 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap28 = new G4LogicalVolume(fSolidHFaceGap28, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap28 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap28, 
       				"FaceGap", 
       				fLogicDetector28, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap29 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap29 = new G4LogicalVolume(fSolidHFaceGap29, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap29 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap29, 
       				"FaceGap", 
       				fLogicDetector29, 
       				false, 
    			    	0);
    			    	
  fSolidHFaceGap30 = new G4Polyhedra("FaceGap", 0.*deg, 360.*deg, 6, 2, zPlane10, rInner, rInner1);
  fLogicFaceGap30 = new G4LogicalVolume(fSolidHFaceGap30, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap30 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap), 
        			fLogicFaceGap30, 
       				"FaceGap", 
       				fLogicDetector30, 
       				false, 
    			    	0);
  
  fLogicGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap2->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap3->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap4->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap5->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap6->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap7->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap8->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap9->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap10->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap11->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap12->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap13->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap14->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap15->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap16->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap17->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap18->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap19->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap20->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap21->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap22->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap23->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap24->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap25->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap26->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap27->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap28->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap29->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap30->SetVisAttributes(fAuxEdgeVisAtt);
  
  fLogicGap->SetVisAttributes(fBlueVisAtt);
  fLogicGap2->SetVisAttributes(fBlueVisAtt);
  fLogicGap3->SetVisAttributes(fBlueVisAtt);
  fLogicGap4->SetVisAttributes(fBlueVisAtt);
  fLogicGap5->SetVisAttributes(fBlueVisAtt);
  fLogicGap6->SetVisAttributes(fBlueVisAtt);
  fLogicGap7->SetVisAttributes(fBlueVisAtt);
  fLogicGap8->SetVisAttributes(fBlueVisAtt);
  fLogicGap9->SetVisAttributes(fBlueVisAtt);
  fLogicGap10->SetVisAttributes(fBlueVisAtt);
  fLogicGap11->SetVisAttributes(fBlueVisAtt);
  fLogicGap12->SetVisAttributes(fBlueVisAtt);
  fLogicGap13->SetVisAttributes(fBlueVisAtt);
  fLogicGap14->SetVisAttributes(fBlueVisAtt);
  fLogicGap15->SetVisAttributes(fBlueVisAtt);
  fLogicGap16->SetVisAttributes(fBlueVisAtt);
  fLogicGap17->SetVisAttributes(fBlueVisAtt);
  fLogicGap18->SetVisAttributes(fBlueVisAtt);
  fLogicGap19->SetVisAttributes(fBlueVisAtt);
  fLogicGap20->SetVisAttributes(fBlueVisAtt);
  fLogicGap21->SetVisAttributes(fBlueVisAtt);
  fLogicGap22->SetVisAttributes(fBlueVisAtt);
  fLogicGap23->SetVisAttributes(fBlueVisAtt);
  fLogicGap24->SetVisAttributes(fBlueVisAtt);
  fLogicGap25->SetVisAttributes(fBlueVisAtt);
  fLogicGap26->SetVisAttributes(fBlueVisAtt);
  fLogicGap27->SetVisAttributes(fBlueVisAtt);
  fLogicGap28->SetVisAttributes(fBlueVisAtt);
  fLogicGap29->SetVisAttributes(fBlueVisAtt);
  fLogicGap30->SetVisAttributes(fBlueVisAtt);
  
  fLogicFaceGap->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap2->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap3->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap4->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap5->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap6->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap7->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap8->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap9->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap10->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap11->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap12->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap13->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap14->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap15->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap16->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap17->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap18->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap19->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap20->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap21->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap22->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap23->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap24->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap25->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap26->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap27->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap28->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap29->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap30->SetVisAttributes(fAuxEdgeVisAtt);
  
  fLogicFaceGap->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap2->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap3->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap4->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap5->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap6->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap7->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap8->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap9->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap10->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap11->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap12->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap13->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap14->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap15->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap16->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap17->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap18->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap19->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap20->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap21->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap22->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap23->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap24->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap25->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap26->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap27->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap28->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap29->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap30->SetVisAttributes(fCyanVisAtt);
  
  fLogicFaceGap1->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicGap1->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap1->SetVisAttributes(fCyanVisAtt);
  fLogicGap1->SetVisAttributes(fCyanVisAtt);
  
  // Aluminum Casing (casing surrounding crystal)
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
    			    	
  fSolidHAlCase2 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase2 = new G4LogicalVolume(fSolidHAlCase2, fAlCaseMaterial, "AlCase");
  fPhysiAlCase2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase2, 
       				"AlCase", 
       				fLogicDetector2, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase3 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase3 = new G4LogicalVolume(fSolidHAlCase3, fAlCaseMaterial, "AlCase");
  fPhysiAlCase3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase3, 
       				"AlCase", 
       				fLogicDetector3, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase4 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase4 = new G4LogicalVolume(fSolidHAlCase4, fAlCaseMaterial, "AlCase");
  fPhysiAlCase4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase4, 
       				"AlCase", 
       				fLogicDetector4, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase5 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase5 = new G4LogicalVolume(fSolidHAlCase5, fAlCaseMaterial, "AlCase");
  fPhysiAlCase5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase5, 
       				"AlCase", 
       				fLogicDetector5, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase6 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase6 = new G4LogicalVolume(fSolidHAlCase6, fAlCaseMaterial, "AlCase");
  fPhysiAlCase6 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase6, 
       				"AlCase", 
       				fLogicDetector6, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase7 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase7 = new G4LogicalVolume(fSolidHAlCase7, fAlCaseMaterial, "AlCase");
  fPhysiAlCase7 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase7, 
       				"AlCase", 
       				fLogicDetector7, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase8 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase8 = new G4LogicalVolume(fSolidHAlCase8, fAlCaseMaterial, "AlCase");
  fPhysiAlCase8 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase8, 
       				"AlCase", 
       				fLogicDetector8, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase9 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase9 = new G4LogicalVolume(fSolidHAlCase9, fAlCaseMaterial, "AlCase");
  fPhysiAlCase9 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase9, 
       				"AlCase", 
       				fLogicDetector9, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase10 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase10 = new G4LogicalVolume(fSolidHAlCase10, fAlCaseMaterial, "AlCase");
  fPhysiAlCase10 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase10, 
       				"AlCase", 
       				fLogicDetector10, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase11 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase11 = new G4LogicalVolume(fSolidHAlCase11, fAlCaseMaterial, "AlCase");
  fPhysiAlCase11 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase11, 
       				"AlCase", 
       				fLogicDetector11, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase12 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase12 = new G4LogicalVolume(fSolidHAlCase12, fAlCaseMaterial, "AlCase");
  fPhysiAlCase12 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase12, 
       				"AlCase", 
       				fLogicDetector12, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase13 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase13 = new G4LogicalVolume(fSolidHAlCase13, fAlCaseMaterial, "AlCase");
  fPhysiAlCase13 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase13, 
       				"AlCase", 
       				fLogicDetector13, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase14 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase14 = new G4LogicalVolume(fSolidHAlCase14, fAlCaseMaterial, "AlCase");
  fPhysiAlCase14 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase14, 
       				"AlCase", 
       				fLogicDetector14, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase15 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase15 = new G4LogicalVolume(fSolidHAlCase15, fAlCaseMaterial, "AlCase");
  fPhysiAlCase15 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase15, 
       				"AlCase", 
       				fLogicDetector15, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase16 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase16 = new G4LogicalVolume(fSolidHAlCase16, fAlCaseMaterial, "AlCase");
  fPhysiAlCase16 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase16, 
       				"AlCase", 
       				fLogicDetector16, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase17 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase17 = new G4LogicalVolume(fSolidHAlCase17, fAlCaseMaterial, "AlCase");
  fPhysiAlCase17 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase17, 
       				"AlCase", 
       				fLogicDetector17, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase18 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase18 = new G4LogicalVolume(fSolidHAlCase18, fAlCaseMaterial, "AlCase");
  fPhysiAlCase18 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase18, 
       				"AlCase", 
       				fLogicDetector18, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase19 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase19 = new G4LogicalVolume(fSolidHAlCase19, fAlCaseMaterial, "AlCase");
  fPhysiAlCase19 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase19, 
       				"AlCase", 
       				fLogicDetector19, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase20 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase20 = new G4LogicalVolume(fSolidHAlCase20, fAlCaseMaterial, "AlCase");
  fPhysiAlCase20 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase20, 
       				"AlCase", 
       				fLogicDetector20, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase21 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase21 = new G4LogicalVolume(fSolidHAlCase21, fAlCaseMaterial, "AlCase");
  fPhysiAlCase21 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase21, 
       				"AlCase", 
       				fLogicDetector21, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase22 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase22 = new G4LogicalVolume(fSolidHAlCase22, fAlCaseMaterial, "AlCase");
  fPhysiAlCase22 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase22, 
       				"AlCase", 
       				fLogicDetector22, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase23 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase23 = new G4LogicalVolume(fSolidHAlCase23, fAlCaseMaterial, "AlCase");
  fPhysiAlCase23 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase23, 
       				"AlCase", 
       				fLogicDetector23, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase24 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase24 = new G4LogicalVolume(fSolidHAlCase24, fAlCaseMaterial, "AlCase");
  fPhysiAlCase24 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase24, 
       				"AlCase", 
       				fLogicDetector24, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase25 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase25 = new G4LogicalVolume(fSolidHAlCase25, fAlCaseMaterial, "AlCase");
  fPhysiAlCase25 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase25, 
       				"AlCase", 
       				fLogicDetector25, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase26 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase26 = new G4LogicalVolume(fSolidHAlCase26, fAlCaseMaterial, "AlCase");
  fPhysiAlCase26 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase26, 
       				"AlCase", 
       				fLogicDetector26, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase27 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase27 = new G4LogicalVolume(fSolidHAlCase27, fAlCaseMaterial, "AlCase");
  fPhysiAlCase27 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase27, 
       				"AlCase", 
       				fLogicDetector27, 
       				false, 
    			    	0);
    			    
  fSolidHAlCase28 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase28 = new G4LogicalVolume(fSolidHAlCase28, fAlCaseMaterial, "AlCase");
  fPhysiAlCase28 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase28, 
       				"AlCase", 
       				fLogicDetector28, 
       				false, 
    			    	0);
    			    
  fSolidHAlCase29 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase29 = new G4LogicalVolume(fSolidHAlCase29, fAlCaseMaterial, "AlCase");
  fPhysiAlCase29 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase29, 
       				"AlCase", 
       				fLogicDetector29, 
       				false, 
    			    	0);
    			    	
  fSolidHAlCase30 = new G4Polyhedra("AlCase", 0.*deg, 360.*deg, 6, 2, zPlane5, rInner2, rOuter2);
  fLogicAlCase30 = new G4LogicalVolume(fSolidHAlCase30, fAlCaseMaterial, "AlCase");
  fPhysiAlCase30 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase), 
        			fLogicAlCase30, 
       				"AlCase", 
       				fLogicDetector30, 
       				false, 
    			    	0);

  // Aluminum Face Casing (casing covering face of crystal)
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
    			    	
  fSolidHFaceAlCase2 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase2 = new G4LogicalVolume(fSolidHFaceAlCase2, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase2, 
       				"FaceAlCase", 
       				fLogicDetector2, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase3 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase3 = new G4LogicalVolume(fSolidHFaceAlCase3, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase3, 
       				"FaceAlCase", 
       				fLogicDetector3, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase4 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase4 = new G4LogicalVolume(fSolidHFaceAlCase4, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase4, 
       				"FaceAlCase", 
       				fLogicDetector4, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase5 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase5 = new G4LogicalVolume(fSolidHFaceAlCase5, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase5, 
       				"FaceAlCase", 
       				fLogicDetector5, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase6 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase6 = new G4LogicalVolume(fSolidHFaceAlCase6, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase6 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase6, 
       				"FaceAlCase", 
       				fLogicDetector6, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase7 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase7 = new G4LogicalVolume(fSolidHFaceAlCase7, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase7 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase7, 
       				"FaceAlCase", 
       				fLogicDetector7, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase8 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase8 = new G4LogicalVolume(fSolidHFaceAlCase8, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase8 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase8, 
       				"FaceAlCase", 
       				fLogicDetector8, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase9 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase9 = new G4LogicalVolume(fSolidHFaceAlCase9, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase9 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase9, 
       				"FaceAlCase", 
       				fLogicDetector9, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase10 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase10 = new G4LogicalVolume(fSolidHFaceAlCase10, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase10 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase10, 
       				"FaceAlCase", 
       				fLogicDetector10, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase11 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase11 = new G4LogicalVolume(fSolidHFaceAlCase11, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase11 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase11, 
       				"FaceAlCase", 
       				fLogicDetector11, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase12 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase12 = new G4LogicalVolume(fSolidHFaceAlCase12, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase12 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase12, 
       				"FaceAlCase", 
       				fLogicDetector12, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase13 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase13 = new G4LogicalVolume(fSolidHFaceAlCase13, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase13 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase13, 
       				"FaceAlCase", 
       				fLogicDetector13, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase14 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase14 = new G4LogicalVolume(fSolidHFaceAlCase14, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase14 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase14, 
       				"FaceAlCase", 
       				fLogicDetector14, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase15 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase15 = new G4LogicalVolume(fSolidHFaceAlCase15, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase15 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase15, 
       				"FaceAlCase", 
       				fLogicDetector15, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase16 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase16 = new G4LogicalVolume(fSolidHFaceAlCase16, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase16 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase16, 
       				"FaceAlCase", 
       				fLogicDetector16, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase17 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase17 = new G4LogicalVolume(fSolidHFaceAlCase17, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase17 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase17, 
       				"FaceAlCase", 
       				fLogicDetector17, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase18 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase18 = new G4LogicalVolume(fSolidHFaceAlCase18, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase18 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase18, 
       				"FaceAlCase", 
       				fLogicDetector18, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase19 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase19 = new G4LogicalVolume(fSolidHFaceAlCase19, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase19 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase19, 
       				"FaceAlCase", 
       				fLogicDetector19, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase20 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase20 = new G4LogicalVolume(fSolidHFaceAlCase20, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase20 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase20, 
       				"FaceAlCase", 
       				fLogicDetector20, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase21 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase21 = new G4LogicalVolume(fSolidHFaceAlCase21, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase21 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase21, 
       				"FaceAlCase", 
       				fLogicDetector21, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase22 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase22 = new G4LogicalVolume(fSolidHFaceAlCase22, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase22 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase22, 
       				"FaceAlCase", 
       				fLogicDetector22, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase23 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase23 = new G4LogicalVolume(fSolidHFaceAlCase23, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase23 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase23, 
       				"FaceAlCase", 
       				fLogicDetector23, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase24 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase24 = new G4LogicalVolume(fSolidHFaceAlCase24, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase24 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase24, 
       				"FaceAlCase", 
       				fLogicDetector24, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase25 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase25 = new G4LogicalVolume(fSolidHFaceAlCase25, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase25 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase25, 
       				"FaceAlCase", 
       				fLogicDetector25, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase26 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase26 = new G4LogicalVolume(fSolidHFaceAlCase26, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase26 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase26, 
       				"FaceAlCase", 
       				fLogicDetector26, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase27 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase27 = new G4LogicalVolume(fSolidHFaceAlCase27, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase27 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase27, 
       				"FaceAlCase", 
       				fLogicDetector27, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase28 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase28 = new G4LogicalVolume(fSolidHFaceAlCase28, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase28 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase28, 
       				"FaceAlCase", 
       				fLogicDetector28, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase29 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase29 = new G4LogicalVolume(fSolidHFaceAlCase29, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase29 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase29, 
       				"FaceAlCase", 
       				fLogicDetector29, 
       				false,
    			    	0);
    			    	
  fSolidHFaceAlCase30 = new G4Polyhedra("FaceAlCase", 0.*deg, 360.*deg, 6, 2, zPlane6, rInner, rOuter1);
  fLogicFaceAlCase30 = new G4LogicalVolume(fSolidHFaceAlCase30, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase30 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase), 
        			fLogicFaceAlCase30, 
       				"FaceAlCase", 
       				fLogicDetector30, 
       				false,
    			    	0);
    			    	
  fLogicAlCase->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase2->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase3->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase4->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase5->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase6->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase7->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase8->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase9->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase10->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase11->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase12->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase13->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase14->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase15->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase16->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase17->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase18->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase19->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase20->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase21->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase22->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase23->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase24->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase25->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase26->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase27->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase28->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase29->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase30->SetVisAttributes(fGreyVisAtt);
   
  fLogicFaceAlCase->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase2->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase3->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase4->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase5->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase6->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase7->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase8->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase9->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase10->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase11->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase12->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase13->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase14->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase15->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase16->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase17->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase18->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase19->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase20->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase21->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase22->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase23->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase24->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase25->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase26->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase27->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase28->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase29->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase30->SetVisAttributes(fGreyVisAtt);
  
  // PMT
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
    			    	
  fSolidPMT2 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT2 = new G4LogicalVolume(fSolidPMT2, fPMTMaterial,"PMT");
  fPhysiPMT2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT2, 
       				"PMT", 
       				fLogicDetector2, 
       				false, 
    			    	0);
    			    	
  fSolidPMT3 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT3 = new G4LogicalVolume(fSolidPMT3, fPMTMaterial,"PMT");
  fPhysiPMT3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT3, 
       				"PMT", 
       				fLogicDetector3, 
       				false, 
    			    	0);
    			    	
  fSolidPMT4 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT4 = new G4LogicalVolume(fSolidPMT4, fPMTMaterial,"PMT");
  fPhysiPMT4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT4, 
       				"PMT", 
       				fLogicDetector4, 
       				false, 
    			    	0);
    			    	
  fSolidPMT5 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT5 = new G4LogicalVolume(fSolidPMT5, fPMTMaterial,"PMT");
  fPhysiPMT5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT5, 
       				"PMT", 
       				fLogicDetector5, 
       				false, 
    			    	0);
    			    	
  fSolidPMT6 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT6 = new G4LogicalVolume(fSolidPMT6, fPMTMaterial,"PMT");
  fPhysiPMT6 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT6, 
       				"PMT", 
       				fLogicDetector6, 
       				false, 
    			    	0);
    			    	
  fSolidPMT7 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT7 = new G4LogicalVolume(fSolidPMT7, fPMTMaterial,"PMT");
  fPhysiPMT7 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT7, 
       				"PMT", 
       				fLogicDetector7, 
       				false, 
    			    	0);
    			    	
  fSolidPMT8 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT8 = new G4LogicalVolume(fSolidPMT8, fPMTMaterial,"PMT");
  fPhysiPMT8 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT8, 
       				"PMT", 
       				fLogicDetector8, 
       				false, 
    			    	0);
    			    	
  fSolidPMT9 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT9 = new G4LogicalVolume(fSolidPMT9, fPMTMaterial,"PMT");
  fPhysiPMT9 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT9, 
       				"PMT", 
       				fLogicDetector9, 
       				false, 
    			    	0);
    			    	
  fSolidPMT10 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT10 = new G4LogicalVolume(fSolidPMT10, fPMTMaterial,"PMT");
  fPhysiPMT10 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT10, 
       				"PMT", 
       				fLogicDetector10, 
       				false, 
    			    	0);
    			    	
  fSolidPMT11 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT11 = new G4LogicalVolume(fSolidPMT11, fPMTMaterial,"PMT");
  fPhysiPMT11 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT11, 
       				"PMT", 
       				fLogicDetector11, 
       				false, 
    			    	0);
    			    	
  fSolidPMT12 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT12 = new G4LogicalVolume(fSolidPMT12, fPMTMaterial,"PMT");
  fPhysiPMT12 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT12, 
       				"PMT", 
       				fLogicDetector12, 
       				false, 
    			    	0);
    			    	
  fSolidPMT13 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT13 = new G4LogicalVolume(fSolidPMT13, fPMTMaterial,"PMT");
  fPhysiPMT13 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT13, 
       				"PMT", 
       				fLogicDetector13, 
       				false, 
    			    	0);
    			    	
  fSolidPMT14 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT14 = new G4LogicalVolume(fSolidPMT14, fPMTMaterial,"PMT");
  fPhysiPMT14 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT14, 
       				"PMT", 
       				fLogicDetector14, 
       				false, 
    			    	0);
    			    	
  fSolidPMT15 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT15 = new G4LogicalVolume(fSolidPMT15, fPMTMaterial,"PMT");
  fPhysiPMT15 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT15, 
       				"PMT", 
       				fLogicDetector15, 
       				false, 
    			    	0);
    			    	
  fSolidPMT16 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT16 = new G4LogicalVolume(fSolidPMT16, fPMTMaterial,"PMT");
  fPhysiPMT16 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT16, 
       				"PMT", 
       				fLogicDetector16, 
       				false, 
    			    	0);
    			    	
  fSolidPMT17 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT17 = new G4LogicalVolume(fSolidPMT17, fPMTMaterial,"PMT");
  fPhysiPMT17 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT17, 
       				"PMT", 
       				fLogicDetector17, 
       				false, 
    			    	0);
    			    	
  fSolidPMT18 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT18 = new G4LogicalVolume(fSolidPMT18, fPMTMaterial,"PMT");
  fPhysiPMT18 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT18, 
       				"PMT", 
       				fLogicDetector18, 
       				false, 
    			    	0);
  
  fSolidPMT19 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT19 = new G4LogicalVolume(fSolidPMT19, fPMTMaterial,"PMT");
  fPhysiPMT19 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT19, 
       				"PMT", 
       				fLogicDetector19, 
       				false, 
    			    	0);
  
  fSolidPMT20 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT20 = new G4LogicalVolume(fSolidPMT20, fPMTMaterial,"PMT");
  fPhysiPMT20 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT20, 
       				"PMT", 
       				fLogicDetector20, 
       				false, 
    			    	0);
    			    	
  fSolidPMT21 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT21 = new G4LogicalVolume(fSolidPMT21, fPMTMaterial,"PMT");
  fPhysiPMT21 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT21, 
       				"PMT", 
       				fLogicDetector21, 
       				false, 
    			    	0);
    			    	
  fSolidPMT22 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT22 = new G4LogicalVolume(fSolidPMT22, fPMTMaterial,"PMT");
  fPhysiPMT22 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT22, 
       				"PMT", 
       				fLogicDetector22, 
       				false, 
    			    	0);
    			    	
  fSolidPMT23 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT23 = new G4LogicalVolume(fSolidPMT23, fPMTMaterial,"PMT");
  fPhysiPMT23 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT23, 
       				"PMT", 
       				fLogicDetector23, 
       				false, 
    			    	0);
    			    	
  fSolidPMT24 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT24 = new G4LogicalVolume(fSolidPMT24, fPMTMaterial,"PMT");
  fPhysiPMT24 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT24, 
       				"PMT", 
       				fLogicDetector24, 
       				false, 
    			    	0);
    			    	
  fSolidPMT25 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT25 = new G4LogicalVolume(fSolidPMT25, fPMTMaterial,"PMT");
  fPhysiPMT25 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT25, 
       				"PMT", 
       				fLogicDetector25, 
       				false, 
    			    	0);
    			    	
  fSolidPMT26 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT26 = new G4LogicalVolume(fSolidPMT26, fPMTMaterial,"PMT");
  fPhysiPMT26 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT26, 
       				"PMT", 
       				fLogicDetector26, 
       				false, 
    			    	0);
    			    	
  fSolidPMT27 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT27 = new G4LogicalVolume(fSolidPMT27, fPMTMaterial,"PMT");
  fPhysiPMT27 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT27, 
       				"PMT", 
       				fLogicDetector27, 
       				false, 
    			    	0);
    			    	
  fSolidPMT28 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT28 = new G4LogicalVolume(fSolidPMT28, fPMTMaterial,"PMT");
  fPhysiPMT28 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT28, 
       				"PMT", 
       				fLogicDetector28, 
       				false, 
    			    	0);
    			    	
  fSolidPMT29 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT29 = new G4LogicalVolume(fSolidPMT29, fPMTMaterial,"PMT");
  fPhysiPMT29 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT29, 
       				"PMT", 
       				fLogicDetector29, 
       				false, 
    			    	0);
    			    	
  fSolidPMT30 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter, 0.5*fPMTLength, 0.*deg, 360.*deg); 	
  fLogicPMT30 = new G4LogicalVolume(fSolidPMT30, fPMTMaterial,"PMT");
  fPhysiPMT30 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT), 
        			fLogicPMT30, 
       				"PMT", 
       				fLogicDetector30, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ = new G4LogicalVolume(fSolidPMTJ, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ, 
       				"Vacuum", 
       				fLogicPMT, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ2 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ2 = new G4LogicalVolume(fSolidPMTJ2, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ2, 
       				"Vacuum", 
       				fLogicPMT2, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ3 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ3 = new G4LogicalVolume(fSolidPMTJ3, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ3, 
       				"Vacuum", 
       				fLogicPMT3, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ4 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ4 = new G4LogicalVolume(fSolidPMTJ4, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ4, 
       				"Vacuum", 
       				fLogicPMT4, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ5 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ5 = new G4LogicalVolume(fSolidPMTJ5, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ5, 
       				"Vacuum", 
       				fLogicPMT5, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ6 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ6 = new G4LogicalVolume(fSolidPMTJ6, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ6 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ6, 
       				"Vacuum", 
       				fLogicPMT6, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ7 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ7 = new G4LogicalVolume(fSolidPMTJ7, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ7 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ7, 
       				"Vacuum", 
       				fLogicPMT7, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ8 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ8 = new G4LogicalVolume(fSolidPMTJ8, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ8 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ8, 
       				"Vacuum", 
       				fLogicPMT8, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ9 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ9 = new G4LogicalVolume(fSolidPMTJ9, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ9 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ9, 
       				"Vacuum", 
       				fLogicPMT9, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ10 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ10 = new G4LogicalVolume(fSolidPMTJ10, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ10 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ10, 
       				"Vacuum", 
       				fLogicPMT10, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ11 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ11 = new G4LogicalVolume(fSolidPMTJ11, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ11 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ11, 
       				"Vacuum", 
       				fLogicPMT11, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ12 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ12 = new G4LogicalVolume(fSolidPMTJ12, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ12 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ12, 
       				"Vacuum", 
       				fLogicPMT12, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ13 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ13 = new G4LogicalVolume(fSolidPMTJ13, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ13 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ13, 
       				"Vacuum", 
       				fLogicPMT13, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ14 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ14 = new G4LogicalVolume(fSolidPMTJ14, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ14 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ14, 
       				"Vacuum", 
       				fLogicPMT14, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ15 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ15 = new G4LogicalVolume(fSolidPMTJ15, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ15 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ15, 
       				"Vacuum", 
       				fLogicPMT15, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ16 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ16 = new G4LogicalVolume(fSolidPMTJ16, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ16 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ16, 
       				"Vacuum", 
       				fLogicPMT16, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ17 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ17 = new G4LogicalVolume(fSolidPMTJ17, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ17 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ17, 
       				"Vacuum", 
       				fLogicPMT17, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ18 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ18 = new G4LogicalVolume(fSolidPMTJ18, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ18 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ18, 
       				"Vacuum", 
       				fLogicPMT18, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ19 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ19 = new G4LogicalVolume(fSolidPMTJ19, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ19 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ19, 
       				"Vacuum", 
       				fLogicPMT19, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ20 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ20 = new G4LogicalVolume(fSolidPMTJ20, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ20 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ20, 
       				"Vacuum", 
       				fLogicPMT20, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ21 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ21 = new G4LogicalVolume(fSolidPMTJ21, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ21 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ21, 
       				"Vacuum", 
       				fLogicPMT21, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ22 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ22 = new G4LogicalVolume(fSolidPMTJ22, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ22 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ22, 
       				"Vacuum", 
       				fLogicPMT22, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ23 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ23 = new G4LogicalVolume(fSolidPMTJ23, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ23 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ23, 
       				"Vacuum", 
       				fLogicPMT23, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ24 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ24 = new G4LogicalVolume(fSolidPMTJ24, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ24 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ24, 
       				"Vacuum", 
       				fLogicPMT24, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ25 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ25 = new G4LogicalVolume(fSolidPMTJ25, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ25 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ25, 
       				"Vacuum", 
       				fLogicPMT25, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ26 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ26 = new G4LogicalVolume(fSolidPMTJ26, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ26 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ26, 
       				"Vacuum", 
       				fLogicPMT26, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ27 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ27 = new G4LogicalVolume(fSolidPMTJ27, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ27 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ27, 
       				"Vacuum", 
       				fLogicPMT27, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ28 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ28 = new G4LogicalVolume(fSolidPMTJ28, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ28 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ28, 
       				"Vacuum", 
       				fLogicPMT28, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ29 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ29 = new G4LogicalVolume(fSolidPMTJ29, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ29 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ29, 
       				"Vacuum", 
       				fLogicPMT29, 
       				false, 
    			    	0);
    			    	
  fSolidPMTJ30 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter, 0.5*fVacuumLength, 0.*deg, 360.*deg);		    	
  fLogicPMTJ30 = new G4LogicalVolume(fSolidPMTJ30, fPMTIntMaterial,"Vacuum");
  fPhysiPMTJ30 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTJ30, 
       				"Vacuum", 
       				fLogicPMT30, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL = new G4LogicalVolume(fSolidPMTL, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL, 
       				"Photocathode", 
       				fLogicPMT, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL2 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL2 = new G4LogicalVolume(fSolidPMTL2, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL2, 
       				"Photocathode", 
       				fLogicPMT2, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL3 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL3 = new G4LogicalVolume(fSolidPMTL3, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL3, 
       				"Photocathode", 
       				fLogicPMT3, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL4 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL4 = new G4LogicalVolume(fSolidPMTL4, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL4, 
       				"Photocathode", 
       				fLogicPMT4, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL5 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL5 = new G4LogicalVolume(fSolidPMTL5, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL5, 
       				"Photocathode", 
       				fLogicPMT5, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL6 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL6 = new G4LogicalVolume(fSolidPMTL6, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL6 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL6, 
       				"Photocathode", 
       				fLogicPMT6, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL7 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL7 = new G4LogicalVolume(fSolidPMTL7, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL7 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL7, 
       				"Photocathode", 
       				fLogicPMT7, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL8 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL8 = new G4LogicalVolume(fSolidPMTL8, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL8 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL8, 
       				"Photocathode", 
       				fLogicPMT8, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL9 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL9 = new G4LogicalVolume(fSolidPMTL9, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL9 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL9, 
       				"Photocathode", 
       				fLogicPMT9, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL10 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL10 = new G4LogicalVolume(fSolidPMTL10, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL10 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL10, 
       				"Photocathode", 
       				fLogicPMT10, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL11 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL11 = new G4LogicalVolume(fSolidPMTL11, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL11 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL11, 
       				"Photocathode", 
       				fLogicPMT11, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL12 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL12 = new G4LogicalVolume(fSolidPMTL12, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL12 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL12, 
       				"Photocathode", 
       				fLogicPMT12, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL13 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL13 = new G4LogicalVolume(fSolidPMTL13, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL13 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL13, 
       				"Photocathode", 
       				fLogicPMT13, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL14 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL14 = new G4LogicalVolume(fSolidPMTL14, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL14 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL14, 
       				"Photocathode", 
       				fLogicPMT14, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL15 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL15 = new G4LogicalVolume(fSolidPMTL15, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL15 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL15, 
       				"Photocathode", 
       				fLogicPMT15, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL16 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL16 = new G4LogicalVolume(fSolidPMTL16, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL16 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL16, 
       				"Photocathode", 
       				fLogicPMT16, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL17 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL17 = new G4LogicalVolume(fSolidPMTL17, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL17 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL17, 
       				"Photocathode", 
       				fLogicPMT17, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL18 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL18 = new G4LogicalVolume(fSolidPMTL18, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL18 = new G4PVPlacement(0, 
        		    G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL18, 
       				"Photocathode", 
       				fLogicPMT18, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL19 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL19 = new G4LogicalVolume(fSolidPMTL19, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL19 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL19, 
       				"Photocathode", 
       				fLogicPMT19, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL20 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL20 = new G4LogicalVolume(fSolidPMTL20, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL20 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL20, 
       				"Photocathode", 
       				fLogicPMT20, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL21 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL21 = new G4LogicalVolume(fSolidPMTL21, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL21 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL21, 
       				"Photocathode", 
       				fLogicPMT21, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL22 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL22 = new G4LogicalVolume(fSolidPMTL22, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL22 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL22, 
       				"Photocathode", 
       				fLogicPMT22, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL23 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL23 = new G4LogicalVolume(fSolidPMTL23, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL23 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL23,
       				"Photocathode", 
       				fLogicPMT23, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL24 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL24 = new G4LogicalVolume(fSolidPMTL24, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL24 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL24, 
       				"Photocathode", 
       				fLogicPMT24, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL25 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL25 = new G4LogicalVolume(fSolidPMTL25, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL25 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL25, 
       				"Photocathode", 
       				fLogicPMT25, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL26 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL26 = new G4LogicalVolume(fSolidPMTL26, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL26 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL26, 
       				"Photocathode", 
       				fLogicPMT26, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL27 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL27 = new G4LogicalVolume(fSolidPMTL27, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL27 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL27, 
       				"Photocathode", 
       				fLogicPMT27, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL28 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL28 = new G4LogicalVolume(fSolidPMTL28, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL28 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL28, 
       				"Photocathode", 
       				fLogicPMT28, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL29 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL29 = new G4LogicalVolume(fSolidPMTL29, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL29 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL29, 
       				"Photocathode", 
       				fLogicPMT29, 
       				false, 
    			    	0);
    			    	
  fSolidPMTL30 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter, 0.5*fPhotoCathodeLength, 0.*deg, 360.*deg);		    	
  fLogicPMTL30 = new G4LogicalVolume(fSolidPMTL30, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTL30 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode), 
        			fLogicPMTL30, 
       				"Photocathode", 
       				fLogicPMT30, 
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
    			    	
  fSolidPMTI = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter1, 0.5*fVacuumLength1, 0.*deg, 360.*deg);			    	
  fLogicPMTI = new G4LogicalVolume(fSolidPMTI, fPMTIntMaterial,"Vacuum");
  fPhysiPMTI = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTI, 
       				"Vacuum", 
       				fLogicPMT1, 
       				false, 
    			    	0);
  
  fSolidPMTK = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter1, 0.5*fPhotoCathodeLength1, 0.*deg, 360.*deg);			    	
  fLogicPMTK = new G4LogicalVolume(fSolidPMTK, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTK = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode1), 
        			fLogicPMTK, 
       				"Photocathode", 
       				fLogicPMT1, 
       				false, 
    			    	0);
    			    	
  fLogicPMTK->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL2->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL3->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL4->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL5->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL6->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL7->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL8->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL9->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL10->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL11->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL12->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL13->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL14->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL15->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL16->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL17->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL18->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL19->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL20->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL21->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL22->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL23->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL24->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL25->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL26->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL27->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL28->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL29->SetVisAttributes(fYellowVisAtt);
  fLogicPMTL30->SetVisAttributes(fYellowVisAtt);
  }
  if (fDetectorGeometry == 5){
  // Timing method arrangement for 5 GRIFFIN LaBr3:Ce detectors.
  
  fSolidDetector1 = new G4Tubs("GRIFFIN LaBr3:Ce Detector 0", 0., 0.5*fTotalDetectorDiameter1, 0.5*fTotalDetectorLength1, 0.*deg, 360.*deg);
  fSolidDetector2 = new G4Tubs("GRIFFIN LaBr3:Ce Detector 1", 0., 0.5*fTotalDetectorDiameter1, 0.5*fTotalDetectorLength1, 0.*deg, 360.*deg);
  fSolidDetector3 = new G4Tubs("GRIFFIN LaBr3:Ce Detector 2", 0., 0.5*fTotalDetectorDiameter1, 0.5*fTotalDetectorLength1, 0.*deg, 360.*deg);
  fSolidDetector4 = new G4Tubs("GRIFFIN LaBr3:Ce Detector 3", 0., 0.5*fTotalDetectorDiameter1, 0.5*fTotalDetectorLength1, 0.*deg, 360.*deg);
  fSolidDetector5 = new G4Tubs("GRIFFIN LaBr3:Ce Detector 4", 0., 0.5*fTotalDetectorDiameter1, 0.5*fTotalDetectorLength1, 0.*deg, 360.*deg);
  
  fLogicDetector1 = new G4LogicalVolume(fSolidDetector1,fDetectorMaterial, "GRIFFIN LaBr3:Ce Detector 0");
  fLogicDetector2 = new G4LogicalVolume(fSolidDetector2,fDetectorMaterial, "GRIFFIN LaBr3:Ce Detector 1");
  fLogicDetector3 = new G4LogicalVolume(fSolidDetector3, fDetectorMaterial, "GRIFFIN LaBr3:Ce Detector 2");
  fLogicDetector4 = new G4LogicalVolume(fSolidDetector4,fDetectorMaterial, "GRIFFIN LaBr3:Ce Detector 3");
  fLogicDetector5 = new G4LogicalVolume(fSolidDetector5, fDetectorMaterial, "GRIFFIN LaBr3:Ce Detector 4");
  
  G4AssemblyVolume* assemblyDetector = new G4AssemblyVolume();
  // 0
  //P.setX(-12*cm); P.setY(0*cm); P.setZ(-10.295*cm);
  P.setX(-2*fTotalDetectorDiameter1 - 2*fAirGap); P.setY(3.175*cm + 0.5*fTotalDetectorDiameter1); P.setZ(fZpos1 - 0.5*fTotalDetectorLength1);
  //P.setX(-2*fTotalDetectorDiameter1 - 2*fAirGap); P.setY(0*cm); P.setZ(-fOutTBoxY - 0.5*fTotalDetectorLength1);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector1,Tr);
  // 1
  //P.setX(-6*cm); P.setY(0*cm); P.setZ(-10.295*cm);
  P.setX(-fTotalDetectorDiameter1 - fAirGap); P.setY(3.175*cm + 0.5*fTotalDetectorDiameter1); P.setZ(fZpos1 - 0.5*fTotalDetectorLength1);
  //P.setX(-fTotalDetectorDiameter1 - fAirGap); P.setY(0*cm); P.setZ(-fOutTBoxY - 0.5*fTotalDetectorLength1);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector2,Tr);
  // 2
  //P.setX(0*cm); P.setY(0*cm); P.setZ(-10.295*cm);
  P.setX(0*cm); P.setY(3.175*cm + 0.5*fTotalDetectorDiameter1); P.setZ(fZpos1 - 0.5*fTotalDetectorLength1);
  //P.setX(0*cm); P.setY(0*cm); P.setZ(-fOutTBoxY - 0.5*fTotalDetectorLength1);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector3,Tr);
  // 3
  //P.setX(6*cm); P.setY(0*cm); P.setZ(-10.295*cm);
  P.setX(fTotalDetectorDiameter1 + fAirGap); P.setY(3.175*cm + 0.5*fTotalDetectorDiameter1); P.setZ(fZpos1 - 0.5*fTotalDetectorLength1);
  //P.setX(fTotalDetectorDiameter1 + fAirGap); P.setY(0*cm); P.setZ(-fOutTBoxY - 0.5*fTotalDetectorLength1);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector4,Tr);
  // 4
  P.setX(12*cm); P.setY(0*cm); P.setZ(-10.295*cm);
  P.setX(2*fTotalDetectorDiameter1 + 2*fAirGap); P.setY(3.175*cm + 0.5*fTotalDetectorDiameter1); P.setZ(fZpos1 - 0.5*fTotalDetectorLength1);
  //P.setX(2*fTotalDetectorDiameter1 + 2*fAirGap); P.setY(0*cm); P.setZ(-fOutTBoxY - 0.5*fTotalDetectorLength1);
  Tr = G4Transform3D(rotm,P);
  assemblyDetector->AddPlacedVolume(fLogicDetector5,Tr);
  // Assembly Placement in the World
  // Rotated 90° along the axes to be perpendicular to the gas target.
  // Rotated 180° along the z-axis the place the array downside up instead of upside down.
  // Coordinate transformation: (X, Y, Z) -> (Z, Y, -X)
  G4ThreeVector WorldP(0,0,0);
  G4RotationMatrix worldrotm = G4RotationMatrix(90.*deg,90.*deg,270.*deg); 
  G4Transform3D WorldTr = G4Transform3D(worldrotm, WorldP);
  assemblyDetector->MakeImprint(fLogicWorld, WorldTr);
  
  fSolidPMTWin1 = new G4Tubs("PMTWin", 0., 0.5*fDetectorDiameter1, 0.5*fPMTWinThickness1, 0.*deg, 360.*deg);
  fLogicPMTWin1 = new G4LogicalVolume(fSolidPMTWin1, fPMTWinMaterial, "PMTWin");
  fPhysiPMTWin1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMTWin1), 
        			fLogicPMTWin1, 
       				"PMTWin", 
       				fLogicDetector1, 
       				false, 
    			    	0);
    			    	
  fSolidPMTWin2 = new G4Tubs("PMTWin", 0., 0.5*fDetectorDiameter1, 0.5*fPMTWinThickness1, 0.*deg, 360.*deg);
  fLogicPMTWin2 = new G4LogicalVolume(fSolidPMTWin2, fPMTWinMaterial, "PMTWin");
  fPhysiPMTWin2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMTWin1), 
        			fLogicPMTWin2, 
       				"PMTWin", 
       				fLogicDetector2, 
       				false, 
    			    	0);
    			    	
  fSolidPMTWin3 = new G4Tubs("PMTWin", 0., 0.5*fDetectorDiameter1, 0.5*fPMTWinThickness1, 0.*deg, 360.*deg);
  fLogicPMTWin3 = new G4LogicalVolume(fSolidPMTWin3, fPMTWinMaterial, "PMTWin");
  fPhysiPMTWin3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMTWin1), 
        			fLogicPMTWin3, 
       				"PMTWin", 
       				fLogicDetector3, 
       				false, 
    			    	0);
    			    	
  fSolidPMTWin4 = new G4Tubs("PMTWin", 0., 0.5*fDetectorDiameter1, 0.5*fPMTWinThickness1, 0.*deg, 360.*deg);
  fLogicPMTWin4 = new G4LogicalVolume(fSolidPMTWin4, fPMTWinMaterial, "PMTWin");
  fPhysiPMTWin4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMTWin1), 
        			fLogicPMTWin4, 
       				"PMTWin", 
       				fLogicDetector4, 
       				false, 
    			    	0);
    			    	
  fSolidPMTWin5 = new G4Tubs("PMTWin", 0., 0.5*fDetectorDiameter1, 0.5*fPMTWinThickness1, 0.*deg, 360.*deg);
  fLogicPMTWin5 = new G4LogicalVolume(fSolidPMTWin5, fPMTWinMaterial, "PMTWin");
  fPhysiPMTWin5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMTWin1), 
        			fLogicPMTWin5, 
       				"PMTWin", 
       				fLogicDetector5, 
       				false, 
    			    	0);
    			    	
  fLogicPMTWin1->SetVisAttributes(fGreenVisAtt);			    	
  fLogicPMTWin2->SetVisAttributes(fGreenVisAtt);
  fLogicPMTWin3->SetVisAttributes(fGreenVisAtt);
  fLogicPMTWin4->SetVisAttributes(fGreenVisAtt);
  fLogicPMTWin5->SetVisAttributes(fGreenVisAtt);
  
  fSolidCrystal1 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter1, 0.5*fDetectorLength1, 0.*deg, 360.*deg);			    	
  fLogicCrystal1 = new G4LogicalVolume(fSolidCrystal1, fCrystalMaterial1, "Crystal");
  fPhysiCrystal1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos1), 
        			fLogicCrystal1, 
       				"Crystal", 
       				fLogicDetector1, 
       				false, 
    			    	0);
    			    	
  fSolidCrystal2 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter1, 0.5*fDetectorLength1, 0.*deg, 360.*deg);			    	
  fLogicCrystal2 = new G4LogicalVolume(fSolidCrystal2, fCrystalMaterial1, "Crystal");
  fPhysiCrystal2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos1), 
        			fLogicCrystal2, 
       				"Crystal", 
       				fLogicDetector2, 
       				false, 
    			    	1);
    			    	
  fSolidCrystal3 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter1, 0.5*fDetectorLength1, 0.*deg, 360.*deg);			    	
  fLogicCrystal3 = new G4LogicalVolume(fSolidCrystal3, fCrystalMaterial1, "Crystal");
  fPhysiCrystal3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos1), 
        			fLogicCrystal3, 
       				"Crystal", 
       				fLogicDetector3, 
       				false, 
    			    	2);
    			    	
  fSolidCrystal4 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter1, 0.5*fDetectorLength1, 0.*deg, 360.*deg);			    	
  fLogicCrystal4 = new G4LogicalVolume(fSolidCrystal4, fCrystalMaterial1, "Crystal");
  fPhysiCrystal4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos1), 
        			fLogicCrystal4, 
       				"Crystal", 
       				fLogicDetector4, 
       				false, 
    			    	3);
    			    	
  fSolidCrystal5 = new G4Tubs("Crystal", 0., 0.5*fDetectorDiameter1, 0.5*fDetectorLength1, 0.*deg, 360.*deg);			    	
  fLogicCrystal5 = new G4LogicalVolume(fSolidCrystal5, fCrystalMaterial1, "Crystal");
  fPhysiCrystal5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZpos1), 
        			fLogicCrystal5, 
       				"Crystal", 
       				fLogicDetector5, 
       				false, 
    			    	4);
    			    	
  fLogicCrystal1->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal2->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal3->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal4->SetVisAttributes(fMagnetaVisAtt);
  fLogicCrystal5->SetVisAttributes(fMagnetaVisAtt);
  
  fSolidFaceGap1 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter1, 0.5*fGapFaceThickness1, 0.*deg, 360.*deg);			    	
  fLogicFaceGap1 = new G4LogicalVolume(fSolidFaceGap1, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap1), 
        			fLogicFaceGap1, 
       				"FaceGap", 
       				fLogicDetector1, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap2 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter1, 0.5*fGapFaceThickness1, 0.*deg, 360.*deg);			    	
  fLogicFaceGap2 = new G4LogicalVolume(fSolidFaceGap2, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap1), 
        			fLogicFaceGap2, 
       				"FaceGap", 
       				fLogicDetector2, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap3 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter1, 0.5*fGapFaceThickness1, 0.*deg, 360.*deg);			    	
  fLogicFaceGap3 = new G4LogicalVolume(fSolidFaceGap3, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap1), 
        			fLogicFaceGap3, 
       				"FaceGap", 
       				fLogicDetector3, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap4 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter1, 0.5*fGapFaceThickness1, 0.*deg, 360.*deg);			    	
  fLogicFaceGap4 = new G4LogicalVolume(fSolidFaceGap4, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap1), 
        			fLogicFaceGap4, 
       				"FaceGap", 
       				fLogicDetector4, 
       				false, 
    			    	0);
    			    	
  fSolidFaceGap5 = new G4Tubs("FaceGap", 0., 0.5*fDetectorDiameter1, 0.5*fGapFaceThickness1, 0.*deg, 360.*deg);			    	
  fLogicFaceGap5 = new G4LogicalVolume(fSolidFaceGap5, fFaceGapMaterial, "FaceGap");
  fPhysiFaceGap5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceGap1), 
        			fLogicFaceGap5, 
       				"FaceGap", 
       				fLogicDetector5, 
       				false, 
    			    	0);
    			    	
  fLogicFaceGap1->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap2->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap3->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap4->SetVisAttributes(fAuxEdgeVisAtt);
  fLogicFaceGap5->SetVisAttributes(fAuxEdgeVisAtt);
    			    
  fLogicFaceGap1->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap2->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap3->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap4->SetVisAttributes(fCyanVisAtt);
  fLogicFaceGap5->SetVisAttributes(fCyanVisAtt);
  
  fSolidAlCase1 = new G4Tubs("AlCase", 0.5*fDetectorDiameter1, 0.5*fAlCaseDiameter1, 0.5*fAlCaseLength1, 0.*deg, 360.*deg); 			    	
  fLogicAlCase1 = new G4LogicalVolume(fSolidAlCase1, fAlCaseMaterial, "AlCase");
  fPhysiAlCase1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase1), 
        			fLogicAlCase1, 
       				"AlCase", 
       				fLogicDetector1, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase2 = new G4Tubs("AlCase", 0.5*fDetectorDiameter1, 0.5*fAlCaseDiameter1, 0.5*fAlCaseLength1, 0.*deg, 360.*deg); 			    	
  fLogicAlCase2 = new G4LogicalVolume(fSolidAlCase2, fAlCaseMaterial, "AlCase");
  fPhysiAlCase2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase1), 
        			fLogicAlCase2, 
       				"AlCase", 
       				fLogicDetector2, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase3 = new G4Tubs("AlCase", 0.5*fDetectorDiameter1, 0.5*fAlCaseDiameter1, 0.5*fAlCaseLength1, 0.*deg, 360.*deg); 			    	
  fLogicAlCase3 = new G4LogicalVolume(fSolidAlCase3, fAlCaseMaterial, "AlCase");
  fPhysiAlCase3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase1), 
        			fLogicAlCase3, 
       				"AlCase", 
       				fLogicDetector3, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase4 = new G4Tubs("AlCase", 0.5*fDetectorDiameter1, 0.5*fAlCaseDiameter1, 0.5*fAlCaseLength1, 0.*deg, 360.*deg); 			    	
  fLogicAlCase4 = new G4LogicalVolume(fSolidAlCase4, fAlCaseMaterial, "AlCase");
  fPhysiAlCase4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase1), 
        			fLogicAlCase4, 
       				"AlCase", 
       				fLogicDetector4, 
       				false, 
    			    	0);
    			    	
  fSolidAlCase5 = new G4Tubs("AlCase", 0.5*fDetectorDiameter1, 0.5*fAlCaseDiameter1, 0.5*fAlCaseLength1, 0.*deg, 360.*deg); 			    	
  fLogicAlCase5 = new G4LogicalVolume(fSolidAlCase5, fAlCaseMaterial, "AlCase");
  fPhysiAlCase5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposAlCase1), 
        			fLogicAlCase5, 
       				"AlCase", 
       				fLogicDetector5, 
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
    			    	
  fSolidFaceAlCase2 = new G4Tubs("FaceAlCase", 0., 0.5*fDetectorDiameter1, 0.5*fAlFaceThickness1, 0.*deg, 360.*deg);  			    	
  fLogicFaceAlCase2 = new G4LogicalVolume(fSolidFaceAlCase2, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase1), 
        			fLogicFaceAlCase2, 
       				"FaceAlCase", 
       				fLogicDetector2, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase3 = new G4Tubs("FaceAlCase", 0., 0.5*fDetectorDiameter1, 0.5*fAlFaceThickness1, 0.*deg, 360.*deg);  			    	
  fLogicFaceAlCase3 = new G4LogicalVolume(fSolidFaceAlCase3, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase1), 
        			fLogicFaceAlCase3, 
       				"FaceAlCase", 
       				fLogicDetector3, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase4 = new G4Tubs("FaceAlCase", 0., 0.5*fDetectorDiameter1, 0.5*fAlFaceThickness1, 0.*deg, 360.*deg);  			    	
  fLogicFaceAlCase4 = new G4LogicalVolume(fSolidFaceAlCase4, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase1), 
        			fLogicFaceAlCase4, 
       				"FaceAlCase", 
       				fLogicDetector4, 
       				false, 
    			    	0);
    			    	
  fSolidFaceAlCase5 = new G4Tubs("FaceAlCase", 0., 0.5*fDetectorDiameter1, 0.5*fAlFaceThickness1, 0.*deg, 360.*deg);  			    	
  fLogicFaceAlCase5 = new G4LogicalVolume(fSolidFaceAlCase5, fFaceAlCaseMaterial, "FaceAlCase");
  fPhysiFaceAlCase5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposFaceAlCase1), 
        			fLogicFaceAlCase5, 
       				"FaceAlCase", 
       				fLogicDetector5, 
       				false, 
    			    	0);
  
  fLogicAlCase1->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase2->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase3->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase4->SetVisAttributes(fGreyVisAtt);
  fLogicAlCase5->SetVisAttributes(fGreyVisAtt);
  
  fLogicFaceAlCase1->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase2->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase3->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase4->SetVisAttributes(fGreyVisAtt);
  fLogicFaceAlCase5->SetVisAttributes(fGreyVisAtt);
  
  fSolidPMT1 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter1, 0.5*fPMTLength1, 0.*deg, 360.*deg); 			    	
  fLogicPMT1 = new G4LogicalVolume(fSolidPMT1, fPMTMaterial,"PMT");
  fPhysiPMT1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT1), 
        			fLogicPMT1, 
       				"PMT", 
       				fLogicDetector1, 
       				false, 
    			    	0);
    			    	
  fSolidPMT2 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter1, 0.5*fPMTLength1, 0.*deg, 360.*deg); 			    	
  fLogicPMT2 = new G4LogicalVolume(fSolidPMT2, fPMTMaterial,"PMT");
  fPhysiPMT2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT1), 
        			fLogicPMT2, 
       				"PMT", 
       				fLogicDetector2, 
       				false, 
    			    	0);
    			    	
  fSolidPMT3 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter1, 0.5*fPMTLength1, 0.*deg, 360.*deg); 			    	
  fLogicPMT3 = new G4LogicalVolume(fSolidPMT3, fPMTMaterial,"PMT");
  fPhysiPMT3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT1), 
        			fLogicPMT3, 
       				"PMT", 
       				fLogicDetector3, 
       				false, 
    			    	0);
    			    	
  fSolidPMT4 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter1, 0.5*fPMTLength1, 0.*deg, 360.*deg); 			    	
  fLogicPMT4 = new G4LogicalVolume(fSolidPMT4, fPMTMaterial,"PMT");
  fPhysiPMT4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT1), 
        			fLogicPMT4, 
       				"PMT", 
       				fLogicDetector4, 
       				false, 
    			    	0);
    			    	
  fSolidPMT5 = new G4Tubs("PMT", 0., 0.5*fPMTDiameter1, 0.5*fPMTLength1, 0.*deg, 360.*deg); 			    	
  fLogicPMT5 = new G4LogicalVolume(fSolidPMT5, fPMTMaterial,"PMT");
  fPhysiPMT5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPMT1), 
        			fLogicPMT5, 
       				"PMT", 
       				fLogicDetector5, 
       				false, 
    			    	0);
    			    				    	
  fSolidPMTI1 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter1, 0.5*fVacuumLength1, 0.*deg, 360.*deg); 			    	
  fLogicPMTI1 = new G4LogicalVolume(fSolidPMTI1, fPMTIntMaterial,"Vacuum");
  fPhysiPMTI1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTI1, 
       				"Vacuum", 
       				fLogicPMT1, 
       				false, 
    			    	0);
    			    	
  fSolidPMTI2 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter1, 0.5*fVacuumLength1, 0.*deg, 360.*deg); 			    	
  fLogicPMTI2 = new G4LogicalVolume(fSolidPMTI2, fPMTIntMaterial,"Vacuum");
  fPhysiPMTI2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTI2, 
       				"Vacuum", 
       				fLogicPMT2, 
       				false, 
    			    	0);
    			    	
  fSolidPMTI3 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter1, 0.5*fVacuumLength1, 0.*deg, 360.*deg); 			    	
  fLogicPMTI3 = new G4LogicalVolume(fSolidPMTI3, fPMTIntMaterial,"Vacuum");
  fPhysiPMTI3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTI3, 
       				"Vacuum", 
       				fLogicPMT3, 
       				false, 
    			    	0);
    			    	
  fSolidPMTI4 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter1, 0.5*fVacuumLength1, 0.*deg, 360.*deg); 			    	
  fLogicPMTI4 = new G4LogicalVolume(fSolidPMTI4, fPMTIntMaterial,"Vacuum");
  fPhysiPMTI4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTI4, 
       				"Vacuum", 
       				fLogicPMT4, 
       				false, 
    			    	0);
    			    	
  fSolidPMTI5 = new G4Tubs("Vacuum", 0., 0.5*fVacuumDiameter1, 0.5*fVacuumLength1, 0.*deg, 360.*deg); 			    	
  fLogicPMTI5 = new G4LogicalVolume(fSolidPMTI5, fPMTIntMaterial,"Vacuum");
  fPhysiPMTI5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,0.), 
        			fLogicPMTI5, 
       				"Vacuum", 
       				fLogicPMT5, 
       				false, 
    			    	0);
    			    	
  fSolidPMTK1 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter1, 0.5*fPhotoCathodeLength1, 0.*deg, 360.*deg);			    	
  fLogicPMTK1 = new G4LogicalVolume(fSolidPMTK1, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTK1 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode1), 
        			fLogicPMTK1, 
       				"Photocathode", 
       				fLogicPMT1, 
       				false, 
    			    	0);
    			    	
  fSolidPMTK2 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter1, 0.5*fPhotoCathodeLength1, 0.*deg, 360.*deg);			    	
  fLogicPMTK2 = new G4LogicalVolume(fSolidPMTK2, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTK2 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode1), 
        			fLogicPMTK2, 
       				"Photocathode", 
       				fLogicPMT2, 
       				false, 
    			    	0);
    			    	
  fSolidPMTK3 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter1, 0.5*fPhotoCathodeLength1, 0.*deg, 360.*deg);			    	
  fLogicPMTK3 = new G4LogicalVolume(fSolidPMTK3, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTK3 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode1), 
        			fLogicPMTK3, 
       				"Photocathode", 
       				fLogicPMT3, 
       				false, 
    			    	0);
    			    	
  fSolidPMTK4 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter1, 0.5*fPhotoCathodeLength1, 0.*deg, 360.*deg);			    	
  fLogicPMTK4 = new G4LogicalVolume(fSolidPMTK4, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTK4 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode1), 
        			fLogicPMTK4, 
       				"Photocathode", 
       				fLogicPMT4, 
       				false, 
    			    	0);
    			    	
  fSolidPMTK5 = new G4Tubs("Photocathode", 0., 0.5*fPhotoCathodeDiameter1, 0.5*fPhotoCathodeLength1, 0.*deg, 360.*deg);			    	
  fLogicPMTK5 = new G4LogicalVolume(fSolidPMTK5, fPhotoCathodeMaterial,"Photocathode");
  fPhysiPMTK5 = new G4PVPlacement(0, 
        			G4ThreeVector(0.,0.,fZposPhotoCathode1), 
        			fLogicPMTK5, 
       				"Photocathode", 
       				fLogicPMT5, 
       				false, 
    			    	0);
  
  fLogicPMT1->SetVisAttributes(fWireFrameVisAtt);
  fLogicPMT2->SetVisAttributes(fWireFrameVisAtt);
  fLogicPMT3->SetVisAttributes(fWireFrameVisAtt);
  fLogicPMT4->SetVisAttributes(fWireFrameVisAtt);
  fLogicPMT5->SetVisAttributes(fWireFrameVisAtt);
 
  fLogicPMTI1->SetVisAttributes(fWireFrameVisAtt);
  fLogicPMTI2->SetVisAttributes(fWireFrameVisAtt);
  fLogicPMTI3->SetVisAttributes(fWireFrameVisAtt);
  fLogicPMTI4->SetVisAttributes(fWireFrameVisAtt);
  fLogicPMTI5->SetVisAttributes(fWireFrameVisAtt);
   			    
  fLogicPMTK1->SetVisAttributes(fYellowVisAtt);
  fLogicPMTK2->SetVisAttributes(fYellowVisAtt);
  fLogicPMTK3->SetVisAttributes(fYellowVisAtt);
  fLogicPMTK4->SetVisAttributes(fYellowVisAtt);
  fLogicPMTK5->SetVisAttributes(fYellowVisAtt);
  }
  if (fDetectorGeometry == 3 || fDetectorGeometry == 4 || fDetectorGeometry == 5){
 
  // DRAGON Outer Gas Target Box
  // Outer Wall 
  G4RotationMatrix rotm  = G4RotationMatrix(0,0,0);     
  rotm.rotateX(90*deg); 
  rotm.rotateY(-90*deg);
  rotm.rotateZ(0*deg);
  // Position target box so z-axis at beam height.
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
  // Position inner gas cell at beam height relative to the inner gas target box volume.     
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
  
  // Entrance Collimator Thin Metal Disk Without Hole
  rotm  = G4RotationMatrix(0,0,0);     
  position = G4ThreeVector(0.,0.,0.1195*cm);
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
  position = G4ThreeVector(5.351*cm,0.,-2.087*cm);
  transform = G4Transform3D(rotm,position); 
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
  pLowNorm = G4ThreeVector(1.0,0.,-sqrt(3));
  pHighNorm = G4ThreeVector(-1.0,0.,sqrt(3));
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
if (fDetectorGeometry == 5){
  G4cout << "\n" << "The geometry is cylindrical. The detectors are stacked in an array for the timing method (5)." << G4endl;
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
  G4cout << "\n" << "Total Array Detector Diameter: " << G4BestUnit(fTotalDetectorDiameter, "Length") << G4endl;
  G4cout << "\n" << "Array Scintillator Diameter: " << G4BestUnit(fDetectorDiameter, "Length") << G4endl;
  G4cout << "\n" << "Total Array Detector Length: " << G4BestUnit(fTotalDetectorLength, "Length") << G4endl;
  G4cout << "\n" << "Array Scintillator Length: " << G4BestUnit(fDetectorLength, "Length") << G4endl;
  G4cout << "\n" << "Gap Thickness: " << 0.1*fGapThickness << " cm" << G4endl;
  G4cout << "\n" << "Al Thickness: " << 0.1*fAlCaseThickness << " cm" << G4endl;
  G4cout << "\n" << "Pb Thickness: " << 0.1*fPbCaseThickness << " cm" << G4endl;
  G4cout << "\n" << "PMT Diameter: " << G4BestUnit(fPMTDiameter, "Length") << G4endl;
  G4cout << "\n" << "PMT Length: " << G4BestUnit(fPMTLength, "Length") << G4endl;
  G4cout << "\n" << "Total LaBr3:Ce Detector Diameter: " << G4BestUnit(fTotalDetectorDiameter1, "Length") << G4endl;
  G4cout << "\n" << "LaBr3:Ce Scintillator Diameter: " << G4BestUnit(fDetectorDiameter1, "Length") << G4endl;
  G4cout << "\n" << "Total LaBr3:Ce Detector Length: " << G4BestUnit(fTotalDetectorLength1, "Length") << G4endl;
  G4cout << "\n" << "LaBr3:Ce Scintillator Length: " << G4BestUnit(fDetectorLength1, "Length") << G4endl;
  G4cout << "\n" << "Gap Thickness: " << 0.1*fGapThickness1 << " cm" << G4endl;
  G4cout << "\n" << "Al Thickness: " << 0.1*fAlCaseThickness1 << " cm" << G4endl;
  G4cout << "\n" << "PMT Diameter: " << G4BestUnit(fPMTDiameter1, "Length") << G4endl;
  G4cout << "\n" << "PMT Length: " << G4BestUnit(fPMTLength1, "Length") << G4endl;
  G4cout << "\n" << "Air Gap: " << G4BestUnit(fAirGap, "Length") << G4endl;
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
  // Search the material by its name.
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fDetectorMaterial != pttoMaterial) {
    fDetectorMaterial = pttoMaterial;                  
    if(fLogicDetector) fLogicDetector->SetMaterial(fDetectorMaterial);
  }
}

void DetectorConstruction::SetCrystalMaterial(G4String materialChoice)
{
  // Search the material by its name.
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fCrystalMaterial != pttoMaterial) {
    fCrystalMaterial = pttoMaterial;                  
    if(fLogicCrystal) fLogicCrystal->SetMaterial(fCrystalMaterial);
  }
}

void DetectorConstruction::SetGapMaterial(G4String materialChoice)
{
  // Search the material by its name.
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fGapMaterial != pttoMaterial) {
    fGapMaterial = pttoMaterial;                  
    if(fLogicGap) fLogicGap->SetMaterial(fGapMaterial);
  }
}

void DetectorConstruction::SetFaceGapMaterial(G4String materialChoice)
{
  // Search the material by its name.
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fFaceGapMaterial != pttoMaterial) {
    fFaceGapMaterial = pttoMaterial;                  
    if(fLogicFaceGap) fLogicFaceGap->SetMaterial(fFaceGapMaterial);
  }
}

void DetectorConstruction::SetAlCaseMaterial(G4String materialChoice)
{
  // Search the material by its name.
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fAlCaseMaterial != pttoMaterial) {
    fAlCaseMaterial = pttoMaterial;                  
    if(fLogicAlCase) fLogicAlCase->SetMaterial(fAlCaseMaterial);
  }
}

void DetectorConstruction::SetFaceAlCaseMaterial(G4String materialChoice)
{
  // Search the material by its name.
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fFaceAlCaseMaterial != pttoMaterial) {
    fFaceAlCaseMaterial = pttoMaterial;                  
    if(fLogicFaceAlCase) fLogicFaceAlCase->SetMaterial(fFaceAlCaseMaterial);
  }
}

void DetectorConstruction::SetPbCaseMaterial(G4String materialChoice)
{
  // Search the material by its name.
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fPbCaseMaterial != pttoMaterial) {
    fPbCaseMaterial = pttoMaterial;                  
    if(fLogicPbCase) fLogicPbCase->SetMaterial(fPbCaseMaterial);
  }
}

void DetectorConstruction::SetPbCollarMaterial(G4String materialChoice)
{
  // Search the material by its name.
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fPbCollarMaterial != pttoMaterial) {
    fPbCollarMaterial = pttoMaterial;                  
    if(fLogicPbCollar) fLogicPbCollar->SetMaterial(fPbCollarMaterial);
  }
}

void DetectorConstruction::SetPMTMaterial(G4String materialChoice)
{
  // Search the material by its name.
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fPMTMaterial != pttoMaterial) {
    fPMTMaterial = pttoMaterial;                  
    if(fLogicPMT) fLogicPMT->SetMaterial(fPMTMaterial);
  }
}

void DetectorConstruction::SetPMTWinMaterial(G4String materialChoice)
{
  // Search the material by its name.
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fPMTWinMaterial != pttoMaterial) {
    fPMTWinMaterial = pttoMaterial;                  
    if(fLogicPMTWin) fLogicPMTWin->SetMaterial(fPMTWinMaterial);
  }
}

void DetectorConstruction::SetAlCaseDiameter1(G4double val)
{
  fAlCaseDiameter1 = val;
}
void DetectorConstruction::SetAlCaseThickness(G4double val)
{
  fAlCaseThickness = val;
}
void DetectorConstruction::SetAirGap(G4double val)
{
  fAirGap = val;
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
