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
/// \file Si_Ion_Chamber_v6/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
// $Id$
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)
{
  fParticleGun = new G4GeneralParticleSource();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the beginning of event
  // 
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//void PrimaryGeneratorAction::GeneratePrimaries2(G4Event* anEvent2)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent2);
//}

//void PrimaryGeneratorAction::GeneratePrimaries3(G4Event* anEvent3)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent3);
//}

//void PrimaryGeneratorAction::GeneratePrimaries4(G4Event* anEvent4)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent4);
//}

//void PrimaryGeneratorAction::GeneratePrimaries5(G4Event* anEvent5)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent5);
//}

//void PrimaryGeneratorAction::GeneratePrimaries6(G4Event* anEvent6)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent6);
//}

//void PrimaryGeneratorAction::GeneratePrimaries7(G4Event* anEvent7)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent7);
//}

//void PrimaryGeneratorAction::GeneratePrimaries8(G4Event* anEvent8)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent8);
//}

//void PrimaryGeneratorAction::GeneratePrimaries9(G4Event* anEvent9)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent9);
//}

//void PrimaryGeneratorAction::GeneratePrimaries10(G4Event* anEvent10)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent10);
//}

//void PrimaryGeneratorAction::GeneratePrimaries11(G4Event* anEvent11)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent11);
//}

//void PrimaryGeneratorAction::GeneratePrimaries12(G4Event* anEvent12)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent12);
//}

//void PrimaryGeneratorAction::GeneratePrimaries13(G4Event* anEvent13)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent13);
//}

//void PrimaryGeneratorAction::GeneratePrimaries14(G4Event* anEvent14)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent14);
//}

//void PrimaryGeneratorAction::GeneratePrimaries15(G4Event* anEvent15)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent15);
//}

//void PrimaryGeneratorAction::GeneratePrimaries16(G4Event* anEvent16)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent16);
//}

//void PrimaryGeneratorAction::GeneratePrimaries17(G4Event* anEvent17)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent17);
//}

//void PrimaryGeneratorAction::GeneratePrimaries18(G4Event* anEvent18)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent18);
//}

//void PrimaryGeneratorAction::GeneratePrimaries19(G4Event* anEvent19)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent19);
//}

//void PrimaryGeneratorAction::GeneratePrimaries20(G4Event* anEvent20)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent20);
//}

//void PrimaryGeneratorAction::GeneratePrimaries21(G4Event* anEvent21)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent21);
//}

//void PrimaryGeneratorAction::GeneratePrimaries22(G4Event* anEvent22)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent22);
//}

//void PrimaryGeneratorAction::GeneratePrimaries23(G4Event* anEvent23)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent23);
//}

//void PrimaryGeneratorAction::GeneratePrimaries24(G4Event* anEvent24)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent24);
//}

//void PrimaryGeneratorAction::GeneratePrimaries25(G4Event* anEvent25)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent25);
//}

//void PrimaryGeneratorAction::GeneratePrimaries26(G4Event* anEvent26)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent26);
//}

//void PrimaryGeneratorAction::GeneratePrimaries27(G4Event* anEvent27)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent27);
//}

//void PrimaryGeneratorAction::GeneratePrimaries28(G4Event* anEvent28)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent28);
//}

//void PrimaryGeneratorAction::GeneratePrimaries29(G4Event* anEvent29)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent29);
//}

//void PrimaryGeneratorAction::GeneratePrimaries30(G4Event* anEvent30)
//{
  ////this function is called at the beginning of event
  //// 
  //fParticleGun->GeneratePrimaryVertex(anEvent30);
//}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
