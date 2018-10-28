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
/// \file Si_Ion_Chamber_v6/include/PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class
//
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4GeneralParticleSource;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();    
  ~PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event*);
  //virtual void GeneratePrimaries2(G4Event*);
  //virtual void GeneratePrimaries3(G4Event*);
  //virtual void GeneratePrimaries4(G4Event*);
  //virtual void GeneratePrimaries5(G4Event*);
  //virtual void GeneratePrimaries6(G4Event*);
  //virtual void GeneratePrimaries7(G4Event*);
  //virtual void GeneratePrimaries8(G4Event*);
  //virtual void GeneratePrimaries9(G4Event*);
  //virtual void GeneratePrimaries10(G4Event*);
  //virtual void GeneratePrimaries11(G4Event*);
  //virtual void GeneratePrimaries12(G4Event*);
  //virtual void GeneratePrimaries13(G4Event*);
  //virtual void GeneratePrimaries14(G4Event*);
  //virtual void GeneratePrimaries15(G4Event*);
  //virtual void GeneratePrimaries16(G4Event*);
  //virtual void GeneratePrimaries17(G4Event*);
  //virtual void GeneratePrimaries18(G4Event*);
  //virtual void GeneratePrimaries19(G4Event*);
  //virtual void GeneratePrimaries20(G4Event*);
  //virtual void GeneratePrimaries21(G4Event*);
  //virtual void GeneratePrimaries22(G4Event*);
  //virtual void GeneratePrimaries23(G4Event*);
  //virtual void GeneratePrimaries24(G4Event*);
  //virtual void GeneratePrimaries25(G4Event*);
  //virtual void GeneratePrimaries26(G4Event*);
  //virtual void GeneratePrimaries27(G4Event*);
  //virtual void GeneratePrimaries28(G4Event*);
  //virtual void GeneratePrimaries29(G4Event*);
  //virtual void GeneratePrimaries30(G4Event*);
  
private:
  G4GeneralParticleSource* fParticleGun;  //pointer a to G4 class
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
