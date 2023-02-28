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
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class RunAction;
class HistoManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
public:
  EventAction(RunAction*, HistoManager*);
  virtual ~EventAction();

  void  BeginOfEventAction(const G4Event*);
  void  EndOfEventAction(const G4Event*);
  
  void AddScint(G4double de) {fEnergyScint += de;};
  void AddResScint(G4double re) {fEnergyResScint += re;}; 
  void AddScint1(G4double de) {fEnergyScint1 += de;};
  void AddResScint1(G4double re) {fEnergyResScint1 += re;};
  void AddScint2(G4double de2) {fEnergyScint2 += de2;};
  void AddResScint2(G4double re2) {fEnergyResScint2 += re2;};
  void AddScint3(G4double de3) {fEnergyScint3 += de3;};
  void AddResScint3(G4double re3) {fEnergyResScint3 += re3;};
  void AddScint4(G4double de4) {fEnergyScint4 += de4;};
  void AddResScint4(G4double re4) {fEnergyResScint4 += re4;};
  void AddScint5(G4double de5) {fEnergyScint5 += de5;};
  void AddResScint5(G4double re5) {fEnergyResScint5 += re5;};
  void AddScint6(G4double de6) {fEnergyScint6 += de6;};
  void AddResScint6(G4double re6) {fEnergyResScint6 += re6;};
  void AddScint7(G4double de7) {fEnergyScint7 += de7;};
  void AddResScint7(G4double re7) {fEnergyResScint7 += re7;};
  void AddScint8(G4double de8) {fEnergyScint8 += de8;};
  void AddResScint8(G4double re8) {fEnergyResScint8 += re8;};
  void AddScint9(G4double de9) {fEnergyScint9 += de9;};
  void AddResScint9(G4double re9) {fEnergyResScint9 += re9;};
  void AddScint10(G4double de10) {fEnergyScint10 += de10;};
  void AddResScint10(G4double re10) {fEnergyResScint10 += re10;};
  void AddScint11(G4double de11) {fEnergyScint11 += de11;};
  void AddResScint11(G4double re11) {fEnergyResScint11 += re11;};
  void AddScint12(G4double de12) {fEnergyScint12 += de12;};
  void AddResScint12(G4double re12) {fEnergyResScint12 += re12;};
  void AddScint13(G4double de13) {fEnergyScint13 += de13;};
  void AddResScint13(G4double re13) {fEnergyResScint13 += re13;};
  void AddScint14(G4double de14) {fEnergyScint14 += de14;};
  void AddResScint14(G4double re14) {fEnergyResScint14 += re14;};
  void AddScint15(G4double de15) {fEnergyScint15 += de15;};
  void AddResScint15(G4double re15) {fEnergyResScint15 += re15;};
  void AddScint16(G4double de16) {fEnergyScint16 += de16;};
  void AddResScint16(G4double re16) {fEnergyResScint16 += re16;};
  void AddScint17(G4double de17) {fEnergyScint17 += de17;};
  void AddResScint17(G4double re17) {fEnergyResScint17 += re17;};
  void AddScint18(G4double de18) {fEnergyScint18 += de18;};
  void AddResScint18(G4double re18) {fEnergyResScint18 += re18;};
  void AddScint19(G4double de19) {fEnergyScint19 += de19;};
  void AddResScint19(G4double re19) {fEnergyResScint19 += re19;};
  void AddScint20(G4double de20) {fEnergyScint20 += de20;};
  void AddResScint20(G4double re20) {fEnergyResScint20 += re20;};
  void AddScint21(G4double de21) {fEnergyScint21 += de21;};
  void AddResScint21(G4double re21) {fEnergyResScint21 += re21;};
  void AddScint22(G4double de22) {fEnergyScint22 += de22;};
  void AddResScint22(G4double re22) {fEnergyResScint22 += re22;};
  void AddScint23(G4double de23) {fEnergyScint23 += de23;};
  void AddResScint23(G4double re23) {fEnergyResScint23 += re23;};
  void AddScint24(G4double de24) {fEnergyScint24 += de24;};
  void AddResScint24(G4double re24) {fEnergyResScint24 += re24;};
  void AddScint25(G4double de25) {fEnergyScint25 += de25;};
  void AddResScint25(G4double re25) {fEnergyResScint25 += re25;};
  void AddScint26(G4double de26) {fEnergyScint26 += de26;};
  void AddResScint26(G4double re26) {fEnergyResScint26 += re26;};
  void AddScint27(G4double de27) {fEnergyScint27 += de27;};
  void AddResScint27(G4double re27) {fEnergyResScint27 += re27;};
  void AddScint28(G4double de28) {fEnergyScint28 += de28;};
  void AddResScint28(G4double re28) {fEnergyResScint28 += re28;};
  void AddScint29(G4double de29) {fEnergyScint29 += de29;};
  void AddResScint29(G4double re29) {fEnergyResScint29 += re29;};
  void AddScint30(G4double de30) {fEnergyScint30 += de30;};
  void AddResScint30(G4double re30) {fEnergyResScint30 += re30;};
  
    
private:
   RunAction*    fRunAct;
   HistoManager* fHistoManager;
      
   G4double  fEnergyScint;
   G4double  fEnergyResScint;
   G4double  fEnergyScint1;
   G4double  fEnergyResScint1;
   G4double  fEnergyScint2;
   G4double  fEnergyResScint2;
   G4double  fEnergyScint3;
   G4double  fEnergyResScint3;
   G4double  fEnergyScint4;
   G4double  fEnergyResScint4;
   G4double  fEnergyScint5;
   G4double  fEnergyResScint5;
   G4double  fEnergyScint6;
   G4double  fEnergyResScint6;
   G4double  fEnergyScint7;
   G4double  fEnergyResScint7;
   G4double  fEnergyScint8;
   G4double  fEnergyResScint8;
   G4double  fEnergyScint9;
   G4double  fEnergyResScint9;
   G4double  fEnergyScint10;
   G4double  fEnergyResScint10;
   G4double  fEnergyScint11;
   G4double  fEnergyResScint11;
   G4double  fEnergyScint12;
   G4double  fEnergyResScint12;
   G4double  fEnergyScint13;
   G4double  fEnergyResScint13;
   G4double  fEnergyScint14;
   G4double  fEnergyResScint14;
   G4double  fEnergyScint15;
   G4double  fEnergyResScint15;
   G4double  fEnergyScint16;
   G4double  fEnergyResScint16;
   G4double  fEnergyScint17;
   G4double  fEnergyResScint17;
   G4double  fEnergyScint18;
   G4double  fEnergyResScint18;
   G4double  fEnergyScint19;
   G4double  fEnergyResScint19;
   G4double  fEnergyScint20;
   G4double  fEnergyResScint20;
   G4double  fEnergyScint21;
   G4double  fEnergyResScint21;
   G4double  fEnergyScint22;
   G4double  fEnergyResScint22;
   G4double  fEnergyScint23;
   G4double  fEnergyResScint23;
   G4double  fEnergyScint24;
   G4double  fEnergyResScint24;
   G4double  fEnergyScint25;
   G4double  fEnergyResScint25;
   G4double  fEnergyScint26;
   G4double  fEnergyResScint26;
   G4double  fEnergyScint27;
   G4double  fEnergyResScint27;
   G4double  fEnergyScint28;
   G4double  fEnergyResScint28;
   G4double  fEnergyScint29;
   G4double  fEnergyResScint29;
   G4double  fEnergyScint30;
   G4double  fEnergyResScint30;
                     
   G4int     fPrintModulo;                             
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
