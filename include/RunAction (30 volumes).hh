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
/// \file analysis/shared/include/RunAction.hh
/// \brief Definition of the RunAction class
//
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class HistoManager;

class RunAction : public G4UserRunAction
{
public:
  RunAction(HistoManager*);
  virtual ~RunAction();

  void BeginOfRunAction(const G4Run*);
  void   EndOfRunAction(const G4Run*);
    
  void fillPerEvent(G4double);
  void fillPerEvent2(G4double);
  void fillPerEvent3(G4double);
  void fillPerEvent4(G4double);
  void fillPerEvent5(G4double);
  void fillPerEvent6(G4double);
  void fillPerEvent7(G4double);
  void fillPerEvent8(G4double);
  void fillPerEvent9(G4double);
  void fillPerEvent10(G4double);
  void fillPerEvent11(G4double);
  void fillPerEvent12(G4double);
  void fillPerEvent13(G4double);
  void fillPerEvent14(G4double);
  void fillPerEvent15(G4double);
  void fillPerEvent16(G4double);
  void fillPerEvent17(G4double);
  void fillPerEvent18(G4double);
  void fillPerEvent19(G4double);
  void fillPerEvent20(G4double);
  void fillPerEvent21(G4double);
  void fillPerEvent22(G4double);
  void fillPerEvent23(G4double);
  void fillPerEvent24(G4double);
  void fillPerEvent25(G4double);
  void fillPerEvent26(G4double);
  void fillPerEvent27(G4double);
  void fillPerEvent28(G4double);
  void fillPerEvent29(G4double);
  void fillPerEvent30(G4double);
  
private:
  HistoManager* fHistoManager;
  G4double fSumELaBr3;
  G4double fSum2ELaBr3;
  G4double fSumELaBr32;
  G4double fSum2ELaBr32;
  G4double fSumELaBr33;
  G4double fSum2ELaBr33;
  G4double fSumELaBr34;
  G4double fSum2ELaBr34;
  G4double fSumELaBr35;
  G4double fSum2ELaBr35;
  G4double fSumELaBr36;
  G4double fSum2ELaBr36;
  G4double fSumELaBr37;
  G4double fSum2ELaBr37;
  G4double fSumELaBr38;
  G4double fSum2ELaBr38;
  G4double fSumELaBr39;
  G4double fSum2ELaBr39;
  G4double fSumELaBr310;
  G4double fSum2ELaBr310;
  G4double fSumELaBr311;
  G4double fSum2ELaBr311;
  G4double fSumELaBr312;
  G4double fSum2ELaBr312;
  G4double fSumELaBr313;
  G4double fSum2ELaBr313;
  G4double fSumELaBr314;
  G4double fSum2ELaBr314;
  G4double fSumELaBr315;
  G4double fSum2ELaBr315;
  G4double fSumELaBr316;
  G4double fSum2ELaBr316;
  G4double fSumELaBr317;
  G4double fSum2ELaBr317;
  G4double fSumELaBr318;
  G4double fSum2ELaBr318;
  G4double fSumELaBr319;
  G4double fSum2ELaBr319;
  G4double fSumELaBr320;
  G4double fSum2ELaBr320;
  G4double fSumELaBr321;
  G4double fSum2ELaBr321;
  G4double fSumELaBr322;
  G4double fSum2ELaBr322;
  G4double fSumELaBr323;
  G4double fSum2ELaBr323;
  G4double fSumELaBr324;
  G4double fSum2ELaBr324;
  G4double fSumELaBr325;
  G4double fSum2ELaBr325;
  G4double fSumELaBr326;
  G4double fSum2ELaBr326;
  G4double fSumELaBr327;
  G4double fSum2ELaBr327;
  G4double fSumELaBr328;
  G4double fSum2ELaBr328;
  G4double fSumELaBr329;
  G4double fSum2ELaBr329;
  G4double fSumELaBr330;
  G4double fSum2ELaBr330;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
