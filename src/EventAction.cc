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

#include "EventAction.hh"

#include "RunAction.hh"
#include "HistoManager.hh"

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* run, HistoManager* histo)
:fRunAct(run),fHistoManager(histo)
{
 fPrintModulo = 1000; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4int evtNb = evt->GetEventID();
  if (evtNb%fPrintModulo == 0) 
    G4cout << "---> Beginning of event: " << evtNb << G4endl;
 
 // initialisation per event
 fEnergyScint = 0.;
 fEnergyResScint = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* )
{
  //accumulates statistic
  //
  fRunAct->fillPerEvent(fEnergyScint);
  fRunAct->fillPerEvent(fEnergyResScint);
  
  //fill histograms

  // do not increment the energy histogram if no energy loss
  if (fEnergyScint > 1*CLHEP::eV)
    fHistoManager->FillHisto(1, fEnergyScint);
  // do not increment the energy histogram if no energy loss
  if (fEnergyResScint > 1*CLHEP::eV)
    fHistoManager->FillHisto(2, fEnergyResScint);
  if ((fEnergyScint > 0.) && (fEnergyScint < 1*CLHEP::eV))
    fHistoManager->FillHisto(3, fEnergyScint);
  

  //fHistoManager->Fill2Histo(2, fEnergyGas,fEnergyDSSSD);

  
  //fill ntuple
  //
//  fHistoManager->FillNtuple(fEnergyAbs, fEnergyGap, fTrackLAbs, fTrackLGap);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
