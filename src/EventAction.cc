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
 fEnergyScint2 = 0.;
 fEnergyResScint2 = 0;
 fEnergyScint3 = 0.;
 fEnergyResScint3 = 0;
 fEnergyScint4 = 0.;
 fEnergyResScint4 = 0;
 fEnergyScint5 = 0.;
 fEnergyResScint5 = 0;
 fEnergyScint6 = 0.;
 fEnergyResScint6 = 0;
 fEnergyScint7 = 0.;
 fEnergyResScint7 = 0;
 fEnergyScint8 = 0.;
 fEnergyResScint8 = 0;
 fEnergyScint9 = 0.;
 fEnergyResScint9 = 0;
 fEnergyScint10 = 0.;
 fEnergyResScint10 = 0;
 fEnergyScint11 = 0.;
 fEnergyResScint11 = 0;
 fEnergyScint12 = 0.;
 fEnergyResScint12 = 0;
 fEnergyScint13 = 0.;
 fEnergyResScint13 = 0;
 fEnergyScint14 = 0.;
 fEnergyResScint14 = 0;
 fEnergyScint15 = 0.;
 fEnergyResScint15 = 0;
 fEnergyScint16 = 0.;
 fEnergyResScint16 = 0;
 fEnergyScint17 = 0.;
 fEnergyResScint17 = 0;
 fEnergyScint18= 0.;
 fEnergyResScint18 = 0;
 fEnergyScint19 = 0.;
 fEnergyResScint19 = 0;
 fEnergyScint20 = 0.;
 fEnergyResScint20 = 0;
 fEnergyScint21 = 0.;
 fEnergyResScint21 = 0;
 fEnergyScint22 = 0.;
 fEnergyResScint22 = 0;
 fEnergyScint23 = 0.;
 fEnergyResScint23 = 0;
 fEnergyScint24 = 0.;
 fEnergyResScint24 = 0;
 fEnergyScint25 = 0.;
 fEnergyResScint25 = 0;
 fEnergyScint26 = 0.;
 fEnergyResScint26 = 0;
 fEnergyScint27 = 0.;
 fEnergyResScint27 = 0;
 fEnergyScint28 = 0.;
 fEnergyResScint28 = 0;
 fEnergyScint29 = 0.;
 fEnergyResScint29 = 0;
 fEnergyScint30 = 0.;
 fEnergyResScint30 = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* )
{
  //accumulates statistic
  //
  fRunAct->fillPerEvent(fEnergyScint);
  fRunAct->fillPerEvent(fEnergyResScint);
  fRunAct->fillPerEvent(fEnergyScint2);
  fRunAct->fillPerEvent(fEnergyResScint2);
  fRunAct->fillPerEvent(fEnergyScint3);
  fRunAct->fillPerEvent(fEnergyResScint3);
  fRunAct->fillPerEvent(fEnergyScint4);
  fRunAct->fillPerEvent(fEnergyResScint4);
  fRunAct->fillPerEvent(fEnergyScint5);
  fRunAct->fillPerEvent(fEnergyResScint5);
  fRunAct->fillPerEvent(fEnergyScint6);
  fRunAct->fillPerEvent(fEnergyResScint6);
  fRunAct->fillPerEvent(fEnergyScint7);
  fRunAct->fillPerEvent(fEnergyResScint7);
  fRunAct->fillPerEvent(fEnergyScint8);
  fRunAct->fillPerEvent(fEnergyResScint8);
  fRunAct->fillPerEvent(fEnergyScint9);
  fRunAct->fillPerEvent(fEnergyResScint9);
  fRunAct->fillPerEvent(fEnergyScint10);
  fRunAct->fillPerEvent(fEnergyResScint10);
  fRunAct->fillPerEvent(fEnergyScint11);
  fRunAct->fillPerEvent(fEnergyResScint11);
  fRunAct->fillPerEvent(fEnergyScint12);
  fRunAct->fillPerEvent(fEnergyResScint12);
  fRunAct->fillPerEvent(fEnergyScint13);
  fRunAct->fillPerEvent(fEnergyResScint13);
  fRunAct->fillPerEvent(fEnergyScint14);
  fRunAct->fillPerEvent(fEnergyResScint14);
  fRunAct->fillPerEvent(fEnergyScint15);
  fRunAct->fillPerEvent(fEnergyResScint15);
  fRunAct->fillPerEvent(fEnergyScint16);
  fRunAct->fillPerEvent(fEnergyResScint16);
  fRunAct->fillPerEvent(fEnergyScint17);
  fRunAct->fillPerEvent(fEnergyResScint17);
  fRunAct->fillPerEvent(fEnergyScint18);
  fRunAct->fillPerEvent(fEnergyResScint18);
  fRunAct->fillPerEvent(fEnergyScint19);
  fRunAct->fillPerEvent(fEnergyResScint19);
  fRunAct->fillPerEvent(fEnergyScint20);
  fRunAct->fillPerEvent(fEnergyResScint20);
  fRunAct->fillPerEvent(fEnergyScint21);
  fRunAct->fillPerEvent(fEnergyResScint21);
  fRunAct->fillPerEvent(fEnergyScint22);
  fRunAct->fillPerEvent(fEnergyResScint22);
  fRunAct->fillPerEvent(fEnergyScint23);
  fRunAct->fillPerEvent(fEnergyResScint23);
  fRunAct->fillPerEvent(fEnergyScint24);
  fRunAct->fillPerEvent(fEnergyResScint24);
  fRunAct->fillPerEvent(fEnergyScint25);
  fRunAct->fillPerEvent(fEnergyResScint25);
  fRunAct->fillPerEvent(fEnergyScint26);
  fRunAct->fillPerEvent(fEnergyResScint26);
  fRunAct->fillPerEvent(fEnergyScint27);
  fRunAct->fillPerEvent(fEnergyResScint27);
  fRunAct->fillPerEvent(fEnergyScint28);
  fRunAct->fillPerEvent(fEnergyResScint28);
  fRunAct->fillPerEvent(fEnergyScint29);
  fRunAct->fillPerEvent(fEnergyResScint29);
  fRunAct->fillPerEvent(fEnergyScint30);
  fRunAct->fillPerEvent(fEnergyResScint30);
  
  //fill histograms

  // do not increment the energy histogram if no energy loss
  if (fEnergyScint > 1*CLHEP::eV)
    fHistoManager->FillHisto(1, fEnergyScint);
  // do not increment the energy histogram if no energy loss
  // if (fEnergyResScint > 1*CLHEP::eV)
  if (fEnergyScint > 1*CLHEP::eV)
    fHistoManager->FillHisto(2, fEnergyResScint);
  //if ((fEnergyScint > 0.) && (fEnergyScint < 1*CLHEP::eV))
    //fHistoManager->FillHisto(3, fEnergyScint);
    
  if (fEnergyScint2 > 1*CLHEP::eV)
    fHistoManager->FillHisto(3, fEnergyScint2);
  if (fEnergyScint2 > 1*CLHEP::eV)
    fHistoManager->FillHisto(4, fEnergyResScint2);
  if (fEnergyScint3 > 1*CLHEP::eV)
    fHistoManager->FillHisto(5, fEnergyScint3);
  if (fEnergyScint3 > 1*CLHEP::eV)
    fHistoManager->FillHisto(6, fEnergyResScint3);
  if (fEnergyScint4 > 1*CLHEP::eV)
    fHistoManager->FillHisto(7, fEnergyScint4);
  if (fEnergyScint4 > 1*CLHEP::eV)
    fHistoManager->FillHisto(8, fEnergyResScint4);
  if (fEnergyScint5 > 1*CLHEP::eV)
    fHistoManager->FillHisto(9, fEnergyScint5);
  if (fEnergyScint5 > 1*CLHEP::eV)
    fHistoManager->FillHisto(10, fEnergyResScint5);
  if (fEnergyScint6 > 1*CLHEP::eV)
    fHistoManager->FillHisto(11, fEnergyScint6);
  if (fEnergyScint6 > 1*CLHEP::eV)
    fHistoManager->FillHisto(12, fEnergyResScint6);
  if (fEnergyScint7 > 1*CLHEP::eV)
    fHistoManager->FillHisto(13, fEnergyScint7);
  if (fEnergyScint7 > 1*CLHEP::eV)
    fHistoManager->FillHisto(14, fEnergyResScint7);
  if (fEnergyScint8 > 1*CLHEP::eV)
    fHistoManager->FillHisto(15, fEnergyScint8);
  if (fEnergyScint8 > 1*CLHEP::eV)
    fHistoManager->FillHisto(16, fEnergyResScint8);
  if (fEnergyScint9 > 1*CLHEP::eV)
    fHistoManager->FillHisto(17, fEnergyScint9);
  if (fEnergyScint9 > 1*CLHEP::eV)
    fHistoManager->FillHisto(18, fEnergyResScint9);
  if (fEnergyScint10 > 1*CLHEP::eV)
    fHistoManager->FillHisto(19, fEnergyScint10);
  if (fEnergyScint10 > 1*CLHEP::eV)
    fHistoManager->FillHisto(20, fEnergyResScint10);
  if (fEnergyScint11 > 1*CLHEP::eV)
    fHistoManager->FillHisto(21, fEnergyScint11);
  if (fEnergyScint11 > 1*CLHEP::eV)
    fHistoManager->FillHisto(22, fEnergyResScint11);
  if (fEnergyScint12 > 1*CLHEP::eV)
    fHistoManager->FillHisto(23, fEnergyScint12);
  if (fEnergyScint12 > 1*CLHEP::eV)
    fHistoManager->FillHisto(24, fEnergyResScint12);
  if (fEnergyScint13 > 1*CLHEP::eV)
    fHistoManager->FillHisto(25, fEnergyScint13);
  if (fEnergyScint13 > 1*CLHEP::eV)
    fHistoManager->FillHisto(26, fEnergyResScint13);
  if (fEnergyScint14 > 1*CLHEP::eV)
    fHistoManager->FillHisto(27, fEnergyScint14);
  if (fEnergyScint14 > 1*CLHEP::eV)
    fHistoManager->FillHisto(28, fEnergyResScint14);
  if (fEnergyScint15 > 1*CLHEP::eV)
    fHistoManager->FillHisto(29, fEnergyScint15);
  if (fEnergyScint15 > 1*CLHEP::eV)
    fHistoManager->FillHisto(30, fEnergyResScint15);
  if (fEnergyScint16 > 1*CLHEP::eV)
    fHistoManager->FillHisto(31, fEnergyScint16);
  if (fEnergyScint16 > 1*CLHEP::eV)
    fHistoManager->FillHisto(32, fEnergyResScint16);
  if (fEnergyScint17 > 1*CLHEP::eV)
    fHistoManager->FillHisto(33, fEnergyScint17);
  if (fEnergyScint17 > 1*CLHEP::eV)
    fHistoManager->FillHisto(34, fEnergyResScint17);
  if (fEnergyScint18 > 1*CLHEP::eV)
    fHistoManager->FillHisto(35, fEnergyScint18);
  if (fEnergyScint18 > 1*CLHEP::eV)
    fHistoManager->FillHisto(36, fEnergyResScint18);
  if (fEnergyScint19 > 1*CLHEP::eV)
    fHistoManager->FillHisto(37, fEnergyScint19);
  if (fEnergyScint19 > 1*CLHEP::eV)
    fHistoManager->FillHisto(38, fEnergyResScint19);
  if (fEnergyScint20 > 1*CLHEP::eV)
    fHistoManager->FillHisto(39, fEnergyScint20);
  if (fEnergyScint20 > 1*CLHEP::eV)
    fHistoManager->FillHisto(40, fEnergyResScint20);
  if (fEnergyScint21 > 1*CLHEP::eV)
    fHistoManager->FillHisto(41, fEnergyScint21);
  if (fEnergyScint21 > 1*CLHEP::eV)
    fHistoManager->FillHisto(42, fEnergyResScint21);
  if (fEnergyScint22 > 1*CLHEP::eV)
    fHistoManager->FillHisto(43, fEnergyScint22);
  if (fEnergyScint22 > 1*CLHEP::eV)
    fHistoManager->FillHisto(44, fEnergyResScint22);
  if (fEnergyScint23 > 1*CLHEP::eV)
    fHistoManager->FillHisto(45, fEnergyScint23);
  if (fEnergyScint23 > 1*CLHEP::eV)
    fHistoManager->FillHisto(46, fEnergyResScint23);
  if (fEnergyScint24 > 1*CLHEP::eV)
    fHistoManager->FillHisto(47, fEnergyScint24);
  if (fEnergyScint24 > 1*CLHEP::eV)
    fHistoManager->FillHisto(48, fEnergyResScint24);
  if (fEnergyScint25 > 1*CLHEP::eV)
    fHistoManager->FillHisto(49, fEnergyScint25);
  if (fEnergyScint25 > 1*CLHEP::eV)
    fHistoManager->FillHisto(50, fEnergyResScint25);
  if (fEnergyScint26 > 1*CLHEP::eV)
    fHistoManager->FillHisto(51, fEnergyScint26);
  if (fEnergyScint26 > 1*CLHEP::eV)
    fHistoManager->FillHisto(52, fEnergyResScint26);
  if (fEnergyScint27 > 1*CLHEP::eV)
    fHistoManager->FillHisto(53, fEnergyScint27);
  if (fEnergyScint27 > 1*CLHEP::eV)
    fHistoManager->FillHisto(54, fEnergyResScint27);
  if (fEnergyScint28 > 1*CLHEP::eV)
    fHistoManager->FillHisto(55, fEnergyScint28);
  if (fEnergyScint28 > 1*CLHEP::eV)
    fHistoManager->FillHisto(56, fEnergyResScint28);
  if (fEnergyScint29 > 1*CLHEP::eV)
    fHistoManager->FillHisto(57, fEnergyScint29);
  if (fEnergyScint29 > 1*CLHEP::eV)
    fHistoManager->FillHisto(58, fEnergyResScint29);
  if (fEnergyScint30 > 1*CLHEP::eV)
    fHistoManager->FillHisto(59, fEnergyScint30);
  if (fEnergyScint30 > 1*CLHEP::eV)
    fHistoManager->FillHisto(60, fEnergyResScint30);
  //if (fEnergyScint1 > 1*CLHEP::eV)
    //fHistoManager->FillHisto(61, fEnergyScint1);
  //if (fEnergyScint1 > 1*CLHEP::eV)
    //fHistoManager->FillHisto(62, fEnergyResScint1);

  //fHistoManager->Fill2Histo(2, fEnergyGas,fEnergyDSSSD);

  
  //fill ntuple
  //
//  fHistoManager->FillNtuple(fEnergyAbs, fEnergyGap, fTrackLAbs, fTrackLGap);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
