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
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"
#include "B1SteppingAction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 10*m, env_sizeZ = 10*m;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;
  G4Material* vacuum =
    new G4Material("Vacuum",      //Name as String
                   1,             //Atomic Number,  in this case we use 1 for hydrogen  
                   1.008*g/mole,  //Mass per Mole "Atomic Weight"  1.008*g/mole for Hydoren 
                   1.e-25*g/cm3,  //Density of Vaccuum  *Cant be Zero, Must be small insted 
                   kStateGas,     //kStateGas for Gas
                   2.73*kelvin,   //Temperature for gas
                   1.e-25*g/cm3); //Pressure for Vaccum

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
//  G4NistManager* man = G4NistManager::Instance();
//  man->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  G4double pressure_of_CO2 = 100./760.;//in atmospheres
  G4Material* CO2 = new G4Material("CO2",0.001836*g/cm3*pressure_of_CO2,2);
  G4Element* C = new G4Element("Carbon","C",6,12.011*g/mole);
  G4Element* O = new G4Element("Oxygen","O",8,15.99*g/mole);
  G4Element* H  = new G4Element("Hydrogen","H", 1,  1.008*g/mole);
  G4Element* D = new G4Element("Deuterium","D",1,2.014102*g/mole);
  G4Element* Cr = new G4Element("Chrome", "Cr", 25,   51.996*g/mole);
  G4Element* Fe = new G4Element("Iron"  , "Fe", 26,   55.845*g/mole);
  G4Element* Co = new G4Element("Cobalt", "Co", 27,   58.933*g/mole);
  G4Element* Ni = new G4Element("Nickel", "Ni", 28,   58.693*g/mole);
  G4Element* W  = new G4Element("Tungsten","W", 74,  183.850*g/mole);
  G4Element* Au = new G4Element("Gold","Au",79,196.97*g/mole);

  G4Element* Al  = new G4Element("Aluminium","Al", 13,  26.98*g/mole);
  G4Material* al = new G4Material("al",2.7*g/cm3,1);
  al->AddElement(Al,1);
  CO2->AddElement(C,1);
  CO2->AddElement(O,2);
  G4double Deuterium_atmosphere = 1.0;//
//  G4Material* D2 = new G4Material("D2",Deuterium_atmosphere*0.000169*g/cm3,1);
  G4Material* D2 = new G4Material("D2",Deuterium_atmosphere*0.00018*g/cm3,1);
  D2->AddElement(D,2);//D2 gas
  G4Material* Gold = new G4Material("Gold",19.32*g/cm3,1);
  Gold->AddElement(Au,1);//Gold
//
  G4Material* Havar =
    new G4Material("Havar", 8.3*g/cm3, 5);
  Havar->AddElement(Cr, 0.1785);
  Havar->AddElement(Fe, 0.1822);
  Havar->AddElement(Co, 0.4452);
  Havar->AddElement(Ni, 0.1310);
  Havar->AddElement(W , 0.0631);

  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        CO2,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
/*
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        CO2,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
*/
  //     
  // Shape 1
  //  
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_Si");
  G4ThreeVector pos1 = G4ThreeVector(-6.22*cm, -3.08*cm+10*cm, 20*cm);
  G4ThreeVector pos2 = G4ThreeVector(6.22*cm, 3.08*cm+10*cm, 20*cm);
  G4ThreeVector pos3 = G4ThreeVector(-6.22*cm, 3.08*cm+10*cm, 20*cm);
  G4ThreeVector pos4 = G4ThreeVector(6.22*cm, -3.08*cm+10*cm, 20*cm);
        
  G4Box* solidShape1 = new G4Box("Si_det",0.5*5*cm,0.5*5*cm,0.5*625*um); 
  G4LogicalVolume* logicShape1 =                         
    new G4LogicalVolume(solidShape1,         //its solid
                        shape1_mat,          //its material
                        "Si_det_log");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    logicShape1,             //its logical volume
                    "Si_det",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  new G4PVPlacement(0,                       //no rotation
                    pos2,                    //at position
                    logicShape1,             //its logical volume
                    "Si_det",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    1,                       //copy number
                    checkOverlaps);          //overlaps checking
  new G4PVPlacement(0,                       //no rotation
                    pos3,                    //at position
                    logicShape1,             //its logical volume
                    "Si_det",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    2,                       //copy number
                    checkOverlaps);          //overlaps checking
  new G4PVPlacement(0,                       //no rotation
                    pos4,                    //at position
                    logicShape1,             //its logical volume
                    "Si_det",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    3,                       //copy number
                    checkOverlaps);          //overlaps checking
  G4Material* CH2 = new G4Material("CH2",0.92*g/cm3,2);
  CH2->AddElement(C,1);
  CH2->AddElement(H,2);
  fScoringVolume=logicWorld;

 //CH2
  G4VSolid* CH2foilSolid = new G4Tubs("CH2foilSolid",0.,5.*cm,28.6*um/2.,0.,360.*deg);
  G4LogicalVolume *fCH2FoilLogical = new G4LogicalVolume(CH2foilSolid,CH2,"CH2foilLogical");
 new G4PVPlacement(0,G4ThreeVector(0.,10.*cm,555.2*mm),fCH2FoilLogical,"CH2foilPhysical",logicWorld,
                   false,0,checkOverlaps);

  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
