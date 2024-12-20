#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "detector.hh"
#include "G4Polyhedra.hh"

#include "G4VisAttributes.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "G4ThreeVector.hh" 
#include "G4RotationMatrix.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
MyDetectorConstruction();
~MyDetectorConstruction();

virtual G4VPhysicalVolume *Construct();

private:
G4LogicalVolume *logicDetector;
virtual void ConstructSDandField();

};
#endif
