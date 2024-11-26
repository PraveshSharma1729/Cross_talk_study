#include "event.hh"
#include "G4SystemOfUnits.hh"
#include <chrono>
#include "Randomize.hh"

MyEventAction::MyEventAction(MyRunAction*)
{
	for (G4int i = 0; i < 192*2; ++i) {
        fEdep[i] = 0.0;
    }
	
}

MyEventAction::~MyEventAction()
{}

void MyEventAction::BeginOfEventAction(const G4Event*)
{
	// Set a new random seed based on the current system time
  	long seeds[2];
  	auto now = std::chrono::high_resolution_clock::now().time_since_epoch();
  	seeds[0] = std::chrono::duration_cast<std::chrono::seconds>(now).count();
  	seeds[1] = std::chrono::duration_cast<std::chrono::nanoseconds>(now).count();
  	CLHEP::HepRandom::setTheSeeds(seeds);
	
	for (G4int i = 0; i < 192*2; ++i) {
        fEdep[i] = 0.0;
        }	
}

void MyEventAction::EndOfEventAction(const G4Event*)
{
G4AnalysisManager *man = G4AnalysisManager::Instance();

  for (G4int i = 0; i <192; ++i) {
       man->FillNtupleDColumn( 0,i, fEdep[i]);
        }
    man->AddNtupleRow(0);

  for (G4int i = 0; i <192; ++i) {
       man->FillNtupleDColumn( 1,i, fEdep[192+i]);
        }
    man->AddNtupleRow(1);
}

void MyEventAction::AddEdep(G4int copyNo, G4double edep) {
        if (copyNo<384){
        fEdep[copyNo] += edep; // Store energy deposition in corresponding array element
        }
    }

