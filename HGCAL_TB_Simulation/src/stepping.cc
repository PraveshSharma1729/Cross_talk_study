#include "stepping.hh"

MySteppingAction::MySteppingAction(MyEventAction *eventAction)
{

fEventAction = eventAction;
}

MySteppingAction::~MySteppingAction()
{}

void MySteppingAction::UserSteppingAction(const G4Step *step)
{
G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

const MyDetectorConstruction *detectorConstruction = static_cast<const MyDetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

G4double edep = step->GetTotalEnergyDeposit();
 // Get the copy number of the volume
    G4int copyNo = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();

for (G4int i=1; i < 3 ;i++){
    // Add the energy deposition to the event action
        fEventAction->AddEdep(copyNo, edep);
    }
}
