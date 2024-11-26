#include "run.hh"
#include <sstream>

MyRunAction::MyRunAction()
{}

MyRunAction::~MyRunAction()
{}

void MyRunAction::BeginOfRunAction(const G4Run*)
{
G4AnalysisManager *man = G4AnalysisManager::Instance();

man->OpenFile("100GeV_Pion_Run_100000_Events_test.root");

man->CreateNtuple("Layer_1","Layer_1");
for(G4int k =0 ;k<192;k++){
    std::stringstream ss;
    ss << "Channel_" << k;
    man->CreateNtupleDColumn(ss.str());
      }
man->FinishNtuple(0);


man->CreateNtuple("Layer_2","Layer_2");
for(G4int k =0 ;k<192;k++){
        std::stringstream ss;
        ss << "Channel_" << k;
        man->CreateNtupleDColumn(ss.str());
}
man->FinishNtuple(1);
}

void MyRunAction::EndOfRunAction(const G4Run*)
{
G4AnalysisManager *man = G4AnalysisManager::Instance();

man->Write();

man->CloseFile();

}

