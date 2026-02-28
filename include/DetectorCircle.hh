#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
//includes for CSV
#include "g4csv_defs.hh"
#include "G4CsvAnalysisManager.hh"


class MySensitiveDetector : public G4VSensitiveDetector {
    public:
        MySensitiveDetector(G4String);
        ~MySensitiveDetector();
        
        virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);

        void ResetTotalDepositedEnergyDetector();
        G4double GetTotalDepositedEnergyDetector();

        void ResetTotalDepositedLightYieldDetector();
        G4double GetTotalDepositedLightYieldDetector();
        G4int GetNPhotonsDet0();
private:
        G4double totalEnergyDepositDetector;
        G4double totalLightYieldDetector;
        G4int nPhotonsDet1;
        G4int nPhotonsDet2;
        G4int nPhotonsDet3;
        G4int nPhotonsDet4;

};

#endif
