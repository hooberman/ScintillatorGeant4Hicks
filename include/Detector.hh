#ifndef DETECTOR_HH
#define DETECTOR_HH

#include <fstream>
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
        static std::ofstream outFile;
private:
        G4double totalEnergyDepositDetector;
        G4double totalLightYieldDetector;

};

#endif
