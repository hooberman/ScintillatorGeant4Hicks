#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"

#include "G4AutoLock.hh"
#include "G4ThreeVector.hh"

#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

// Thread-safe file writing (important if you ever enable MT)
namespace {
  G4Mutex gHitFileMutex = G4MUTEX_INITIALIZER;

  double WrapPhi0To2Pi(double phi) {
    const double twopi = 2.0 * pi;
    while (phi < 0.0) phi += twopi;
    while (phi >= twopi) phi -= twopi;
    return phi;
  }

  class MuonKinematicsSD : public G4VSensitiveDetector {
   public:
    explicit MuonKinematicsSD(const G4String& name, const G4String& outFile)
      : G4VSensitiveDetector(name),
        m_outFile(outFile) {}

    ~MuonKinematicsSD() override = default;

    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override {
      if (!step) return false;

      G4Track* track = step->GetTrack();
      if (!track) return false;

      // Only record muons
      const int pdg = track->GetDefinition()->GetPDGEncoding();
      if (pdg != 13 && pdg != -13) return false;

      // Record only on entry to the sensitive volume (avoid multiple steps)
      const auto pre = step->GetPreStepPoint();
      if (!pre) return false;
      if (pre->GetStepStatus() != fGeomBoundary) return false;

      const G4ThreeVector dir = track->GetMomentumDirection();
      const double cz = std::clamp(dir.z(), -1.0, 1.0);

      const double theta = std::acos(cz);                  // radians
      double phi = std::atan2(dir.y(), dir.x());           // radians
      phi = WrapPhi0To2Pi(phi);

      const double E_GeV = track->GetTotalEnergy() / GeV;  // total energy

      const int eventID = G4RunManager::GetRunManager()
                            ->GetCurrentEvent()
                            ->GetEventID();
      const int trackID = track->GetTrackID();

      {
        G4AutoLock lock(&gHitFileMutex);
        std::ofstream out(m_outFile, std::ios::out | std::ios::app);
        out << std::fixed << std::setprecision(10)
            << eventID << " "
            << trackID << " "
            << (theta / deg) << " "
            << (phi / deg) << " "
            << E_GeV
            << "\n";
      }

      return true;
    }

   private:
    G4String m_outFile;
  };
}

MyDetectorConstruction::MyDetectorConstruction() {}

MyDetectorConstruction::~MyDetectorConstruction() {}

G4VPhysicalVolume* MyDetectorConstruction::Construct() {

  // ------------------------------------------------------------
  // Materials (no optical properties needed)
  // ------------------------------------------------------------
  auto nist = G4NistManager::Instance();

  // World material: air is fine, or use vacuum if you want
  G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");

  // Detector material: arbitrary solid (plastic scintillator-like) is fine
  // since we are not simulating optical photons.
  G4Material* detMat = nist->FindOrBuildMaterial("G4_POLYSTYRENE");

  // ------------------------------------------------------------
  // World volume
  // ------------------------------------------------------------
  G4Box* solidWorld = new G4Box("solidWorld", 10*m, 10*m, 10*m);
  G4LogicalVolume* logicWorld =
      new G4LogicalVolume(solidWorld, worldMat, "logicWorld");

  G4VPhysicalVolume* physWorld =
      new G4PVPlacement(
          nullptr,
          G4ThreeVector(),
          logicWorld,
          "physWorld",
          nullptr,
          false,
          0,
          true
      );

  // ------------------------------------------------------------
  // Simplified detector: single cylinder
  // Radius = 5 inches, height = 1 meter
  // ------------------------------------------------------------
  const G4double detRadius = 6.0 * cm;     // 12.7 cm
  const G4double detHalfZ  = 0.5 * m;      // height 1 m

  G4Tubs* solidDet =
      new G4Tubs("solidDetector", 0.0, detRadius, detHalfZ, 0.0, 360.0*deg);

  G4LogicalVolume* logicDet =
      new G4LogicalVolume(solidDet, detMat, "logicDetector");

  new G4PVPlacement(
      nullptr,
      G4ThreeVector(0, 0, 0),
      logicDet,
      "physDetector",
      logicWorld,
      false,
      0,
      true
  );

  // ------------------------------------------------------------
  // Sensitive detector: record muon theta, phi, energy on entry
  // ------------------------------------------------------------
  auto sdManager = G4SDManager::GetSDMpointer();

  const G4String sdName = "MuonKinematicsSD";
  const G4String outFile = "muon_detector_hits_theta_phi_E.txt";

  if (!sdManager->FindSensitiveDetector(sdName, false)) {
    auto sd = new MuonKinematicsSD(sdName, outFile);
    sdManager->AddNewDetector(sd);
    logicDet->SetSensitiveDetector(sd);
  } else {
    logicDet->SetSensitiveDetector(sdManager->FindSensitiveDetector(sdName));
  }

  return physWorld;
}
