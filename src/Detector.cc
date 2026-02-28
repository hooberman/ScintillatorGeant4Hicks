#include "Detector.hh"
std::ofstream MySensitiveDetector::outFile;

MySensitiveDetector::MySensitiveDetector(G4String name)
  : G4VSensitiveDetector(name), totalEnergyDepositDetector(0.0), totalLightYieldDetector(0.0) {

}

MySensitiveDetector::~MySensitiveDetector() {}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory * )
{
    G4Track* track = aStep->GetTrack();

    // Only optical photons
    if(track->GetParticleDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
        return false;

    // Kill after detection
    track->SetTrackStatus(fStopAndKill);

    // Photon energy (eV)
    G4double energy_eV = aStep->GetTrack()->GetKineticEnergy() / eV;

    // Momentum direction
    G4ThreeVector mom = track->GetMomentumDirection();

    // Identify which SiPM (via copy number)
    const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();
    G4int copyNo = touch->GetCopyNumber();

    // Tile index (1–25), SiPM index (1–4)
    G4int tile = (copyNo / 4) + 1;
    G4int sipm = (copyNo % 4) + 1;

    /*
    // Surface normal depending on SiPM
    G4ThreeVector normal;
    switch (sipm-1) {
        case 0: normal = G4ThreeVector( 1, 0, 0 ); break; // +x
        case 1: normal = G4ThreeVector( 0, 1, 0 ); break; // +y
        case 2: normal = G4ThreeVector(-1, 0, 0 ); break; // -x
        case 3: normal = G4ThreeVector( 0,-1, 0 ); break; // -y
    }

    // Incidence angle (degrees)
    //double angle_deg =
    //  std::acos( mom.unit().dot(normal) ) * 180.0 / CLHEP::pi;

    // Surface normal depending on SiPM
    
    // Incidence angle (degrees)
    // Proper incidence angle in [0, 90] degrees
    double cosTheta = mom.unit().dot(normal);
    cosTheta = std::fabs(cosTheta);            // make it 0..1
    if (cosTheta > 1.0) cosTheta = 1.0;        // numeric safety

    double angle_deg = std::acos(cosTheta) * 180.0 / CLHEP::pi;
    */
    
    // Surface normal for cylindrical side: radial direction in x–y plane
    G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector radial(pos.x(), pos.y(), 0.0);

    // Guard against pathological case at the exact center
    G4ThreeVector normal = radial.unit();
    if (radial.mag2() == 0.0) {
        // fallback (should basically never happen)
        normal = G4ThreeVector(1.0, 0.0, 0.0);
    }

    // Incidence angle (degrees), folded into [0, 90]
    double cosTheta = mom.unit().dot(normal);
    cosTheta = std::fabs(cosTheta);          // 0..1
    if (cosTheta > 1.0) cosTheta = 1.0;      // numeric safety

    double angle_deg = std::acos(cosTheta) * 180.0 / CLHEP::pi;












    
    //double angle_deg = std::acos( mom.unit().dot(normal) ) * 180.0 / CLHEP::pi;

 
    // // Incidence angle in degrees, always in [0, 90]
    // G4double cosInc = mom.unit().dot(normal);

    // // Take magnitude – we don't care about which side of the surface
    // cosInc = std::fabs(cosInc);

    // // Clamp to valid range to avoid numerical issues
    // if (cosInc > 1.0) cosInc = 1.0;
    // if (cosInc < -1.0) cosInc = -1.0;

    // G4double angle_deg = std::acos(cosInc) * 180.0 / CLHEP::pi;

    
    // --- Write to text file ---
    //static std::ofstream outfile("sipm_hits.txt", std::ios::app);
    outFile
        << std::setprecision(3)   // 4 significant figures
        << tile  << "     "   // 5 spaces
        << sipm  << "     "
        << energy_eV << "     "
        << angle_deg
        << "\n";       // 5 blank lines after each hit

    return true;
}





G4double MySensitiveDetector::GetTotalDepositedEnergyDetector() 
{
    // G4cout << "GET TOTAL ENERGY METHOD WAS CALLED: "<< totalEnergyDepositDetector << G4endl;
    return totalEnergyDepositDetector;
}

void MySensitiveDetector::ResetTotalDepositedEnergyDetector() 
{
    totalEnergyDepositDetector = 0.0;
}

G4double MySensitiveDetector::GetTotalDepositedLightYieldDetector() 
{
    // G4cout << "GET TOTAL LIGHT YIELD METHOD WAS CALLED: "<< totalLightYieldDetector << G4endl;
    return totalLightYieldDetector;
}

//G4int MySensitiveDetector::GetNPhotonsDet0(){
//  return nPhotonsDet0;
//}

void MySensitiveDetector::ResetTotalDepositedLightYieldDetector() 
{
    totalLightYieldDetector = 0.0;
}

