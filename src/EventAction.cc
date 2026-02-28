#include "EventAction.hh"
#include <vector>
#include <algorithm>
#include <cmath>

extern int gNumberOfEvents;
extern std::string gLabel;

MyEventAction::MyEventAction(MySteppingAction* steppingAction, MySensitiveDetector* sensitiveDetector)
    : fSteppingAction(steppingAction), fSensitiveDetector(sensitiveDetector)
{}

MyEventAction::~MyEventAction() {}

// Called at the beginning of each event
void MyEventAction::BeginOfEventAction(const G4Event* event)
{
    // Reset the total optical photon energy at the beginning of each event
    
    fSteppingAction->ResetTotalOpticalPhotonEnergy();
    fSteppingAction->ResetTotalLightYield();
    fSteppingAction->ResetTotalDepositedEnergy();

    fSensitiveDetector->ResetTotalDepositedEnergyDetector();
    fSensitiveDetector->ResetTotalDepositedLightYieldDetector();

    // Open (and overwrite) the hits file for this event
    //MySensitiveDetector::outFile.open("output/sipm_hits.txt",
    //                                  std::ios::out | std::ios::app);

    // std::string filename;
    
    // if (gNumberOfEvents > 0) {
    //   filename = "output/sipm_hits_" + std::to_string(gNumberOfEvents) + "events.txt";
    // } else {
    //   filename = "output/sipm_hits.txt";
    // }

    std::string filename = "output/sipm_hits";

    if (!gLabel.empty()) {
      filename += "_" + gLabel;
    }
    
    if (gNumberOfEvents > 0) {
      filename += "_" + std::to_string(gNumberOfEvents) + "events";
    }

    filename += ".txt";

    MySensitiveDetector::outFile.open(filename, std::ios::out | std::ios::app);

    // Compute and write muon info for this event
    G4double mu_polar_deg   = 0.0;
    G4double mu_azimuth_deg = 0.0;
    G4double mu_z_mid       = 0.0;

    ComputeMuonInfo(event, mu_polar_deg, mu_azimuth_deg, mu_z_mid);

    //if( mu_z_mid > -39.0 ){
    MySensitiveDetector::outFile << "MuonPolarAngle[deg]  MuonPhiAngle[deg]   MuonZAtMidPoint[cm]\n";
    
    // First line of the file for this event:
    // polar angle [deg], azimuthal angle [deg], z midpoint [cm]
    MySensitiveDetector::outFile
      << mu_polar_deg   << " "
      << mu_azimuth_deg << " "
      << (mu_z_mid / cm) << "\n\n";
    
    MySensitiveDetector::outFile << "Tile  SiPM  PhotonEnergy[eV]  PhotonIncidenceAngle[degrees]\n";

    //}
}

// Called at the end of each event
void MyEventAction::EndOfEventAction(const G4Event* event)
{
    // ::::::::::::::::::::::: Values from Stepping Action :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    // Get and print the total optical photon energy accumulated during the event
    G4double totalEnergy = fSteppingAction->GetTotalOpticalPhotonEnergy();
    G4cout << "Total scintillation photon energy for this event: "
        << totalEnergy / MeV << " MeV" << G4endl; // divide or multiply by MeV?
     
    // Get and print the total light yield during the event
    G4double totalYield = fSteppingAction->GetTotalLightYield();
    
    G4cout << "Total light yield for this event: "
            << totalYield << " scintillation photons were emitted." <<G4endl;
    
    // Check results with Birk's Law Prediction:
    G4double totalDeposited = fSteppingAction->GetTotalDepositedEnergy();

    // Convert the total light yield to the light yield per unit energy deposited:
     G4cout << "Light yield per MeV for this event: "
            << totalYield/totalDeposited << G4endl;

    // Get and print total muon deposited energy
     G4cout << "Total energy deposited by the muon for this event: " 
           << totalDeposited / MeV << "MeV" << G4endl;


     //G4cout << "# PHOTONS STRIKING SIPM0 " << fSensitiveDetector->GetNPhotonsDet0() << G4endl;
     
    // ::::::::::::::::::::::: Values from Sensitive Detector :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    G4double sipmTotalEnergy = fSensitiveDetector->GetTotalDepositedEnergyDetector();
    G4double sipmTotalYield = fSensitiveDetector->GetTotalDepositedLightYieldDetector();

    //G4cout << "TOTAL PHOTONS STRIKING SIPMS       " << sipmTotalYield << G4endl;
    //G4cout << "TOTAL PHOTON ENERGY STRIKING SIPMS " << sipmTotalEnergy << G4endl;

    // Add 5 blank lines after this event and close file
    MySensitiveDetector::outFile << "\n\n\n\n\n";
    MySensitiveDetector::outFile.close();

    
    // G4cout << "ENERGY: " << sipmTotalEnergy << "MeV" << G4endl;

    /* This was How I filled the ROOT output.root histograms ::::::::::::::::::::::::::
    // Store in ROOT Histogram:
 
    G4AnalysisManager* manager = G4AnalysisManager::Instance();

    manager->FillH1(0, sipmTotalYield);
    manager->FillH1(1, sipmTotalEnergy);
    manager->FillH1(2, totalEnergy);
    // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    */ 

    // ::::::::::::::::::::::: Store Deposited Energy in ROOT Ntuple :::::::::::::::::::::::::::::::::::::::::::::::::::

/*
    G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

    G4AnalysisManager *manager = G4AnalysisManager::Instance();
    manager->FillNtupleIColumn(0, evt);
    manager->FillNtupleDColumn(1, totalEnergyDepositedInDetector);
    manager->FillNtupleDColumn(2, totalYieldDepositedInDetector);
    manager->AddNtupleRow(0);
*/


    // I may have learned why geant4 doesn't implement SiPMs. Scintillation photons do not deposit energy in sensitive detectors
    // in the same way that charged particles do. This is because it doesn't interact with the material via excitation or ionization (?)
    // which is the cause of energy deposition.
    // I may need to track the number of scintillation photons (light yield) and the energy of each particle. 

    // So instead of implementing a sensitive detector which is similar to a SiPM, I need to know the amount of scintillation photons 
    // which hit the detector, the energies of these photons (use kinetic energy because massless particle means KE = E),
    // and the quantum efficiency of the SiPM. 

}


void MyEventAction::ComputeMuonInfo(const G4Event* event,
                                    G4double& polar_deg,
                                    G4double& azimuth_deg,
                                    G4double& z_mid)
{
    // Get primary muon’s vertex and direction
    const G4PrimaryVertex* vtx = event->GetPrimaryVertex(0);
    const G4PrimaryParticle* primary = vtx->GetPrimary(0);

    G4ThreeVector pos0 = vtx->GetPosition();                  // starting position
    G4ThreeVector dir  = primary->GetMomentumDirection().unit(); // direction (unit)

    // 1) Polar angle wrt +z (top/bottom faces normal)
    polar_deg = std::acos(dir.z()) * 180.0 / CLHEP::pi;

    // 2) Azimuthal angle in x–y plane
    azimuth_deg = std::atan2(dir.y(), dir.x()) * 180.0 / CLHEP::pi;
    if (azimuth_deg < 0.0) azimuth_deg += 360.0;

    // 3) z midpoint between entry and exit through bounding cylinder
    //    Cylinder is aligned with z, radius R, total length L, centered at z=0.
    const G4double R      = 3.75*cm;   // scintillator disk radius
    const G4double L      = 49.0*cm;   // total detector length (25 disks + gaps)
    const G4double halfL  = 0.5*L;

    const G4double ux = dir.x();
    const G4double uy = dir.y();
    const G4double uz = dir.z();
    const G4double px = pos0.x();
    const G4double py = pos0.y();
    const G4double pz = pos0.z();

    std::vector<G4double> ts;

    // Intersections with cylindrical side (x^2 + y^2 = R^2)
    G4double A = ux*ux + uy*uy;
    G4double B = 2.0*(px*ux + py*uy);
    G4double C = px*px + py*py - R*R;

    if (A > 0.0) {
        G4double disc = B*B - 4.0*A*C;
        if (disc >= 0.0) {
            G4double sqrtDisc = std::sqrt(disc);
            G4double t1 = (-B - sqrtDisc)/(2.0*A);
            G4double t2 = (-B + sqrtDisc)/(2.0*A);

            for (G4double t : {t1, t2}) {
                G4double z = pz + t*uz;
                if (z >= -halfL && z <= halfL)
                    ts.push_back(t);
            }
        }
    }

    // Intersections with top/bottom planes z = ±halfL
    if (std::fabs(uz) > 0.0) {
        for (G4double zPlane : { -halfL, halfL }) {
            G4double t = (zPlane - pz)/uz;
            G4double x = px + t*ux;
            G4double y = py + t*uy;
            if (x*x + y*y <= R*R)
                ts.push_back(t);
        }
    }

    if (ts.size() >= 2) {
        std::sort(ts.begin(), ts.end());
        G4double tEnter = ts.front();
        G4double tExit  = ts.back();

        G4double zEnter = pz + tEnter*uz;
        G4double zExit  = pz + tExit *uz;

        z_mid = 0.5*(zEnter + zExit);
    } else {
        // Fallback: if something goes weird, just use starting z
        z_mid = pz;
    }
}
