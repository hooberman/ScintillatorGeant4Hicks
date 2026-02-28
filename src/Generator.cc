#include "Generator.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"

#include <cmath>
#include <algorithm>

MyPrimaryGenerator::MyPrimaryGenerator() {
    fParticleGun = new G4ParticleGun(1);
}

MyPrimaryGenerator::~MyPrimaryGenerator() {
    delete fParticleGun;
}

// Return true if the ray (x0,y0,z0) + t*(ux,uy,uz), t>=0
// intersects a finite cylinder of radius R and z-range [zmin, zmax],
// with a "close" margin applied by inflating radius and z-extent.
static bool RayHitsCylinderWithMargin(
    G4double x0, G4double y0, G4double z0,
    G4double ux, G4double uy, G4double uz,
    G4double R, G4double zmin, G4double zmax,
    G4double margin
) {
    const G4double Rtol  = R + margin;
    const G4double zminT = zmin - margin;
    const G4double zmaxT = zmax + margin;

    // Helper: check point at parameter t
    auto withinCaps = [&](G4double t) -> bool {
        if (t < 0.0) return false;
        const G4double x = x0 + t*ux;
        const G4double y = y0 + t*uy;
        const G4double z = z0 + t*uz;
        const G4double r2 = x*x + y*y;
        return (r2 <= Rtol*Rtol) && (z >= zminT) && (z <= zmaxT);
    };

    // 1) Side intersection: solve (x0+ux t)^2 + (y0+uy t)^2 = Rtol^2
    const G4double a = ux*ux + uy*uy;
    const G4double b = 2.0*(x0*ux + y0*uy);
    const G4double c = (x0*x0 + y0*y0) - Rtol*Rtol;

    if (a > 0.0) {
        const G4double disc = b*b - 4.0*a*c;
        if (disc >= 0.0) {
            const G4double sdisc = std::sqrt(disc);
            const G4double t1 = (-b - sdisc) / (2.0*a);
            const G4double t2 = (-b + sdisc) / (2.0*a);

            // Check z at the side-hit parameters
            if (t1 >= 0.0) {
                const G4double z = z0 + t1*uz;
                if (z >= zminT && z <= zmaxT) return true;
            }
            if (t2 >= 0.0) {
                const G4double z = z0 + t2*uz;
                if (z >= zminT && z <= zmaxT) return true;
            }
        }
    } else {
        // Ray is parallel to cylinder axis (ux=uy=0): radial distance is constant.
        const G4double r2 = x0*x0 + y0*y0;
        if (r2 <= Rtol*Rtol) {
            // Does z(t) for t>=0 overlap [zminT, zmaxT]?
            if (uz > 0.0) {
                // z increases with t
                if (z0 <= zmaxT) {
                    const G4double t_enter = (zminT - z0) / uz;
                    const G4double t_exit  = (zmaxT - z0) / uz;
                    if (t_exit >= 0.0 && t_enter <= t_exit) return true;
                }
            } else if (uz < 0.0) {
                // z decreases with t (not expected for your upward-going muons, but safe)
                if (z0 >= zminT) {
                    const G4double t_enter = (zmaxT - z0) / uz;
                    const G4double t_exit  = (zminT - z0) / uz;
                    if (t_exit >= 0.0 && t_enter <= t_exit) return true;
                }
            } else {
                // uz==0: ray is stationary in z; then it "hits" only if z0 inside slab
                if (z0 >= zminT && z0 <= zmaxT) return true;
            }
        }
    }

    // 2) Endcap intersection: intersect with planes z=zminT and z=zmaxT
    if (uz != 0.0) {
        const G4double t_zmin = (zminT - z0) / uz;
        const G4double t_zmax = (zmaxT - z0) / uz;
        if (withinCaps(t_zmin)) return true;
        if (withinCaps(t_zmax)) return true;
    }

    return false;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent) {

    //------------------------------------------------------------
    // Detector bounding cylinder (matches your stack envelope)
    //------------------------------------------------------------
    const G4double detR = 3.75*cm;   // scintRadius = 3.75 cm
    const G4double detZ = 49.0*cm;   // totalLength = 49 cm for 25 disks (1 cm thick) + 24 gaps (1 cm)
    const G4double zmin = -0.5*detZ;
    const G4double zmax =  0.5*detZ;

    //------------------------------------------------------------
    // Step 1) fixed start plane at z = -400 mm
    //------------------------------------------------------------
    const G4double z0_plane = -400.0*mm;  // NOTE: must be inside World volume to track!

    //------------------------------------------------------------
    // Select particle (mu+)
    //------------------------------------------------------------
    G4ParticleDefinition *mu =
        G4ParticleTable::GetParticleTable()->FindParticle("mu+");

    //------------------------------------------------------------
    // Step 3) theta distribution for I(θ) ∝ cos^2 θ (per sr),
    //         including solid-angle factor: p(θ) dθ ∝ cos^2θ sinθ dθ.
    //------------------------------------------------------------
    const G4double theta_max = 0.95 * (90.0*deg); // per your requirement
    const G4double umin = std::cos(theta_max);    // u = cos(theta) in [umin, 1]

    //------------------------------------------------------------
    // Step 5) "comes close" margin to avoid edge effects
    //------------------------------------------------------------
    const G4double margin = 10.0*mm; // inflate radius and z-range slightly

    {
      const bool test = RayHitsCylinderWithMargin(
						  0.0, 0.0, z0_plane,
						  0.0, 0.0, 1.0,
						  detR, zmin, zmax,
						  margin
						  );
      G4cout << "Sanity (vertical center ray) hits? " << test << G4endl;
    }
    
    //------------------------------------------------------------
    // Step 2) choose (x,y) uniformly on a plane large enough to cover
    //         all rays up to theta_max that could hit the detector.
    //
    // Use worst-case Δz from plane to the far end of the detector.
    //------------------------------------------------------------
    const G4double dz_to_top = (zmax - z0_plane); // zmax > z0_plane, positive
    const G4double r_plane_max = (detR + margin) + dz_to_top * std::tan(theta_max);

    // We'll keep sampling until the ray passes the intersection test.
    //for (int tries = 0; tries < 10000000000000; ++tries) {
    while(true){
    
      static int printed = 0;
      if (printed < 5) {
	G4cout << "z0_plane = " << z0_plane/mm << " mm\n";
	G4cout << "detR = " << detR/mm << " mm, detZ = " << detZ/mm << " mm\n";
	G4cout << "zmin,zmax = " << zmin/mm << "," << zmax/mm << " mm\n";
	G4cout << "theta_max = " << theta_max/deg << " deg\n";
	G4cout << "dz_to_top = " << dz_to_top/mm << " mm\n";
	G4cout << "r_plane_max = " << r_plane_max/mm << " mm\n";
	printed++;
      }


        //------------------------------------------------------------
        // Step 2) x,y uniform on [-r_plane_max, +r_plane_max]
        //------------------------------------------------------------
        const G4double x0 = (2.0*G4UniformRand() - 1.0) * r_plane_max;
        const G4double y0 = (2.0*G4UniformRand() - 1.0) * r_plane_max;
        const G4double z0 = z0_plane;

        //------------------------------------------------------------
        // Step 3) sample u=cos(theta) with p(u) ∝ u^2 on [umin, 1]
        // If X is uniform in [umin^3, 1], then u = X^(1/3).
        //------------------------------------------------------------
        const G4double X     = umin*umin*umin + (1.0 - umin*umin*umin) * G4UniformRand();
        const G4double ucost = std::cbrt(X);
        const G4double theta = std::acos(ucost);

        //------------------------------------------------------------
        // Step 4) phi uniform 0..360 deg
        //------------------------------------------------------------
        const G4double phiMu = 2.0 * M_PI * G4UniformRand();

        //------------------------------------------------------------
        // Direction unit vector (upward-going; uz>0 since theta in [0, theta_max<90°])
        //------------------------------------------------------------
        const G4double ux = std::sin(theta) * std::cos(phiMu);
        const G4double uy = std::sin(theta) * std::sin(phiMu);
        const G4double uz = std::cos(theta);

        //------------------------------------------------------------
        // Step 5) keep event if the ray intersects (or comes close to) the detector
        //------------------------------------------------------------
        const bool keep = RayHitsCylinderWithMargin(
            x0, y0, z0,
            ux, uy, uz,
            detR, zmin, zmax,
            margin
        );

        if (!keep) {
            // Step 6) reject and go back to Step 2
            continue;
        }

        //------------------------------------------------------------
        // Step 6) fire the muon
        //------------------------------------------------------------
        fParticleGun->SetParticleDefinition(mu);
        fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux, uy, uz));
        fParticleGun->SetParticleMomentum(100.0 * GeV);

        fParticleGun->GeneratePrimaryVertex(anEvent);
        break;
    }

    //G4Exception("MyPrimaryGenerator::GeneratePrimaries", "GEN001", FatalException,
    //		"Failed to find an accepted muon after 1e6 tries; check r_plane_max/theta_max.");
    
}
