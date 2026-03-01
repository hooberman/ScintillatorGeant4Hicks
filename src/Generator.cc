#include "Generator.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4Event.hh"
#include "G4ios.hh"

#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>

// ------------------------------------------------------------
// Your existing ray–cylinder intersection with margin
// ------------------------------------------------------------
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

    // If ray is almost parallel to cylinder axis, handle separately
    const G4double a = ux*ux + uy*uy;

    // Intersect with infinite cylinder x^2 + y^2 = Rtol^2
    if (a > 1e-20) {
        const G4double b = 2.0*(x0*ux + y0*uy);
        const G4double c = x0*x0 + y0*y0 - Rtol*Rtol;
        const G4double disc = b*b - 4.0*a*c;
        if (disc >= 0.0) {
            const G4double sdisc = std::sqrt(disc);
            const G4double t1 = (-b - sdisc) / (2.0*a);
            const G4double t2 = (-b + sdisc) / (2.0*a);

            // Check z-range at intersections
            const G4double z1 = z0 + t1*uz;
            const G4double z2 = z0 + t2*uz;

            if (t1 >= 0.0 && z1 >= zminT && z1 <= zmaxT) return true;
            if (t2 >= 0.0 && z2 >= zminT && z2 <= zmaxT) return true;
        }
    } else {
        // ux ~ uy ~ 0, ray parallel to z axis: check if within radius already
        const G4double r2 = x0*x0 + y0*y0;
        if (r2 <= Rtol*Rtol) {
            // Then it hits if it passes through z-range at any t>=0
            // Compute t for entering zmin/zmax
            if (std::abs(uz) < 1e-20) return false; // degenerate
            const G4double t_zmin = (zminT - z0)/uz;
            const G4double t_zmax = (zmaxT - z0)/uz;
            if (withinCaps(t_zmin) || withinCaps(t_zmax)) return true;
        }
        return false;
    }

    // Also check intersection with the endcaps (planes z=zminT and z=zmaxT)
    if (std::abs(uz) > 1e-20) {
        const G4double t_zmin = (zminT - z0)/uz;
        const G4double t_zmax = (zmaxT - z0)/uz;
        if (withinCaps(t_zmin)) return true;
        if (withinCaps(t_zmax)) return true;
    }

    return false;
}

// ------------------------------------------------------------
// Gaisser sea-level differential muon flux (shape only)
// d^2N/(dE dΩ) ∝ E^{-2.7} * [ 1/(1+1.1 E cosθ / 115) + 0.054/(1+1.1 E cosθ / 850) ]
// Units of the standard prefactor are not needed for sampling (we sample shape).
// ------------------------------------------------------------
static double Gaisser_d2NdEdOmega(double E_GeV, double cosTheta) {
    // Guard against weird cosTheta
    cosTheta = std::clamp(cosTheta, 0.0, 1.0);

    // Avoid E<=0
    if (E_GeV <= 0.0) return 0.0;

    const double E = E_GeV;
    const double gamma = 2.7;

    const double term1 = 1.0 / (1.0 + 1.1 * E * cosTheta / 115.0);
    const double term2 = 0.054 / (1.0 + 1.1 * E * cosTheta / 850.0);

    return std::pow(E, -gamma) * (term1 + term2);
}

// ------------------------------------------------------------
// Sample correlated (theta, E) from f(theta,E) ∝ (d^2N/dE dΩ)*sinθ
// with theta ∈ [0, theta_max], E ∈ [Emin, Emax], and flat phi.
// Implementation: rejection sampling in (u=cosθ, y=lnE) with a precomputed fmax.
// Proposal:
//   u uniform in [umin, 1]
//   y uniform in [lnEmin, lnEmax] => E = exp(y) (log-uniform in energy)
// Target weight (up to constant factors):
//   w(u,E) ∝ (d^2N/dE dΩ)(E,u) * sinθ * Jacobians
// Since:
//   theta from u: dθ = -du/sinθ
//   and y=lnE: dE = E dy
// With proposal uniform in u and y, the acceptance ratio can use:
//   w(u,E) ∝ (d^2N/dE dΩ)(E,u) * sinθ * (dΩ factor already in sinθ dθ dφ)
// but sampling variables are (u,y). The correct unnormalized density in (u,y) is:
//   f(u,y) ∝ (d^2N/dE dΩ)(E,u) * sinθ * |dθ/du| * |dE/dy|
//          ∝ (d^2N/dE dΩ)(E,u) * sinθ * (1/sinθ) * E
//          ∝ (d^2N/dE dΩ)(E,u) * E
// Nice cancellation: no explicit sinθ needed in the acceptance when sampling uniform u.
// That’s the point of choosing u=cosθ.
// ------------------------------------------------------------
struct ThetaEnergySampler {
    double theta_max_rad;
    double umin;
    double lnEmin;
    double lnEmax;
    double fmax; // maximum of f(u,y) ∝ (d2N/dEdΩ)*E over the domain

    explicit ThetaEnergySampler(double theta_max, double Emin_GeV, double Emax_GeV)
        : theta_max_rad(theta_max),
          umin(std::cos(theta_max)),
          lnEmin(std::log(Emin_GeV)),
          lnEmax(std::log(Emax_GeV)),
          fmax(0.0)
    {
        // Precompute a conservative fmax on a grid (cheap, done once).
        // f(u,E) ∝ Gaisser(E,u)*E
        const int Nu = 200;
        const int Ny = 200;

        double maxv = 0.0;
        for (int iu = 0; iu < Nu; ++iu) {
            const double fu = (iu + 0.5) / Nu;
            const double u = umin + (1.0 - umin) * fu;

            for (int iy = 0; iy < Ny; ++iy) {
                const double fy = (iy + 0.5) / Ny;
                const double y = lnEmin + (lnEmax - lnEmin) * fy;
                const double E = std::exp(y);

                //const double v = Gaisser_d2NdEdOmega(E, u) * E;
		const double v = Gaisser_d2NdEdOmega(E, u) * E * u; // extra cos(theta) for plane-crossing flux
                if (v > maxv) maxv = v;
            }
        }

        // Inflate a bit to be safe against grid miss
        fmax = 1.2 * maxv;

        // If something went wrong, avoid division by zero later
        if (!(fmax > 0.0)) {
            fmax = 1.0;
        }
    }

    void Sample(double& theta_out, double& E_GeV_out) const {
        while (true) {
            const double u = umin + (1.0 - umin) * G4UniformRand();        // uniform in cosθ
            const double y = lnEmin + (lnEmax - lnEmin) * G4UniformRand(); // uniform in lnE
            const double E = std::exp(y);

            //const double f = Gaisser_d2NdEdOmega(E, u) * E; // target in (u, lnE)
	    const double f = Gaisser_d2NdEdOmega(E, u) * E * u;
 
            const double accept = f / fmax;
            if (G4UniformRand() < accept) {
                theta_out = std::acos(std::clamp(u, 0.0, 1.0));
                E_GeV_out = E;
                return;
            }
        }
    }
};

// ------------------------------------------------------------
// File logging helpers (append mode).
// Writes: theta(deg) energy(GeV)
// ------------------------------------------------------------
static void AppendThetaE(const char* fname, double theta_rad, double E_GeV) {
    static const int prec = 10;
    std::ofstream out(fname, std::ios::out | std::ios::app);
    out << std::fixed << std::setprecision(prec)
        << (theta_rad / deg) << " " << E_GeV << "\n";
}

// ------------------------------------------------------------
// Primary generator class
// ------------------------------------------------------------
MyPrimaryGenerator::MyPrimaryGenerator() {
    fParticleGun = new G4ParticleGun(1);
}

MyPrimaryGenerator::~MyPrimaryGenerator() {
    delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event* anEvent) {

    //------------------------------------------------------------
    // Detector bounding cylinder (matches your stack envelope)
    //------------------------------------------------------------
    const G4double detR = 3.75*cm;   // scintRadius = 3.75 cm
    const G4double detZ = 49.0*cm;   // totalLength = 49 cm
    const G4double zmin = -0.5*detZ;
    const G4double zmax =  0.5*detZ;

    //------------------------------------------------------------
    // Fixed start plane at z = -400 mm (keep your current geometry)
    //------------------------------------------------------------
    const G4double z0_plane = -400.0*mm;  // must be inside World volume

    //------------------------------------------------------------
    // "Comes close" margin to avoid edge effects
    //------------------------------------------------------------
    const G4double margin = 10.0*mm;

    //------------------------------------------------------------
    // Muon particle (mu+)
    //------------------------------------------------------------
    G4ParticleDefinition* mu =
        G4ParticleTable::GetParticleTable()->FindParticle("mu+");

    //------------------------------------------------------------
    // Domain for correlated (theta, E) sampling
    //------------------------------------------------------------
    const G4double theta_max = 0.95 * (90.0*deg); // keep your theta cap
    const double Emin_GeV = 1.0;    // adjust if you want (e.g., 0.5, 2.0, etc.)
    const double Emax_GeV = 1000.0; // adjust if you want

    // Build sampler once (static) with your chosen ranges
    static ThetaEnergySampler sampler(theta_max, Emin_GeV, Emax_GeV);

    //------------------------------------------------------------
    // Plane size: large enough to cover all rays up to theta_max
    //------------------------------------------------------------
    const G4double dz_to_top = (zmax - z0_plane);
    const G4double r_plane_max = (detR + margin) + dz_to_top * std::tan(theta_max);

    //------------------------------------------------------------
    // Keep sampling until ray passes intersection test
    //------------------------------------------------------------
    while (true) {

        // Step: x,y uniform on plane square [-r_plane_max, +r_plane_max]
        const G4double x0 = (2.0*G4UniformRand() - 1.0) * r_plane_max;
        const G4double y0 = (2.0*G4UniformRand() - 1.0) * r_plane_max;
        const G4double z0 = z0_plane;

        // Sample correlated theta and E from Gaisser parameterization
        double theta = 0.0;
        double E_GeV = 0.0;
        sampler.Sample(theta, E_GeV);

        // phi uniform 0..2π
        const G4double phiMu = 2.0 * M_PI * G4UniformRand();

        // Direction unit vector (keep your upward-going convention)
        const G4double ux = std::sin(theta) * std::cos(phiMu);
        const G4double uy = std::sin(theta) * std::sin(phiMu);
        const G4double uz = std::cos(theta);

        // Log every generated (theta,E) trial, regardless of keep
        AppendThetaE("muon_all_theta_E.txt", theta, E_GeV);

        // Keep event if ray intersects (or comes close to) detector cylinder
        const bool keep = RayHitsCylinderWithMargin(
            x0, y0, z0,
            ux, uy, uz,
            detR, zmin, zmax,
            margin
        );

        if (!keep) {
            continue;
        }

        // Log only kept events
        AppendThetaE("muon_kept_theta_E.txt", theta, E_GeV);

        // Fire muon with the sampled energy
        fParticleGun->SetParticleDefinition(mu);
        fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux, uy, uz));

        // Interpret E_GeV as TOTAL energy in GeV (as in Gaisser formula).
        // Geant4 SetParticleEnergy expects kinetic energy.
        const G4double m = mu->GetPDGMass();           // energy units (MeV by default internally but with units it's fine)
        const G4double Etot = E_GeV * GeV;
        G4double Ekin = Etot - m;
        if (Ekin < 0.0) Ekin = 0.0;

        fParticleGun->SetParticleEnergy(Ekin);

        fParticleGun->GeneratePrimaryVertex(anEvent);
        break;
    }
}
