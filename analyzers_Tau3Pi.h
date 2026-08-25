///////////////////////////////////////////////////////////////////////////////////
///// Future Circular Collider Software - Tau-lepton Mass Peak Reconstruction /////
/////          Helper file - C++ custom functions for the analysis            /////
/////         Base file provided by Willy Weber and Aurelien Martens          /////
///////////////////////////////////////////////////////////////////////////////////

#ifndef FCCANA_ADDITIONAL_ANALYZERS_H
#define FCCANA_ADDITIONAL_ANALYZERS_H

// C++ standard library
#include <cmath>
#include <random>
#include <chrono>
#include <utility>
#include <array>

// ROOT
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "ROOT/RVec.hxx"

// EDM4hep
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/VertexData.h"
#include "edm4hep/MCParticleData.h"

// FCCAnalyses
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/ReconstructedParticle2Track.h"
#include "FCCAnalyses/VertexFitterSimple.h"

namespace FCCAnalyses :: VertexAnalysis {
    const double PION_MASS = 0.13957;  // GeV
    const double Z_BOSON_MASS = 91.188;  // GeV
    const double E_TAU = Z_BOSON_MASS/2; // GeV, NOMINAL ENERGY!
    const double TAU_MASS = 1.77693;
    const double c = 299792458; // m/s

    ///// Computes the mass of the tau at the MC level /////  
    std::pair<double,double> TauMassMC(const double Et, const double M,
                                       const double rx, const double ry, const double rz,
                                       const double px, const double py, const double pz) {
        double p3pi = std::sqrt(px*px + py*py + pz*pz);
        double M2 = M*M;
        double E3pi = std::sqrt(p3pi*p3pi + M2);
        double r = std::sqrt(rx*rx + ry*ry + rz*rz);

        double cb = (rx*px + ry*py + rz*pz)/(p3pi*r);
        cb = std::clamp(cb, -1.0, 1.0);

        double A = 2*E3pi*Et;
        double B = M2;
        double C = 2*p3pi*p3pi*cb*cb;
        double D = 2*p3pi*cb*std::sqrt(C/2 + Et*Et - A + B);

        return {std::sqrt(A - B - C - D), std::acos(cb)}; //returns both the mass and the angle (beta)
    }

    ///// Computes the reconstructed tau mass /////  
    std::pair<double,double> TauMass(const VertexingUtils::FCCAnalysesVertex& vertex,
                                     const float vx, const float vy, const float vz,
                                     const float BS_x, const float BS_y, const float BS_z) {
        TLorentzVector tau;
        ROOT::VecOps::RVec<TVector3> momenta = vertex.updated_track_momentum_at_vertex;
        int n = momenta.size();
        if (n!=3)
          cout << n << endl;
        for (int ileg=0; ileg < n; ileg++) {
          TVector3 track_momentum = momenta[ileg];
          TLorentzVector leg;
          leg.SetXYZM(track_momentum[0], track_momentum[1], track_momentum[2], PION_MASS) ;
          tau += leg;
        }

        const double Et = Z_BOSON_MASS/2; // GeV

        const double rx = vx - BS_x;
        const double ry = vy - BS_y;
        const double rz = vz - BS_z;

        double px = tau.Px();
        double py = tau.Py();
        double pz = tau.Pz();

        double E3pi = tau.E();
        double p3pi = tau.P();
        double r = std::sqrt(rx*rx + ry*ry + rz*rz);

        double cb = (rx*px + ry*py + rz*pz)/(p3pi*r);
        cb = std::clamp(cb, -1.0, 1.0);

        double A = 2*E3pi*Et;
        double B = tau.M2();
        double C = 2*p3pi*p3pi*cb*cb;
        double D = 2*p3pi*cb*std::sqrt(C/2 + Et*Et - A + B);

        return {std::sqrt(A - B - C - D), std::acos(cb)};
    }

    ///// Computes the reconstructed tau mass, but using the MC values for 4-momentum of the 3pi. /////
    ///// Used to study the impact of the 3pi reconstruction on the final value of tau mass       /////
    std::pair<double,double> TauMass(const float px, const float py, const float pz, // MC values
                                     const double m3pi_MC, // MC value
                                     const float vx, const float vy, const float vz,
                                     const float BS_x, const float BS_y, const float BS_z) {
        ROOT::Math::PxPyPzMVector tau;
        tau.SetCoordinates(px, py, pz, m3pi_MC);

        const double Et = Z_BOSON_MASS/2; // GeV

        const double rx = vx - BS_x;
        const double ry = vy - BS_y;
        const double rz = vz - BS_z;

        double E3pi = tau.E();
        double p3pi = tau.P();
        double r = std::sqrt(rx*rx + ry*ry + rz*rz);

        double cb = (rx*px + ry*py + rz*pz)/(p3pi*r);
        cb = std::clamp(cb, -1.0, 1.0);

        double A = 2*E3pi*Et;
        double B = tau.M2();
        double C = 2*p3pi*p3pi*cb*cb;
        double D = 2*p3pi*cb*std::sqrt(C/2 + Et*Et - A + B);

        return {std::sqrt(A - B - C - D), std::acos(cb)};
    }

    ///// Computes the 3pi 3-momentum at vertex after the vertex fit. /////
    ///// Used to extract the components of the 3-pion system         /////
    TVector3 TauVertexMomentum(const FCCAnalyses::VertexingUtils::FCCAnalysesVertex& v) {
        TVector3 p;
        for (const auto& p_tr : v.updated_track_momentum_at_vertex) {
          p += p_tr;
        }
        return p;
    }

    ///// Estimator of the primary vertex (PV) constrained by the beam spot resolution.        /////
    ///// Check that the sigmas are properly passed in mm (they are typically reported in um), /////
    ///// since v1 and v2 are in mm                                                            /////
    inline edm4hep::Vector3f EstimatePV(const float v1_x, const float v1_y, const float v1_z,
                                        const float v2_x, const float v2_y, const float v2_z,
                                        const float sig_x, const float sig_y, const float sig_z) {   
        // direction vector between the two decay vertices
        edm4hep::Vector3f d = {v2_x - v1_x, v2_y - v1_y, v2_z - v1_z};

        // weights set by the beam spot constraints
        float wx = 1.0/(sig_x*sig_x);
        float wy = 1.0/(sig_y*sig_y);
        float wz = 1.0/(sig_z*sig_z);

        // numerator and denominator from optimization (derivative of chi^2)
        float num = wx*v1_x*d.x + wy*v1_y*d.y + wz*v1_z*d.z;
        float den = wx*d.x*d.x + wy*d.y*d.y + wz*d.z*d.z;

        float t = -num/den; // solve

        return {v1_x + t*d.x, v1_y + t*d.y, v1_z + t*d.z};
    }

    ///// Computes the SUM of the two taus lifetimes' in an event, given both SV /////
    double TauLifetime(const float SV1_x, const float SV1_y, const float SV1_z,
                       const float SV2_x, const float SV2_y, const float SV2_z) {
        const double lambda_tau = std::sqrt( (SV1_x - SV2_x)*(SV1_x - SV2_x) + (SV1_y - SV2_y)*(SV1_y - SV2_y) + (SV1_z - SV2_z)*(SV1_z - SV2_z) ) / 1000; // flight distance of the tau in mm

        return ( ( (lambda_tau * TAU_MASS) / std::sqrt( E_TAU*E_TAU - TAU_MASS*TAU_MASS) ) / c) * 1e15; // lifetime in fs
    }
    
    ///// Computes the post-vertex-fit mass of the 3-pion system /////
    double tau3pi_vertex_mass(const VertexingUtils::FCCAnalysesVertex& vertex) {
        TLorentzVector tau;

        ROOT::VecOps::RVec<TVector3> momenta = vertex.updated_track_momentum_at_vertex;
        int n = momenta.size();
        
        if (n!=3)
          cout << n << endl;
        
        for (int ileg=0; ileg < n; ileg++) {
          TVector3 track_momentum = momenta[ileg];
          TLorentzVector leg;
          leg.SetXYZM(track_momentum[0], track_momentum[1], track_momentum[2], PION_MASS) ;
          tau += leg;
        }

        return tau.M();
    }

    ///// Computes the pre-vertex-fit mass of the 3-pion system /////
    double tau3pi_raw_mass(const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& legs) {
        TLorentzVector tau;
        
        int n = legs.size();
        //cout << n << endl;
        
        for (int ileg=0; ileg < 3; ileg++) {
          TLorentzVector leg;
          leg.SetXYZM(legs[ileg].momentum.x, legs[ileg].momentum.y, legs[ileg].momentum.z, PION_MASS);
          tau += leg;
        }

        return tau.M();
    }

    ///// Computes the mass of the 3-pion system at the MC level /////
    double tau3pi_MC_mass(const ROOT::VecOps::RVec<int>& indices, const ROOT::VecOps::RVec<edm4hep::MCParticleData>& parts) {
        ROOT::Math::PxPyPzMVector tau;
        
        for (int ileg=0;ileg<3;ileg++) {
         ROOT::Math::PxPyPzMVector leg;
         edm4hep::MCParticleData part = parts.at(indices[ileg+1]);
         leg.SetCoordinates(part.momentum.x, part.momentum.y, part.momentum.z, part.mass);
         tau += leg;
        }
        
        return tau.M();
    }
}

#endif /* FCCANA_ADDITIONAL_ANALYZERS_H */