//////////////////////////////////////////////////////////////////////////////////////
///// Future Circular Collider Software - Tau-lepton Mass Peak Reconstruction    /////
///// This is a catalog of several different functions available to extend and   /////
///// complement the study of the tau mass reconstruction and provide the tools  /////
///// for analysis beyond to the ones available from the original rootfile. This /////    
///// file should be included in the "further_analysis.C" macro                  /////
//////////////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <iostream>

//////////////////////////////////////////////////////////////////////////////////////
// Constants declaration
//////////////////////////////////////////////////////////////////////////////////////

//const double PION_MASS = 0.13957;  // GeV
const double E_TAU = 91.188/2; // GeV
//const double E_TAU = 90.930509/2; // average MC sqrt(s) of the taus, GeV 
const double E_TAU_M = E_TAU - 1e-4; // (M_Z/2) minus uncertainty
const double E_TAU_P = E_TAU + 1e-4; // (M_Z/2) plus uncertainty

const double TAU_MASS = 1.77693; // +- 9e-5 GeV, from PDG average
const double TAU_MASS_M = 1.77693 - 9e-5; // m_tau minus uncertainty
const double TAU_MASS_P = 1.77693 + 9e-5; // m_tau plus uncertainty

const double c = 299792458; // m/s

//////////////////////////////////////////////////////////////////////////////////////

/* 1. The original mass function; used in the steering file. Can be modified to use a
   different value of E_TAU or to return the quantity cos(beta) */
double Mass(float SV_x, float SV_y, float SV_z,
            float SV_x_CC, float SV_y_CC, float SV_z_CC,
            double m3pi,
            double px, double py, double pz)
{
    double p3pi = std::sqrt(px*px + py*py + pz*pz);    
    double M2 = m3pi*m3pi;    
    double E3pi = std::sqrt(p3pi*p3pi + M2);
    double r = std::sqrt( (SV_x - SV_x_CC)*(SV_x - SV_x_CC) + (SV_y - SV_y_CC)*(SV_y - SV_y_CC) + (SV_z - SV_z_CC)*(SV_z - SV_z_CC) );
      
    double cb = ( (SV_x - SV_x_CC)*px + (SV_y - SV_y_CC)*py + (SV_z - SV_z_CC)*pz)/(p3pi*r);
    cb = std::clamp(cb, -1.0, 1.0);
    if (cb < 0) {
        cb = -cb;
    }
      
    double A = 2*E3pi*E_TAU;
    double B = M2;
    double C = 2*p3pi*p3pi*cb*cb;
    double D = 2*p3pi*cb*std::sqrt(C/2 + E_TAU*E_TAU - A + B);

    return std::sqrt(A - B - C - D);
    //return cb;
}

/* 2. Numerated "Mass" functions just change the data types of the arguments.
    This mass function has the scheme: float, double, double, double.
    Different types of arguments might be needed in order to study combinations
    of truth and reco values, pre- and post-fit values, etc. The data types
    may be consulted within the rootfile employed for the new analysis */
double Mass2(float SV_x, float SV_y, float SV_z,
             double SV_x_CC, double SV_y_CC, double SV_z_CC,
             double m3pi,
             double px, double py, double pz)
{
    double p3pi = std::sqrt(px*px + py*py + pz*pz);    
    double M2 = m3pi*m3pi;    
    double E3pi = std::sqrt(p3pi*p3pi + M2);
    double r = std::sqrt( (SV_x - SV_x_CC)*(SV_x - SV_x_CC) + (SV_y - SV_y_CC)*(SV_y - SV_y_CC) + (SV_z - SV_z_CC)*(SV_z - SV_z_CC) );
      
    double cb = ( (SV_x - SV_x_CC)*px + (SV_y - SV_y_CC)*py + (SV_z - SV_z_CC)*pz)/(p3pi*r);
    cb = std::clamp(cb, -1.0, 1.0);
      
    double A = 2*E3pi*E_TAU;
    double B = M2;
    double C = 2*p3pi*p3pi*cb*cb;
    double D = 2*p3pi*cb*std::sqrt(C/2 + E_TAU*E_TAU - A + B);

    return std::sqrt(A - B - C - D);
}

/* 3. Mass function with the scheme: float, float, double, float */
double Mass3(float SV_x, float SV_y, float SV_z,
             float SV_x_CC, float SV_y_CC, float SV_z_CC,
             double m3pi,
             float px, float py, float pz)
{
    double p3pi = std::sqrt(px*px + py*py + pz*pz);    
    double M2 = m3pi*m3pi;    
    double E3pi = std::sqrt(p3pi*p3pi + M2);
    double r = std::sqrt( (SV_x - SV_x_CC)*(SV_x - SV_x_CC) + (SV_y - SV_y_CC)*(SV_y - SV_y_CC) + (SV_z - SV_z_CC)*(SV_z - SV_z_CC) );
      
    double cb = ( (SV_x - SV_x_CC)*px + (SV_y - SV_y_CC)*py + (SV_z - SV_z_CC)*pz)/(p3pi*r);
    cb = std::clamp(cb, -1.0, 1.0);
    if (cb < 0) {
        cb = -cb;
    }
      
    double A = 2*E3pi*E_TAU;
    double B = M2;
    double C = 2*p3pi*p3pi*cb*cb;
    double D = 2*p3pi*cb*std::sqrt(C/2 + E_TAU*E_TAU - A + B);

    return std::sqrt(A - B - C - D);
    //return cb;
}

/* 4. Mass function with the scheme: double, double, double, double */
double Mass4(double SV_x, double SV_y, double SV_z,
             double SV_x_CC, double SV_y_CC, double SV_z_CC,
             double m3pi,
             double px, double py, double pz)
{
    double p3pi = std::sqrt(px*px + py*py + pz*pz);    
    double M2 = m3pi*m3pi;    
    double E3pi = std::sqrt(p3pi*p3pi + M2);
    double r = std::sqrt( (SV_x - SV_x_CC)*(SV_x - SV_x_CC) + (SV_y - SV_y_CC)*(SV_y - SV_y_CC) + (SV_z - SV_z_CC)*(SV_z - SV_z_CC) );
      
    double cb = ( (SV_x - SV_x_CC)*px + (SV_y - SV_y_CC)*py + (SV_z - SV_z_CC)*pz)/(p3pi*r);
    cb = std::clamp(cb, -1.0, 1.0);
    if (cb < 0) {
        cb = -cb;
    }
      
    double A = 2*E3pi*E_TAU;
    double B = M2;
    double C = 2*p3pi*p3pi*cb*cb;
    double D = 2*p3pi*cb*std::sqrt(C/2 + E_TAU*E_TAU - A + B);

    return std::sqrt(A - B - C - D);
}

/* 5. Mass function with the scheme: double, double, double, float */
double Mass5(double SV_x, double SV_y, double SV_z,
             double SV_x_CC, double SV_y_CC, double SV_z_CC,
             double m3pi,
             float px, float py, float pz)
{
    double p3pi = std::sqrt(px*px + py*py + pz*pz);    
    double M2 = m3pi*m3pi;    
    double E3pi = std::sqrt(p3pi*p3pi + M2);
    double r = std::sqrt( (SV_x - SV_x_CC)*(SV_x - SV_x_CC) + (SV_y - SV_y_CC)*(SV_y - SV_y_CC) + (SV_z - SV_z_CC)*(SV_z - SV_z_CC) );
      
    double cb = ( (SV_x - SV_x_CC)*px + (SV_y - SV_y_CC)*py + (SV_z - SV_z_CC)*pz)/(p3pi*r);
    cb = std::clamp(cb, -1.0, 1.0);
    if (cb < 0) {
        cb = -cb;
    }
      
    double A = 2*E3pi*E_TAU;
    double B = M2;
    double C = 2*p3pi*p3pi*cb*cb;
    double D = 2*p3pi*cb*std::sqrt(C/2 + E_TAU*E_TAU - A + B);

    return std::sqrt(A - B - C - D);
}

/* 6. Mass function with E_TAU = 45.5939 GeV. Used to determine one systematic */
double Mass_m(float SV_x, float SV_y, float SV_z,
              float SV_x_CC, float SV_y_CC, float SV_z_CC,
              double m3pi,
              double px, double py, double pz)
{
    double p3pi = std::sqrt(px*px + py*py + pz*pz);    
    double M2 = m3pi*m3pi;    
    double E3pi = std::sqrt(p3pi*p3pi + M2);
    double r = std::sqrt( (SV_x - SV_x_CC)*(SV_x - SV_x_CC) + (SV_y - SV_y_CC)*(SV_y - SV_y_CC) + (SV_z - SV_z_CC)*(SV_z - SV_z_CC) );
      
    double cb = ( (SV_x - SV_x_CC)*px + (SV_y - SV_y_CC)*py + (SV_z - SV_z_CC)*pz)/(p3pi*r);
    cb = std::clamp(cb, -1.0, 1.0);
      
    double A = 2*E3pi*E_TAU_M;
    double B = M2;
    double C = 2*p3pi*p3pi*cb*cb;
    double D = 2*p3pi*cb*std::sqrt(C/2 + E_TAU_M*E_TAU_M - A + B);

    return std::sqrt(A - B - C - D);
}

/* 7. Mass function with E_TAU = 45.5941 GeV. Used to determine one systematic */
double Mass_p(float SV_x, float SV_y, float SV_z,
              float SV_x_CC, float SV_y_CC, float SV_z_CC,
              double m3pi,
              double px, double py, double pz)
{
    double p3pi = std::sqrt(px*px + py*py + pz*pz);    
    double M2 = m3pi*m3pi;    
    double E3pi = std::sqrt(p3pi*p3pi + M2);
    double r = std::sqrt( (SV_x - SV_x_CC)*(SV_x - SV_x_CC) + (SV_y - SV_y_CC)*(SV_y - SV_y_CC) + (SV_z - SV_z_CC)*(SV_z - SV_z_CC) );
      
    double cb = ( (SV_x - SV_x_CC)*px + (SV_y - SV_y_CC)*py + (SV_z - SV_z_CC)*pz)/(p3pi*r);
    cb = std::clamp(cb, -1.0, 1.0);
      
    double A = 2*E3pi*E_TAU_P;
    double B = M2;
    double C = 2*p3pi*p3pi*cb*cb;
    double D = 2*p3pi*cb*std::sqrt(C/2 + E_TAU_P*E_TAU_P - A + B);

    return std::sqrt(A - B - C - D);
}

/* 8. Mass function which takes cos(beta) as an argument insted of the vertices */
double Mass_cos(double m3pi, double cosBeta,
                double px, double py, double pz)
{
    double p3pi = std::sqrt(px*px + py*py + pz*pz);    
    double M2 = m3pi*m3pi;    
    double E3pi = std::sqrt(p3pi*p3pi + M2);
      
    double A = 2*E3pi*E_TAU;
    double B = M2;
    double C = 2*p3pi*p3pi*cosBeta*cosBeta;
    double D = 2*p3pi*cosBeta*std::sqrt(C/2 + E_TAU*E_TAU - A + B);

    return std::sqrt(A - B - C - D);
}

/* 9. Mass function which takes the tau energy as an explicit argument.
      Scheme: double, double, double, float, float. Used with truth energy */
double MassWithEnergy_MC(double SV_x, double SV_y, double SV_z,
                         double SV_x_CC, double SV_y_CC, double SV_z_CC,
                         double m3pi,
                         float E_t,
                         float px, float py, float pz)
{
    double p3pi = std::sqrt(px*px + py*py + pz*pz);    
    double M2 = m3pi*m3pi;    
    double E3pi = std::sqrt(p3pi*p3pi + M2);
    double r = std::sqrt( (SV_x - SV_x_CC)*(SV_x - SV_x_CC) + (SV_y - SV_y_CC)*(SV_y - SV_y_CC) + (SV_z - SV_z_CC)*(SV_z - SV_z_CC) );
      
    double cb = ( (SV_x - SV_x_CC)*px + (SV_y - SV_y_CC)*py + (SV_z - SV_z_CC)*pz)/(p3pi*r);
    cb = std::clamp(cb, -1.0, 1.0);
    if (cb < 0) {
        cb = -cb;
    }

    double A = 2*E3pi*E_t;
    double B = M2;
    double C = 2*p3pi*p3pi*cb*cb;
    double D = 2*p3pi*cb*std::sqrt(C/2 + E_t*E_t - A + B);

    return std::sqrt(A - B - C - D);
}

/* 10. Mass function which takes the tau energy as an explicit argument.
       Scheme: float, float, double, float, double */
double MassWithEnergy(float SV_x, float SV_y, float SV_z,
                      float SV_x_CC, float SV_y_CC, float SV_z_CC,
                      double m3pi,
                      float E_t,
                      double px, double py, double pz)
{
    double p3pi = std::sqrt(px*px + py*py + pz*pz);    
    double M2 = m3pi*m3pi;    
    double E3pi = std::sqrt(p3pi*p3pi + M2);
    double r = std::sqrt( (SV_x - SV_x_CC)*(SV_x - SV_x_CC) + (SV_y - SV_y_CC)*(SV_y - SV_y_CC) + (SV_z - SV_z_CC)*(SV_z - SV_z_CC) );
      
    double cb = ( (SV_x - SV_x_CC)*px + (SV_y - SV_y_CC)*py + (SV_z - SV_z_CC)*pz)/(p3pi*r);
    cb = std::clamp(cb, -1.0, 1.0);
    if (cb < 0) {
        cb = -cb;
    }
      
    double A = 2*E3pi*E_t;
    double B = M2;
    double C = 2*p3pi*p3pi*cb*cb;
    double D = 2*p3pi*cb*std::sqrt(C/2 + E_t*E_t - A + B);

    return std::sqrt(A - B - C - D);
}

/* 11. Computes the kinematic variable cos(beta). Scheme: float, float, double */
double cos_beta(float SV_x, float SV_y, float SV_z,
                float SV_x_CC, float SV_y_CC, float SV_z_CC,
                double px, double py, double pz)
{
    double p3pi = std::sqrt(px*px + py*py + pz*pz);
    double r = std::sqrt( (SV_x - SV_x_CC)*(SV_x - SV_x_CC) + (SV_y - SV_y_CC)*(SV_y - SV_y_CC) + (SV_z - SV_z_CC)*(SV_z - SV_z_CC) );
      
    double cb = ( (SV_x - SV_x_CC)*px + (SV_y - SV_y_CC)*py + (SV_z - SV_z_CC)*pz)/(p3pi*r);
    cb = std::clamp(cb, -1.0, 1.0);
    if (cb < 0) {
        cb = -cb;
    }
    
    return cb;
}

/* 12. Computes the kinematic variable cos(beta). Scheme: double, double, float */
double cos_beta2(double SV_x, double SV_y, double SV_z,
                 double SV_x_CC, double SV_y_CC, double SV_z_CC,
                 float px, float py, float pz)
{
    double p3pi = std::sqrt(px*px + py*py + pz*pz);
    double r = std::sqrt( (SV_x - SV_x_CC)*(SV_x - SV_x_CC) + (SV_y - SV_y_CC)*(SV_y - SV_y_CC) + (SV_z - SV_z_CC)*(SV_z - SV_z_CC) );
      
    double cb = ( (SV_x - SV_x_CC)*px + (SV_y - SV_y_CC)*py + (SV_z - SV_z_CC)*pz)/(p3pi*r);
    cb = std::clamp(cb, -1.0, 1.0);
    if (cb < 0) {
        cb = -cb;
    }
    
    return cb;
}

/* 13. Computes the polar angle */
double cos_theta(const double px, const double py, const double pz)
{
    return pz/std::sqrt(px * px + py * py + pz * pz);
}

/* 14. Computes the lifetime of one tau given the SV and an estimate of the PV */
double Lifetime(const float PV_x, const float PV_y, const float PV_z,
                const float SV_x, const float SV_y, const float SV_z)
{
    const double lambda_tau = std::sqrt( (SV_x - PV_x)*(SV_x - PV_x) + (SV_y - PV_y)*(SV_y - PV_y) + (SV_z - PV_z)*(SV_z - PV_z) ) / 1000;
    
    return ( ( (lambda_tau * TAU_MASS) / std::sqrt( E_TAU*E_TAU - TAU_MASS*TAU_MASS) ) / c) * 1e15; // lifetime in fs
}

/* 15. Compute the sum of both taus' lifetimes. No PV estimate is required */
double SumLifetimes(const float SV1_x, const float SV1_y, const float SV1_z,
                    const float SV2_x, const float SV2_y, const float SV2_z)
{
    const double lambda_tau = std::sqrt( (SV1_x - SV2_x)*(SV1_x - SV2_x) + (SV1_y - SV2_y)*(SV1_y - SV2_y) + (SV1_z - SV2_z)*(SV1_z - SV2_z) ) / 1000;
    
    return ( ( (lambda_tau * TAU_MASS) / std::sqrt( E_TAU*E_TAU - TAU_MASS*TAU_MASS) ) / c) * 1e15; // sum of lifetimes in fs
}

/* 16. Cosine between two vectors. Is used to compute the acolinearity between 
       the reconstructed vectors of the two 3-pion systems in an event */
double DotProductUni(const double ux, const double uy, const double uz,
                     const double vx, const double vy, const double vz)
{
    double abs_u = std::sqrt(ux*ux + uy*uy + uz*uz);
    double abs_v = std::sqrt(vx*vx + vy*vy + vz*vz);
    return (ux*vx + uy*vy + uz*vz)/(abs_u*abs_v);
}

/* 17. Introduces corrections to the nominal value used for the tau energy, E_TAU.
       Might be used, for example, to substract the energy radiated through photons
       via ISR/FSR in an event. Takes the energy to be substracted as the argument */
float Etau(const float e)
{
    return ((E_TAU * 2) - e)/2;
}

///////////////////////////////////////////////////////////////////////////////////////
// The purpose of the functions in the block below was trying to smear the truth values
// of the vertex using the modeled resolution of the vertex detector, but didn't
// produced the expected result, so it wasn't studied anymore
///////////////////////////////////////////////////////////////////////////////////////

struct Mix3Gauss {
  double s1, s2, s3;
  double f1, f2; // Weights: f1, f2, and (1 - f1 - f2)
};

/* 18. Build the detector resolution as a sum of three Gaussians with shared mu, different sigma */
double SampleMixture(TRandom3 &rng, const Mix3Gauss &p)
{
  const double u = rng.Rndm();
  if (u < p.f1) {
    return rng.Gaus(0.0, p.s1);
  } else if (u < p.f1 + p.f2) {
    return rng.Gaus(0.0, p.s2);
  }
  return rng.Gaus(0.0, p.s3);
}

constexpr double K = 1.0; // Constant to improve or inflate vertex resolution
constexpr Mix3Gauss kX{K*4.429e-3, K*1.746e-2, K*4.028e-2, 1.914e-1, 5.202e-1}; // Resolution in x
constexpr Mix3Gauss kY{K*4.431e-3, K*1.743e-2, K*4.023e-2, 1.917e-1, 5.194e-1}; // Resolution in y
constexpr Mix3Gauss kZ{K*6.467e-3, K*2.837e-2, K*1.326e-1, 3.314e-1, 5.619e-1}; // Resolution in z

TRandom3 rng(0);

/* 19. Smears the truth-level vertex positions given the resolutions for each direction */
double smear_vtx(double r, const Mix3Gauss& k) {return r + SampleMixture(rng, k);}

///////////////////////////////////////////////////////////////////////////////////////