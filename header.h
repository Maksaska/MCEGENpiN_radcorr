#include <vector>
#include <iostream>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <ctime>
#include <chrono>
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;

extern double E0, R, L, W_min, W_max, Q2_min, Q2_max, delta;
extern vector<vector<double>> data, data_interp;
extern int h, N;
extern double seed_;
extern string path, source, source_interp;
extern bool channel, method, histogram, rad_corr;
extern double m_p, m_n, m_pip, m_pi0, m_e;
extern vector<double> values_rad;

extern double Q2_line, W_line;

extern TH2* h1;
extern TH1F* h3;
extern TH1F* h4;
extern TH2* h5;

extern int Total_Number, Accepted_Number;

/* -------- kinematics.cpp -------- */

double fRand(const double& fMin, const double& fMax); /*   random function from [fMin,fMax] area   */
double q(const double& W); /*   Pion momentum   */
double k(const double& W); /*   Photon equivalent energy   */
double k_mod(const double& W, const double& Q2); /*   Virtual photon momentum   */
double E1(const double& W, const double& Q2); /*   init. electron energy in cm   */
double E2(const double& W); /*   final electron energy in cm   */
double Eq(const double& W, const double& Q2); /*   virtual photon energy in cm   */
double lambda_q(const double& W, const double& Q2);
double cos_1(const double& W, const double& Q2); /*   cos(theta)_1 - for polar angle of init. electron   */
double cos_2(const double& W, const double& Q2); /*   cos(theta)_2 - for polar angle of final electron   */
double sin_1(const double& W, const double& Q2); /*   sin(theta)_1 - for polar angle of init. electron   */
double sin_2(const double& W, const double& Q2); /*   sin(theta)_2 - for polar angle of final electron   */

/* -------- general.cpp -------- */

void input_check(int argc, char* argv[]); /*   This function checks if there any option were passed in main()   */
void Reading(string Path,vector<vector<double>>&V); /*   Reads the data from csv   */
double P(const int& der,const int& n, const double& theta); /*   Legendre polynomials   */
vector<complex<double>> Finder(const double& W,  const double& Q2); /*   it searches for multipoles EMS   */
void generate_particle(const int& k); /*
                * This function generates the event
                * It creates kin. parameters from given E0 [W_min, W_max] && [Q2_min, Q2_max]
                * Gets the cross-section
                * Creates the Particle
                * Creates the Event
                * Writes the event in output file
                                    */

vector<double> Coefficients_lin(const double& x1, const double& y1, const double& x2, const double& y2); /*   Coef. for linear interp.   */
double Linear(const double& W,  const double& Q2, const double& theta,  const double& phi); /*   linear interpolation   */
double Section(const double& W,  const double& Q2, const double& theta,  const double& phi); /*
                                 Cross_section evaluation from Helicity ampl.
                                                         */

vector<complex<double>> Helicity_amplitudes(const double& W,  const double& Q2, const double& theta);/*   Helicity ampl. eval.   */
double Section_int(const double& W, const double& Q2, const double& E_beam); /*   Integral cross section dS/dOmegadE from dataset   */
double Section_interp_int(const double& W, const double& Q2, const double& E_beam); /*   Integral cross section dS/dOmegadE for random W, Q2   */
double Spence(const double& x); /*   Spence function   */
double delta_r(const double& W, const double& Q2);
double ts(const double& W, const double& Q2, const double& iter);
double tp(const double& W, const double& Q2, const double& iter);
double Section_interp_Q2_extra(const double& W,  const double& Q2, const double& theta,  const double& phi); /*    Q2 extrapolation with W in [1.08, 2.0] for MAID data   */
double Section_interp_int_Q2_extra(const double& W, const double& Q2); /*  Q2 extrapolation with W in [1.08, 2.0] for R3 */
double Section_interp_int_Q2_extra_R2(const double& W, const double& Q2, const double& E0); /*  Q2 extrapolation with W in [1.08, 2.0] for R2 */
double R1(const double& W, const double& Q2); /*   soft region for RC   */
double R2(const double& W, const double& Q2); /*   Hard radiation with init. elec.   */
double R3(const double& W, const double& Q2); /*   Hard radiation with final elec.   */
double Gen_omega_init(const double& W, const double& Q2); /*   Rad. photon energy generation from init. elec.   */
double Gen_omega_fin(const double& W, const double& Q2); /*   Rad. photon energy generation from fin. elec.   */

double Normal(const double& mean, const double& sigma); /*   Normal distro   */
