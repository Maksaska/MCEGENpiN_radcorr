#include "header.h"

double fRand(const double& fMin, const double& fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double q(const double& W)
{
	if(channel){return sqrt((W*W + m_pi0*m_pi0 - m_p*m_p)*(W*W + m_pi0*m_pi0 - m_p*m_p)/(4*W*W) - m_pi0*m_pi0);}
	else{return sqrt((W*W + m_pip*m_pip - m_n*m_n)*(W*W + m_pip*m_pip - m_n*m_n)/(4*W*W) - m_pip*m_pip);}	
}

double k(const double& W)
{
	return (W*W - m_p*m_p)/(2*W);
}

double k_mod(const double& W, const double& Q2)
{
	return sqrt(Q2 + pow(W*W - m_p*m_p - Q2, 2)/(4*W*W));
}

double Eq(const double& W, const double& Q2) /*   virtual photon energy in cm   */
{
	return (W*W - Q2 - m_p*m_p)/(2*W);
}

double lambda_q(const double& W, const double& Q2)
{
	return 4*W*W*(Eq(W, Q2)*Eq(W, Q2) + Q2);
}

double E1(const double& W, const double& Q2) /*   init. electron energy in cm   */
{
	return (2*E0*m_p - Q2)/(2*W);
}

double E2(const double& W) /*   final electron energy in cm   */
{
	return (2*E0*m_p - W*W + m_p*m_p)/(2*W);
}

double cos_1(const double& W, const double& Q2) /*   cos(theta)_1 - for polar angle of init. electron   */
{
	return ((2*E0*m_p - Q2)*(W*W - Q2 - m_p*m_p) + 2*Q2*W*W)/(2*W*sqrt(E1(W, Q2)*E1(W, Q2) - m_e*m_e)*sqrt(lambda_q(W, Q2)));
}

double cos_2(const double& W, const double& Q2) /*   cos(theta)_2 - for polar angle of final electron   */
{
	return ((2*E0*m_p - W*W + m_p*m_p)*(W*W - Q2 - m_p*m_p) - 2*Q2*W*W)/(2*W*sqrt(E2(W)*E2(W) - m_e*m_e)*sqrt(lambda_q(W, Q2)));
}

double sin_1(const double& W, const double& Q2) /*   sin(theta)_1 - for polar angle of init. electron   */
{
	return sqrt(1 - cos_1(W, Q2)*cos_1(W, Q2));
}

double sin_2(const double& W, const double& Q2) /*   sin(theta)_2 - for polar angle of final electron   */
{
	return sqrt(1 - cos_2(W, Q2)*cos_2(W, Q2));
}


