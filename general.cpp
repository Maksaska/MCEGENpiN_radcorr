#include "header.h"
#include "event.h"
#include <getopt.h>
#include <fstream>

using namespace std;

double E0(6.5), R(1), L(10), W_min(1.08), W_max(2.0), Q2_min(0.05), Q2_max(5.0);
int h(0), N(1000000); double seed_(0);
bool channel(true), method(true), histogram(false), rad_corr(false);
string path = "MCEGENpiN_radcorr.dat";
string source = "pi0p.csv";
string source_interp = "pi0p_int.csv";
vector<vector<double>> data, data_interp;
double m_p(0.93827), m_n(0.93957), m_pip(0.13957), m_pi0(0.13498), m_e(0.000511), delta(0.01); 
vector<double> values_rad{9999, 1, 9999, 1};

void add_hadrons(const double& W_, const double& Q2_, const double& theta, const double& phi, Event& bob); /*   Adds hadron system in cm into Event object   */
void add_electron(const double& W, const double& Q2, Event& bob); /*   Adds outgoing electron in cm into Event object   */
void add_rad_photon(const double& Erad, const double& ang1, const double& ang2, Event& bob); /*   Adds outgoing rad photon in cm into Event object   */

void Reading(string Path, vector<vector<double>>& V)
{
	string line; stringstream ss;
	ifstream File;
	double dub; int w(1); vector<double> Numbers;

	File.open(Path,fstream::in | fstream::out | fstream::app);

	if (!File.is_open())
	{
		cout << "Can't open " << Path << " !" << endl;
	}
	else
	{	
		while(!File.eof())
		{		
			if (w == 1)
			{
				getline(File,line);
				w = 0;
			}

			getline(File,line); 	
			ss << line;

			if(File.eof())
			{
				break;
			}

			while(ss>>dub)
			{							
				Numbers.push_back(dub);
				while(ss.peek() == ',')
            			ss.ignore();
			}
	
			V.push_back(Numbers);
			
			Numbers.clear();

			ss.clear();
		}
	}
	File.close();
}

void input_check(int argc, char* argv[])
{
	const char* short_options = "nwe:r:l:N:h:z:x:c:v:pmas:"; int rez; int option_index;
	
	const struct option long_options[] = {
						{"beam_energy", required_argument, NULL, 'e'},
	       				{"target_R", required_argument, NULL, 'r'},
	        				{"target_L", required_argument, NULL, 'l'},
	        				{"W_min", required_argument, NULL, 'z'},
	        				{"W_max", required_argument, NULL, 'x'},
	        				{"Q2_min", required_argument, NULL, 'c'},
	        				{"Q2_max", required_argument, NULL, 'v'},
	        				{"hist", no_argument, NULL, 'p'},
	        				{"RC", no_argument, NULL, 'm'}, 
	        				{"trig", required_argument, NULL, 'N'}, 
	        				{"docker", no_argument, NULL, 'a'}, 
	        				{"seed", required_argument, NULL, 's'}, 
	        				{NULL, 0, NULL, 0}
								};
	while ((rez=getopt_long(argc, argv, short_options, long_options, &option_index)) != -1)
	{ 
		switch(rez)
		{
			case 'n': 
			{
				channel = false;
				path = "MCEGENpiN_radcorr.dat"; 
				source = "pin.csv";
				source_interp = "pin_int.csv";
				break;
			};			
			case 'w': 
			{
				method = false;
				break;
			};
			case 'p': 
			{
				histogram = true;
				break;
			};	
			case 'm': 
			{
				rad_corr = true;
				break;
			};
			case 'e': {
				E0 = atof(optarg);
				break;
			};
			case 'r': {
				R = atof(optarg);
				break;
			};
			case 'l': {
				L = atof(optarg);
				break;
			};
			case 's': {
				seed_ = atof(optarg);
				break;
			};
			case 'N': {
				N = atoi(optarg);
				break;
			};
			case 'h': {
				h = atoi(optarg);
				break;
			};
			case 'z': {
				W_min = atof(optarg);
				if(W_min < 1.08){W_min = 1.08;}
				break;
			};
			case 'x': {
				W_max = atof(optarg);
				if(W_max > 2.0){W_max = 2.0;}
				break;
			};
			case 'c': {
				Q2_min = atof(optarg);
				if(Q2_min < 0.05){Q2_min = 0.05;}
				break;
			};
			case 'v': {
				Q2_max = atof(optarg);
				if(Q2_max > 5.0){Q2_max = 5.0;}
				break;
			};
			case 'a': {
				E0 = 6.5; R = 1; L = 10; W_min = 1.08; W_max = 2.0; Q2_min = 0.05; Q2_max = 5.0;
				h = 0; N = 1000000; 
				channel = true; method = true; histogram = false; rad_corr = false;
				break;
			};	
			case '?': default: {
				cerr << "Unkhown option" << endl;
				break;
			};
		};
	};
	
	cout << " ------------------------------------------------------------------- " << endl;
	cout << "| Monte Carlo event generator for exclusive pion electroproduction  | \n| with radiative corrections              \"MCEGENpiN_radcorr V7e\"   |       \n|                                                                   |\n|     Authors: Davydov M. - MSU, Physics dep.                       |\n|              Isupov E.  - MSU, SINP                               |\n|                                                                   |\n| https://github.com/Maksaska/pi0p-pin-generator                    |\n ------------------------------------------------------------------- " << endl;
	
	cout << endl;	
	
	cout << "Beam energy is E = " << E0 << " GeV with";
	if(h != 0){cout << " polarization h = " << h << endl;}
	else{cout << " no polarization" << endl;}
	cout << "Chosen kinematic area of (W, Q^2) values is W: " << W_min << " - " << W_max << " GeV , Q2: " << Q2_min << " - " << Q2_max << " GeV^2" << endl;
	if(channel){cout << "Channel: pi0p with decay" << endl;}
	else{cout << "Channel: pi+n" << endl;}
	cout << "Uniform distribution with weights" << endl;
	cout << "Number of events: " << N << endl;
	if(histogram){cout << "Histograms  will be created\n" << endl;}
	if(rad_corr){cout << "Radiative corrections: Enabled\n" << endl;}
	else{cout << "Radiative corrections: Disabled\n" << endl;}	
	
	Reading(source, data);	
	Reading(source_interp, data_interp);
}

double P(const int& der,const int& n, const double& theta) 
{
	double Value(0);

	switch(n)
	{
	
		case 0: if(der == 0)
			{
				Value = 1;
			}
			else
			{
				Value = 0;
			}
			break;

		case 1: if(der == 0){Value = cos(theta);} 
			if(der == 1){Value = 1;}
			if(der == 2){Value = 0;}
			break;

		case 2: if(der == 0){Value = 0.5*(3*pow(cos(theta),2) - 1);} 
			if(der == 1){Value = 3*cos(theta);}
			if(der == 2){Value = 3;}
			break;

		case 3: if(der == 0){Value = 0.5*(5*pow(cos(theta),3) - 3*cos(theta));} 
			if(der == 1){Value = 7.5*pow(cos(theta),2)-1.5;}
			if(der == 2){Value = 15*cos(theta);}
			break;

		case 4: if(der == 0){Value = 0.125*(35*pow(cos(theta),4) - 30*pow(cos(theta),2) + 3);} 
			if(der == 1){Value = 17.5*pow(cos(theta),3) - 7.5*cos(theta);}
			if(der == 2){Value = 52.5*pow(cos(theta),2) - 7.5;}
			break;

		case 5: if(der == 0){Value = 0.125*(63*pow(cos(theta),5) - 70*pow(cos(theta),3) + 15*cos(theta));} 
			if(der == 1){Value = 315*pow(cos(theta),4)/8 - 210*pow(cos(theta),2)/8 + 15/8;}
			if(der == 2){Value = 315*pow(cos(theta),3)/2 - 210*cos(theta)/4;}
			break;

		case 6: if(der == 0){Value = 231*pow(cos(theta),6)/16 - 315*pow(cos(theta),4)/16 + 105*pow(cos(theta),2)/16 - 5/16;} 
			if(der == 1){Value = 3*231*pow(cos(theta),5)/8 - 315*pow(cos(theta),3)/4 + 105*cos(theta)/8  ;}
			if(der == 2){Value = 15*231*pow(cos(theta),4)/8 - 945*pow(cos(theta),2)/4 + 105/8;}
			break;

		default: Value = 0; break;
	}
	
	return Value;
}

vector<complex<double>> Finder(const double& W,  const double& Q2)
{
	vector<complex<double>> mult;
	complex<double> buff; 
	
	for(long unsigned int i = 0; i < data.size(); i++)
	{
		if(data[i][0] == W and data[i][1] == Q2)
		{
			for(long unsigned int j = 26; j < data[i].size()-1; j++)
			{
				if(j == 38 or j == 50 or j == 62)
				{
					buff.real(0);
					buff.imag(0);
					mult.push_back(buff);
				}
				buff.real(data[i][j]);
				buff.imag(data[i][j+1]);
				mult.push_back(buff); j++;
			}
			
			buff.real(0);
			buff.imag(0);
			mult.push_back(buff);
			
			for(int j = 2; j < 25; j++)
			{
				if(j == 14)
				{
					buff.real(0);
					buff.imag(0);
					mult.push_back(buff);
				}
				buff.real(data[i][j]);
				buff.imag(data[i][j+1]);
				mult.push_back(buff); j++;
			}
			
			buff.real(0);
			buff.imag(0);
			mult.push_back(buff);
			
			return mult;	
		}
	}

	return mult;
}

vector<complex<double>> Helicity_amplitudes(const double& W,  const double& Q2, const double& theta)
{
	vector<complex<double>> AMP, H; double l;
	complex<double> buff; buff = 0;
	
	AMP = Finder(W, Q2);
	
	for(int i = 0; i < 6; i++) // H1
	{
		l = double(i);
		buff += (AMP[l] - AMP[l+14] - AMP[l+8] - AMP[l+22])*(P(2,i,theta) - P(2,i+1,theta));
	} 
	buff = buff*sin(theta)*cos(theta/2)/sqrt(2);
	H.push_back(buff); buff = 0;
	
	for(int i = 0; i < 6; i++) // H2
	{
		l = double(i);
		buff += ((l+2)*AMP[l] + l*AMP[l+14] + l*AMP[l+8] - (l+2)*AMP[l+22])*(P(1,i,theta) - P(1,i+1,theta));
	} 
	buff = buff*cos(theta/2)/sqrt(2);
	H.push_back(buff); buff = 0;
	
	for(int i = 0; i < 6; i++) // H3
	{
		l = double(i);
		buff += (AMP[l] - AMP[l+14] + AMP[l+8] + AMP[l+22])*(P(2,i,theta) + P(2,i+1,theta));
	} 
	buff = buff*sin(theta)*sin(theta/2)/sqrt(2);
	H.push_back(buff); buff = 0;
	
	for(int i = 0; i < 6; i++) // H4
	{
		l = double(i);
		buff += ((l+2)*AMP[l] + l*AMP[l+14] - l*AMP[l+8] + (l+2)*AMP[l+22])*(P(1,i,theta) + P(1,i+1,theta));
	} 
	buff = buff*sin(theta/2)/sqrt(2);
	H.push_back(buff);buff = 0;
	
	for(int i = 0; i < 6; i++) // H5
	{
		l = double(i);
		buff += (l+1)*(AMP[l+28] + AMP[l+36])*(P(1,i,theta) - P(1,i+1,theta));
	} 
	buff = buff*cos(theta/2)*sqrt(Q2)/k_mod(W,Q2);
	H.push_back(buff);buff = 0;
	
	for(int i = 0; i < 6; i++) // H6
	{
		l = double(i);
		buff += (l+1)*(AMP[l+28] - AMP[l+36])*(P(1,i,theta) + P(1,i+1,theta));
	} 
	buff = buff*sin(theta/2)*sqrt(Q2)/k_mod(W,Q2);
	H.push_back(buff); 
	
	AMP.clear();
	return H;
}

double Section(const double& W,  const double& Q2, const double& theta,  const double& phi)
{
	double S_t, S_l, S_tt, S_lt, S_lt_pr, S, Gamma_flux, eps, nu, L(0.14817); 
	vector<complex<double>> H;
	
	nu =  (W*W + Q2 - m_p*m_p)/(2*m_p);
	eps = 1/(1 + 2*(nu*nu + Q2)/(4*(E0 - nu)*E0 - Q2));		
	Gamma_flux = (E0 - nu)*(W*W - m_p*m_p)/(137*4*M_PI*M_PI*m_p*E0*(1 - eps)*Q2);
	
	H = Helicity_amplitudes(W, Q2, theta);
	
	S_t = pow(L, 2)*q(W)*(abs(H[0])*abs(H[0]) + abs(H[1])*abs(H[1]) + abs(H[2])*abs(H[2]) + abs(H[3])*abs(H[3]))/(2*k(W));
	S_l = pow(L, 2)*q(W)*(abs(H[4])*abs(H[4]) + abs(H[5])*abs(H[5]))/k(W);
	S_tt = pow(L, 2)*q(W)*real(H[2]*conj(H[1]) - H[3]*conj(H[0]))/k(W);
	S_lt = -pow(L, 2)*q(W)*real((H[0] - H[3])*conj(H[4]) + (H[1] + H[2])*conj(H[5]))/(sqrt(2)*k(W));
	S_lt_pr = -pow(L, 2)*q(W)*imag((H[0] - H[3])*conj(H[4]) + (H[1] + H[2])*conj(H[5]))/(sqrt(2)*k(W));

	S = S_t + eps*S_l + eps*S_tt*cos(2*phi) + sqrt(2*eps*(1+eps))*S_lt*cos(phi) + double(h)*S_lt_pr*sin(phi); 

	S = Gamma_flux*S; 
	
	return S;
}

vector<double> Coefficients_lin(const double& x1, const double& y1, const double& x2, const double& y2)
{
	vector<double> Abc(2);

	Abc[0] = (y1 - y2)/(x1 - x2);
	Abc[1] = (x1*y2 - y1*x2)/(x1 - x2);

	return Abc;
}

double Linear(const double& W,  const double& Q2, const double& theta,  const double& phi)
{
	double S; vector<double> ab;
	double W_min_d, W_max_d, Q2_min_d, Q2_max_d, S1, S2, S3, S4, S_1, S_2;

	W_min_d = floor(W*100)/100;
	W_max_d = ceil(W*100)/100;
	Q2_min_d = floor(Q2*20)/20;
	Q2_max_d = ceil(Q2*20)/20;

	S1 = Section(W_min_d, Q2_min_d, theta, phi);
	S2 = Section(W_min_d, Q2_max_d, theta, phi); 
	S3 = Section(W_max_d, Q2_min_d, theta, phi);
	S4 = Section(W_max_d, Q2_max_d, theta, phi); 
	
	if(S1 == S3){S_1 = S1;}
	else
	{
		ab = Coefficients_lin(W_min, S1, W_max, S3);
		S_1 = ab[0]*W + ab[1]; ab.clear();	
	}  
	
	if(S2 == S4){S_2 = S2;}
	else
	{
		ab = Coefficients_lin(W_min, S2, W_max, S4);
		S_2 = ab[0]*W + ab[1]; ab.clear();	
	}
		
	if(S_1 == S_2){return S_1;}
	else
	{
		ab = Coefficients_lin(Q2_min, S_1, Q2_max, S_2);
		S = ab[0]*Q2 + ab[1]; ab.clear(); 	
	}

	return S;
}

double Section_int(const double& W, const double& Q2, const double& E_ini)
{
	double nu, eps, Gamma_flux;
	
	nu =  (W*W + Q2 - m_p*m_p)/(2*m_p);
	eps = 1/(1 + 2*(nu*nu + Q2)/(4*(E_ini - nu)*E_ini - Q2));		
	Gamma_flux = (E_ini - nu)*(W*W - m_p*m_p)/(137*4*M_PI*M_PI*m_p*E_ini*(1 - eps)*Q2);
	
	if(std::isnan(Gamma_flux) or std::isnan(eps)){return 0;}
	
	for(auto i:data_interp)
	{
		if(W == i[0] and Q2 == i[1])
		{
			return Gamma_flux*(i[2] + eps*i[3]);			
		}
	}	
	return 0;
}

double Section_interp_int(const double& W, const double& Q2, const double& E_ini)
{
	double S; vector<double> ab;
	double W_min_d, W_max_d, Q2_min_d, Q2_max_d, S1, S2, S3, S4, S_1, S_2;

	W_min_d = floor(W*100)/100;
	W_max_d = ceil(W*100)/100;
	Q2_min_d = floor(Q2*20)/20;
	Q2_max_d = ceil(Q2*20)/20;

	S1 = Section_int(W_min_d, Q2_min_d, E_ini); 
	S2 = Section_int(W_min_d, Q2_max_d, E_ini); 
	S3 = Section_int(W_max_d, Q2_min_d, E_ini);  
	S4 = Section_int(W_max_d, Q2_max_d, E_ini); 
	
	if(S1 == 0 or S2 == 0 or S3 == 0 or S4 == 0)
	{
		return 0;
	}
	
	if(S1 == S3){S_1 = S1;}
	else
	{
		ab = Coefficients_lin(W_min, S1, W_max, S3);
		S_1 = ab[0]*W + ab[1]; ab.clear(); 	
	}  
	
	if(S2 == S4){S_2 = S2;}
	else
	{
		ab = Coefficients_lin(W_min, S2, W_max, S4);
		S_2 = ab[0]*W + ab[1]; ab.clear();	
	}
		
	if(S_1 == S_2){return S_1;}
	else
	{
		ab = Coefficients_lin(Q2_min, S_1, Q2_max, S_2);
		S = ab[0]*Q2 + ab[1]; ab.clear(); 	
	}
	
	if(std::isnan(S))
	{
		return S_2;
	}
	
	return S;
}

double Section_interp_int_Q2_extra(const double& W, const double& Q2)
{
	double C1, C2, C3, y1, y2, y3, Q2_1(4.9), Q2_2(4.95), Q2_3(5.0);
	
	y1 = Section_interp_int(W, Q2_1, E0);
	y2 = Section_interp_int(W, Q2_2, E0);
	y3 = Section_interp_int(W, Q2_3, E0);
	
	C3 = ((y3 - y1)*(1/Q2_2 - 1/Q2_1)/(1/Q2_3 - 1/Q2_1) - (y2 - y1))/((1/(Q2_3*Q2_3) - 1/(Q2_1*Q2_1))*(1/Q2_2 - 1/Q2_1)/(1/Q2_3 - 1/Q2_1) - (1/(Q2_2*Q2_2) - 1/(Q2_1*Q2_1)));
	C2 = ((y2 - y1) - C3*(1/(Q2_2*Q2_2) - 1/(Q2_1*Q2_1)))/(1/Q2_2 - 1/Q2_1);
	C1 = y1 - C2/Q2_1 - C3/(Q2_1*Q2_1);
	
	return C1 + C2/Q2 + C3/(Q2*Q2);
}

double Spence(const double& x)
{
	double result(0), S1(1), S2;
	
	for(double iter = x/10; iter <= x; iter += x/10)
	{
		S2 = -log(abs(1 - iter))/iter;
		result += (S1 + S2)*x/20;
		S1 = S2;
	}
	
	return result;
}

double delta_r(const double& W, const double& Q2)
{
	TLorentzVector s, p;
	double nu =  (W*W + Q2 - m_p*m_p)/(2*m_p);
	double E_out = E0 - nu;
	
	return -(28/9 - 13*log(Q2/(m_e*m_e))/6 + (log(E0/delta) + log(E_out/delta))*(log(Q2/(m_e*m_e)) - 1) - Spence(-nu/E_out) - Spence(nu/E0))/(M_PI*137);
}

double ts(const double& W, const double& Q2, const double& iter)
{
	return (0.5*(1 + pow(iter/E0, 2))*log(Q2/(m_e*m_e)) - iter/E0)/(M_PI*137);
}

double tp(const double& W, const double& Q2, const double& iter)
{
	double nu =  (W*W + Q2 - m_p*m_p)/(2*m_p);
	double E_out = E0 - nu;
	
	return (0.5*(1 + pow(E_out/iter, 2))*log(Q2/(m_e*m_e)) - E_out/iter)/(M_PI*137);
}

double R1(const double& W, const double& Q2)
{
	return Section_interp_int(W, Q2, E0)*exp(delta_r(W, Q2));
}

double R2(const double& W, const double& Q2) 
{
	double result(0), S1, S2, W_, Es_min;
	double nu =  (W*W + Q2 - m_p*m_p)/(2*m_p);
	
	if(channel){Es_min = (m_pi0*m_pi0 + 2*m_p*m_pi0 + 2*m_p*(E0 - nu))/(2*m_p - Q2/E0);}
	else{Es_min = (m_pip*m_pip + 2*m_n*m_pip + 2*m_p*(E0 - nu))/(2*m_p - Q2/E0);}

	W_ = sqrt(m_p*m_p - Q2*Es_min/E0 + 2*(Es_min - E0 + nu)*m_p);
	
	if(std::isnan(W_))
	{
		S1 = 0;
	} else if(W_ < 1.08)
	{
		S1 = 0;
	} else 
	{	
		if(W_ > 2.0){W_ = 2.0;}
		S1 = Section_interp_int(W_, Q2*Es_min/E0, Es_min)*ts(W, Q2, Es_min)/(E0 - Es_min);	
	}

	for(double iter = Es_min + (E0 - delta - Es_min)/10; iter <= E0 - delta; iter += (E0 - delta - Es_min)/10)
	{
		W_ = sqrt(m_p*m_p - Q2*iter/E0 + 2*(iter - E0 + nu)*m_p);
		if(std::isnan(W_))
		{
			S2 = 0;
			result += (S1 + S2)*(E0 - delta - Es_min)/20;
			S1 = S2; continue;
		}
		
		if(W_ >= 1.08)
		{
			if(W_ > 2.0){W_ = 2.0;}
			S2 = Section_interp_int(W_, Q2*iter/E0, iter)*ts(W, Q2, iter)/(E0 - iter); 
		} else
		{
			S2 = 0;
		}
	
		result += (S1 + S2)*(E0 - delta - Es_min)/20;
		S1 = S2;
	}
	
	return result;
}

double R3(const double& W, const double& Q2) 
{
	double result(0), S1, S2, W_, Ep_max;
	double nu =  (W*W + Q2 - m_p*m_p)/(2*m_p);
	
	if(channel){Ep_max = (-m_pi0*m_pi0 - 2*m_p*m_pi0 + 2*m_p*E0)/(2*m_p + Q2/(E0 - nu));}
	else{Ep_max = (-m_pip*m_pip - 2*m_n*m_pip + 2*m_p*E0)/(2*m_p + Q2/(E0 - nu));}

	W_ = sqrt(m_p*m_p - Q2*(E0 - nu + delta)/(E0 - nu) + 2*(nu - delta)*m_p);
	
	if(std::isnan(W_))
	{
		S1 = 0;
	} else if(W_ < 1.08)
	{
		S1 = 0;
	}else
	{
		if(W_ > 2.0){W_ = 2.0;}
if(Q2*(E0 - nu + delta)/(E0 - nu) > 5.0){S2 = Section_interp_int_Q2_extra(W_, Q2*(E0 - nu + delta)/(E0 - nu))*tp(W, Q2, E0 - nu + delta)/delta;}
		else{S1 = Section_interp_int(W_, Q2*(E0 - nu + delta)/(E0 - nu), E0)*tp(W, Q2, E0 - nu - delta)/delta;}			
	}	

	for(double iter = E0 - nu + delta + (Ep_max - E0 + nu - delta)/10; iter <= Ep_max; iter += (Ep_max - E0 + nu - delta)/10)
	{
		W_ = sqrt(m_p*m_p - Q2*iter/(E0 - nu) + 2*(E0 - iter)*m_p);
		
		if(std::isnan(W_))
		{
			S2 = 0;
			result += (S1 + S2)*(Ep_max - E0 + nu - delta)/20;
			S1 = S2; continue;
		}
		
		if(W_ > 1.08)
		{
			if(W_ > 2.0){W_ = 2.0;}
			if(Q2*iter/(E0 - nu) > 5.0){S2 = Section_interp_int_Q2_extra(W_, Q2*iter/(E0 - nu))*tp(W, Q2, iter)/(iter - E0 + nu);}
			else{S2 = Section_interp_int(W_, Q2*iter/(E0 - nu), E0)*tp(W, Q2, iter)/(iter - E0 + nu);}		
		} else
		{
			S2 = 0;
		}
		
		result += (S1 + S2)*(Ep_max - E0 + nu - delta)/20;
		S1 = S2;
	}
	
	return result;
}

double Gen_omega_init(const double& W, const double& Q2)
{
	double Ep = E0 - (W*W + Q2 - m_p*m_p)/(2*m_p);
	double omega_max, rho(1.0), Uni(2.0), omega;
		
	if(channel)
	{
		omega_max = E0 - (m_pi0*m_pi0 + 2*m_p*m_pi0 + 2*m_p*Ep)/(2*m_p - Q2/E0);
	} else
	{
		omega_max = E0 - (m_pip*m_pip + 2*m_p*m_pip + 2*m_p*Ep)/(2*m_p - Q2/E0);
	}	
	
	omega = fRand(delta, omega_max);
	rho = delta*(0.5*(1 + pow(1 - omega/E0, 2))*log(Q2/(m_e*m_e)) - 1 + omega/E0)/(omega*(0.5*(1 + pow(1 - delta/E0, 2))*log(Q2/(m_e*m_e)) - 1 + delta/E0));
	Uni = fRand(0.0, 1.0);	
	
	if(values_rad[0] > omega_max)
	{
		values_rad[0] = fRand(delta, omega_max);
		values_rad[1] = delta*(0.5*(1 + pow(1 - values_rad[0]/E0, 2))*log(Q2/(m_e*m_e)) - 1 + values_rad[0]/E0)/(values_rad[0]*(0.5*(1 + pow(1 - delta/E0, 2))*log(Q2/(m_e*m_e)) - 1 + delta/E0));
	}
	
	if(Uni < rho/values_rad[1])
	{
		values_rad[0] = omega;
		values_rad[1] = rho;
		return omega;
	}
	else
	{
		return values_rad[0];
	}
}

double Gen_omega_fin(const double& W, const double& Q2)
{
	double nu = (W*W + Q2 - m_p*m_p)/(2*m_p);
	double Ep = E0 - nu;
	double omega_max, rho(1.0), Uni(2.0), omega;
	
	if(channel)
	{
		omega_max = (2*m_p*E0 - 2*m_p*m_pi0 - m_pi0*m_pi0)/(2*m_p + Q2/Ep) - Ep;
	} else
	{
		omega_max = (2*m_p*E0 - 2*m_p*m_pip - m_pip*m_pip)/(2*m_p + Q2/Ep) - Ep;
	}	
	
	omega = fRand(delta, omega_max);
	rho = delta*(0.5*(1 + pow((E0 - nu)/(E0 - nu + omega), 2))*log(Q2/(m_e*m_e)) - (E0 - nu)/(E0 - nu + omega))/(omega*(0.5*(1 + pow((E0 - nu)/(E0 - nu + delta), 2))*log(Q2/(m_e*m_e)) - (E0 - nu)/(E0 - nu + delta)));
	Uni = fRand(0.0, 1.0);
	
	if(values_rad[2] > omega_max)
	{
		values_rad[2] = fRand(delta, omega_max);
		values_rad[3] = delta*(0.5*(1 + pow((E0 - nu)/(E0 - nu + values_rad[2]), 2))*log(Q2/(m_e*m_e)) - (E0 - nu)/(E0 - nu + values_rad[2]))/(values_rad[2]*(0.5*(1 + pow((E0 - nu)/(E0 - nu + delta), 2))*log(Q2/(m_e*m_e)) - (E0 - nu)/(E0 - nu + delta)));
	
	}
	
	if(Uni < rho/values_rad[3])
	{
		values_rad[2] = omega;
		values_rad[3] = rho;
		return omega;
	}
	else
	{
		return values_rad[2];
	}
	
	
	return omega;
}

void add_hadrons(const double& W_, const double& Q2_, const double& theta, const double& phi, Event& bob)
{
	double Epi, Ep, p, ang1, ang2; Particle proxy;
	
	if(channel)
	{
		Epi = (W_*W_ + m_pi0*m_pi0 - m_p*m_p)/(2*W_);
		Ep = (W_*W_ + m_p*m_p - m_pi0*m_pi0)/(2*W_);
		p = sqrt(Epi*Epi - m_pi0*m_pi0);
		
		proxy.mass = m_p; 
		proxy.id = 2212;
		(proxy.p).SetPxPyPzE(-p*cos(phi)*sin(theta), -p*sin(phi)*sin(theta), -p*cos(theta), Ep);
		bob.add_particle(proxy);
		
		TVector3 beta;
		ang1 = fRand(0, M_PI);
		ang2 = fRand(0, 2*M_PI);
		
		beta.SetXYZ( 0., 0., p/Epi);
		
		proxy.mass = 0; 
		proxy.id = 22;
		(proxy.p).SetPxPyPzE(m_pi0*cos(ang2)*sin(ang1)/2 ,m_pi0*sin(ang2)*sin(ang1)/2 ,m_pi0*cos(ang1)/2 ,m_pi0/2);
		(proxy.p).Boost(beta);
		(proxy.p).RotateY(theta);
		(proxy.p).RotateZ(phi);
		bob.add_particle(proxy);

		(proxy.p).SetPxPyPzE(-m_pi0*cos(ang2)*sin(ang1)/2 ,-m_pi0*sin(ang2)*sin(ang1)/2 ,-m_pi0*cos(ang1)/2 ,m_pi0/2);
		(proxy.p).Boost(beta);
		(proxy.p).RotateY(theta);
		(proxy.p).RotateZ(phi);
		bob.add_particle(proxy);	
	}
	else
	{
		Epi = (W_*W_ + m_pip*m_pip - m_n*m_n)/(2*W_);
		Ep = (W_*W_ + m_n*m_n - m_pip*m_pip)/(2*W_);
		p = sqrt(Epi*Epi - m_pip*m_pip);
		
		proxy.mass = m_n; 
		proxy.id = 2112;
		(proxy.p).SetPxPyPzE(-p*cos(phi)*sin(theta), -p*sin(phi)*sin(theta), -p*cos(theta), Ep);
		bob.add_particle(proxy);
		
		proxy.mass = m_pip; 
		proxy.id = 211;
		(proxy.p).SetPxPyPzE(p*cos(phi)*sin(theta), p*sin(phi)*sin(theta), p*cos(theta), Epi);
		bob.add_particle(proxy);	
	}
}

void add_electron(const double& W, const double& Q2, Event& bob)
{
	Particle proxy;
	proxy.mass = m_e; 
	proxy.id = 11;
	(proxy.p).SetPxPyPzE(sqrt(E2(W)*E2(W) - m_e*m_e)*sin_2(W, Q2), 0, sqrt(E2(W)*E2(W) - m_e*m_e)*cos_2(W, Q2), E2(W));
	bob.add_particle(proxy);
}

void add_rad_photon(const double& Erad, const double& ang1, const double& ang2, Event& bob)
{
	Particle proxy;
	proxy.mass = 0; 
	proxy.id = 22;
	(proxy.p).SetPxPyPzE(Erad*cos(ang2)*sin(ang1), Erad*sin(ang2)*sin(ang1), Erad*cos(ang1), Erad);
	bob.add_particle(proxy);
}

void generate_particle(const int& k)
{ 
	double W, Q2, theta, phi, S, factor(1), nu, S_1;
	double Epi, Ep, p;
	double r1, r2, r3;
	double Erad(0.0), Q2_, W_;
	double ang1, ang2;
	double fff; int count(1); int ghost(0);
	
	if(N - k < 3)
	{
		count = N % 3;
	}
	
	Event bob; 
	
	S = nan("");
	
	while(std::isnan(S))
	{
		W = fRand(W_min, W_max); 
		Q2 = fRand(Q2_min, Q2_max); 
		theta = fRand(0, M_PI);
		phi = fRand(0, 2*M_PI);		
		S = Linear(W, Q2, theta, phi);	
	}
	
	W_ = W;
	Q2_ = Q2;
	nu = (W*W + Q2 - m_p*m_p)/(2*m_p);
	
	if(rad_corr)
	{
		r1 = R1(W, Q2); 
		r2 = R2(W, Q2); 
		r3 = R3(W, Q2); 
		
		factor = (r1 + r2 + r3)/Section_interp_int(W, Q2, E0); 
		S *= factor;
		
		bob.set_beam(E0, h);
		
		bob.set_cm_system(); // r1
		bob.set_coordinates(R, L); S_1 = S*r1/(r1 + r2 + r3);
		bob.set_section(S_1);
		add_hadrons(W_, Q2_, theta, phi, bob);
		add_electron(W, Q2, bob);
		bob.cm_to_lab(W, Q2);
		bob.print_lund(path);
		
		if(histogram)
		{
			h1 -> Fill(W_, Q2_, S_1);
			h3 -> Fill(W_, S_1);
			h4 -> Fill(Q2_, S_1); fff = cos(theta);
			h5 -> Fill(phi, fff, S_1);	
		}
		bob.clear_event();
		
		if(N - k < 3)
		{
			count--;
		}
		
		if(count > 0)
		{
			bob.set_cm_system(); // r2
			bob.set_coordinates(R, L); S_1 = S*r2/(r1 + r2 + r3);
			bob.set_section(S_1);			
			W_ = 0;
			
			while((W_ < m_n + m_pip or std::isnan(W_)) and ghost < 10) 
			{
				Erad = Gen_omega_init(W, Q2); 
				Q2_ = Q2*(E0 - Erad)/E0;
				W_ = sqrt(m_p*m_p - Q2_ + 2*(- Erad + nu)*m_p); ghost++;		
			}
			
			if(ghost == 10)
			{
				W_ = W;
				Q2_ = Q2; S_1 = S; bob.set_section(S_1);
			}
					
			add_hadrons(W_, Q2_, theta, phi, bob);
			add_electron(W, Q2, bob);
			bob.cm_to_lab(W, Q2);
			
			if(ghost < 10)
			{
				ang1 = fRand(0, M_PI*sqrt(m_e/E0)/2);
				ang2 = fRand(0, 2*M_PI);
				add_rad_photon(Erad, ang1, ang2, bob);
			}
			
			bob.print_lund(path);
			
			if(histogram)
			{
				h1 -> Fill(W_, Q2_, S_1);
				h3 -> Fill(W_, S_1);
				h4 -> Fill(Q2_, S_1); fff = cos(theta);
				h5 -> Fill(phi, fff, S_1);	
			}
			bob.clear_event(); ghost = 0;
			
			if(N - k < 3)
			{
				count--;
			}
		}		
		
		if(count > 0)
		{
			bob.set_cm_system(); // r3
			bob.set_coordinates(R, L); S_1 = S*r3/(r1 + r2 + r3);
			bob.set_section(S_1);
			
			W_ = 0; 
			
			while((W_ < m_n + m_pip or std::isnan(W_)) and ghost < 10) 
			{
				Erad = Gen_omega_fin(W, Q2);
				Q2_ = Q2*(E0 - nu - Erad)/(E0 - nu);
				W_ = sqrt(m_p*m_p - Q2_ + 2*(nu + Erad)*m_p); ghost++;		
			}
			
			if(ghost == 10)
			{
				W_ = W;
				Q2_ = Q2; S_1 = S; bob.set_section(S_1);
			}
					
			add_hadrons(W_, Q2_, theta, phi, bob);
			add_electron(W, Q2, bob);
			bob.cm_to_lab(W, Q2);
			
			if(ghost < 10)
			{
				ang1 = acos(1 - Q2/(2*E0*(E0 - nu))) + fRand(-M_PI*sqrt(m_e/(E0 - nu))/2, M_PI*sqrt(m_e/(E0 - nu))/2);
				ang2 = fRand(0, 2*M_PI);
				add_rad_photon(Erad, ang1, ang2, bob);			
			}

			bob.print_lund(path);
			
			if(histogram)
			{
				h1 -> Fill(W_, Q2_, S_1);
				h3 -> Fill(W_, S_1);
				h4 -> Fill(Q2_, S_1); fff = cos(theta);
				h5 -> Fill(phi, fff, S_1);	
			}
		}				
	} else
	{
		bob.set_beam(E0, h);		
		bob.set_cm_system(); // Born
		bob.set_coordinates(R, L); 
		bob.set_section(S);
		add_hadrons(W, Q2, theta, phi, bob);
		add_electron(W, Q2, bob);
		bob.cm_to_lab(W, Q2);
		bob.print_lund(path);
		
		if(histogram)
		{
			h1 -> Fill(W, Q2, S);
			h3 -> Fill(W, S);
			h4 -> Fill(Q2, S); fff = cos(theta);
			h5 -> Fill(phi, fff, S);	
		}
	}
}





