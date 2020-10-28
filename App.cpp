#include <iostream> 
#include <fstream>
#include <vector>
#include <string.h>
#include <sstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include <complex>
#include <ctime>

using namespace std;

auto c1 = new TCanvas("c1", "Histogram", 1280, 1080); 
  
TH2* h1 = new TH2F("h1", "Histogram (W,Q^{2})_2d", 186, 1.08, 2, 440, 0, 10);
TH1F* h3 = new TH1F("h3", "Histogram W", 186, 1.08, 2);
TH1F* h4 = new TH1F("h4", "Histogram Q^{2}", 440, 0, 10);

vector<double> Settings(5); vector<int> Settings_mode(5); int polarization(0), weight_mode(0); 
double length(0), Radius_c(0), Q2_degree_extr(0);	



void Reading(string Path,vector<vector<double>>& V2, vector<string>& V3)
{
	string lane, line, word; stringstream ss;
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
				getline(File,lane);
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
	
			V2.push_back(Numbers);
			
			Numbers.clear();

			ss.clear();


		}
	}

	string str = ",";
	str.insert(0,lane);
	lane = str;
	string delim(",");
	size_t prev = 0;
	size_t next;
	size_t delta = delim.length();

	while(( next = lane.find( delim, prev ) ) != string::npos)
	{
		string tmp = lane.substr( prev, next-prev );
   		V3.push_back( lane.substr( prev, next-prev ) );
   		prev = next + delta;
	}
	File.close();
}

double fRand(const double& fMin, const double& fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void PrintVectorS(vector<string>& V3)
{
	for (vector<string>::iterator it=V3.begin();it!=V3.end();it++)
	{
		cout << *it << "\t\t";
	}
	cout << endl;
}

double mmin(double& x)
{
	return (x < 1) ? x:1;
}

bool attempt(double& p)
{
	double index = fRand(0, 1);
	return (index < p) ? true : false;
}

void Print(vector<complex<double>>& V3)
{
	for (vector<complex<double>>::iterator it=V3.begin();it!=V3.end();it++)
	{
		if(imag(*it) < 0)
		{
			cout << real(*it) << imag(*it) <<  "i\t\t";
		}
		else
		{
			cout << real(*it) << " + " << imag(*it) <<  "i\t\t";
		}
		
	}
	cout << endl;
}


void PrintVectorD(const vector<double>& V3)
{
	for (const auto& i : V3)
	{
		cout << i << "\t"; 
	}
	cout << endl;
	cout << "\n\n\n\n\nCol-vo:" << V3.size() << endl; 
}

void PrintBiggy(const vector<vector<double>>& V)
{
	int y(0);
	for ( int i = 0; i < V.size(); i++)
	{
		for( int j = 0; j < V[i].size(); j++)
		{
			cout << V[i][j] << "\t";
		}
		cout << endl; y++;
	}
	
	cout << "Col-vo strok: " << y << endl;
}

double P(const int& der,const int& n, double& theta) //legendre polynomials P(cos)n
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

double Sections(vector<double>& info,const int& t, double& W, double& Q2,const int& phi,const double& E0)
{
	double theta(t*M_PI/180), mp(0.93827), mpi(0.13498), nu; //radians 
	vector<complex<double>> Mp, Mm, Ep, Em, Lp, Lm; complex<double> Value; double Re, Im, eps, S, Gamma_flux(0);
	complex<double> F1(0,0), F2(0,0), F3(0,0), F4(0,0), F5(0,0), F6(0,0);
	double Rt, Rl, Rtl, Rtt, Rtl2, Rtt2;

	if(Settings_mode[3] == 2)
	{
		mpi = 0.13957; mp = 0.93957;
	} 

	double C,CC, L(0.14817); 

	double Epi = (W*W + mpi*mpi - mp*mp)/(2*W);
	double Ppi = sqrt(Epi*Epi - mpi*mpi);
	
	nu =  (W*W + Q2 - mp*mp)/(2*mp);
	
	C = 2*W*Ppi/sqrt((pow((W*W - mp*mp), 2) + Q2*(2*mp*mp + Q2))); CC = sqrt(Q2)/nu;
	for(int l = 0; l < 6; l++)
	{
		Re = info[50 + 2*l]; Im = info[51 + 2*l]; Value = complex<double>(Re,Im); Mp.push_back(Value);
		Re = info[62 + 2*l]; Im = info[63 + 2*l]; Value = complex<double>(Re,Im); Mm.push_back(Value);
		Re = info[26 + 2*l]; Im = info[27 + 2*l]; Value = complex<double>(Re,Im); Ep.push_back(Value);
		Re = info[38 + 2*l]; Im = info[39 + 2*l]; Value = complex<double>(Re,Im); Em.push_back(Value);
		Re = info[2 + 2*l]; Im = info[3 + 2*l]; Value = complex<double>(Re,Im); Lp.push_back(Value);
		Re = info[14 + 2*l]; Im = info[15 + 2*l]; Value = complex<double>(Re,Im); Lm.push_back(Value);
	}

	for(int i = 0; i < 6; i++)
	{
		double l = double(i);
		F1 = F1 + (l*Mp[l] + Ep[l])*P(1, l+1, theta) + ((l+1)*Mm[l] + Em[l])*P(1, l-1, theta); 
		F2 = F2 + ((l+1)*Mp[l] + l*Mm[l])*P(1, l, theta);
		F3 = F3 + (Ep[l] - Mp[l])*P(2, l+1, theta) + (Em[l] + Mm[l])*P(2, l-1, theta);
		F4 = F4 + (Mp[l] - Ep[l] - Mm[l] - Em[l])*P(2, l, theta);
		F5 = F5 + (l+1)*Lp[l]*P(1, l+1, theta) - l*Lm[l]*P(1, l-1, theta);
		F6 = F6 + (l*Lm[l] - (l+1)*Lp[l])*P(1, l, theta); }
		
		Rt = pow(L, 2)*(abs(F1)*abs(F1) + abs(F2)*abs(F2) + 0.5*pow(sin(theta), 2)*(abs(F3)*abs(F3) + abs(F4)*abs(F4)) - real(2*cos(theta)*conj(F1)*F2 - pow(sin(theta),2)*(conj(F1)*F4 + conj(F2)*F3 + cos(theta)*conj(F3)*F4)));
		Rl = CC*CC*pow(L, 2)*(abs(F5)*abs(F5) + abs(F6)*abs(F6) + 2*cos(theta)*real(conj(F5)*F6));
		Rtl = CC*pow(L, 2)*(-sin(theta)*real((conj(F2) + conj(F3) + cos(theta)*conj(F4))*F5 + (conj(F1) + conj(F4) + cos(theta)*conj(F3))*F6));
	Rtt = pow(L, 2)*pow(sin(theta), 2)*(0.5*(abs(F3)*abs(F3) + abs(F4)*abs(F4)) + real(conj(F1)*F4 + conj(F2)*F3 + cos(theta)*conj(F3)*F4));

		Rtl2 = CC*pow(L, 2)*(-sin(theta)*imag((conj(F2) + conj(F3) + cos(theta)*conj(F4))*F5 + (conj(F1) + conj(F4) + cos(theta)*conj(F3))*F6));

	eps = 1/(1 + 2*(nu*nu + Q2)/(4*(E0 - nu)*E0 - Q2));

	if(Q2 > 5)
	{
	S = pow(sqrt(5), Q2_degree_extr)*(C*Rt + eps*C*Rl + sqrt(2*eps*(1 + eps))*Rtl*cos(phi*M_PI/180)*C + eps*C*Rtt*cos(phi*M_PI/90) + polarization*C*sqrt(2*eps*(1 - eps))*Rtl2*sin(phi*M_PI/180))/pow(sqrt(Q2), Q2_degree_extr);
	} else 
	{
		S = C*Rt + eps*C*Rl + sqrt(2*eps*(1 + eps))*Rtl*cos(phi*M_PI/180)*C + eps*C*Rtt*cos(phi*M_PI/90) + polarization*C*sqrt(2*eps*(1 - eps))*Rtl2*sin(phi*M_PI/180);
	}

	Gamma_flux = W*(W*W - mp*mp)/(137*4*M_PI*mp*mp*E0*E0*(1 - eps)*Q2);

	S = Gamma_flux*S;

	Mp.clear(); Mm.clear(); Ep.clear(); Em.clear(); Lp.clear(); Lm.clear(); 

	return S;
}

vector<double> Coefficients_lin(const double& x1, const double& y1, const double& x2, const double& y2)
{
	vector<double> Abc(2);

	Abc[0] = (y1 - y2)/(x1 - x2); 
	Abc[1] = (x1*y2 - y1*x2)/(x1 - x2); 

	return Abc;
}

double Linear(vector<vector<double>>& V, const double& W,  const double& Q2, const double& theta,  const double& phi, const double& E0)
{
	double S; vector<double> ab;
	double W_min, W_max, Q2_min, Q2_max, S1, S2, S3, S4, S_1, S_2;
	vector<double> buff1, buff2, buff3, buff4;

	W_min = floor(W*100)/100;
	W_max = ceil(W*100)/100;
	Q2_min = floor(Q2*20)/20;
	Q2_max = ceil(Q2*20)/20;

	int i1, i2, i3, i4;			
	
	i1 = (W_min - 1.08)*100 + Q2_min*20*93;
	i2 = (W_min - 1.08)*100 + Q2_max*20*93;
	i3 = (W_max - 1.08)*100 + Q2_min*20*93;
	i4 = (W_max - 1.08)*100 + Q2_max*20*93;		

	for(int j = 0; j < V[i1].size(); j++)
	{
		buff1.push_back(V[i1][j]);
	}

	for(int j = 0; j < V[i2].size(); j++)
	{
		buff2.push_back(V[i2][j]);
	}

	for(int j = 0; j < V[i3].size(); j++)
	{
		buff3.push_back(V[i3][j]);
	}

	for(int j = 0; j < V[i4].size(); j++)
	{
		buff4.push_back(V[i4][j]);
	}

	S1 = Sections(buff1, theta, W_min, Q2_min, phi, E0);
	S2 = Sections(buff2, theta, W_min, Q2_max, phi, E0);
	S3 = Sections(buff3, theta, W_max, Q2_min, phi, E0);
	S4 = Sections(buff4, theta, W_max, Q2_max, phi, E0); 

	buff1.clear();
	buff2.clear();
	buff3.clear();
	buff4.clear();

	ab = Coefficients_lin(W_min, S1, W_max, S3);
	S_1 = ab[0]*W + ab[1]; ab.clear();
	ab = Coefficients_lin(W_min, S2, W_max, S4);
	S_2 = ab[0]*W + ab[1]; ab.clear();
	ab = Coefficients_lin(Q2_min, S_1, Q2_max, S_2);
	S = ab[0]*Q2 + ab[1]; ab.clear();	 
	
//S = ((S4 - S3 - S2 + S1)*(Q2 - Q2_min)/(Q2_max - Q2_min) + S3 - S1)*(W - W_min)/(W_max - W_min) + (S2 - S1)*(Q2 - Q2_min)/(Q2_max - Q2_min) + S1;

	return S;
}

double Section_transition(vector<vector<double>>& V, double& W, double& Q2, const double& theta,  const double& phi, const double& E0)
{
	double Si; int index; vector<double> buff;

	index = (W - 1.08)*100 + Q2*20*93;

	for(int j = 0; j < V[index].size(); j++)
	{
		buff.push_back(V[index][j]);
	}

	Si = Sections(buff, theta, W, Q2, phi, E0);

	buff.clear();

	return Si;
}

vector<double> Coefficients(const double& x1, const double& y1, const double& x2, const double& y2, const double& x3, const double& y3)
{
	vector<double> Abc(3); double det;

	det = (x2 - x3)*(x1*x1 - x1*(x2 + x3) + x2*x3);

	Abc[0] = (y1*(x2 - x3) + y2*(x3 - x1) + y3*(x1 - x2))/det; 
	Abc[1] = (y1*(x3*x3 - x2*x2) + y2*(x1*x1 - x3*x3) + y3*(x2*x2 - x1*x1))/det;  
	Abc[2] = (y1*x2*x3*(x2 - x3) + y2*x1*x3*(x3 - x1) + y3*x1*x2*(x1 - x2))/det; 

	return Abc;
}

double Quadratic(vector<vector<double>>& V, const double& W,  const double& Q2, const double& theta,  const double& phi, const double& E0)
{
	double S, W_1, W_2, W_3, Q2_1, Q2_2, Q2_3, S1, S2, S3;
	vector<double> S_i(9); vector<double> abc;

	W_1 = floor(W*100)/100;
	W_2 = ceil(W*100)/100;
	Q2_1 = floor(Q2*20)/20;
	Q2_2 = ceil(Q2*20)/20;

	if(Q2_2 == 11)
	{
		Q2_3 = 11;
		Q2_2 -= 0.05;		
		Q2_1 -= 0.05;
	} else {Q2_3 = Q2_2 + 0.05;}

	if(W_2 == 2)
	{
		W_3 = 2;
		W_2 -= 0.01;		
		W_1 -= 0.01;
	} else {W_3 = W_2 + 0.01;} 

	S_i[0] = Section_transition(V, W_1, Q2_1, theta, phi, E0);
	S_i[1] = Section_transition(V, W_2, Q2_1, theta, phi, E0);
	S_i[2] = Section_transition(V, W_3, Q2_1, theta, phi, E0);

	S_i[3] = Section_transition(V, W_1, Q2_2, theta, phi, E0);
	S_i[4] = Section_transition(V, W_2, Q2_2, theta, phi, E0);
	S_i[5] = Section_transition(V, W_3, Q2_2, theta, phi, E0);

	S_i[6] = Section_transition(V, W_1, Q2_3, theta, phi, E0);
	S_i[7] = Section_transition(V, W_2, Q2_3, theta, phi, E0);
	S_i[8] = Section_transition(V, W_3, Q2_3, theta, phi, E0); 

	abc = Coefficients(W_1, S_i[0], W_2, S_i[1], W_3, S_i[2]);
	S1 = abc[0]*W*W + abc[1]*W + abc[2]; abc.clear();

	abc = Coefficients(W_1, S_i[3], W_2, S_i[4], W_3, S_i[5]);
	S2 = abc[0]*W*W + abc[1]*W + abc[2]; abc.clear();

	abc = Coefficients(W_1, S_i[6], W_2, S_i[7], W_3, S_i[8]);
	S3 = abc[0]*W*W + abc[1]*W + abc[2]; abc.clear();

	abc = Coefficients(Q2_1, S1, Q2_2, S2, Q2_3, S3);
	S = abc[0]*Q2*Q2 + abc[1]*Q2 + abc[2]; abc.clear();	
	
	S_i.clear(); 

	return S;
}

void All(vector<vector<double>>& V, const int& FileNumber)
{
	double W, Q2, theta, phi, S(0), S1, P(0), Ep, Epi, p, mp(0.93827), mpi(0.13498), nu, ang1, ang2, weight(0), z, x, y,W_1, Q2_1, theta_1, phi_1;
	int v(0); bool trigger(true), decision; double ratio;

	double E0, W_min_d, W_max_d, Q2_min_d, Q2_max_d; int N, Hist, FileNumberAll, decay_m;
	E0 = Settings[0];		N = Settings_mode[0];
	W_min_d = Settings[1];		Hist = Settings_mode[1];
	W_max_d = Settings[2];		FileNumberAll = Settings_mode[2];
	Q2_min_d = Settings[3];		decay_m = Settings_mode[3];
	Q2_max_d = Settings[4];
	
	
	TLorentzVector e, adron, meson, gamma1, gamma2; 

	TVector3 beta;

	char FileName[100];

	if(Settings_mode[3] == 1 or Settings_mode[3] == 0)
	{
		sprintf(FileName,"pi0p_W_%g_%g_Q2_%g_%g_(%i).lund", W_min_d, W_max_d, Q2_min_d, Q2_max_d, FileNumber);
	} else 
	{
		sprintf(FileName,"pin_W_%g_%g_Q2_%g_%g_(%i).lund", W_min_d, W_max_d, Q2_min_d, Q2_max_d, FileNumber);
		mpi = 0.13957; mp = 0.93957;
	}	

	cout << "\n\tThe name of output file will be " << FileName << endl;

	ofstream File;

	File.open(FileName);

 

	while(v < N)
	{
		if(trigger)
		{
			theta = fRand(0, 180); 
			phi = fRand(0, 360); 
			W = fRand(W_min_d, W_max_d);
			Q2 = fRand(Q2_min_d, Q2_max_d);

			if(Settings_mode[4] == 0)
			{
				S = Linear(V, W, Q2, theta, phi, E0);
			} else if(Settings_mode[4] == 1)
			{
				S = Quadratic(V, W, Q2, theta, phi, E0);
			}
		}

		if(weight_mode == 0)
		{
			P = 0.00000000001*(rand() % 100000000000); 
		} else if(weight_mode == 1) 
		{
			weight = S; P = 0;
		} else if(weight_mode == 2)
		{
			trigger = false;
			theta_1 = fRand(0, 180); 
			phi_1 = fRand(0, 360); 
			W_1 = fRand(W_min_d, W_max_d);
			Q2_1 = fRand(Q2_min_d, Q2_max_d);

			if(Settings_mode[4] == 0)
			{
				S1 = Linear(V, W_1, Q2_1, theta_1, phi_1, E0);
			} else if(Settings_mode[4] == 1)
			{
				S1 = Quadratic(V, W_1, Q2_1, theta_1, phi_1, E0);
			}

			ratio = S1/S;

			P = mmin(ratio);

			decision = attempt(P);

			if(decision)
			{
				theta = theta_1; 
				phi = phi_1; 
				W = W_1;
				Q2 = Q2_1;
				S = S1;
			}
		}

		if((S/60) >= P) //
		{
			v++; 

			if(length != 0)
			{
				z = 0.5*fRand(-length, length);
			} else { z = 0;}

			if(Radius_c != 0)
			{
				x = Radius_c*2; y = Radius_c*2;
				while(x*x + y*y > Radius_c*Radius_c)
				{
					x = fRand(-Radius_c, Radius_c);
					y = fRand(-Radius_c, Radius_c);
				}
			} else { x = 0; y = 0;}
			
			if(Hist == 1)
			{
				h1 -> Fill(W, Q2);
				h3 -> Fill(W);
				h4 -> Fill(Q2);
			}
	
			if((100*v)%N == 0)
			{
				cout << 100*v/N << "%\tFile " << FileNumber << " of " << FileNumberAll << endl;
			}			

			Epi = (W*W + mpi*mpi - mp*mp)/(2*W);
			Ep = (W*W + mp*mp - mpi*mpi)/(2*W);
			p = sqrt(Epi*Epi - mpi*mpi);

			nu =  (W*W + Q2 - mp*mp)/(2*mp);

			ang1 = fRand(0, 180); 
			ang2 = fRand(0, 360);
			
	meson.SetPxPyPzE(p*cos(phi*M_PI/180)*sin(theta*M_PI/180) ,p*sin(phi*M_PI/180)*sin(theta*M_PI/180) ,p*cos(theta*M_PI/180) , Epi);
	adron.SetPxPyPzE(-p*cos(phi*M_PI/180)*sin(theta*M_PI/180) ,-p*sin(phi*M_PI/180)*sin(theta*M_PI/180) ,-p*cos(theta*M_PI/180) , Ep);
	gamma1.SetPxPyPzE(mpi*cos(ang2*M_PI/180)*sin(ang1*M_PI/180)/2 ,mpi*sin(ang2*M_PI/180)*sin(ang1*M_PI/180)/2 ,mpi*cos(ang1*M_PI/180)/2 ,mpi/2);
gamma2.SetPxPyPzE(-mpi*cos(ang2*M_PI/180)*sin(ang1*M_PI/180)/2 ,-mpi*sin(ang2*M_PI/180)*sin(ang1*M_PI/180)/2 ,-mpi*cos(ang1*M_PI/180)/2 ,mpi/2);

			beta.SetXYZ( 0., 0., p/Epi);
			gamma1.Boost(beta);
			gamma2.Boost(beta);

			gamma1.RotateY(theta*M_PI/180);
			gamma2.RotateY(theta*M_PI/180);
			gamma1.RotateZ(phi*M_PI/180); 
			gamma2.RotateZ(phi*M_PI/180);			

			beta.SetXYZ( 0., 0., sqrt(nu*nu + Q2)/(nu + mp));

			adron.Boost(beta);
			meson.Boost(beta);
			gamma1.Boost(beta);
			gamma2.Boost(beta);

      			adron.RotateY(acos((Q2 + 2*E0*nu)/(2*E0*sqrt(nu*nu + Q2))));
			meson.RotateY(acos((Q2 + 2*E0*nu)/(2*E0*sqrt(nu*nu + Q2))));
			gamma1.RotateY(acos((Q2 + 2*E0*nu)/(2*E0*sqrt(nu*nu + Q2))));
			gamma2.RotateY(acos((Q2 + 2*E0*nu)/(2*E0*sqrt(nu*nu + Q2))));

			e.SetPxPyPzE( (E0 - nu)*sqrt(1 - pow(1 - Q2/(2*E0*(E0 - nu)),2)), 0, (E0 - nu)*(1 - Q2/(2*E0*(E0 - nu))), E0 - nu);

			ang2 = fRand(0, 360);
	
			adron.RotateZ(ang2*M_PI/180);
			meson.RotateZ(ang2*M_PI/180);
			gamma1.RotateZ(ang2*M_PI/180);
			gamma2.RotateZ(ang2*M_PI/180);
			e.RotateZ(ang2*M_PI/180);

			if(decay_m == 1)
			{
	File << 4 << " " << 0 << " " << 0 << " " << 0 << " " << polarization << " " << 11 << " " << E0 << " " << 2212 << " " << 0 << " " << weight << endl;
File << 0 << " " << 0 << " " << 0 << " " << 11 << " " << 0 << " " << 0 << " " <<  e.Px() << " " << e.Py() << " " << e.Pz() << " " << e.E() << " " << 0.0005 << " " << x << " " << y << " " << z << endl;
File << 0 << " " << 0 << " " << 0 << " " << 2212 << " " << 0 << " " << 0 << " " <<  adron.Px() << " " << adron.Py() << " " << adron.Pz() << " " << adron.E() << " " << mp << " " << x << " " << y << " " << z << endl;
File << 0 << " " << 0 << " " << 0 << " " << 22 << " " << 0 << " " << 0 << " " <<  gamma1.Px() << " " << gamma1.Py() << " " << gamma1.Pz() << " " << gamma1.E() << " " << 0 << " " << x << " " << y << " " << z << endl;
File << 0 << " " << 0 << " " << 0 << " " << 22 << " " << 0 << " " << 0 << " " <<  gamma2.Px() << " " << gamma2.Py() << " " << gamma2.Pz() << " " << gamma2.E() << " " << 0 << " " << x << " " << y << " " << z << endl;
			} else if(decay_m == 0)
			{
	File << 3 << " " << 0 << " " << 0 << " " << 0 << " " << polarization << " " << 11 << " " << E0 << " " << 2212 << " " << 0 << " " << weight << endl;
			
File << 0 << " " << 0 << " " << 0 << " " << 11 << " " << 0 << " " << 0 << " " <<  e.Px() << " " << e.Py() << " " << e.Pz() << " " << e.E() << " " << 0.0005 << " " << x << " " << y << " " << z << endl;
File << 0 << " " << 0 << " " << 0 << " " << 2212 << " " << 0 << " " << 0 << " " <<  adron.Px() << " " << adron.Py() << " " << adron.Pz() << " " << adron.E() << " " << mp << " " << x << " " << y << " " << z << endl;
File << 0 << " " << 0 << " " << 0 << " " << 111 << " " << 0 << " " << 0 << " " <<  meson.Px() << " " << meson.Py() << " " << meson.Pz() << " " << meson.E() << " " << mpi << " " << x << " " << y << " " << z << endl;
			} else if(decay_m == 2)
			{
	File << 3 << " " << 0 << " " << 0 << " " << 0 << " " << polarization << " " << 11 << " " << E0 << " " << 2212 << " " << 0 << " " << weight << endl;
			
File << 0 << " " << 0 << " " << 0 << " " << 11 << " " << 0 << " " << 0 << " " <<  e.Px() << " " << e.Py() << " " << e.Pz() << " " << e.E() << " " << 0.0005 << " " << x << " " << y << " " << z << endl;
File << 0 << " " << 0 << " " << 0 << " " << 2112 << " " << 0 << " " << 0 << " " <<  adron.Px() << " " << adron.Py() << " " << adron.Pz() << " " << adron.E() << " " << mp << " " << x << " " << y << " " << z << endl;
File << 0 << " " << 0 << " " << 0 << " " << 211 << " " << 0 << " " << 0 << " " <<  meson.Px() << " " << meson.Py() << " " << meson.Pz() << " " << meson.E() << " " << mpi << " " << x << " " << y << " " << z << endl;
			}			
		}
	}

	cout << "\tCompleted 100%" << endl;

	File.close();	
}
	
int main(int argc, char **argv)
{
	vector<vector<double>> Biggy; 
	vector<string> VecShap; 

	cout << " ------------------------------------------------------------------- " << endl;
	cout << "| Welcome to event builder for Pi0p and Pi+n channels of meson      | \n| electroproduction reaction!                                       |       \n|                                                                   |\n|     Authors: Davydov M. - MSU, Physics dep.                       |\n|              Isupov E.  - MSU, SINP                               |\n|                                                   Version 6.0     |\n| https://github.com/Maksaska/pi0p-pin-generator                    |\n ------------------------------------------------------------------- " << endl;
	
	cout << endl;	

	for(double& i : Settings) 
	{
		cin >> i;
	}

	cin >> polarization;

	cin >> Q2_degree_extr;

	cin >> length;

	cin >> Radius_c;

	for(int& i : Settings_mode) 
	{
		cin >> i;
	}

	cin >> weight_mode;

	string FileName = (Settings_mode[3] == 0 or Settings_mode[3] == 1) ? "pi0p_e.csv":"pin_e.csv";

	c1 -> Divide(2,2);

	if(Settings[0] < 0)
	{
		cout << "The beam energy is below zero!" << endl;
		return 0;
	}
 
	if(Settings[1] < 1.08 or Settings[1] > 2)
	{
		cout << "W_min is wrong!\nChoose the other one!\n" << endl;
		cout << "Choose the kinematic area of (W, Q^2) values in the range W: 1.08 - 2.0 GeV , Q2: 0 - 10 GeV^2" << endl;
		return 0;
	}
 
	if(Settings[2] < 1.08 or Settings[2] > 2)
	{
		cout << "W_max is wrong!\nChoose the other one!\n" << endl;
		cout << "Choose the kinematic area of (W, Q^2) values in the range W: 1.08 - 2.0 GeV , Q2: 0 - 10 GeV^2" << endl;

		if(Settings[2] < Settings[1])
		{
			cout << "Choose the bigger one!    W_max>" << Settings[1] << endl;
		}
		return 0;
	}

	if(Settings[2] < Settings[1])
	{
		cout << "W_max is wrong!\nChoose the bigger one!    W_max>" << Settings[1] << endl;
		return 0;
	}
 
	if(Settings[3] < 0 or Settings[3] > 10)
	{
		cout << "Q2_min is wrong!\nChoose the other one!\n" << endl;
		cout << "Choose the kinematic area of (W, Q^2) values in the range W: 1.08 - 2.0 GeV , Q2: 0 - 10 GeV^2" << endl;
		return 0;
	}	
 
	if(Settings[4] < 0 or Settings[4] > 10)
	{
		cout << "Q2_max is wrong!\nChoose the other one!\n" << endl;
		cout << "Choose the kinematic area of (W, Q^2) values in the range W: 1.08 - 2.0 GeV , Q2: 0 - 10 GeV^2" << endl;

		if(Settings[4] < Settings[3])
		{
			cout << "Choose the bigger one!    Q2_max >" << Settings[3] << endl;
		}
		return 0;
	}

	if(Settings[4] < Settings[3])
	{
		cout << "Choose the bigger one!    Q2_max >" << Settings[3] << endl;
		return 0;
	}

	if(Settings_mode[0] < 0)
	{
		cout << "The number of generated particles can't be negative!" << endl;
		return 0;
	}	

	cout << "The beam energy is E = " << Settings[0] << " GeV" << endl;
	cout << "Choosen kinematic area of (W, Q^2) values is W: " << Settings[1] << " - " << Settings[2] << " GeV , Q2: " << Settings[3] << " - " << Settings[4] << " GeV^2" << endl;
	cout << "Total amount of files: " << Settings_mode[2] << "\nHistogram need: ";
	if(Settings_mode[1] == 1){cout << "Yes" << endl;} else {cout << "No" << endl;}
	cout << "Interpolation mode: ";
	if(Settings_mode[4] == 1){cout << "Biquadratic" << endl;} else {cout << "Bilinear" << endl;}
	if(Settings_mode[3] == 1 or Settings_mode[3] == 0)
	{
		cout << "Channel: pi0p" << endl;
		cout << "Pi^0 meson decay: ";
		if(Settings_mode[3] == 1){cout << "Yes" << endl;} else {cout << "No" << endl;}
	} else 
	{
		cout << "Channel: pin" << endl;
	}	
	
	cout << "Stand by...\n" << endl;	

	Reading(FileName,Biggy,VecShap);

	srand(time(NULL));

	for(int yy = 1; yy <= Settings_mode[2]; yy++)
	{
		All(Biggy, yy);
	}

	if(Settings_mode[1] == 1)
	{
		c1->cd(1);
		h1->Draw("COL");

		c1->cd(2);
		h1->Draw("SURF2");

		c1->cd(3);
		h3->Draw();

		c1->cd(4);
		h4->Draw();

		char FileName1[100];

		sprintf(FileName1,"Histograms_reg%i_degree_%g_W_%g_%g_Q2_%g_%g.jpeg", Settings_mode[3], Q2_degree_extr, Settings[1], Settings[2], Settings[3], Settings[4]);

		c1 -> Print(FileName1);
	}

	Biggy.clear(); 
	VecShap.clear(); 
	Settings.clear(); Settings_mode.clear();

	return 0;
}
