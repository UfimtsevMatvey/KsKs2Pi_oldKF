#include <iostream>
#include <fstream>
#include <string>
using namespace std;
void NZero(int, int, double, double, double, double, double*, double*);
double error_div(double, double, double, double, double, double);
void Sigma_osb()
{
	double Enew[20];
	const int N_point = 5;
	double N0_30, N0_40;
	double N0_e_30, N0_e_40;
	double N1_30, N1_40;
	double N2_30, N2_40;
	double a_30, a_40;
	double a_e_30, a_e_40;
	double b_30, b_40;
	double b_e_30, b_e_40;
	double Eff_30, Eff_40;
	double Eff_e_30, Eff_e_40;
	double NULL_V;
	double S_30, S_e_30;
	double S_40, S_e_40;
	double LU, LUEr, EnL;
	double LE[21], L[21];
	double E0[] = { 987.5,978.,975.,962.5,955.,951.1,950.,945.,936.,925.,912.5,900.
					,887.5,875,862.5,850.
					,837.5,825.,812.5,800. };
	double E[] = { 2000., 1900., 1800.,1700.,1600. };
	std::ifstream lum;
	lum.open("lum2019.dat", ios::in);
	std::ifstream Feff;
	Feff.open("Good_MC_efficiency.dat", ios::in);
	std::ifstream B_file;
	B_file.open("Coefficient_B_and_Nevents.dat", ios::in);
	std::ifstream A_file;
	A_file.open("Coefficient_A.dat", ios::in);
	std::ifstream B_mc;
	B_mc.open("Number_events_mc.dat", ios::in);
	std::ofstream Sigma;
	Sigma.open("Sigma_Obs.dat",ios::out);
	while (lum >> EnL >> LU >> LUEr)
	{
		for(int i = 0; i < 20; i++){
			if(EnL == E0[i]){
				Enew[i] = E0[i];
				L[i] = LU;
				LE[i] = LUEr;
			}
		}
	}
	L[20] = 0;
	LE[20] = 0;
	double LR[5];
	double LER[5];
	for(int i = 0; i < 5; i++){
		for(int j = 0; j < 5; j++){
			LR[j] = LR[j] + L[4*i+j];
			LER[j] = LER[j] + LE[4*i+j];
		}
	}
	for(int i = 0; i < N_point; i++){
		Feff >> NULL_V >> Eff_30 >> Eff_e_30 >> Eff_40 >> Eff_e_40;
		B_file >> NULL_V >> b_30 >> b_e_30 >> b_40 >> b_e_40 >> N1_30 >> N1_40 >> N2_30 >> N2_40;
		A_file >> NULL_V >> NULL_V >> NULL_V >> NULL_V >> a_30 >> a_40 >> a_e_30 >> a_e_40 >> NULL_V;
		//B_mc >> NULL_V >> NULL_V >> NULL_V >> NULL_V >>NULL_V >> NULL_V >> NULL_V >> NULL_V >>b_30 >> b_e_30 >> b_40 >> b_e_40;
		cout << N1_30 << '\t'<< N1_40 << '\t'<< N2_30 << '\t'<< N2_40<<'\n';
		
		NZero(N1_30, N2_30, a_30, a_e_30, b_30, b_e_30,&N0_30, &N0_e_30);
		NZero(N1_40, N2_40, a_40, a_e_40, b_40, b_e_40,&N0_40, &N0_e_40);
		cout << a_30 << '\t' << a_40 << '\n';
		cout << b_30 << '\t' << b_40 << '\n';
		cout << N0_40 << '\t' << N0_30 << '\n';
		cout << L[i*4] << '\n';
		S_30 = (N0_30/Eff_30)/LR[i];
		S_40 = (N0_40/Eff_40)/LR[i];
		//S_e_30 = error_div(Eff_30, N0_30, N0_e_30, Eff_e_30)/LR[i];
		//S_e_40 = error_div(Eff_40, N0_40, N0_e_40, Eff_e_40)/LR[i];
		S_e_30 = error_div(Eff_30, N0_30, N0_e_30, Eff_e_30,LR[i], LER[i]);
		S_e_40 = error_div(Eff_40, N0_40, N0_e_40, Eff_e_40,LR[i], LER[i]);
		Sigma << E[i]/2 << '\t' << S_30 << '\t' << S_e_30<< '\t' << S_40 << '\t' << S_e_40 << '\t' << '\n';
	}
	B_mc.close();
	lum.close();
	A_file.close();
	Sigma.close();
	B_file.close();
	Feff.close();
	
}

void NZero(int N1, int N2, double a, double a_e, double b, double b_e,double* N0, double* Error)
{
	double R1 = 0.0, R2 = 0.0;
	double N1_e = sqrt(N1);
	double N2_e = sqrt(N2);
	double c = b - a;
	R1 = (N1*b - N2)/c;
	R2 = pow(b*N1_e/c,2.)+pow(N2_e/c,2.) + 
		pow(b_e*(N1/c - R1/c),2.0) + pow(a_e*R1/c,2.);
	R2 = sqrt(R2);
	(*N0) = R1;
	(*Error) = R2;
}
double error_div(double y, double x, double dx, double dy, double z, double dz)
{
	double R = pow(dx/(y*z),2) + pow(x*dy/(y*y*z),2) + pow(x*dz/(z*z*y),2);
	return sqrt(R);
}