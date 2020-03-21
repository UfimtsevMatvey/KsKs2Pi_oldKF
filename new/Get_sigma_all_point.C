#include <iostream>
#include <fstream>
#include <string>
using namespace std;
void NZero(int, int, double, double, double, double, double*, double*);
double error_div(double, double, double, double, double, double);
void Get_sigma_all_point()
{
	double Enew[20];
	const int N_point = 20;
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
	L[20] = 0;
	LE[20] = 0;

    double b_summ_30 = 0, b_summ_40 = 0;
    double a_summ_30 = 0, a_summ_40 = 0;

    double b_e_summ_30 = 0, b_e_summ_40 = 0;
    double a_e_summ_30 = 0, a_e_summ_40 = 0;

    double MC_ef_fit1[20];
	double MC_ef_fit2[20];
	
	double MC_ef_fit1_errors[20];
	double MC_ef_fit2_errors[20];

	double E0[] = { 987.5,978.,975.,962.5,955.,951.1,950.,945.,936.,925.,912.5,900.
					,887.5,875,862.5,850.
					,837.5,825.,812.5,800. };
	double E[] = { 2000., 1900., 1800.,1700.,1600. };
	std::ifstream lum;
	lum.open("lum2019.dat", ios::in);
	std::ifstream NumbE;
	NumbE.open("Number_events.dat", ios::in);
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
    TGraphErrors* gr1 = new TGraphErrors("Good_MC_efficiency.dat", "%lg %lg %lg %*lg %*lg");
	TGraphErrors* gr2 = new TGraphErrors("Good_MC_efficiency.dat", "%lg %*lg %*lg %lg %lg");
	TFitResultPtr r1 = gr1->Fit("pol1","S");
	TFitResultPtr r2 = gr2->Fit("pol1","S");
    Double_t chi2   = r1->Chi2();
	Double_t par0   = r1->Value(0);
	Double_t par1   = r1->Value(1);
	Double_t err0   = r1->Error(0);
	Double_t err1   = r1->Error(1);

    for(int u = 0; u < 20; u++){
		MC_ef_fit1[u] = par0 + par1*E0[u]*2;
		MC_ef_fit1_errors[u] = err0 + err1*E0[u]*2;
	}
	chi2   = r2->Chi2();
	par0   = r2->Value(0);
	par1   = r2->Value(1);
	err0   = r2->Error(0);
	err1   = r2->Error(1);
	for(int u = 0; u < 20; u++){
		MC_ef_fit2[u] = par0 + par1*E0[u]*2;
		MC_ef_fit2_errors[u] = err0 + err1*E0[u]*2;
	}

    for (int i = 0; i < 5; i++)
    {
        B_file >> NULL_V >> b_30 >> b_e_30 >> b_40 >> b_e_40 >> N1_30 >> N1_40 >> N2_30 >> N2_40;
		A_file >> NULL_V >> NULL_V >> NULL_V >> NULL_V >> a_30 >> a_40 >> a_e_30 >> a_e_40 >> NULL_V;
        
        b_summ_30 = b_summ_30 + b_30;
        b_summ_40 = b_summ_40 + b_40;
        a_summ_30 = a_summ_30 + a_30;
        a_summ_40 = a_summ_40 + a_40;

        b_e_summ_30 = b_e_summ_30 + b_e_30;
        b_e_summ_40 = b_e_summ_40 + b_e_40;
        a_e_summ_30 = a_e_summ_30 + a_e_30;
        a_e_summ_40 = a_e_summ_40 + a_e_40;
    }
    b_summ_30 = b_summ_30/5.;
    b_summ_40 = b_summ_40/5.;
    a_summ_30 = a_summ_30/5.;
    a_summ_40 = a_summ_40/5.;
    
    b_e_summ_30 = 0;
    b_e_summ_40 = 0;
    a_e_summ_30 = 0;
    a_e_summ_30 = 0;
    /*b_e_summ_30 = b_e_summ_30/5;
    b_e_summ_40 = b_e_summ_40/5;
    a_e_summ_30 = a_e_summ_30/5;
    a_e_summ_30 = a_e_summ_40/5;*/

	for(int i = 0; i < N_point; i++){
		//B_mc >> NULL_V >> NULL_V >> NULL_V >> NULL_V >>NULL_V >> NULL_V >> NULL_V >> NULL_V >>b_30 >> b_e_30 >> b_40 >> b_e_40;
        NumbE >>  NULL_V >> NULL_V >> NULL_V >> NULL_V >> N1_30 >> N2_30 >> N1_40 >> N2_40 >> NULL_V >> NULL_V >> NULL_V >> NULL_V;
		cout << N1_30 << '\t'<< N2_30 << '\t'<< N1_40 << '\t'<< N2_40<<'\n';
		
		NZero(N1_30, N2_30, a_summ_30, a_e_summ_30, b_summ_30, b_e_summ_30,&N0_30, &N0_e_30);
		NZero(N1_40, N2_40, a_summ_40, a_e_summ_40, b_summ_40, b_e_summ_40,&N0_40, &N0_e_40);
		cout << a_summ_30 << '\t' << a_summ_40 << '\n';
		cout << b_summ_30 << '\t' << b_summ_40 << '\n';
		cout << N0_40 << '\t' << N0_30 << '\n';
		cout << L[i] << '\n';
		S_30 = (N0_30/MC_ef_fit1[i])/L[i];
		S_40 = (N0_40/MC_ef_fit2[i])/L[i];
		//S_e_30 = error_div(Eff_30, N0_30, N0_e_30, Eff_e_30)/LR[i];
		//S_e_40 = error_div(Eff_40, N0_40, N0_e_40, Eff_e_40)/LR[i];
		S_e_30 = error_div(MC_ef_fit1[i], N0_30, N0_e_30, MC_ef_fit1_errors[i],L[i], LE[i]);
		S_e_40 = error_div(MC_ef_fit2[i], N0_40, N0_e_40, MC_ef_fit2_errors[i],L[i], LE[i]);
		Sigma << E0[i] << '\t' << S_30 << '\t' << S_e_30<< '\t' << S_40 << '\t' << S_e_40 << '\t' << '\n';
	}
	B_mc.close();
	lum.close();
	A_file.close();
	Sigma.close();
	B_file.close();
	NumbE.close();
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