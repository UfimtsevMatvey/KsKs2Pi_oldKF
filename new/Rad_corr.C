const double alpha = 1.0/137.02;
const double m_e = 0.511;
const double Ks_mass = 497.614;
const double Pi_mass_Q = 139.57;
const double Pi_mass_0 = 0.279;
#define Mass_cut 637.184
double x_av;
double sigma_x;
double n;
double bound;
double a0;
double a1;
double a2;
double a3;
double a4;
double a5;
using namespace std;
double F(double x, double E)
{
    double R = 0.0;
    //double E = sqrt(s)/2.0;
	double s = 4*E*E;
    double L = log(s/(m_e*m_e));
    double B = 2 * alpha * (L - 1)/M_PI;
    double F1 = 0.0, F2 = 0.0;
    F1 = B * pow(x,B - 1)*(1 + (alpha/M_PI)*(M_PI*M_PI/3 - 0.5) + 0.75*B-B*B/24.0*(L/3.0 + 2*M_PI*M_PI-37.0/4.0)) - 
        B*(1-0.5*x) + 0.125*B*B*(4*(2-x)*log(1/x) + (1+3*(1-x)*(1-x)/x)*log(1/(1-x)) - 6 + x);
    if(x > 2 * m_e/E)
        F2 = pow(x-2*m_e/E,B)/(6*x)*pow(log(s*x*x/(m_e*m_e)) - 5.0/3.0,2) * (2-2*x+x*x + B/3.0*(log(s*x*x/(m_e*m_e)) - 5.0/3.0)) + 
            0.5*L*L*(2.0/3.0 * (1-pow(1-x,3))/(1-x) - (2-x)*log(1.0/(1-x)) + 0.5*x);
    else 
        F2 = 0.0;
    R = F1 + F2*alpha*alpha/(M_PI*M_PI);
    return R;
}
double pol_2(double E)
{
    double x = E;//sqrt(s)/2.0;
    if(E < a2) return 0;
    // return a0+a1*(x-Mass_cut)+a2*pow(x-Mass_cut,2); /*a0*exp(-(x-a1)*(x-a1)/(2*a2*a2));*///a0 + a1*x + a2*x*x;
    return a0*(1-exp(-(x-a2)*(x-a2)/(a1*a1)))*(a3+a4*x);
    //return a1/(1+exp(-a0*(x-637.184)))*((x - 637.184)+a2*pow(x-637.184,2)+ a3*pow(x-637.184,3));
    //return a0*(1-exp(-pow(x-637.184,2)/(a1*a1)))*(1 + a2*(x-637.184)+ a3*pow(x-637.184,2));
    //return a0*(1-exp(-pow(x-637.184,2)/pow(a1,2)))*(1 + a2*pow(x-637.184,4)+ a3*pow(x-637.184,2));
    //return a0*(1-exp(-pow(x-637.184,2)/pow(a1,2)));
}
double Simpson_methood(double (*Func)(double x, double E), double (*Sigma)(double p), double x_min, double x_max, int N, double E)
{
    double R = 0.0;
    double s = 4*E*E;
    double h = (x_max - x_min)/N;
    
    double* F_arr = (double*)calloc(N, sizeof(double));
    for(int i = 0; i < N; i++){
      *(F_arr + i) = Func(x_min + i*h,E)*Sigma(E*sqrt(1-x_min - i*h));
		//cout << "F = "<< Func(x_min + i*h,E)*Sigma(E*sqrt(1-x_min - i*h))<<"\n";
    }

    for (int i = 0; i < N; i++) { 
        if (i == 0 || i == N) 
            R += *(F_arr + i); 
        else if (i % 2 != 0) 
            R += 4 * (*(F_arr + i)); 
        else
            R += 2 * (*(F_arr + i)); 
		//cout << "R = " << R <<'\n';
    } 
    R = R * (h / 3); 
    free(F_arr);
    return R;
}
double S_trans(double t, double x_max, double x_min)
{
	return 0.5*(x_max + x_min) + 0.5*(x_max - x_min)*tanh(t);
}
double Simpson_methood_special(double (*Func)(double x, double E), double (*Sigma)(double p), double x_min, double x_max, int N, double E)
{
    double R = 0.0;
    double s = 4*E*E;
	double t = 0.0;
	double T = 19.0;
	double h = 2*T/N;
    double* F_arr = (double*)calloc(N, sizeof(double));
    for(int i = 0; i < N; i++){
        *(F_arr + i) = 0.5*(x_max - x_min)*Func(S_trans(-T + i*h,x_max,x_min),E)*Sigma(E*sqrt(1.0-S_trans(-T + i*h,x_max,x_min)))/(cosh(-T + i*h)*cosh(-T + i*h));
		//cout << "F = "<<0.5*(x_max - x_min)*Func(S_trans(-T + i*h,x_max,x_min),E)*Sigma(E*sqrt(1.0-S_trans(-T + i*h,x_max,x_min)))/(cosh(-T + i*h)*cosh(-T + i*h))<<"\n";
	}
	
    for (int i = 0; i < N; i++) { 
        if (i == 0 || i == N) 
            R += *(F_arr + i); 
        else if (i % 2 != 0) 
            R += 4 * (*(F_arr + i)); 
        else
            R += 2 * (*(F_arr + i)); 
		//cout << "R = " << R <<'\n';
    } 
    R = R * (h / 3); 
    free(F_arr);
    return R;
}
void Rad_corr()
{
	double N = 1000;
	std::ifstream Sigma_file;
	Sigma_file.open("Sigma_Obs.dat", ios::in);
	
    std::ofstream Sigma_rad_file;
	Sigma_rad_file.open("Sigma.dat", ios::out);
	double epsilon = 1e-4;
	double delta_azero = 1e-4;
	double R0[] = {0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9};
	//double R0[] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
	double R1[] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
	double R_Sub = 0.0;
	double R_Sub2 = 0.0;
    double Energy[] = {987.5,978.,975.,962.5,955.,951.1,950.,945.,936.,925.,912.5,900.
					,887.5,875,862.5,850.
					,837.5,825.,812.5,800.};
    double Sigma_50[20];
	double Sigma_100[20];
	double Sigma_200[20];
	double Sigma_temp[20];
	double Sigma_50_error[20];
	double Sigma_100_error[20];
	double Sigma_200_error[20];
    double delta = 1e-7;
	double Cor_zero = 0;
	double Sigma0[20];
	double Sigma1[20];
    double B = 0.0;
	double TT = 0.0;
	int i = 0;
	//TGraphErrors* gr1;
	//TFitResultPtr r1;
	while (Sigma_file >> Energy[i] >> Sigma_200[i] >> Sigma_200_error[i] >> Sigma_100[i] >> Sigma_100_error[i])
		i++;
	Sigma_file.close();
	//double h = -(1e-15- 1e-12)/N;
	double T = 19.0;
	double h = 2*T/N;
	double* F_arr = (double*)calloc(N, sizeof(double));
	double* x_arr = (double*)calloc(N, sizeof(double));
	double temp = 0.0;
	B = 2 * alpha * (log(4*900*900/(m_e*m_e)) - 1)/M_PI;
	a0 = 0.49;
	a1 = 951.0;
	a2 = 60.0;
    for(int i = 0; i < N; i++){
		*(x_arr + i) = -T + i*h;
        //*(F_arr + i) = F(1e-15 + i*h,900)/pow(1e-15 + i*h,B-1);
		temp = 950*sqrt(1.0 - S_trans(-T + i*h,delta_azero,0));
		TT = tanh(-T);//S_trans(-T + i*h,delta_azero,0);
		cout << "F(delta) = " << TT << '\n';
		if(temp < Mass_cut) 
			*(F_arr + i) = 0;
		else
			*(F_arr + i) = 0.5*(delta_azero)*F(S_trans(-T + i*h,delta_azero,0),950)*pol_2(temp)/(cosh(-T + i*h)*cosh(-T + i*h));
		cout << *(F_arr + i) << "\n";
	}
	//TGraph* gr;
	TGraph* gr1 = new TGraph(N,x_arr,F_arr);// Chi2 < 100
	//gr1->Draw("apl");
	free(x_arr);
	free(F_arr);
	for(int j = 0; j< 20; j++)
		Sigma0[j] = Sigma_100[j]/R0[j];
	//TF1 *f1 = new TF1("f1","[0]*(x-637.184)+[1]*(x-637.184)*(x-637.184) + [2]*pow(x-637.184,3)",Mass_cut,1000);
	/*TF1 *f2 = new TF1("f2","exp(-[0]*pow(x-637.184,-1))*([1] + x*[2])",Mass_cut,1000);
	TF1 *f3 = new TF1("f3","[1]/(1+exp(-[0]*(x-637.184)))*((x - 637.184)+[2]*pow(x-637.184,2) + [3]*pow(x-637.184,3))",Mass_cut,1000);
	TF1 *f4 = new TF1("f4","[0]*(1-exp(-pow(x-637.184,2)/pow([1],2)))*(1 + [2]*pow(x-637.184,4) + [3]*pow(x-637.184,2)))",Mass_cut,1000);/*/
	TF1 *f5 = new TF1("f5","[0]*(1-exp(-pow(x-[2],2)/pow([1],2)))*([3]+[4]*x)",Mass_cut,1000);
    //double N = 0.0, x_av = 0.0, sigma_x = 0.0, n = 0.0, bound = 0.0;
    //Find quadric parametrs fitting
	do
	{
		cout << "Sigma0[2] = " << Sigma0[2] << '\n';
		
		TGraphErrors* gr = new TGraphErrors(20,Energy,Sigma0,0,Sigma_100_error);// Chi2 < 100
		//TFitResultPtr r1 = gr->Fit([1]*(x-600)+[2]*pow(x-600,2) + [3]*pow(x-600,3),"S");	
		  //TFitResultPtr r1 = gr->Fit("pol2","S");
		
		f5->SetParameters(0.6,100,700,0.1,0.1);
		f5->SetParLimits(0,0.05,1.0);
		f5->SetParLimits(1,50.0,500.0);
		TFitResultPtr r1 = gr->Fit("f5","S");
		//	gr->Draw("apl");
	        
	        
		a0 = r1->Value(0);
		a1 = r1->Value(1);
		a2 = r1->Value(2);
	    a3 = r1->Value(3);
		a4 = r1->Value(4);
		//a5 = r1->Value(5);
		gr->Draw();
	        f5->Draw("same");
	       
		cout << "a2 = " << a2 <<'\t';
		R_Sub = 0.0;
		for(int j = 0; j < 20; j++){
		  
		  Cor_zero = Simpson_methood_special(F,pol_2,0,delta_azero,10000,Energy[j]);
		  //B = 2 * alpha * (log(4*Energy[j]*Energy[j]/(m_e*m_e)) - 1)/M_PI;
		  //cout << "betta = " << B << '\n';
		  //Cor_zero = F(delta_azero*0.5,Energy[j])*pol_2(Energy[j]*sqrt(1-delta_azero*0.5))*pow(0.5*delta_azero,-B + 1)*pow(delta_azero,B)/B * 0.5;
		  cout << "Cor = " << Cor_zero << '\n';
		  Sigma_temp[j] = Simpson_methood(F,pol_2,delta_azero,1,10000,Energy[j]);
		  Sigma_temp[j] = Sigma_temp[j] + Cor_zero;
		  cout << "Sigma temp = "<< Sigma_temp[j] <<"Sigma0 = "<<Sigma_100[j] << '\n';

		  //R1[j] = Sigma_temp[j]/Sigma0[j];
		  R1[j] = Sigma_temp[j]/pol_2(Energy[j]);
		  //			Sigma1[j] = Sigma_temp[j]/R0[j];
		  
		  //R1[j] = Sigma_100[j]/Sigma1[j];
		  cout <<"R0 = "<< R0[j]<<'\t' <<"R1 = " <<R1[j] << '\t' << "j = " << j << '\t' <<'\n';
		  
		  //R0[j] = Sigma_temp[j]/Sigma_100[j];
		  //cout << "R_sub =" << R_Sub <<"\n";
		  R_Sub = R_Sub + (R1[j] - R0[j])*(R1[j] - R0[j]);
		  
		  
		  Sigma0[j] = Sigma_100[j]/R0[j];
		  R0[j] = R1[j];
		}
	}while(R_Sub > epsilon);
	
	for(int y = 0; y < 20; y++){
	  Sigma_rad_file << 2*Energy[y] << '\t' << Sigma_100[y]/R1[y] <<'\t'<< Sigma_100_error[y]/R1[y] <<'\t' <<R1[y] <<'\n';
	  cout << "E = " << 2*Energy[y] << "'\t'(1 + R) = " << R1[y] << "'\t'Sigma = " << Sigma_100[y]/R1[y] << "'\t' Sigma0 = "<<Sigma_100[y] <<'\n';
	}
	Sigma_rad_file.close();
}
