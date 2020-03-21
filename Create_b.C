#include <iostream>
#include <fstream>
#include <string>
using namespace std;
double error_div(int N1, int N2);
void Create_b()
{
	double En[] = { 987.5,978.,975.,962.5,955.,951.1,950.,945.,936.,925.,912.5,900.
					,887.5,875,862.5,850.
					,837.5,825.,812.5,800. };
	double En_mc[] = {2000, 1900, 1800, 1700, 1600};
	double B_30;
	double B_40;
	double Be_30;
	double Be_40;
	double N1_30, N1_40, N2_30, N2_40;
	double _N1_30, _N1_40, _N2_30, _N2_40;
	double S_N1_30, S_N1_40, S_N2_30, S_N2_40;
	double SN1_30, SN1_40, SN2_30, SN2_40;
	SN1_30 = 0;
	SN1_40 = 0;
	SN2_30 = 0;
	SN2_40 = 0;
	S_N1_30 = 0;
	S_N1_40 = 0;
	S_N2_30 = 0;
	S_N2_40 = 0;
	double NULL_V;
	std::ifstream real;
	real.open("Number_events.dat", ios::in);
	std::ofstream B_file;
	B_file.open("Coefficient_B_and_Nevents.dat", ios::out);
	int j = 0;
	for(int i = 0; i < 21; i++){
		
		if(((i%4 == 0)&&(i != 0))||(i == 20)){
			B_30 = SN2_30/SN1_30;
			Be_30 = error_div(SN1_30,SN2_30);
			
			B_40 = SN2_40/SN1_40;
			Be_40 = error_div(SN1_40,SN2_40);
			B_file << En_mc[j] << '\t' << B_30 << '\t' << Be_30 << '\t' << B_40 << '\t' << Be_40 << '\t' << S_N1_30 
							<< '\t' << S_N1_40 <<'\t' << S_N2_30 << '\t' << S_N2_40 << '\n';
			SN1_30 = 0;
			SN1_40 = 0;
			SN2_30 = 0;
			SN2_40 = 0;
			
			S_N1_30 = 0;
			S_N1_40 = 0;
			S_N2_30 = 0;
			S_N2_40 = 0;
			j++;
		}
		if(i != 20){
			real >> N1_30 >> N1_40>> N2_30 >> N2_40 >> _N1_30 >> _N2_30 >> _N1_40 >> _N2_40>> NULL_V>> NULL_V>> NULL_V>> NULL_V;
			SN1_30 = SN1_30 + N1_30;
			SN1_40 = SN1_40 + N1_40;
			SN2_30 = SN2_30 + N2_30;
			SN2_40 = SN2_40 + N2_40;
			
			S_N1_30 = S_N1_30 + _N1_30;
			S_N1_40 = S_N1_40 + _N1_40;
			S_N2_30 = S_N2_30 + _N2_30;
			S_N2_40 = S_N2_40 + _N2_40;
		}
	}
	B_file.close();
	real.close();
}
double error_div(int N1, int N2)
{
	double R = (double)N2/((double)N1*N1) + (double)N2*N2/((double)N1*N1*N1);
	return sqrt(R);
}