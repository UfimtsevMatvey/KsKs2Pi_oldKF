#include <iostream>
#include <fstream>
#include <string>
using namespace std;
void Create_good_file()
{
	double A,B,C,D;
	int i = 0;
	double E[] = { 2000., 1900., 1800.,1700.,1600. };
	std::ofstream proces;
	proces.open("Good_MC_efficiency.dat", ios::out);
	std::ifstream sourse;
	sourse.open("efficiency_mc_only_central_sq.dat", ios::in);
	
	while(sourse >> A >> B >> C >> D){
		proces << E[i] << '\t' << A << '\t' << B << '\t' << C << '\t' << D << '\t' << '\n';
		i++;
	}
	sourse.close();
	proces.close();
}					