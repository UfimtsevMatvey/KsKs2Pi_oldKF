#define KSKS2pi_cxx
#include "KSKS2pi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile.h"
#include "TFile.h"
#include "TMatrixTSym.h"
#include "TMath.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TBenchmark.h"
#include <iostream>
#include <string>
#include "Cmd3KF.C"
#include "CovPPhiTheta2PxPyPz.C"
#include "CovPtPhiTheta2PxPyPz.C"
#include <iostream>
#include <fstream>
#include "badrunfunction.C"
bool Cond_2KS(double, double, double,double);
bool Cond_2KS_out(double,double,double,double);
double error_a(int, int);
double error_b(int, int);
void Auto_get(double, double, double*, double*);
void NZero(int N1, int N2, double a, double a_e, double b, double b_e,double* N0, double* Error);
void KSKS2pi::Loop(std::string histFileName)
{
//   In a ROOT session, you can do:
//      Root > .L KSKS2pi.C
//      Root > KSKS2pi t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TFile *top = new TFile(histFileName.c_str(), "recreate");

    // minv, mism, ptot for 5-6 good tracks

   KFParticle CurParticle;
   vector<KFParticle> InParticles;
   vector<KFParticle> OutParticles;
   double KsMass = 497.614;
   double KMass = 493.677;
   double PiMass = 139.57;
   double NormChi2 = 0;
   int N_e = 0;
   int N_b = 0;
   int N_e2 = 0;
   int N_b2 = 0;
   int N_e3 = 0;
   int N_b3 = 0;
   int N_e20 = 0,N_b20 = 0, N_b30 = 0, N_e30 = 0,N_e40 = 0,N_b40 = 0, N_e60 = 0, N_b60 = 0, N_e50 = 0, N_b50 = 0, N_e100 = 0, N_b100 = 0, N_e200 = 0, N_b200 = 0; 
   int N_sig20 = 0, N_sig30 = 0, N_sig50 = 0,N_sig40 = 0, N_sig60 = 0, N_sig100 = 0, N_sig200 = 0;
   int N_sig2 = 0;
   int N_sig3 = 0;
   int N_sig = 0;
   int N_ev = 0;
   
   double Nb1_40 = 0, Nb2_40 = 0;
   double Nb1_30 = 0, Nb2_30 = 0;
   
   double Na01_40 = 0, Na0_40 = 0;
   double Na01_30 = 0, Na0_30 = 0;

   double Nb1_corr2 = 0, Nb2_corr2 = 0;
   double Nb1_corr1 = 0, Nb2_corr1 = 0;
   
   double Na01_corr2 = 0, Na0_corr2 = 0;
   double Na01_corr1 = 0, Na0_corr1 = 0;

   double N1_40 = 0, N2_40 = 0;
   double N1_30 = 0, N2_30 = 0;
   
   double error1 = 0.0;
   double error2 = 0.0;
   double error3 = 0.0;
   
   double error20 = 0.0;
   double error30 = 0.0;
   double error40 = 0.0;
   double error50 = 0.0;
   
   double error60 = 0.0;
   double error100 = 0.0;
   double error200 = 0.0;
    // KSKS
   /****Ufimtsev_2020****/
 
   

   //************************
   
   //*************************
   
   //*************************
   TH1D* Chi2_SQ_1  = new TH1D("Chi2_SQ_15","Chi2_SQ_15",100,0.0,100.0);
   TH1D* Chi2_8SQ_1  = new TH1D("Chi2_8SQ_15","Chi2_8SQ_15",100,0.0,100.0);
   
   TH1D* Chi2_SQ_2  = new TH1D("Chi2_SQ_30","Chi2_SQ_30",100,0.0,100.0);
   TH1D* Chi2_8SQ_2  = new TH1D("Chi2_8SQ_30","Chi2_8SQ_30",100,0.0,100.0);
   
   TH1D* Chi2_SQ_3  = new TH1D("Chi2_SQ_40","Chi2_SQ_40",100,0.0,100.0);
   TH1D* Chi2_8SQ_3  = new TH1D("Chi2_8SQ_40","Chi2_8SQ_40",100,0.0,100.0);
   //-----------------------------
   TH1D* LenKs = new TH1D("LenKs","",100,0,5);
   TH1D* Len2Pi = new TH1D("Len2Pi","",100,0,5);
   TH1D* Minv2Pi = new TH1D("Minv2Pi","",1024, 250, 1000);
   TH1D* MinvKs2Pi = new TH1D("MinvKs2Pi","",1024, 800,1800);
   TH1D* MinvKsPi = new TH1D("MinvKsPi","",1024,400,1200);
   
   
   //*******************************
   //*******************************
   TH2D* KsvsKS_without_cut = new TH2D("KsvsKS_without_cut","KsvsKS_without_cut",100,400.0,600.0,100,400.0,600.0);


   /****Ufimtsev_2020****/
   
	/*Ufimtsev add*/
   
   TH1D* KsMomIN = new TH1D("KsMomIN","KsMomIN",100.0,0.0,600);
   TH1D* KsMomOUT = new TH1D("KsMomOUT","KsMomOUT",100.0,0.0,600);
   TH1D* KsMomMOD= new TH1D("KsMomMOD","KsMomMOD",100.0,0.0,600);

   TH1D* PipMomIN = new TH1D("PipMomIN","PipMomIN",100.0,0.0,600);
   TH1D* PipMomOUT = new TH1D("PipMomOUT","PipMomOUT",100.0,0.0,600);
   TH1D* PipMomMOD= new TH1D("PipMomMOD","PipMomMOD",100.0,0.0,600);

   TH1D* KsMomIN_50 = new TH1D("KsMomIN","KsMomIN",100.0,0.0,600);
   TH1D* KsMomOUT_50 = new TH1D("KsMomOUT","KsMomOUT",100.0,0.0,600);

   TH1D* PipMomIN_50 = new TH1D("PipMomIN","PipMomIN",100.0,0.0,600);
   TH1D* PipMomOUT_50 = new TH1D("PipMomOUT","PipMomOUT",100.0,0.0,600);

   TH1D* Chiks1 = new TH1D("ChiSqCur_Ks1","ChiSqCur1",100,0.0,5000);
   TH1D* Chiks2 = new TH1D("ChiSqCur_Ks2","ChiSqCur2",250,0.0,5);
   TH1D* Chisq2 = new TH1D("Chisq2","Chisq2",100,0.0,100.0);

   TH1D* Chisq21 = new TH1D("Chisq21","Chisq21",200,0.0,200.0);

   TH1D* Inv_out_mass1 = new TH1D("Inv_out_mass1","Inv_out_mass1",100,600.0,1600.0);
   TH1D* Inv_in_mass1 = new TH1D("Inv_in_mass1","Inv_in_mass1",100,600.0,1600.0);
   
   TH1D* delta_mom_ks_in = new TH1D("delta_mom_ks_in","delta_mom_ks_in",100.0,-100.0,100.0);
   TH1D* delta_mom_ks_out = new TH1D("delta_mom_ks_out","delta_mom_ks_out",100.0,-100.0,100.0);

   TH1D* delta_in = new TH1D("delta_in","delta_in",200,-100.0,100.0);

    TH1D* delta_out = new TH1D("delta_out","delta_out",200,-100.0,100.0);

    TH1D* Inv_in_massmodel = new TH1D("Inv_in_massmodel","Inv_in_massmodel",100.0,600.0,1600.0);
	TH1D* Inv_in_massmodelall= new TH1D("Inv_in_massmodelall","Inv_in_massmodelall",100.0,600.0,1600.0);
    TH1D* Inv_massPipKs1 = new TH1D("Inv_massPipKs1","Inv_massPipKs1",100,600.0,1300.0);
	
	TH2D* KsPimVsKsPip = new TH2D("KsPimVsKsPip","KsPimVsKsPip",200,400.,1600.,200,400.,1600.);
	
    TH1D* Inv_NKs = new TH1D("InvNKs","InvNKs",100,300,700);
    
    TH1D* Inv_out_mass2 = new TH1D("Inv_out_mass2","Inv_out_mass2",100,600.0,1600.0);
    TH1D* Inv_in_mass2 = new TH1D("Inv_in_mass2","Inv_in_mass2",100,600.0,1600.0);
	
  	TH1D* Inv_outKsKs = new TH1D("Inv_outKsKs","Inv_outKsKs",300,500.0,1700.0);
	TH1D* Inv_inKsKs = new TH1D("Inv_inKsKs","Inv_inpKsKs",300,500.0,1700.0);
	
	/*-----------------------*/
    TH1D* hdPhiKS= new TH1D("hdPhiKS","dPhiKS",640,-3.2,3.2);  
    TH1D* hdTheKS= new TH1D("hdTheKS","dTheKS",640,-3.2,3.2);
    TH1D* hdPtKS= new TH1D("hdPtKS","dPtKS",400,-100,100);
    TH2D* hmks1vsmks2 = new TH2D("hmks1vsmks2","mKS1 vs mKS2",200,400.,600.,200,400.,600.);
    TH2D* hmks1vsmks2cut = new TH2D("hmks1vsmks2cut","mKS1 vs mKS2",200,400.,600.,200,400.,600.);
    TH2D* hmks1vsmks2cutBkg = new TH2D("hmks1vsmks2cutBkg","mKS1 vs mKS2 for background",200,400.,600.,200,400.,600.);
    TH2D* hmks1vsmks2all = new TH2D("hmks1vsmks2all","mKS1 vs mKS2 for all candidates",200,400.,600.,200,400.,600.);
    TH1D* hmkscut = new TH1D("hmkscut","mKS1 and mKS2",200,400.,600.);
    TH1D* hmkscutBkg = new TH1D("hmkscutBkg","mKS1 and mKS2 for bkg",200,400.,600.);

    TH1D* hRhonotks = new TH1D("hRhonotks","Distance to beam tracks not from KS",200,0,2.);
    TH1D* hRhopiks = new TH1D("hRhopiks","Distance to beam tracks from KS",300,0,6.);
    TH2D* hRmaxnotks = new TH2D("hRmaxnotks","Distance to beam from not KS tracks",200,0,2.0,200,0.,2.);
    TH2D* hRvtxks = new TH2D("hRvtxks","Distance to beam from KS vertex vs KS mss",300,0,6.0,200,400.,600.);
    TH1D* hRvtxks1D = new TH1D("hRvtxks1D","Distance to beam from KS verteses",120,0,6.0);
    TH1D* hRvtxks1Dbkg = new TH1D("hRvtxks1Dbkg","Distance to beam from KS verteses for Chi2 bkg",120,0,6.0);

    TH1D* hChi2min= new TH1D("hChi2min","Chi2min",100,0,10000);

    TH2D* hEtotv6Ptot= new TH2D("hEtotv6Ptot","Total energy - 2ebeam for 6 good tracks vs. Ptot",300,0.,600.,500,-700.,300.);
    TH2D* hEtotv5Ptot= new TH2D("hEtotv5Ptot","Total energy - 2*ebeamfor 5 good tracks vs. Ptot",300,0.,600.,500,-700.,300.);

    TH2D* hPtot2pi= new TH2D("hPtot2pi","Ptot1 vs Ptot2 for pions not from KS",300,0.,600.,300,0.,600.);

    TH1D* hThMis5Ngood= new TH1D("hThMis5Ngood","Missing polar angle  for 5 good tracks",62,0.,3.2);
    TH1D* hThMis5NgoodBkg= new TH1D("hThMis5NgoodBkg","Missing polar angle  for 5 good tracks - bkg",62,0.,3.2);
    TH1D* hThMis6Ngood= new TH1D("hThMis6Ngood","Missing polar angle  for 6 good tracks",62,0.,3.2);

    TH1D* hTheta5Ngood= new TH1D("hTheta5Ngood","Polar angle  for 5 good tracks",62,0.,3.2);
    TH1D* hTheta6Ngood= new TH1D("hTheta6Ngood","Polar angle  for 6 good tracks",62,0.,3.2);

    TH1D* hTheta5NgoodBkg= new TH1D("hTheta5NgoodBkg","Polar angle  for 5 good tracks bkg",62,0.,3.2);
    TH1D* hTheta6NgoodBkg= new TH1D("hTheta6NgoodBkg","Polar angle  for 6 good tracks bkg",62,0.,3.2);

    TH1D* hEtotv6Ngood= new TH1D("hEtotv6Ngood","Total energy - 2ebeam for 6 good tracks",200,-700.,300.);
    TH1D* hEtotv6NgoodBkg= new TH1D("hEtotv6NgoodBkg","Total energy - 2ebeam for 6 good tracks Chi2 Bkg",200,-700.,300.);
    TH1D* hEtotv5Ngood= new TH1D("hEtotv5Ngood","Total energy - 2*ebeamfor 5 good tracks",100,-700.,300.);
    TH1D* hEtotv5NgoodBkg= new TH1D("hEtotv5NgoodBkg","Total energy - 2*ebeamfor 5 good tracks for Chi2 bkg",100,-700.,300.);


    TH1D* hMinvKSpip= new TH1D("hMinvKSpip","Minv for KS + pip",80,600.,1400.);
    TH1D* hMinvKSpim= new TH1D("hMinvKSpim","Minv for KS + pim",80,600.,1400.);
    TH1D* hMinvKS2pi= new TH1D("hMinvKS2pi","Minv for KS + pim +pip",100,800.,1800.);
    TH1D* hMinvKS2piBkg= new TH1D("hMinvKS2piBkg","Minv for KS + pim +pip - background",100,800.,1800.);

    TH2D* hKSpipvsKSpim= new TH2D("hKSpipvsKSpim","Minv for KSpip vs KSpim",80,600.,1400.,80,600.,1400.);


    TH1D* hMinv2pinotks= new TH1D("hMinv2pinotks","Minv for pi+pi- not from KS",100,200.,1200.);
    TH1D* hMinv2pinotksBkg= new TH1D("hMinv2pinotksBkg","Minv for pi+pi- not from KS - background",100,200.,1200.);

    TH2D* hGTPvdEdX= new TH2D("hGTPvdEdX","Good Track Momentum vs dEdX",1500,0,1500., 1000, 0, 20000);
    TH2D* hGTPivdEdX= new TH2D("hGTPivdEdX","Good Pi Momentum vs dEdX",1500,0,1500., 1000, 0, 20000);

   Long64_t nbytes = 0, nb = 0;

   TLorentzVector PGood[10];
   double Norma_10 = 1.0;
   //M cout<<"Number of events to process = "<<nentries<<endl;
   int N_all_all = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     N_ev++;
     //cout << "ksErr \t " << kserr[0][0][0] << "\t"<<kserr[0][1][1]<<"\t"<<kserr[0][2][2]<<endl;
	N_all_all = N_all_all + 1;
     //if((jentry%100000)==0) cout<<jentry<<endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


      double pi = 3.1415927;

      double cutChi2r = 20.;
      double cutChi2z = 20.;
      double cutChi2ndf = 15.;
      int cutNhit = 5; //for 7 pi use 6
      //int cutNhit = 10;
      double cutRmin = 0.35;
      double cutKsl = 0.0;
      double cutKsalgn = 0.85;
      double cutZtrack = 15.; 
      double cutPtot = ebeam;
      double cutPmin = 40;
      //double cutPmin = 70;
      double cutEtot = 180.;
      double cutMindPhi = 0.0;
      double cutMaxdPhi = 3.15;
      double cutPmis = 180.;
      double cutEph = 20.;
      double cutAsy = 0.95;
      //double cutRad = 15;
      double cutRad = 13; //used for 7pi
      double cutRmax = 0.1;
      double cutEneu = 300;
      double cutChiKS = 25.;
      //double cutChiKS = 30.;// see modification of bkg subtraction

      double cutEn = 70.;

      double cutdEdXmax = 1500;
      double cutdEdXmin = -2000;

      int Ngood = 0;
      int charge = 0;
      int NgoodK = 0;
      int NgoodPi = 0;
      double PavGood = 0;
      double Pav0Good = 0;
      double Z0Good = 0;
      double rmax = 0.;


      double massPi0 = 134.98;
      double massEta = 547.862;
      double massKS = 497.67;
      double massPiC = 139.67;

      const int ntmax=10;
      int IndGood[ntmax]{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
      int IndRaw[ntmax]{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
      TLorentzVector PGood[10];
      TLorentzVector Pbeam(0.,0.,0.,2*ebeam); 
      TLorentzVector PFitPi[ntmax];
      TLorentzVector PFitK[ntmax];
      TLorentzVector PFitPiK[ntmax];
      double GTdedx[ntmax];
      double GTclus[ntmax];
      double Pidedx[ntmax];
      double Kdedx[ntmax];
      int flagK[ntmax];
      int chrg[ntmax];
      double Rad[ntmax];

      TLorentzVector PsumPi(0.,0.,0.,0.);
      TLorentzVector PsumnotksPi(0.,0.,0.,0.);
      TLorentzVector PsumPiK(0.,0.,0.,0.);
      TLorentzVector PPi0(0.,0.,0.,0.);
      TLorentzVector PEta(0.,0.,0.,0.);

      TLorentzVector Pks1(0.,0.,0.,0.);
      TLorentzVector Pks2(0.,0.,0.,0.);

      TLorentzVector Pkspip1(0.,0.,0.,0.);
      TLorentzVector Pkspim1(0.,0.,0.,0.);
      TLorentzVector Pkspip2(0.,0.,0.,0.);
      TLorentzVector Pkspim2(0.,0.,0.,0.);

      TLorentzVector Pks2pi1(0.,0.,0.,0.);
      TLorentzVector Pks2pi2(0.,0.,0.,0.);

      if(nt>10) cout<<"error in size of nt "<<nt<<endl;
      // cycle over all tracks
      for(int i = 0; i<nt && i<ntmax; ++i){

	//skip runs with new trigger config
	
	/****Ufimtsev_2020****/
	//if(badrun(runnum)) continue;
	/****Ufimtsev_2020****/
	
	if((runnum>16566 && runnum <16581))continue;

//skip runs with bad trigger and DC eficiency

	if(runnum==42476 || (runnum>=42480 && runnum <=42499) || runnum==42525)continue;// E=960 MeV

	if(runnum==43080 || runnum==43081)continue;// E=990 MeV

	if(runnum==40800 || runnum==40932 ||  runnum==40933 || runnum==40934 || runnum==40980 ||
	    runnum==41378 || runnum==41379 ||
	    (runnum>41406 && runnum <41733)||  
	    runnum==41737 || runnum==41738 ||
	    runnum==41747 || runnum==41748 ||
	    runnum==41752 || runnum==41753 ||  runnum==41754 || runnum==41755 || runnum==41766
	   )continue;//E=936

	if(runnum==41812 || runnum==41813 ||  runnum==41814 || 
	    (runnum>41822 && runnum <41842)||  
	    runnum==41865 || runnum==41869 || runnum==41870 || runnum==41871 ||
	    runnum==41875 || runnum==41887 ||
	    runnum==41937 || runnum==41942 ||  runnum==41948 || runnum==41972 || runnum==41973 || runnum==41975
	   )continue;//E=942

	if(runnum==43753 || runnum==43754 ||  runnum==43755 ||  runnum==43756 || runnum==43818 ||  runnum==43841 ||  runnum==43876 || 
	    runnum==43891 ||  runnum==43892 ||  runnum==43913 || runnum==43966 ||  runnum==43970  
	   )continue;//E=938.9

	if(runnum==44052 || runnum==44053 ||  runnum==44209 ||  runnum==44211 || runnum==44244 ||  runnum==44246 ||  runnum==44247 || 
	    runnum==44248   
	   )continue;//E=939.6

	if(runnum==44268 || runnum==44269 ||  runnum==44332 ||  runnum==44387 
	   )continue;//E=940.2

	if(runnum==44476 || runnum==44554 ||  runnum==44555 ||  runnum==44623 || runnum==44624 
	   )continue;//E=938.3

	if(runnum==44702 || runnum==44703 ||  runnum==44747 ||  
	    (runnum>44809 && runnum <44828)  
	   )continue;//E=937.5_1

	if(runnum==45002 || runnum==45003 || runnum==45008  || (runnum>=45011 && runnum<=45018) || 
	    runnum==45021 || runnum==45027 || runnum==45047 ||  runnum==45048 || runnum==45064 ||  
	    (runnum>45136 && runnum <45153) 
	   )continue;//E=940.8

	if(runnum==45212 
	   )continue;//E=937.5_2

	// select good tracks
	//	cout<<"Ptot "<<tptotv[i]<<" Z "<<tz[i]<<" Chi2 "<<tchi2r[i]<<" Nhit "<<tnhit[i]<<" Rmin "<<trho[i]<<endl;

		if(tptot[i]< cutPtot 
		   //&& fabs(trho[i])<cutRmin 
		   && fabs(tz[i])<cutZtrack
 		   && tnhit[i]>cutNhit
		   && tchi2ndf[i]<cutChi2ndf
		   //&& tchi2z[i]<cutChi2z 
		   && tptot[i]>cutPmin
		   ){
	  
	  
	  double Dz = 20.-fabs(tz[i]);
	  double The = tan(tthv[i]);
	  if(tthv[i]>3.1415927/2.)The = tan(3.1415927-tthv[i]);
	  Rad[i] = Dz*The;
	  if(Rad[i]>35.)Rad[i]=35;
	  //if(Rad[i]<cutRad)continue;// skip "short" tracks
	  

	  //Pidedx[Ngood] =  300000./(tptotv[i]-20.)+700.+0.9*tptotv[i]-tdedx[i]; //for 2011,2012 data
	  Pidedx[Ngood] =  260000./(tptotv[i]-20.)+700.+0.9*tptotv[i]-tdedx[i];//for 2017 data
	  Kdedx[Ngood] =  800000./(tptotv[i]-10.)+300.+0.9*tptotv[i]-tdedx[i];

	  
	  PFitPi[Ngood] = TLorentzVector(tptot[i]*TMath::Cos(tphi[i])*TMath::Sin(tth[i]),
			   tptot[i]*TMath::Sin(tphi[i])*TMath::Sin(tth[i]),
			   tptot[i]*TMath::Cos(tth[i]),
			   TMath::Hypot(tptot[i],139.67)); 
	  
	  PFitK[Ngood] = TLorentzVector(tptotv[i]*TMath::Cos(tphiv[i])*TMath::Sin(tthv[i]),
			  tptotv[i]*TMath::Sin(tphiv[i])*TMath::Sin(tthv[i]),
			  tptotv[i]*TMath::Cos(tthv[i]),
			  TMath::Hypot(tptotv[i],497.6)); 

	  flagK[Ngood] = 0;
	  double mass = 139.67;
	  if((Pidedx[Ngood]<-2000 || Kdedx[Ngood]<600) && Kdedx[Ngood]>-2500){
	    //	  if(tdedx[i]>200000/(tptot[i]-20)+1450+0.5*tptot[i]){
	    mass = 497.6;
	    flagK[Ngood] = 1;
	    NgoodK ++;
	  }
	  PFitPiK[Ngood] = TLorentzVector(tptotv[i]*TMath::Cos(tphiv[i])*TMath::Sin(tthv[i]),
			   tptotv[i]*TMath::Sin(tphiv[i])*TMath::Sin(tthv[i]),
			   tptotv[i]*TMath::Cos(tthv[i]),
			   TMath::Hypot(tptotv[i],mass));


	  Pav0Good += tptot[i];
	  PavGood += tptotv[i];
	  charge += tcharge[i];
	  chrg[Ngood] = tcharge[i];
	  Z0Good += tz[i];
	  GTdedx[Ngood] = tdedx[i];
	  GTclus[Ngood] = ten[i];
	  //hRmin->Fill(trho[i]);
	  //if(fabs(trho[i])>rmax)rmax = fabs(trho[i]);
	  
	  if(Pidedx[Ngood]>cutdEdXmin && Pidedx[Ngood]<cutdEdXmax)NgoodPi++;
	  
	  IndGood[Ngood] = i;
	  IndRaw[i] = Ngood;
	  
	  //PsumPi += PFitPi[Ngood];
	  PsumPiK += PFitPiK[Ngood];

	  Ngood++;
	}
      }

      if(nks==2){}// cout<<"Indexes  "<<ksvind[0][0]<<" "<<ksvind[0][1]<<"  "<<ksvind[1][0]<<" "<<ksvind[1][1]<<endl;

      double ecurr = emeas;
      if(emeas == 0) ecurr = ebeam;


      if(nks <2 || Ngood<5 || Ngood>6 )continue;


      // check de/dx

      for(int i=0;i<Ngood;++i){

	hGTPvdEdX->Fill(PFitPi[i].P(),tdedx[IndGood[i]]);
	if(NgoodPi==Ngood)hGTPivdEdX->Fill(PFitPi[i].P(),tdedx[IndGood[i]]);
      }

      //      if(Ngood != NgoodPi)continue;
      
      //look fot two best KS's

      double Chi2min = 10000000;
      double InvKsPi[4];

      double Chi1 = 0.0, Chi2 = 0.0;
      double Invin1=0.0, Invin2 =0.0;
      double Invout1=0.0, Invout2 =0.0;
      int Pi_indep[10];
      int Indep_Pi = 0;
      double InvKsPiout[3];
      double Count_KChi[5];
      double Count_Kin[5];
      double Count_Kout[5];
      double Inv_mass_Kin[20][5];
      double Inv_mass_Kout[20][5];
      double ChiSqks1Cur[20][5];//-
      double ChiSqks2Cur[5 + 4 + 3 +2 + 1];
      int nn =0;
      int indks1=-1;
      int indks2=-1;
      double InvMassin = 0.0, InvMassout = 0.0;
      int Npi10 = 0;
      int Npi20 = 0;

      int indpiks1[2]{-1,-1};
      int indpiks2[2]{-1,-1};
      int indpinotks[2]{-1,-1};
      int NU = 0, NU1 = 0;//-
      // cout<<" New cycle "<<endl;
      
      /*Ufimtev code*/
      


      int NI = 0;
      double Kmass = 892.0;
      int KsMin = 0, Ks = 0;
      int ChargI = 0;
      int TT = 0;
      double min = 10000000.0;
      int ksBest = -1;
      int ksBest2 = -1;
      int Npi1u = 0;
      int Npi2u = 0;

      int indpi1u[2]{-1,-1};
      int indpi2u[2]{-1,-1};
      TLorentzVector modelKs[nsim];
      TLorentzVector modelKs2(0.0,0.0,0.0,0.0);
      TLorentzVector modelPi[nsim];
      TLorentzVector modelindPi2(0.0,0.0,0.0,0.0);
	
      TLorentzVector tVextor1(0.0,0.0,0.0,0.0);
      TLorentzVector tVector2(0.0,0.0,0.0,0.0);

      TLorentzVector Pipmodel;
      TLorentzVector Pimmodel;
      TLorentzVector Pipin;
      TLorentzVector Pimin;

      TLorentzVector Pipout;
      TLorentzVector Pimout;

      TLorentzVector Ksm;
      TLorentzVector Ksrin;
      TLorentzVector Ksrout;
	
      double Ks1Pi1m = 0;//(modelKs[1]+modelPi[1]).M();
      double Ks1Pi2m = 0;//(modelKs[1]+modelPi[2]).M();
      double Ks1Pi1r = 0;//(PFitPi[Pi_indep[0]]+PGood[1]).M();
      double Ks1Pi2r = 0;//(PFitPi[Pi_indep[1]]+PGood[1]).M();

      double Ks2Pi1m = 0;//(modelKs[2]+modelPi[1]).M();
      double Ks2Pi2m = 0;//(modelKs[2]+modelPi[2]).M();
      double Ks2Pi1r = 0;//(PFitPi[Pi_indep[0]]+PGood[3]).M();
      double Ks2Pi2r = 0;//(PFitPi[Pi_indep[1]]+PGood[3]).M();
      int Dec_event = 0;
      min = 10000000.0;
      for(int k = 0; k < nks; k++){
		for(int j = 0; j < nks;j++){	
		
		if(kstype[k]==0 && kstype[j]==0
			&& ksalign[k]>cutKsalgn && ksalign[j]>cutKsalgn
			){
			Npi1u = 0;
			Npi2u = 0;	
			for(int i=0;i<Ngood;++i){
				if((abs(tphi[IndGood[i]] - kspiphi[k][0])<0.15 &&
				abs(tth[IndGood[i]] - kspith[k][0])<0.15 &&
				abs(tptot[IndGood[i]] - kspipt[k][0])<15.) ||
				(abs(tphi[IndGood[i]] - kspiphi[k][1])<0.15 &&
				abs(tth[IndGood[i]] - kspith[k][1])<0.15 &&
				abs(tptot[IndGood[i]] - kspipt[k][1])<15.)
				){
	      if(Npi1u<2){
		indpi1u[Npi1u] = i;
		Npi1u++;
	      }
	    }
	  }

	   
       for(int i=0;i<Ngood;++i){
	    if((abs(tphi[IndGood[i]] - kspiphi[j][0])<0.15 &&
		abs(tth[IndGood[i]] - kspith[j][0])<0.15 &&
		abs(tptot[IndGood[i]] - kspipt[j][0])<15.) ||
	       (abs(tphi[IndGood[i]] - kspiphi[j][1])<0.15 &&
		abs(tth[IndGood[i]] - kspith[j][1])<0.15 &&
		abs(tptot[IndGood[i]] - kspipt[j][1])<15.)
	       ){
	      if(Npi2u<2 && i != indpi1u[0] && i != indpi1u[1]){
		indpi2u[Npi2u] = i;
		Npi2u++;
	      }
	    }
	  }
	   if(Npi1u == 2 && Npi2u ==2
	     ){
	  if(k != j){
	  if((ksminv[j]-KsMass)*(ksminv[j]-KsMass)+(ksminv[k]-KsMass)*(ksminv[k]-KsMass)<min)
	      {
		min =(ksminv[j]-KsMass)*(ksminv[j]-KsMass)+(ksminv[k]-KsMass)*(ksminv[k]-KsMass);
		ksBest2 = j;
		ksBest = k;
		indpiks1[0] = indpi1u[0];
		indpiks1[1] = indpi1u[1];
		indpiks2[0] = indpi2u[0];
		indpiks2[1] = indpi2u[1];
	      }
	  }
	   }
	  }
	}
	}
    if(ksBest<0 || ksBest2<0)continue;
       int counter = 0;
	
	for(int i=0;i<Ngood;++i){

	    
	    if(i != indpiks1[0] && i != indpiks1[1] 
	       && i != indpiks2[0] && i != indpiks2[1]
	       && counter == 0
	       ) {
	      counter++; 
	      Pi_indep[0] = i;
	    }

	    if(i != indpiks1[0] && i != indpiks1[1] 
	       && i != indpiks2[0] && i != indpiks2[1]
	       && i != Pi_indep[0] 
	       && counter == 1
	       ) {
	      counter++; 
	      Pi_indep[1] = i;
	    }
	}
	
	
	if(counter != 2) continue;
	if(Ngood != 6) continue;
	if((kstype[ksBest] == 1)||(kstype[ksBest2] == 1)) continue;
	NI = Pi_indep[0];	
	
	PGood[0] = TLorentzVector(tptot[NI]*sin(tth[NI])*cos(tphi[NI]),tptot[NI]*sin(tth[NI])*sin(tphi[NI]),
	tptot[NI]*cos(tth[NI]),TMath::Sqrt(tptot[NI]*tptot[NI]+139.57*139.57));
	nn = ksBest;
	PGood[1] = TLorentzVector(ksptot[nn]*sin(ksth[nn])*cos(ksphi[nn]),ksptot[nn]*sin(ksth[nn])*sin(ksphi[nn]),
	ksptot[nn]*cos(ksth[nn]),TMath::Sqrt(ksptot[nn]*ksptot[nn]+KsMass*KsMass));
					
	NI = Pi_indep[1];	
	PGood[2]  = TLorentzVector(tptot[NI]*sin(tth[NI])*cos(tphi[NI]),tptot[NI]*sin(tth[NI])*sin(tphi[NI]),
	tptot[NI]*cos(tth[NI]),TMath::Sqrt(tptot[NI]*tptot[NI]+139.57*139.57));
	nn = ksBest2;
	PGood[3] = TLorentzVector(ksptot[nn]*sin(ksth[nn])*cos(ksphi[nn]),ksptot[nn]*sin(ksth[nn])*sin(ksphi[nn]),
	ksptot[nn]*cos(ksth[nn]),TMath::Sqrt(ksptot[nn]*ksptot[nn]+KsMass*KsMass));
	
	/*------*/
	InParticles.clear();
	OutParticles.clear();

	
	CurParticle.P = PFitPi[Pi_indep[0]];
	CurParticle.Cov = GetTrErrorMatrix(PFitPi[Pi_indep[0]],terr[Pi_indep[0]]);
	InParticles.push_back(CurParticle);
	CurParticle.P = PFitPi[Pi_indep[1]];
	CurParticle.Cov = GetTrErrorMatrix(PFitPi[Pi_indep[1]],terr[Pi_indep[1]]);
	InParticles.push_back(CurParticle);
	CurParticle.P = PGood[1];
       	CurParticle.Cov = GetTrErrorMatrix(PGood[1],kserr[ksBest]);

	double len_1 = kslen[ksBest];
	double len_2 = kslen[ksBest2];
		double C_cut = 0;
	//M cout << "ksErr \t " << kserr[ksBest][0][0] << "\t"<<kserr[ksBest][1][1]<<"\t"<<kserr[ksBest][2][2]<<endl;
	InParticles.push_back(CurParticle);
       	CurParticle.P = PGood[3];
	CurParticle.Cov = GetTrErrorMatrix(PGood[3],kserr[ksBest2]);
	InParticles.push_back(CurParticle);
	
	Chi2 = Cmd3KF(ebeam,InParticles,OutParticles);
	double Chi2_cut = 26;
	double Delta_Chi2 = 1;
	Auto_get(Chi2, Chi2_cut,&N1_40,&N2_40);
	Auto_get(Chi2, Delta_Chi2 + Chi2_cut,&N1_30,&N2_30);
	

	Chisq2->Fill(Chi2);
	Chisq2->Fill(Chi2);
	N_e++;
	/********************/
	if(tcharge[Pi_indep[1]]>0)	KsPimVsKsPip->Fill((PFitPi[Pi_indep[1]]+PGood[1]).M(),(PFitPi[Pi_indep[0]]+PGood[3]).M());
	else	KsPimVsKsPip->Fill((PFitPi[Pi_indep[0]]+PGood[1]).M(),(PFitPi[Pi_indep[1]]+PGood[3]).M());
	/******************************************************/
	// 211: Pi+
	//-211: Pi-
	// 310: Ks
	int counterPi = -1, counterKs = -1;
	int modelPinotKs[2] = {-1,-1};
	int flag1 = 0;
	int flag2 = 0;
	int YY[2];
	int tN = 0;
	int ttN = 0;
	int tNN = 0;
	double dx = 0.0, dy = 0.0,dz = 0.0, dr = 0.0;
	double mindr = 100000000.0;

	Ks1Pi1m = (modelKs[1]+modelPi[1]).M();
	Ks1Pi2m = (modelKs[1]+modelPi[2]).M();
	Ks1Pi1r = (PFitPi[Pi_indep[0]]+PGood[1]).M();
	Ks1Pi2r = (PFitPi[Pi_indep[1]]+PGood[1]).M();

	Ks2Pi1m = (modelKs[2]+modelPi[1]).M();
	Ks2Pi2m = (modelKs[2]+modelPi[2]).M();
	Ks2Pi1r = (PFitPi[Pi_indep[0]]+PGood[3]).M();
	Ks2Pi2r = (PFitPi[Pi_indep[1]]+PGood[3]).M();
	
	double KsPipm = 0.0;
	double KsPIpr = 0.0;
	counterPi = 0;
	counterKs = 0;
	
	if(tcharge[Pi_indep[0]] > 0){
	  Pipin = PFitPi[Pi_indep[0]];
	  Pimin = PFitPi[Pi_indep[1]];
	  
	  Pipout = OutParticles[0].P;
	  Pimout = OutParticles[1].P;
	}
	else{
	  Pipin = PFitPi[Pi_indep[1]];
	  Pimin = PFitPi[Pi_indep[0]];
	  
	  Pipout = OutParticles[1].P;
	  Pimout = OutParticles[0].P;
	}
	bool flag = true;
	if(abs((OutParticles[2].P+Pipout).M()-892) >abs((OutParticles[3].P+Pipout).M()-892)){
	  Ksrin = InParticles[3].P;
	  Ksrout = OutParticles[3].P;
	  flag = true;
	}
	else{
	  Ksrin = InParticles[2].P;
	  Ksrout = OutParticles[2].P;
	  flag = false;
	}
	/************/
	double massKS_cut = 27.0;
	bool R = true;
	bool R2 = true;
	R = Cond_2KS(KsMass,massKS_cut,ksminv[ksBest],ksminv[ksBest2]);
	R2 = Cond_2KS_out(KsMass,massKS_cut,ksminv[ksBest],ksminv[ksBest2]);
	if(R){
	  Chi2_SQ_3->Fill(Chi2);
	  if(nsim != 0){
	    C_cut = Delta_Chi2 + Chi2_cut;
	    Auto_get(Chi2, C_cut,&Na0_corr1,&Na01_corr1);
	    C_cut = Chi2_cut;
	    Auto_get(Chi2, C_cut,&Na0_corr2,&Na01_corr2);
	  }
	}
	if(R2){
	  Chi2_8SQ_3->Fill(Chi2);
	  C_cut = Delta_Chi2 + Chi2_cut;
	  Auto_get(Chi2, C_cut,&Nb1_corr1,&Nb2_corr1);
	  C_cut = Chi2_cut;
	  Auto_get(Chi2, C_cut,&Nb1_corr2,&Nb2_corr2);
	}
	//********************************
	massKS_cut = 25.0;
	R = Cond_2KS(KsMass,massKS_cut,ksminv[ksBest],ksminv[ksBest2]);
	R2 = Cond_2KS_out(KsMass,massKS_cut,ksminv[ksBest],ksminv[ksBest2]);
	nn = ksBest;
	double ksP = (OutParticles[2].P).Rho(), ksTh = (OutParticles[2].P).Theta(), ksPhi = (OutParticles[2].P).Phi();
	TLorentzVector new_methood_4vKs1 = TLorentzVector(ksP*sin(ksTh)*cos(ksPhi),ksP*sin(ksTh)*sin(ksPhi),
	ksP*cos(ksTh),TMath::Sqrt(ksP*ksP+KsMass*KsMass));
        ksP = (OutParticles[3].P).Rho(); ksTh = (OutParticles[3].P).Theta(); ksPhi = (OutParticles[3].P).Phi();
	TLorentzVector new_methood_4vKs2 = TLorentzVector(ksP*sin(ksTh)*cos(ksPhi),ksP*sin(ksTh)*sin(ksPhi),
	ksP*cos(ksTh),TMath::Sqrt(ksP*ksP+KsMass*KsMass))
	if(R){
	  Chi2_SQ_2->Fill(Chi2);
	  if(nsim != 0){
	    C_cut = Delta_Chi2 + Chi2_cut;
	    Auto_get(Chi2, C_cut,&Na0_40,&Na01_40);
	    C_cut = Chi2_cut;
	    Auto_get(Chi2, C_cut,&Na0_30,&Na01_30);
	  }
	  if(Chi2 < Chi2_cut){
	    LenKs->Fill(len_1);
	    LenKs->Fill(len_2);
	    
	    Minv2Pi->Fill((Pipout + Pimout).M());
	    MinvKs2Pi->Fill((Pipout + Pimout + new_methood_4vKs1).M());
	    MinvKs2Pi->Fill((Pipout + Pimout + new_methood_4vKs2).M());
	    if(flag)
	      MinvKsPi->Fill((Pipout + new_methood_4vKs2).M());
	    else
	      MinvKsPi->Fill((Pipout + new_methood_4vKs1).M());
	  }
	  
	}
	if(R2){
	  Chi2_8SQ_2->Fill(Chi2);
	  C_cut = Delta_Chi2 + Chi2_cut;
	  Auto_get(Chi2, C_cut,&Nb1_40,&Nb2_40);
	  C_cut = Chi2_cut;
	  Auto_get(Chi2, C_cut,&Nb1_30,&Nb2_30);
	} 
	//*********************************
	KsvsKS_without_cut->Fill(ksminv[ksBest],ksminv[ksBest2]);
	//***********Get Histogram Len Ks, Len Pi, m(Pi-Pi+), m(Ks 2Pi)
	//****Ufimtsev_2020****//
	if(Chi2 > 100) continue;

	Inv_in_mass2->Fill((Ksrin + Pipin).M());
	KsMomIN->Fill(Ksrin.Rho());
	PipMomIN->Fill(Pipin.Rho());

	Inv_out_mass2->Fill((Ksrout + Pipout).M());
	KsMomOUT->Fill(Ksrout.Rho());
	PipMomOUT->Fill(Pipout.Rho());

	
	cout << "nsim = "<< nsim << '\n';
	if (nsim == 0) continue;

	for(int yy = 0; yy < nsim; yy++){
	  
	  
 	  if(((simtype[yy] == 211)||(simtype[yy] == -211))){
	    if(simorig[yy] != 310){ 
	      counterPi++; 
	      modelPi[counterPi] = TLorentzVector(1000*simmom[yy]*TMath::Cos(simphi[yy])*TMath::Sin(simtheta[yy]),
							1000*simmom[yy]*TMath::Sin(simphi[yy])*TMath::Sin(simtheta[yy]),
							1000*simmom[yy]*TMath::Cos(simtheta[yy]),
							TMath::Hypot(1000*simmom[yy],139.57));
			YY[counterPi] = yy;
	    }
	  }

	  if(simtype[yy] == 310){
	    counterKs++;
	    modelKs[counterKs] = TLorentzVector(1000*simmom[yy]*TMath::Cos(simphi[yy])*TMath::Sin(simtheta[yy]),
							1000*simmom[yy]*TMath::Sin(simphi[yy])*TMath::Sin(simtheta[yy]),
							1000*simmom[yy]*TMath::Cos(simtheta[yy]),
							TMath::Hypot(1000*simmom[yy],massKS));
	    cout << "type = " << yy << " "<<simtype[yy]<<"Mom = "<<simmom[yy] << endl;
							}
	
	}
	
	 Ks1Pi1m = (modelKs[1]+modelPi[1]).M();
	 Ks1Pi2m = (modelKs[1]+modelPi[2]).M();
	 Ks1Pi1r = (PFitPi[Pi_indep[0]]+PGood[1]).M();
	 Ks1Pi2r = (PFitPi[Pi_indep[1]]+PGood[1]).M();

	 Ks2Pi1m = (modelKs[2]+modelPi[1]).M();
	 Ks2Pi2m = (modelKs[2]+modelPi[2]).M();
	 Ks2Pi1r = (PFitPi[Pi_indep[0]]+PGood[3]).M();
	 Ks2Pi2r = (PFitPi[Pi_indep[1]]+PGood[3]).M();
	 
	if(simtype[YY[1]] > 0){
	  Pipmodel = modelPi[1];
	  Pimmodel = modelPi[2];
	}
	else{
	  Pipmodel = modelPi[2];
	  Pimmodel = modelPi[1];
	}

	if(abs((modelKs[1]+Pipmodel).M()-892) >abs((modelKs[2]+Pipmodel).M()-892))
	  Ksm = modelKs[2];
	else
	  Ksm = modelKs[1];
 	
	

	Inv_in_massmodelall->Fill((Ksm+Pipmodel).M());
	KsMomMOD->Fill(Ksm.Rho());
	PipMomMOD->Fill(Pipmodel.Rho());



	delta_mom_ks_in->Fill(Ksm.Rho()-Ksrin.Rho());
	delta_mom_ks_out->Fill(Ksm.Rho()-Ksrout.Rho());
	
	delta_in->Fill((Ksm+Pipmodel).M() - (Ksrin + Pipin).M());
	delta_out->Fill((Ksm+Pipmodel).M() - (Ksrout + Pipout).M());
	      // end of loop		       
   }
	double C_b_40 = Nb2_40/Nb1_40;
	double C_b_30 = Nb2_30/Nb1_30;

	double C_b_corr1 = Nb2_corr1/Nb1_corr1;
	double C_b_corr2 = Nb2_corr2/Nb1_corr2;

	double C_a_40 = 0;
	double C_a_30 = 0;
	double a_e = 0., b_e = 0.;
	double N0 = 0., N0_error = 0.;
	double MC_eficiency_30 = 0., MC_eficiency_40 = 0.;
	double MC_eficiency_30_error = 0., MC_eficiency_40_error = 0.;
	ofstream myfile;
	if(nsim != 0){
		C_a_40 = Na01_40/Na0_40;
		C_a_30 = Na01_30/Na0_30;
		ofstream file_sim;
		file_sim.open("Coefficient_A.dat", ios::app);
		file_sim << Na01_40 << '\t' << Na0_40 << '\t' << Na01_30 << '\t' << Na0_30 << '\t'
			 << C_a_30 << '\t' << C_a_40<<'\t'<< error_a(Na0_30,Na01_30)<<'\t'<<error_a(Na0_40,Na01_40) <<'\t' <<N_all_all<< '\n';
		file_sim.close();
		
		
		a_e = error_a(Na0_30,Na01_30);
		b_e = error_a(Nb1_30,Nb2_30);
		
		NZero(N1_30,N2_30,C_a_30,a_e,C_b_30,b_e,&N0,&N0_error);
		
		MC_eficiency_30 = N0/N_all_all;
		MC_eficiency_30_error = N0_error/N_all_all;
		a_e = error_a(Na0_40,Na01_40);
		b_e = error_a(Nb1_40,Nb2_40);
		NZero(N1_40,N2_40,C_a_40,a_e,C_b_40,b_e,&N0,&N0_error);
		MC_eficiency_40 = N0/N_all_all;
		MC_eficiency_40_error = N0_error/N_all_all;
		myfile.open ("efficiency_mc_only_central_sq.dat", ios::app);
		myfile << MC_eficiency_30 << '\t' << MC_eficiency_30_error <<'\t'<< MC_eficiency_40 << '\t' << MC_eficiency_40_error << '\n';
	}		
	if(nsim == 0){
	  ofstream file1;
	  file1.open("Number_events.dat",ios::app);
	  file1  << Nb1_40 << '\t'<< Nb2_40 <<'\t' << Nb1_30 << '\t'<< Nb2_30 << '\t'<< N1_30 << '\t'
		 << N2_30 << '\t' << N1_40 << '\t' << N2_40 << '\t' << C_b_30 << '\t' << error_b(Nb1_30,Nb2_30) 
		 << '\t' << C_b_40 << '\t' << error_b(Nb1_40,Nb2_40) << '\n';
	  file1.close();
	  ofstream file2;
	  file2.open("Number_event_corr.dat",ios::app);
	  file2 << Nb1_corr2 << '\t'<< Nb2_corr2 <<'\t' << Nb1_corr1 << '\t'<< Nb2_corr1 << '\t'<< N1_30 << '\t'
		 << N2_30 << '\t' << N1_40 << '\t' << N2_40 << '\t' << C_b_corr1 << '\t' << error_b(Nb1_corr1,Nb2_corr1) 
		 << '\t' << C_b_corr2 << '\t' << error_b(Nb1_corr2,Nb2_corr2) << '\n';
	  file2.close();
	}
	
   N_sig20 = N_e20 - N_b20;
   N_sig30 = N_e30 - N_b30;
   N_sig40 = N_e40 - N_b40;
   N_sig60 = N_e60 - N_b60;
   N_sig50 = N_e50 - N_b50;
   N_sig100 = N_e100 - N_b100;
   N_sig200 = N_e200 - N_b200;
   error20 = sqrt(N_e20 + N_b20);
   error30 = sqrt(N_e30 + N_b30);
   error40 = sqrt(N_e40 + N_b40);
   error60 = sqrt(N_e60 + N_b60);
   error50 = sqrt(N_e50 + N_b50);
   error100 = sqrt(N_e100 + N_b100);
   error200 = sqrt(N_e200 + N_b200);
   
   if(nsim == 0){
	myfile.open ("efficiency_real_only_central_sq.dat", ios::app);
	myfile << N_sig20 <<'\t'<<N_sig30<<'\t'<<N_sig50<<'\t'<<N_sig60 <<'\t'<<N_sig100<<'\t'<<N_sig200<< '\t'  << error20  <<'\t'<<error30 << '\t' << error50<<
		'\t' << error60<< '\t' << error100<< '\t' << error200<<'\t'<<N_sig40 <<error40<<"\n";
   }
   
   myfile.close();
   top->Write();
   top->Save();
	
}
bool Cond_2KS(double KsMass,double massKS_cut,double ksBest,double ksBest2)
{
  bool R = (ksBest > KsMass - massKS_cut) && (ksBest < KsMass + massKS_cut) 
    && (ksBest2 > KsMass - massKS_cut) && (ksBest2 < KsMass + massKS_cut);
  return R;
}
bool Cond_2KS_out(double KsMass, double massKS_cut, double ksBest, double ksBest2)
{
  double massKS_out = 3 * massKS_cut;
  //if(massKS_out > 100) massKS_out = 100;
  bool R = true;
  bool T = Cond_2KS(KsMass, massKS_cut,ksBest, ksBest2);
  if(not T)
    R = (ksBest > KsMass - massKS_out) && (ksBest < KsMass + massKS_out) 
      && (ksBest2 > KsMass - massKS_out) && (ksBest2 < KsMass + massKS_out);
  else 
    R = false;
  return R;
}
double error_a(int N1, int N2)
{
	double R = (double)N2/((double)N1*N1) + (double)N2*N2/((double)N1*N1*N1);
	return sqrt(R);
}

double error_b(int N1, int N2)
{
	double R = (double)N2/((double)N1*N1) + (double)N2*N2/((double)N1*N1*N1);
	return sqrt(R);
}

void Auto_get(double A, double B, double* C, double* D)
{
	if(A < B)
		(*C)++;
	if((A > B)&&(A < 2*B))
		(*D)++;
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
