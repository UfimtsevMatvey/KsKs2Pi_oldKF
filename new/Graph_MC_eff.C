void Graph_MC_eff(){
   gROOT->SetStyle("ATLAS");
   TCanvas *c1 = new TCanvas("c1","Rad Cor");
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0);
  TGraphErrors *gr5 = new TGraphErrors("Sigma.dat","%lg %*lg %*lg %lg");
  gr5->SetMarkerStyle(24);
  // gr5->Draw("ap");
  gr5->GetYaxis()->SetTitle("Rad corection #sigma(e^{+}e^{-}#rightarrow K_{s} K_{s} #pi^{+} #pi^{-})");
  gr5->GetYaxis()->SetTitleOffset(1);
  gr5->GetXaxis()->SetTitle("Enegry C.M MeV");
  gr5->Draw("ap");
  /* TGraphErrors* gr1 = new TGraphErrors("efficiency_MC_new_4", "%lg %lg %*lg %*lg %lg %*lg %*lg");
  TGraphErrors* gr2 = new TGraphErrors("efficiency_MC_new_4", "%lg %*lg %lg %*lg %*lg %lg %*lg");
  TGraphErrors* gr3 = new TGraphErrors("efficiency_MC_new_4", "%lg %*lg %*lg %lg %*lg %*lg %lg");
  gr2->GetYaxis()->SetRangeUser(0.055,0.085);
  gr1->SetLineColor(4);
  gr2->SetLineColor(5);
  gr3->SetLineColor(6);
  gr2->Draw("apl");
  gr3->Draw("pl same");
  gr1->Draw("pl same");*/
  TCanvas *c2 = new TCanvas("c2","Crossection");
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0);
  // TGraph* gr4 = new TGraph("efficiency_real_new_2", "%lg %lg %*lg %*lg");
  //TGraph* gr5 = new TGraph("efficiency_real_new_2", "%lg %*lg %lg %*lg");
  //TGraph* gr6 = new TGraph("efficiency_real_new_2", "%lg %*lg %*lg %lg");
  //TGraphErrors *gr4 = new TGraphErrors("Sigma_2019","%lg %lg %*lg %*lg %lg %*lg %*lg");
  //TGraphErrors *gr5 = new TGraphErrors("Sigma_2019","%lg %*lg %lg %*lg %*lg %lg %*lg");
  //TGraphErrors *gr6 = new TGraphErrors("Sigma_2019","%lg %*lg %*lg %lg %*lg %*lg %lg");
  TGraphErrors *gr6 = new TGraphErrors("Sigma.dat","%lg %lg %lg %*lg");
  
  // gr5->GetYaxis()->SetRangeUser(-1.0,110.0);
 // gr4->SetMarkerStyle(16);
  //gr4->SetMarkerColor(6);
 // gr5->SetMarkerStyle(20);
  
 //gr6->SetMarkerColor(4);
  gr6->SetMarkerStyle(24);
 // gr5->Draw("ap");
 gr6->GetYaxis()->SetTitle("#sigma(e^{+}e^{-}#rightarrow K_{s} K_{s} #pi^{+} #pi^{-})");
 gr6->GetYaxis()->SetTitleOffset(1);
 gr6->GetXaxis()->SetTitle("Enegry C.M MeV");
 gr6->Draw("ap");
 // gr4->Draw("p same");
  //gr4->Draw("apl");
}
