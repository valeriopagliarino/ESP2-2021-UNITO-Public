/*******************************************
* Esperimentazioni II - 04/05/2016

  Esempio lettura da file

  per eseguire:
  root -l 'letturaFile.C++("lente_biconvessa.txt",40,200,300)'
  legge il file lente_biconvessa.txt e riempie istogrammi con 40 bin e intervallo asse x = 200-300

* R. Bellan - riccardo.bellan@unito.it
* S. Argiro'- stefano.argiro@unito.it updates
*******************************************/


#include <iostream>
#include <fstream>
#include <string>
#include <algorithm> 

#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TF1.h>
#include <TMath.h>

using namespace std;

void letturaFile (const std::string &nome_file, const double & nbins = 40, const double& minb = 200, const double& maxb = 300) {
  ifstream myfile;
  myfile.open(nome_file.c_str());


  if(!myfile.is_open()) {
    cout << "File " << nome_file << " not found." << endl;
    exit(1);
  }
  
  string line;
  // Skip two comment lines 
  getline (myfile,line);
  getline (myfile,line);
  

  
  TH1F *hLente_serieA = new TH1F("Lente_serieA","lente serie A",nbins, minb, maxb);
  TH1F *hLente_serieB = new TH1F("Lente_serieB","lente serie B",nbins, minb, maxb);
  
  char serie = 'X';

 
  double a,b;
 
  while (myfile>>a>>b) {
  
    hLente_serieA->Fill(a);
    hLente_serieB->Fill(b);

  }  
  myfile.close();

  float sx = 1;

  TCanvas *cX = new TCanvas("x","x",200,10,600,400);
  cX->cd();

  hLente_serieA->GetYaxis()->SetRangeUser(0,(hLente_serieA->GetMaximum() > hLente_serieB->GetMaximum() ?
					     hLente_serieA->GetMaximum() :
					     hLente_serieB->GetMaximum())*1.10);
  hLente_serieA->SetStats(kFALSE);
  hLente_serieA->GetXaxis()->SetTitle("x [mm]");
  hLente_serieA->GetYaxis()->SetTitle("Conteggi");
  

  hLente_serieA->SetLineColor(2);
  hLente_serieA->SetLineWidth(2);
  hLente_serieA->Draw();
  hLente_serieA->Fit("gaus","ME"); 

  TF1 *fitA = hLente_serieA->GetFunction("gaus");
  fitA->SetLineColor(1);
  //fitA->Draw("same");

  TF1 *fitB = new TF1("gausB","gaus",200,300);
  hLente_serieB->SetLineColor(4);
  hLente_serieB->SetLineStyle(9);
  hLente_serieB->Draw("same");
  hLente_serieB->Fit(fitB,"ME0");

  fitB->SetLineColor(1);
  fitB->Draw("same");

  TLegend *legenda = new TLegend(0.15,0.7,0.4,0.85);
  legenda->AddEntry(hLente_serieA,"serie A");
  legenda->AddEntry(hLente_serieB,"serie B");
  legenda->SetFillColor(0);
  legenda->Draw();
}
