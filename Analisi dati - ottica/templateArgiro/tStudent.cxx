/*******************************************
* Esperimentazioni II - 04/05/2016

  Mini-simulazione misure ripetute, distribuzione t-Student
  Situazione:
  - N esperimenti che misurano la quanità x
  - ogni esperimento fa M misure. 

  Per M piccoli (<40) si osserva che la compatibilità tra valore vero e stimato è descritta dalla statistica di Student.
  Per M grandi la statistica diventa normale.  
  Uso: Testare l'ipotesi nulla che il valore vero della media della popolazione sia mu_0

  per eseguire (esempio):
  root -l 'tStudent.C++(50000,5,200,-5,5)'
  genera 5000 esperimenti con 5 misurazioni ciascuno



* R. Bellan - riccardo.bellan@unito.it
*******************************************/


#include <iostream>
#include <TCanvas.h>
#include <TH1F.h>
#include <TAxis.h>
#include <TF1.h>

#include <TRandom3.h>

using namespace std;


void tStudent(int numeroEsperimenti, int neventi, const double & nbins = 40, const double& minb = -10, const double& maxb = 10) {

  // valor vero e sigma
  double valor_vero  = 235;
  double sigma = 5;

  // --------------------------------------------------------- //

  // Generatore di eventi casuali (http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/ARTICLES/mt.pdf)
  TRandom3 *random = new TRandom3();

  double variabili_t[numeroEsperimenti];

  // Referenza calcolo varianza in un passaggio: http://www.jstor.org/stable/1266577?seq=1#page_scan_tab_contents
  for (int esp = 0; esp < numeroEsperimenti; ++esp){
    double delta = 0.;
    double media = 0.;
    double M2    = 0.;
    for (int i=0; i<neventi; ++i){
      double d = random->Gaus(valor_vero,sigma);
      delta = d - media;
      media += delta/(i+1);
      M2    += delta*(d-media);
    }
    double std = neventi > 1 ? sqrt(M2/(neventi-1)) : -9.;
    variabili_t[esp] = (media - valor_vero)/(std/sqrt(neventi));
    // Se anziché usare la deviazione standard campionaria si usasse la deviazione standard vera, allora t si distribuirebbe secondo la statistica di
    // Gauss anche per poche misure. Per verificarlo, provare a usare la riga qui sotto anziché quella sopra.
    // variabili_t[esp] = (media - valor_vero)/(sigma/sqrt(neventi));
  }

  TCanvas *ct = new TCanvas("t","t",200,10,600,400);
  ct->cd();

  // Istanza dell'oggetto istogramma (https://root.cern.ch/doc/master/classTH1F.html)
  TH1F *htStudent = new TH1F("sStudent","tStudent",nbins, minb, maxb);
  htStudent->SetStats(kFALSE);
  htStudent->GetXaxis()->SetTitle("t");
  htStudent->GetYaxis()->SetTitle("Conteggi");
  
  for(int esp=0; esp<numeroEsperimenti;++esp)
    htStudent->Fill(variabili_t[esp]);
  

  htStudent->Draw();
  htStudent->Fit("gaus","ME");
  TF1 *fit = htStudent->GetFunction("gaus");
  cout << "Chi^2:" << fit->GetChisquare() << ", number of DoF: " << fit->GetNDF() << " (Probability: " << fit->GetProb() << ")." << endl;
}

