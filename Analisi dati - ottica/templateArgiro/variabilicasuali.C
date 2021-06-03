/*******************************************
* Esperimentazioni II - 04/05/2016

  Mini-simulazione misure ripetute
  Situazione:
  - Quantità X da stimare.
  - Due modi indipendenti (con errori sistematici differenti) di misurare la stessa quantità,
    ad esempio mediante due tecniche diverse, due diversi esperimenti, due diversi operatori.
   
  La gestione dell opzioni non è elegante, ma solo funzionale all'esposizione della lezione
  opt = 0 (default) --> solo esperimento A
  opt = 1 --> solo esperimento A + fit gaussiano
  opt = 2 --> esperimento A e B sullo stesso canvas
  opt = 3 --> esperimento A e B sullo stesso canvas + fit gaussiani
  opt = 4 --> esperimento A e B sullo stesso canvas + fit gaussiani + test statistici 
  opt = 5 --> esperimento A e B sullo stesso canvas + fit gaussiani + test statistici + combinazione misure

  per eseguire (esempi):
  root -l 'variabilicasuali.C++(40)'
  genera 40 eventi e visualizza solo l'esprimento A
  root -l 'variabilicasuali.C++(100,3,40,200,300)'
  genera 100 eventi per esperimento, visualizza su un istrogramma con 40 bin con asse x che va da 200 a 300 entrambi i risultati con relativo fit

* R. Bellan - riccardo.bellan@unito.it
*******************************************/

#include <iostream>

#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TF1.h>
#include <TMath.h>

#include <TRandom3.h>

using namespace std;

double gaussWithPrecision(TRandom3 *random, const double &mean, const double &sigma, const double &precision){
  double raw = random->Gaus(mean,sigma);
  double rounded = round(raw/precision) * precision;
  return rounded;
}

void variabilicasuali (int neventi, int opt = 0, const double & nbins = 40, const double& minb = 200, const double& maxb = 300) {

  // Media e sigma della primo esperimento
  double meanGausA = 235;
  double sigmaGausA = 5;

  // Media e sigma della secondo esperimento
  // una differente media simula la presenza di un errore sistematico in una (o entrambe) le misure
  // una differente sigma simula la differente risposta dell'esperimento (e.g., due esperimentatori
  // possono avere abilità differenti nell'effettuare la misura, anche se utilizzano lo stesso strumento di misura, 
  // quindi le misure effettuate dai due esperimenti avranno differenti dispersioni.
  //double meanGausB = 270;
  //double sigmaGausB = 3;

  double meanGausB = 239;
  double sigmaGausB = 6;

  
  // Sensibilità dello strumento impiegato. Per semplicità la sensibilità dello strumento
  // è assunta essere uguale nei due esperimenti. La generalizzazione a due differenti risoluzioni è banale.
  double sx = 0.01;
  
  // --------------------------------------------------------- //

  // Generatore di eventi casuali (http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/ARTICLES/mt.pdf)
  TRandom3 *random = new TRandom3();

  TCanvas *cX = new TCanvas("x","x",200,10,600,400);
  cX->cd();

  // Istanza dell'oggetto istogramma (https://root.cern.ch/doc/master/classTH1F.html)
  TH1F *hCampioneA = new TH1F("CampioneA","Campione A",nbins, minb, maxb);
  hCampioneA->SetStats(kFALSE);
  hCampioneA->GetXaxis()->SetTitle("x [xyz]");
  hCampioneA->GetYaxis()->SetTitle("Conteggi");

  // Riempiamo l'istogramma con eventi casuali generati a partire da una gaussiana
  // con media=meanGausA e sigma=sigmaGausA. I valori sono discretizzati secondo la 
  // risoluzione sx dello strumento
  for (int i=0; i<neventi; ++i)
    hCampioneA->Fill(gaussWithPrecision(random,meanGausA,sigmaGausA,sx));

  // Secondo set di dati
  TH1F *hCampioneB = new TH1F("CampioneB","Campione B",nbins, minb, maxb);
  for (int i=0; i<neventi; ++i)
    hCampioneB->Fill(gaussWithPrecision(random,meanGausB,sigmaGausB,sx));

 
  hCampioneA->SetLineColor(2);
  hCampioneA->SetLineWidth(2);
  hCampioneA->Draw();


  // Ricava e stampa un po' di quantità inerenti l'istogramma A. Da notare che GetRMS^2 è la varianza con
  // bias, quindi bisogna correggere per N/(N-1) se il valore di aspettazione non è noto a priori
  double stdA = hCampioneA->GetRMS()*sqrt(hCampioneA->GetEntries()/(hCampioneA->GetEntries()-1));
  double nA   = hCampioneA->GetEntries();
  cout << "\n\n\nPopolazione A - media: " << hCampioneA->GetMean() << ", devizione standard popolazione: " <<  stdA << ", deviazione standard media: " << stdA/sqrt(nA) << endl;

  // La gestione dell opzioni non è elegante, ma solo funzionale all'esposizione della lezione
  if(opt == 0) return;
  TF1 *fitA = 0;
  if(opt != 2){
    hCampioneA->Fit("gaus","ME");
    fitA = hCampioneA->GetFunction("gaus");
    cout << "Chi^2:" << fitA->GetChisquare() << ", number of DoF: " << fitA->GetNDF() << " (Probability: " << fitA->GetProb() << ")." << endl;
    fitA->SetLineColor(1);
  }

  if(opt == 1) return;
  hCampioneB->Draw("same");
  // Ricava e stampa un po' di quantità inerenti l'istogramma B. Da notare che GetRMS^2 è la varianza con
  // bias, quindi bisogna correggere per N/(N-1) se il valore di aspettazione non è noto a priori
  double stdB = hCampioneB->GetRMS()*sqrt(hCampioneB->GetEntries()/(hCampioneB->GetEntries()-1));
  double nB  = hCampioneB->GetEntries();
  cout << "Popolazione B - media: " << hCampioneB->GetMean() << ", devizione standard popolazione: " <<  stdB << ", deviazione standard media: " << stdB/sqrt(nB) << endl;

  // Per visualizzare i due fit sullo stesso canvas bisogna istanziare a parte la seconda funzione di fit
  TF1 *fitB = new TF1("gausB","gaus",minb,maxb);
  hCampioneB->SetLineColor(4);
  hCampioneB->SetLineStyle(9);
  if(opt != 2){
    hCampioneB->Fit(fitB,"ME0"); // opzione "0" necessaria per non disegnare il fit che altrimenti sovrascriverebbe il canvas
    cout << "Chi^2:" << fitB->GetChisquare() << ", number of DoF: " << fitB->GetNDF() << " (Probability: " << fitB->GetProb() << ")." << endl;
    fitB->SetLineColor(1);
    fitB->Draw("same");
  }

  // Per non confodersi, disegnamo una legenda
  TLegend *legenda = new TLegend(0.15,0.7,0.35,0.85);
  legenda->AddEntry(hCampioneA,"esperimento A");
  legenda->AddEntry(hCampioneB,"esperimento B");
  legenda->SetFillColor(0);
  legenda->Draw();
  
  // Accorgimento per avere un canvas che contiene completamente entrambi gli istogrammi
  hCampioneA->GetYaxis()->SetRangeUser(0,(hCampioneA->GetMaximum() > hCampioneB->GetMaximum() ? hCampioneA->GetMaximum() : hCampioneB->GetMaximum())*1.10);

  if(opt == 2 || opt == 3) return;

  // Controlliamo la compatibilità dei due set di misure

  // p_value = int [x,inf] f(t|H0) dt : assumendo H0, probabilita' di trovare un valore piu' grande di quello osservato
  
  // Test gaussiano
  cout << "\n\n--- Test Normale ---" << endl;
  cout << "Popolazione A - media: " << hCampioneA->GetMean() << " sigma popolazione: " << stdA << " incertezza sulla media (sigma media): " << stdA/sqrt(nA) << endl;
  cout << "Popolazione B - media: " << hCampioneB->GetMean() << " sigma popolazione: " << stdB << " incertezza sulla media (sigma media): " << stdB/sqrt(nB) << endl;
  double diffMedieGauss  = hCampioneA->GetMean()-hCampioneB->GetMean();
  double sdiffMedieGauss = sqrt(pow(stdA,2)/nA + pow(stdB,2)/nB);
  double z = fabs(diffMedieGauss)/sdiffMedieGauss;
  cout << "Scarto: " << diffMedieGauss << " +- " << sdiffMedieGauss << endl;  
  cout << "z = " << z << ". Probabilità di avere discrepanza maggiore di quella osservata: " << TMath::Erfc(z/sqrt(2.)) << endl << endl;
  // Nota. Erfc(x) = 2/sqrt(pi) Int [x,inf] exp(-t^2) dt. Quindi x = z/sqrt(2).
  // gaussian : exp(-1/2* ((x-mu)/sigma)^2)
  // erf(a/(sigma * sqrt(2)) esprime la probabilità che  una singola misura si trovi fra -a e +a. (mu =0)

  
  // Test di Student
  cout << "--- Test di Student ---" << endl;
 
 
  double diffMedie = hCampioneA->GetMean() - hCampioneB->GetMean();
  double ndof = nA+nB-2;
  double sxAxB    = sqrt( (1./nA + 1./nB) * ((nA-1)*stdA*stdA + (nB-1)*stdB*stdB)/ndof );
  double tStudent = fabs(diffMedie)/sxAxB;

  cout << "Scarto: " << diffMedie << " +- " << sxAxB << endl; 
  cout << "t-Student: " << tStudent << " ndof: " << ndof  << ". Probabilità di avere discrepanza maggiore di quella osservata: " << 2*(1-TMath::StudentI(tStudent,ndof)) << endl;
  // Sotto la null hypothesis che le medie siano uguali, osservo un t piu' grande nel p*100% dei casi.
  // Nota. StudentI(x,ndof) integra da -inf a x 
  
  // Welch t-test
  cout << "\n\n--- Test di Welch ---" << endl;
  double sxAxBW = sqrt(stdA*stdA/nA + stdB*stdB/nB);
  double ndofW  = pow(sxAxBW,4)/( pow(stdA*stdA/nA,2)/(nA-1) + pow(stdB*stdB/nB,2)/(nB-1));
  double welch  =   fabs(diffMedie)/sxAxBW;
  cout << "Welch: " << welch  << " ndof: " << ndofW  << ". Probabilità di avere discrepanza maggiore di quella osservata: " << 2*(1-TMath::StudentI(welch,ndofW)) << endl;

  if(opt < 5) return;


  // Media pesata e semidifferenza come incertezza sistematica
  cout <<"\n\nMedia pesata e semidifferenza come incertezza" << endl;
  double wA = 1./(fitA->GetParError(1)*fitA->GetParError(1));
  double wB = 1./(fitB->GetParError(1)*fitB->GetParError(1));
  double smediaPesata = 1/sqrt(wA+wB);
  double mediaPesata = (fitA->GetParameter(1)*wA + fitB->GetParameter(1)*wB)/(wA+wB);
  double semidif  = fabs(fitA->GetParameter(1)-fitB->GetParameter(1))/2.;
  cout << "Misura: " << mediaPesata << " +- " << smediaPesata << " (stat.) +- " << semidif << " (sist.)" << endl;


  // Birge ratio (Chi2 a la PDG)
  // http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf pagina 10
  // https://ia601302.us.archive.org/35/items/numericalcompari8124tayl/numericalcompari8124tayl.pdf
  cout << "\nBirge ratio (a la PDG)" << endl;
  double chi2 = pow(mediaPesata - fitA->GetParameter(1),2)*wA + pow(mediaPesata - fitB->GetParameter(1),2)*wB;
  double fscala = sqrt(chi2);
  cout << "Chi2: " << chi2 << " fattore di scala " << fscala << endl;
  cout << "Misura: " << mediaPesata << " +- " << fscala*smediaPesata << endl;


  /* Altre possibilità di combinare le misure
  // Envelope
  cout << "\nEnvelope" << endl;

  double vals[4] =  {fitA->GetParameter(1)+fitA->GetParError(1),fitA->GetParameter(1)-fitA->GetParError(1),
		     fitB->GetParameter(1)+fitB->GetParError(1),fitB->GetParameter(1)-fitB->GetParError(1)};
    

  //double min = std::min( fitA->GetParameter(1), fitB->GetParameter(1));
  double max = *std::max_element(vals,vals+4);
  double min = *std::min_element(vals,vals+4);

  double media     = (max+min)/2.;
  double inviluppo = (max-min)/2.;
  cout << "Misura: " << media << " +- " << inviluppo << endl;


  // Envelope + media pesata
  cout << "\nEnvelope con media pesata" << endl;
  cout << "Misura: " << mediaPesata << " + " << max-mediaPesata << " - " << mediaPesata-min  << endl;
  */


}
