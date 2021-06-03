#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include "statlib.cxx"
#include "csvdata.cxx"
#include "espToolkit.cxx"
using namespace std;

statlib st;

//Percorsi Valerio
const char rootfilepath[256]      = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper5.root";
const char outfilepath[256]       = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper5_results.txt";
const char IVledCsvp[256]          = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 4-5/I_V_LED.csv";


//Percorsi Federica
//const char rootfilepath[256]      = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper5.root";
//const char outfilepath[256]       = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper5_results.txt";
//const char IVledCsvp[256]          = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 4-5/I_V_LED.csv";

//Percorsi Filippo




void esper5_1()
{
    cerr << "\nStart esper1.cxx data analysis";
    //Avviamo il timer che misura il tempo di calcolo impiegato dalla CPU per completare l'analisi
    auto start = std::chrono::system_clock::now();

    //Apriamo un nuovo file di output in formato ROOT (.root) dove salvare i plots
    TFile out_file(rootfilepath, "RECREATE");
    
    //Apriamo un nuovo file di testo in cui salvare l'output dell'analisi
    fstream printout;
    printout.open(outfilepath, std::fstream::out);
    printout << setprecision(4);
    

    //Stampiamo l'intestazione 
    printout << endl;
    printout << "  +---------------------------------------------------------------+ \n";
    printout << "  |  Oreglia, Sibilla, Pagliarino - Corso B - C.d.L. in Fisica    | \n";
    printout << "  |           Universita' degli Studi di Torino -  ESP2           | \n";
    printout << "  +---------------------------------------------------------------+ \n";
    printout << "  |             ANALISI DATI  - C++11 + CERN ROOT 6               | \n";
    printout << "  +---------------------------------------------------------------+ \n";
    printout << "  | Experiment:  Ottica - esperienza 5                            | \n";
    printout << "  | Date:        13/12/3030                                       | \n";
    printout << "  | Revision:    2.1.0                                            | \n";
    printout << "  | Description: Effetto Fotoelettrico                            | \n";
    printout << "  +---------------------------------------------------------------+ \n\n";


    //Apriamo il file CSV esistente in cui è stata salvata la tabella con i dati dello spettroscopio
    csvdata IVledCsv(IVledCsvp);
    IVledCsv.setDelimiters(';');

    //---------------------------- IMPORTAZIONE DATI SPERIMENTALI --------------------------------

    const double c = 299792458/1.0003;           //m/s 

    double lam1 = 639.96e-9;   //nm     
    double lam2 = 475.35e-9;   //nm    
    double lam3 = 425.57e-9;   //nm    
    double lam4 = 549.21e-9;   //nm    
    double lam5 = 605.05e-9;   //nm     
    double lam[5] = {lam1, lam2, lam3, lam4, lam5};

    double nu[6];
    nu[0] = (c) / (lam1);     //GHz 
    nu[1] = (c) / (lam2);     //GHz 
    nu[2] = (c) / (lam3);     //GHz 
    nu[3] = (c) / (lam4);     //GHz 
    nu[4] = (c) / (lam5);     //GHz    

    double lam_err[5] = {13.54e-9/2., 26.19e-9/2., 77.00e-9 /2., 40.11e-9 / 2., 14.39e-9 / 2.};
    double nu_err[5];



    for (int i = 0; i < 5; i++)
    {
        nu_err[i] = c / (lam[i]*lam[i]) * lam_err[i];
        printout << "\n" << nu_err[i] << "  " << nu[i];
    }

    int N_mis = 12;

    double V_lam1[N_mis];    //V
    double V_lam2[N_mis];    //V
    double V_lam3[N_mis];    //V
    double V_lam4[N_mis];    //V
    double V_lam5[N_mis];    //V

    double V_l1_err[N_mis];     //V
    double V_l2_err[N_mis];     //V
    double V_l3_err[N_mis];     //V
    double V_l4_err[N_mis];     //V
    double V_l5_err[N_mis];     //V

    double I_lam1[N_mis];    //A
    double I_lam2[N_mis];    //A
    double I_lam3[N_mis];    //A
    double I_lam4[N_mis];    //A
    double I_lam5[N_mis];    //A
  
    double I_l1_err[N_mis];     //A
    double I_l2_err[N_mis];     //A
    double I_l3_err[N_mis];     //A
    double I_l4_err[N_mis];     //A
    double I_l5_err[N_mis];     //A

    /*
   

    for (int i = 2; i < N_mis + 2; i++)
    {
        V_lam1[i - 2] = IVledCsv.getDouble(i, 0 )*1e-3;
        V_l1_err[i - 2]  = IVledCsv.getDouble(i, 1 )*1e-3;
        I_lam1[i - 2] = IVledCsv.getDouble(i, 2 )*1e-9;
        I_l1_err[i - 2]  = IVledCsv.getDouble(i, 3 )*1e-9;     
        V_lam2[i - 2] = IVledCsv.getDouble(i, 4 )*1e-3;
        V_l2_err[i - 2]  = IVledCsv.getDouble(i, 5 )*1e-3;
        I_lam2[i - 2] = IVledCsv.getDouble(i, 6 )*1e-9;
        I_l2_err[i - 2]  = IVledCsv.getDouble(i, 7 )*1e-9;
        V_lam3[i - 2] = IVledCsv.getDouble(i, 8 )*1e-3;
        V_l3_err[i - 2]  = IVledCsv.getDouble(i, 9 )*1e-3;
        I_lam3[i - 2] = IVledCsv.getDouble(i, 10)*1e-9;
        I_l3_err[i - 2]  = IVledCsv.getDouble(i, 11)*1e-9;
        V_lam4[i - 2] = IVledCsv.getDouble(i, 12)*1e-3;
        V_l4_err[i - 2]  = IVledCsv.getDouble(i, 13)*1e-3;
        I_lam4[i - 2] = IVledCsv.getDouble(i, 14)*1e-9;
        I_l4_err[i - 2]  = IVledCsv.getDouble(i, 15)*1e-9;
        V_lam5[i - 2] = IVledCsv.getDouble(i, 16)*1e-3;
        V_l5_err[i - 2]  = IVledCsv.getDouble(i, 17)*1e-3;
        I_lam5[i - 2] = IVledCsv.getDouble(i, 18)*1e-9;
        I_l5_err[i - 2]  = IVledCsv.getDouble(i, 19)*1e-9;

    
    }    
    */

    //------------------------------- ANALISI DATI SPERIMENTALI ----------------------------------
    
    /*TGraphErrors * t[5];
    t[0] =  new TGraphErrors(N_mis, V_lam1, I_lam1, V_l1_err, I_l1_err);
    t[1] =  new TGraphErrors(N_mis, V_lam2, I_lam2, V_l2_err, I_l2_err);
    t[2] =  new TGraphErrors(N_mis, V_lam3, I_lam3, V_l3_err, I_l3_err);
    t[3] =  new TGraphErrors(N_mis, V_lam4, I_lam4, V_l4_err, I_l4_err);
    t[4] =  new TGraphErrors(N_mis, V_lam5, I_lam5, V_l5_err, I_l5_err);
    for (int j = 0; j < 5; ++j)
    {
        t[j]->SetName((std::string("t") + std::to_string(j)).c_str());
    }*/

    /*
    
    TGraphErrors * IV_1= new TGraphErrors(N_mis, V_lam1, I_lam1, V_l1_err, I_l1_err);
    TGraphErrors * IV_2= new TGraphErrors(N_mis, V_lam2, I_lam2, V_l2_err, I_l2_err);
    TGraphErrors * IV_3= new TGraphErrors(N_mis, V_lam3, I_lam3, V_l3_err, I_l3_err);
    TGraphErrors * IV_4= new TGraphErrors(N_mis, V_lam4, I_lam4, V_l4_err, I_l4_err);
    TGraphErrors * IV_5= new TGraphErrors(N_mis, V_lam5, I_lam5, V_l5_err, I_l5_err);

    TCanvas * cIV = new TCanvas("cIV","Curve I(V)",10,3,600,400);
    cIV->cd();
    cIV->SetTitle("Curve I(V)");

    TMultiGraph  *mg  = new TMultiGraph();
    */
    /*for (int i = 0; i < 5; i++)
    {
        t[i]->Sort();
        mg->Add(t[i]);
    }*/
    /*

    IV_1->Sort();
    IV_2->Sort();
    IV_3->Sort();
    IV_4->Sort();
    IV_5->Sort();

    IV_1->Write("IV1");
    IV_2->Write("IV2");
    IV_3->Write("IV3");
    IV_4->Write("IV4");
    IV_5->Write("IV5");
    
    mg->Add(IV_1);
    mg->Add(IV_2);
    mg->Add(IV_3);
    mg->Add(IV_4);
    mg->Add(IV_5);
   

    mg->GetXaxis()->SetTitle("V controcampo [V]");
    mg->GetYaxis()->SetTitle("Fotocorrente [A]");
    mg ->Draw();
    mg->Write("Multigraph");

    cIV->Draw();
    */
   
 
    /*//force drawing of canvas to generate the fit TPaveStats
    cIV->Update();
    auto stats1 = (TPaveStats*)IV_1->GetListOfFunctions()->FindObject("stats");
    auto stats2 = (TPaveStats*)IV_2->GetListOfFunctions()->FindObject("stats");
    auto stats3 = (TPaveStats*)IV_3->GetListOfFunctions()->FindObject("stats");
    auto stats4 = (TPaveStats*)IV_4->GetListOfFunctions()->FindObject("stats");
    auto stats5 = (TPaveStats*)IV_5->GetListOfFunctions()->FindObject("stats");
    auto stats6 = (TPaveStats*)IV_6->GetListOfFunctions()->FindObject("stats");
    stats1->SetTextColor(kBlue);
    stats2->SetTextColor(kRed);
    cIV->Modified();*/

    /*cIV->SetGrid();
    cIV->Draw();
    cIV->Write("cIV");*/

    //---------------------------------- ricerca dei Vc0 -----------------------------------------

    double Vc0[5] = {0.6369, 1.2973, 1.456, 1.0157, 0.7716};
    double Vc0_err[5] = {0.026, 0.08, 0.076, 0.10, 0.14};

    TCanvas * fitVc0 = new TCanvas("fitVc0","Potenziali d'arresto",10,3,600,400);
    fitVc0->cd();
    fitVc0->SetTitle("Potenziali d'arresto");
     
    /*TF1 * a[5];

    a[0] = new TF1("a1", "[0]*exp([1]*(x-[2]))-[0]", -0.78, -0,41);
    a[1] = new TF1("a2", "[0]*exp([1]*(x-[2]))-[0]", -1.47, -0.96);
    a[2] = new TF1("a3", "[0]*exp([1]*(x-[2]))-[0]", -1.5, -1.1);
    a[3] = new TF1("a4", "[0]*exp([1]*(x-[2]))-[0]", 2.4, 6.8);
    a[4] = new TF1("a5", "[0]*exp([1]*(x-[2]))-[0]", 2.4, 6.8);

    for (int j = 0; j < 5; ++j)
    {
        printout << "\n\n----- a[" << j << "] -----";
        a[j]->SetParameter(0, 1e-10);
        a[j]->SetParameter(1, 10);
        a[j]->SetParameter(2, -1);
        return;
        t[j]->Fit(a[j], "r");
        Vc0[j]     =  abs(a[j]->GetParameter(2));
        Vc0_err[j] =  a[j]->GetParError(2);
    

        t[j]->Write((std::string("I(V)") + std::to_string(j)).c_str());
    }*/
    
    /*TF1 * schottky1    = new TF1("schottky1", "[0]*(exp([1]*(x-[2]))-1)",-0.78,-0.41); //  [0]*(exp([1]*(x-[2]))-1)
    schottky1->SetParameter(0, 1e-11); 
    schottky1->SetParameter(1, 17); 
    schottky1->SetParameter(2, -0.6); 
    IV_1   -> Fit("schottky1","R");
    Vc0[0]     =  abs(schottky1->GetParameter(2));
    Vc0_err[0] =  schottky1->GetParError(2);
    IV_1->Draw();

    
    TF1 * schottky2    = new TF1("schottky2", "[0]*(exp([1]*(x-[2]))-1)",-1.47,-0.96); //  [0]*(exp([1]*(x-[2]))-1)
    schottky2->SetParameter(0, 1e-9); 
    schottky2->SetParameter(1, 10); 
    schottky2->SetParameter(2, -1.3); 
    IV_2->Fit("schottky2","R");
    Vc0[1]     =  abs(schottky2->GetParameter(2));
    Vc0_err[1] =  schottky2->GetParError(2);
    IV_2->Draw();

    TF1 * schottky3    = new TF1("schottky3", "[0]*(exp([1]*(x-[2]))-1)",-1.5,-1.1); //  [0]*(exp([1]*(x-[2]))-1)
    schottky3->SetParameter(0, 1e-9); 
    schottky3->SetParameter(1, 0); 
    schottky3->SetParameter(2, -1); 
    IV_3->Fit("schottky3","R");
    Vc0[2]     =  abs(schottky3->GetParameter(2));
    Vc0_err[2] =  schottky3->GetParError(2);
    IV_3->Draw();

    TF1 * schottky4    = new TF1("schottky4", "[0]*(exp([1]*(x-[2]))-1)",-1.2,-0.7); //  [0]*(exp([1]*(x-[2]))-1)
    schottky4->SetParameter(0, 1e-9); 
    schottky4->SetParameter(1, 0); 
    schottky4->SetParameter(2, -0.5); 
    IV_4->Fit("schottky4","R");
    IV_4->Draw();

    TF1 * schottky5    = new TF1("schottky5", "[0]*(exp([1]*(x-[2]))-1)",-0.8,-0.3); //  [0]*(exp([1]*(x-[2]))-1)
    schottky5->SetParameter(0, 1e-9); 
    schottky5->SetParameter(1, 0); 
    schottky5->SetParameter(2, -0.5); 
    IV_5->Fit("schottky5","R");
    IV_5->Draw();*/

    

    printout << "\n\n";
    printout << "\n\nVc0:\n";
    printout << "-----------------------------";
    for (int i = 0; i < 5; ++i)
    {
        printout << "\n = " << Vc0[i] << "  pm  " << Vc0_err[i];  
    }

    printout << "\n\nnu:\n";
    printout << "-----------------------------";

    for (int i = 0; i < 5; ++i)
    {
        printout << "\n = " << nu[i] << "  pm  " << nu_err[i];  
    }

    fitVc0->Draw();

    //--------------------------------- COSTANTE DI PLANK -----------------------------------------

    TCanvas * cVlam = new TCanvas("cVlam","Retta Vc0(f)",10,3,600,400);
    cVlam->cd();
    cVlam->SetTitle("Vc0(f)");


    TGraphErrors * V_lam = new TGraphErrors(5, nu, Vc0 , nu_err, Vc0_err); 
    V_lam->Sort();

    TF1 * retta   = new TF1("retta", "[0]*x+[1]",0,5e12);
    retta->SetParameter(0, 1e-15); 
    //retta->SetParameter(1, 1e10); 
    V_lam->Fit(retta);
    V_lam->Draw();

    cVlam->SetGrid();
    V_lam->GetXaxis()->SetTitle("Frequenza [Hz]");
    V_lam->GetYaxis()->SetTitle("Potenziale d'arresto [V]");
    cVlam->Draw();
    cVlam->Write("cVlam");

    double e = 1.602176634e-19; //C

    printout << "\n\n";
    printout << "\n\nFit retta Vc0(nu):\n";
    printout << "-----------------------------";
    printout << "\nchi^2 = " << retta->GetChisquare();    
    printout << "\nndf   = " << retta->GetNDF();           
    printout << "\np-val = " << retta->GetProb();       
    printout << "\np1    = " << retta->GetParameter(0)   << " +/- " <<  retta->GetParError(0) << "J*s / C";  
    printout << "\nh     = " << retta->GetParameter(0) * e   << " +/- " <<  retta->GetParError(0) * e << "J*s";
    printout << "\nq     = " << retta->GetParameter(1)   << " +/- " <<  retta->GetParError(1) << "V";
    double h     = retta->GetParameter(0) * e;
    double h_err = retta->GetParError(0)  * e;


    double h_ref = 6.62607015e-34;
    double z     = (h - h_ref) / h_err;
    double pval  = st.pvalZtwotailed(z);

    printout << "\n\n\nTest Z compatibilità con valore di riferimento:\nZ = " << z << "\np = " << pval << "\nref =" << h_ref << " J/Hz";

    /*
    root [1] z = (6.567e-34 - 6.62607015e-34) / 1.491e-34
    (double) -0.0396178
    root [2] st.pvalZtwotailed(z)
    (double) 0.968398
    root [3] 
    https://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf
    */


    //--------------- SCRITTURA DEL TEMPO MACCHINA E CHIUSURA DEI FILES DI OUTPUT ----------------

    printout << "\n\n";
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    printout << "\nFinished computation at " << std::ctime(&end_time)
          << "elapsed time: " << elapsed_seconds.count() << "s\n";

    out_file.Close(); 
    printout.close();
    cout << "\n\n";

    std::ifstream f(outfilepath);
    if (f.is_open())
    {
        std::cout << f.rdbuf();
    }

/*
    //----------------- PREPARAZIONE DEL PACCHETTO PER LA CONSEGNA ANALISI ONLINE ----------------
    std::string cmd = std::string("zip " + archive + " " +
                                  " " + rootfilepath + " " + outfilepath + 
                                  " " + file1 + " " + file2 + " " + file3 + " " +file4);
    system(cmd.c_str());
*/

}