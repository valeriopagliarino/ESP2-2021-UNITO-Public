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
const char rootfilepath[256]      = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper1.root";
const char outfilepath[256]       = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper1_results.txt";
const char HgLine[256]            = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 1/E1_HgLines.csv";
const char nGlass[256]            = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 1/E1_nGlass.csv";


//Percorsi Federica
//const char rootfilepath[256]      = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper1.root";
//const char outfilepath[256]       = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper1_results.txt";
//const char HgLine[256]            = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 1/E1_HgLines.csv";
//const char nGlass[256]            = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 1/E1_nGlass.csv";

//Percorsi Filippo


void esper1_demo()
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
    printout << "  | Experiment:  Ottica - esperienza 1                            | \n";
    printout << "  | Date:        19/04/2021                                       | \n";
    printout << "  | Revision:    2.1.0                                            | \n";
    printout << "  | Description: Spettroscopia                                    | \n";
    printout << "  +---------------------------------------------------------------+ \n\n";


    //Apriamo il file CSV esistente in cui è stata salvata la tabella con i dati dello spettroscopio
    csvdata HgLineCsv(HgLine);
    HgLineCsv.setDelimiters(';');
    csvdata nGlassCsv(nGlass);
    nGlassCsv.setDelimiters(';');

    //---------------------------- IMPORTAZIONE DATI SPERIMENTALI --------------------------------
    
    double alpha        = 45.;           //deg   Angolo al vertice del prisma
    double grating      = 1. / 0.00057; //m     Distanza tra le righe del reticolo (passo)
    double lampPower    = 0;            //W     Potenza della lampada spettrale (mercurio)

    int N_HeLines       = 6;
    int N_nGlass        = N_HeLines;

    //Incertezze:
    double  errLetturaAngolare = 0.5;   //deg
    double  errLambdaDig       = 0.0;   //nm
    double  errGrating         = 0.0;

    //N_HeLines
    double  lineID_He[N_HeLines];
    double  sx4[N_HeLines]; //m
    double  sx3[N_HeLines]; //m
    double  sx2[N_HeLines]; //m
    double  sx1[N_HeLines]; //m
    double  dx1[N_HeLines]; //m
    double  dx2[N_HeLines]; //m
    double  dx3[N_HeLines]; //m
    double  dx4[N_HeLines]; //m
    double  ldm[N_HeLines]; //nm

    //N_nGlass
    double  lineID_Glass[N_nGlass];
    double  lambda[N_nGlass];           //nm
    double  minDevAngle[N_nGlass];      //deg
    double  minDevAngle_err[N_nGlass];  //deg

    for (int i = 1; i < N_HeLines + 1; i++)
    {
        lineID_He[i - 1] = HgLineCsv.getDouble(i, 0);
        sx4[i - 1] = abs(HgLineCsv.getDouble(i, 2 ));
        sx3[i - 1] = abs(HgLineCsv.getDouble(i, 3 ));
        sx2[i - 1] = abs(HgLineCsv.getDouble(i, 4 ));
        sx1[i - 1] = abs(HgLineCsv.getDouble(i, 5 ));
        dx1[i - 1] = abs(HgLineCsv.getDouble(i, 6 ));
        dx2[i - 1] = abs(HgLineCsv.getDouble(i, 7 ));
        dx3[i - 1] = abs(HgLineCsv.getDouble(i, 8 ));
        dx4[i - 1] = abs(HgLineCsv.getDouble(i, 9 ));
        ldm[i - 1] = abs(HgLineCsv.getDouble(i, 10));
    }

    for (int i = 1; i < N_nGlass + 1; ++i)
    {  
        lineID_Glass[i-1]    = nGlassCsv.getDouble(i, 0);
        lambda[i-1]          = nGlassCsv.getDouble(i, 2);
        minDevAngle[i-1]     = nGlassCsv.getDouble(i, 3);
        minDevAngle_err[i-1] = nGlassCsv.getDouble(i, 4);
    }


    //------------------------------- ANALISI DATI SPERIMENTALI ----------------------------------


    //----------------------------------- Righe del mercurio -------------------------------------

    printout << "\n\n--------- Righe del mercurio ---------\n\n";
    TCanvas * C1 = new TCanvas("C1","Linee del Mercurio",10,3,600,400);
    C1->cd();
    C1->SetTitle("Linee del Mercurio");

    auto * g1 = new TGraph(1);
    g1->SetPoint(1, 0.,0.);
    g1->GetXaxis()->SetLimits(-50, 50); //Fondo scala plot delle righe del mercurio
    g1->GetYaxis()->SetLimits(0, 1.);
    g1->Draw();
    C1->SetGrid();
    g1->GetXaxis()->SetTitle("angolo in gradi");
    g1->SetTitle("Linee del mercurio");
    TLine * line1;
    TLine * line2;
    TLine * line3;
    TLine * line4;
    TLine * line5;
    TLine * line6;
    TLine * line7;
    TLine * line8;

    for (int i = 0; i < N_HeLines; ++i)
    {
        line1 = new TLine(-1. * sx4[i],0.,-1. * sx4[i],1.);
        line2 = new TLine(-1. * sx3[i],0.,-1. * sx3[i],1.);
        line3 = new TLine(-1. * sx2[i],0.,-1. * sx2[i],1.);
        line4 = new TLine(-1. * sx1[i],0.,-1. * sx1[i],1.);
        line5 = new TLine(dx1[i],0.,dx1[i],1.);
        line6 = new TLine(dx2[i],0.,dx2[i],1.);
        line7 = new TLine(dx3[i],0.,dx3[i],1.);
        line8 = new TLine(dx4[i],0.,dx4[i],1.);

        line1->SetLineColor(i+1);
        line2->SetLineColor(i+1);
        line3->SetLineColor(i+1);
        line4->SetLineColor(i+1);
        line5->SetLineColor(i+1);
        line6->SetLineColor(i+1);
        line7->SetLineColor(i+1);
        line8->SetLineColor(i+1);

        line1->Draw();
        line2->Draw();
        line3->Draw();
        line4->Draw();
        line5->Draw();
        line6->Draw();
        line7->Draw();
        line8->Draw();
    }

    C1->Draw();
    C1->Write("Righe mercurio");

    double lambda_r[N_HeLines];

    for (int i = 0; i < N_HeLines; ++i)
    {
        double lambda_ord1 = sin(abs(dx1[i] + sx1[i]) / 2.) / (1. * (1./grating));
        double lambda_ord2 = sin(abs(dx2[i] + sx2[i]) / 2.) / (2. * (1./grating));
        double lambda_ord3 = sin(abs(dx3[i] + sx3[i]) / 2.) / (3. * (1./grating));
        double lambda_ord4 = sin(abs(dx4[i] + sx4[i]) / 2.) / (4. * (1./grating));
        lambda_r[i] = (1./2.) * (abs(lambda_ord1) + abs(lambda_ord2));// + lambda_ord3 + lambda_ord4);
        printout << "\n" << i+1 << "] [m]  ordine 1 (" << lambda_ord1 << ")  ordine 2 (" << lambda_ord2 << ")  ordine 3 (" << lambda_ord3 << ")   ordine 4 (" << lambda_ord4 << ")";
        printout << "\n         media = " << lambda_r[i];
    }

    TCanvas * C2 = new TCanvas("C2","Linee del Mercurio",10,3,600,400);
    C2->cd();
    C2->SetTitle("Linee del Mercurio");

    auto * g2 = new TGraph(1);
    g2->SetPoint(1, 0.,0.);
    g2->GetXaxis()->SetLimits(1e-9, 1500e-9); //Fondo scala plot delle righe del mercurio
    g2->GetYaxis()->SetLimits(0, 1.);
    g2->Draw();
    C2->SetGrid();
    g2->GetXaxis()->SetTitle("Lunghezza d'onda");
    g2->SetTitle("Linee del mercurio");
    TLine * line1b;

    for (int i = 0; i < N_HeLines; ++i)
    {
        line1b = new TLine(lambda_r[i],0.,lambda_r[i],1.);
        line1b->Draw();
    }
    C2->Draw();
    C2->Write("Lambda mercurio");


    //-------------------------MISURA DI n(lambda) NEL PRISMA IN VETRO----------------------------

    double glass_indices[N_nGlass];

    printout << "\n\nIndice di rifrazione del vetro:";
    
    for (int i = 0; i < N_nGlass; ++i)
    {
        glass_indices[i] = sin((minDevAngle[i] + alpha) / 2.) / (sin(alpha / 2.));
        printout << "\n" << i + 1 << "] n = " << glass_indices[i];
    }
    
    auto * t1 = new TGraphErrors(N_nGlass);
    for (int i = 0; i < N_nGlass; ++i)
    {
        t1->SetPoint(i, lambda_r[i], glass_indices[i]);
        //Errore da propagare...
    }

    TCanvas * C3 = new TCanvas("C3","Indice di rifrazione del vetro",10,3,600,400);
    C3->cd();
    C3->SetTitle("Indice di rifrazione del vetro");
    t1->SetTitle("Indice di rifrazione del vetro");
    t1->GetXaxis()->SetTitle("Lunghezza d'onda");
    t1->GetYaxis()->SetTitle("Indice di rifrazione");
    C3->SetGrid();
    t1->Draw();
    C3->Draw();
    t1->Write("n_glass");
    C3->Write("n_glass_canvas");

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


//PROCEDURA:
/*
1) Montare l'apparato (forse) come specificato nel "LD Physics Leaflets".
2) Controllare che le linee a dx e a sx siano in posizione simmetrica 
    (ruotare il piatto della semidifferenza della distanza angolare tra due massmi)
3) Assicurarsi che il fascio sia perpendicolare al reticolo di diffrazione
4) Segnare le misure sul file CSV per tutte le linee nei due ordini, a destra e a sinistra

*/