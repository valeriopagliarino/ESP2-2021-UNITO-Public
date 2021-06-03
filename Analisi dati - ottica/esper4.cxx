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
const char rootfilepath[256]      = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper4.root";
const char outfilepath[256]       = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper4_results.txt";
const char biconvessa[256]        = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 4-5/E4_biconvessa.csv";
const char biconcava[256]         = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 4-5/E4_biconcava.csv";


//Percorsi Federica
//const char rootfilepath[256]      = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper4.root";
//const char outfilepath[256]       = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper4_results.txt";
//const char biconvessa[256]        = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 4-5/E4_biconvessa.csv";
//const char biconcava[256]         = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 4-5/E4_biconcava.csv";

///Percorsi Filippo
const char rootfilepath[256]      = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper4.root";
const char outfilepath[256]       = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper4_results.txt";
const char biconvessa[256]            = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 4-5/E4_biconvessa.csv";
const char biconcava[256]             = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 4-5/E4_biconcava.csv";


void esper4()
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
    printout << "  | Experiment:  Ottica - esperienza 4                            | \n";
    printout << "  | Date:        13/12/3030                                       | \n";
    printout << "  | Revision:    2.1.0                                            | \n";
    printout << "  | Description: Lunghezze focali delle lenti                     | \n";
    printout << "  +---------------------------------------------------------------+ \n\n";


    //Apriamo il file CSV esistente in cui è stata salvata la tabella con i dati dello spettroscopio
    csvdata biconvessaCsv(biconvessa);
    biconvessaCsv.setDelimiters(';');
    csvdata biconcavaCsv(biconcava);
    biconcavaCsv.setDelimiters(';');

    //---------------------------- IMPORTAZIONE DATI SPERIMENTALI --------------------------------

    int N_p = 1;
    int N_q = 30;
    int N   = N_p * N_q;

    TH1D *  h_cv[N_p]; //Lente montaggio standard
    TH1D *  k_cv[N_p]; //Lente montaggio invertito
    double  p_cv[N_p];
    double  p_cv_err[N_p];

    TH1D *  h_co[N_p]; //Lente montaggio standard
    TH1D *  k_co[N_p]; //Lente montaggio invertito
    double  p_co[N_p];
    double  p_co_err[N_p];

    double zeroCondensatore = 281.e-3;
    double spessoreFinestra = 3.94e-3;

    double zeroCondensatore_err = 3.;
    double spessoreFinestra_err = 0.;

    double p_offset = 1.e-2;
    double q_offset = 0.;

    double ref0 = zeroCondensatore - spessoreFinestra;

    for (int i = 1; i < N+1; ++i)
    {
        const char * nameh = (std::string("cv_h") + std::to_string(i)).c_str();
        const char * namek = (std::string("cv_k") + std::to_string(i)).c_str();
        h_cv[i-1] = new TH1D(nameh, nameh, 50, 0., 1);
        k_cv[i-1] = new TH1D(namek, namek, 50, 0., 1);

        p_cv[i-1]     = abs(biconvessaCsv.getDouble(i, 0) - ref0) + p_offset;
        p_cv_err[i-1] = abs(biconvessaCsv.getDouble(i, 1));

        for (int j = 0; i < N_q; ++i)
        {
            h_cv[i-1]->Fill(abs(biconvessaCsv.getDouble(i+j, 2)) - ref0);
            k_cv[i-1]->Fill(abs(biconvessaCsv.getDouble(i+j, 3)) - ref0);
        }
    }

    for (int i = 0; i < N; ++i)
    {
        const char * nameh = (std::string("co_h") + std::to_string(i)).c_str();
        const char * namek = (std::string("co_k") + std::to_string(i)).c_str();
        h_co[i] = new TH1D(nameh, nameh, 50, 0., 1);
        k_co[i] = new TH1D(namek, namek, 50, 0., 1);

        p_co[i]     = abs(biconcavaCsv.getDouble(i, 0));
        p_co_err[i] = abs(biconcavaCsv.getDouble(i, 1));

        for (int j = 0; i < N_q; ++i)
        {
            h_co[i-1]->Fill(abs(biconcavaCsv.getDouble(i+j, 2)));
            k_co[i-1]->Fill(abs(biconcavaCsv.getDouble(i+j, 3)));
        }
    }

    h_cv[0]->Draw();
    return;

    //------------------------------- ANALISI DATI SPERIMENTALI ----------------------------------

    //Medie
    double cv_means_h[N_p];
    double co_means_h[N_p];
    double cv_means_k[N_p];
    double co_means_k[N_p];

    //Deviazioni standard DELLA MEDIA
    double cv_stdev_h[N_p];
    double co_stdev_h[N_p];
    double cv_stdev_k[N_p];
    double co_stdev_k[N_p];


    for (int i = 0; i < N_p; ++i)
    {
        cv_means_h[i] = h_cv[i]->GetMean();
        co_means_h[i] = h_co[i]->GetMean();
        cv_stdev_h[i] = h_cv[i]->GetMeanError();
        co_stdev_h[i] = h_co[i]->GetMeanError();
        cv_means_k[i] = k_cv[i]->GetMean();
        co_means_k[i] = k_co[i]->GetMean();
        cv_stdev_k[i] = k_cv[i]->GetMeanError();
        co_stdev_k[i] = k_co[i]->GetMeanError();
    }
    
    

    auto * ConvessaDir = new TGraphErrors(N_p);
    auto * ConvessaInv = new TGraphErrors(N_p);
    auto * BConcavaDir = new TGraphErrors(N_p);
    auto * BConcavaInv = new TGraphErrors(N_p);

    for (int i = 0; i < N_p; ++i)
    {
        ConvessaDir->SetPoint(i, 1./(p_cv[i]), 1./(co_means_h[i])); // 1/p vs 1/q
        ConvessaInv->SetPoint(i, 1./(p_cv[i]), 1./(co_means_k[i]));
        BConcavaDir->SetPoint(i, 1./(p_co[i]), 1./(cv_means_h[i]));
        BConcavaInv->SetPoint(i, 1./(p_co[i]), 1./(cv_means_k[i]));
/*
        ConvessaDir->SetPointError(i, ..., ... ); // 1/p vs 1/q
        ConvessaInv->SetPointError(i, ..., ... );
        BConcavaDir->SetPointError(i, ..., ... );
        BConcavaInv->SetPointError(i, ..., ... );
*/
    }

    /*
        FUNZIONE DI FIT:

        1/p + 1/q = 1/f
        1/q = 1/f - 1/p
        y   = 1/f - x
    */

    auto * f1 = new TF1("f1", "(1./[0])-x", 0., 1.);
    auto * f2 = new TF1("f1", "(1./[0])-x", 0., 1.);
    auto * f3 = new TF1("f1", "(1./[0])-x", 0., 1.);
    auto * f4 = new TF1("f1", "(1./[0])-x", 0., 1.);

    f1->SetParameter(0, 1.);
    f2->SetParameter(0, 1.);
    f3->SetParameter(0, 1.);
    f4->SetParameter(0, 1.);

    ConvessaDir->Fit("f1");
    ConvessaInv->Fit("f2");
    BConcavaDir->Fit("f3");
    BConcavaInv->Fit("f4");

    double f_cv_dir = f1->GetParameter(0);
    double f_cv_inv = f2->GetParameter(0);
    double f_co_dir = f3->GetParameter(0);
    double f_co_inv = f4->GetParameter(0);

    double f_cv_dir_err = f1->GetParError(0);
    double f_cv_inv_err = f2->GetParError(0);
    double f_co_dir_err = f3->GetParError(0);
    double f_co_inv_err = f4->GetParError(0);

    ConvessaDir->Write("ConvessaDir");
    ConvessaInv->Write("ConvessaInv");
    BConcavaDir->Write("BConcavaDir");
    BConcavaInv->Write("BConcavaInv");

    f1->Write("fit_ConvessaDir");
    f2->Write("fit_ConvessaInv");
    f3->Write("fit_BConcavaDir");
    f4->Write("fit_BConcavaInv");

    printout << "\n\n Compatibilità tra misurazioni con lente BICONVESSA diritta e ruotata\n";
    st.quickZtwotailed(f_cv_dir, f_cv_inv, f_cv_dir_err, f_cv_inv_err, printout);
    printout << "\n\n Compatibilità tra misurazioni con lente BICONCAVA diritta e ruotata\n";
    st.quickZtwotailed(f_co_dir, f_co_inv, f_co_dir_err, f_co_inv_err, printout);
/*
    double f_cv      = 0.5*(f_cv_dir + f_cv_inv);
    double f_co_res  = 0.5*(f_co_dir + f_co_inv);
    double f_cv_err  = ...
    double f_cor_err = ...

    printout << "\n Lunghezza Focale lente convergente, biconvessa, f = " << f_cv << " +/- " << f_cv_err;

    //La lente divergente, biconcava, in realtà non è stata misurata da sola, ma come sistema di lenti
    //  insieme a quella biconvessa. Pertanto utilizziamo la composizione di lenti sottili per trovare
    //  la sua lunghezza focale.

    // 1/f_res = 1/f_conv + 1/f_conc
    // 1/f_conc = 1/f_res - 1/f_conv

    double f_conc       = 1./(1./f_co_res - 1./f_cv);
    double f_conc_err   = ...

    printout << "\n\n Lunghezza Focale lente divergente, biconcava, f = " << f_conc << " +/- " << f_conc_err;
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


}