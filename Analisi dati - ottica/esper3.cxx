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
#include <cmath>
using namespace std;

statlib st;

//Percorsi Valerio
const char rootfilepath[256]     = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper3.root";
const char outfilepath[256]      = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper3_results.txt";
const char zeri[256]             = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/Azzeramento_polarimetro.csv";
const char saccarosio[256]       = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/Misure_ripetute_saccarosio.csv";
const char fruttosio[256]        = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/Misure_ripetute_fruttosio.csv";
const char malus[256]            = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/Malus.csv";

//Percorsi Federica
//const char rootfilepath[256]     = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper1.root";
//const char outfilepath[256]      = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper1_results.txt";
//const char zeri[256]             = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/Azzeramento_polarimetro.csv";
//const char saccarosio[256]       = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/Misure_ripetute_saccarosio.csv";
//const char fruttosio[256]        = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/Misure_ripetute_fruttosio.csv";
//const char malus[256]            = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/Malus.csv";

//Percorsi Filippo:
//const char rootfilepath[256]     = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper3.root";
//const char outfilepath[256]      = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper3_results.txt";
//const char zeri[256]             = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/Azzeramento_polarimetro.csv";
//const char saccarosio[256]       = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/Misure_ripetute_saccarosio.csv";
//const char fruttosio[256]        = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/Misure_ripetute_fruttosio.csv";
//const char malus[256]            = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/Malus.csv";

void esper3()
{
    cerr << "\nStart esper3.cxx data analysis";
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
    printout << "  | Experiment:  Ottica - esperienza 3                            | \n";
    printout << "  | Date:        13/12/3030                                       | \n";
    printout << "  | Revision:    2.1.0                                            | \n";
    printout << "  | Description: Polarimetro                                      | \n";
    printout << "  +---------------------------------------------------------------+ \n\n";


    //Apriamo il file CSV esistente in cui è stata salvata la tabella con i dati dello spettroscopio
    csvdata csvZeri(zeri);
    csvZeri.setDelimiters(';');
    csvdata csvSacc(saccarosio);
    csvSacc.setDelimiters(';');
    csvdata csvFrutt(fruttosio);
    csvFrutt.setDelimiters(';');
    csvdata csvMalus(malus);
    csvMalus.setDelimiters(';');
    //csv.setDelimiters(';');

    //---------------------------- IMPORTAZIONE DATI SPERIMENTALI --------------------------------

    double zero_min[12];          //sensibilità nonio 0,1 grado
    double mis_sacc_20[12];       //sensibilità nonio 0,1 grado
    double mis_sacc_err_20[12];   //sensibilità nonio 0,1 grado
    double mis_sacc_10[12];       //sensibilità nonio 0,1 grado
    double mis_sacc_err_10[12];   //sensibilità nonio 0,1 grado
    double mis_frut[12];          //sensibilità nonio 0,1 grado
    double mis_frut_err[12];      //sensibilità nonio 0,1 grado

  


    double prov_sacc_20 = 0.203;  //m
    double prov_sacc_10 = 0.105;  //m

    double pot_rot_sacc = 663.7; //deg*cm^3/(g*m)
    double pot_rot_frut = -920;  //deg*cm^3/(g*m)

    for (int i = 1; i < 13; i++)
    {
        zero_min[i-1]        =   csvZeri.getDouble(i, 0 );
        mis_sacc_20[i-1]     =   csvSacc.getDouble(i, 2 );
        mis_sacc_err_20[i-1] =   csvSacc.getDouble(i, 3 );
        mis_sacc_10[i-1]     =   csvSacc.getDouble(i, 6 );
        mis_sacc_err_10[i-1] =   csvSacc.getDouble(i, 7 );
        mis_frut[i-1]        =   csvFrutt.getDouble(i, 2 );
        mis_frut_err[i-1]    =   csvFrutt.getDouble(i, 3 );
        
    }

    double alpha[36];
    double alpha_err[36];
    double I[36];
    double I_err[36];

    for (int i = 2; i < 38; i++)
    {
        alpha[i-2]       = csvMalus.getDouble(i, 0 );
        alpha_err[i-2]   = csvMalus.getDouble(i, 1 );
        I[i-2]           = csvMalus.getDouble(i, 2 );
        I_err[i-2]       = csvMalus.getDouble(i, 3 );
    }


    //------------------------------- ANALISI DATI SPERIMENTALI ----------------------------------
    
    TGraphErrors * I_alpha = new TGraphErrors(36, alpha, I, alpha_err, I_err);

    auto * malus_fit = new TF1("malus_fit", "[0]*cos((x*pi/180)+[1])*cos((x*pi/180)+[1])+[2]");
    malus_fit->SetParameter(0,5);
    malus_fit->SetParameter(1,2);

    I_alpha->Fit(malus_fit);
    
    printout << "\n\n\n\nFit legge di Malus";
    printout << "\n---------------------------------------------"; 
    printout << "\nFit chi square = " << malus_fit->GetChisquare() << "  ndof = " << malus_fit->GetNDF() << "   chi/ndf = " << malus_fit->GetChisquare() / malus_fit->GetNDF(); 
    printout << "\np-value = " << malus_fit->GetProb();

    TCanvas * laser = new TCanvas("laser","Verifica della legge di Malus",10,3,600,400);
    laser->cd();
    I_alpha->SetTitle("Verifica della legge di Malus");
    I_alpha->GetXaxis()->SetTitle("angolo [deg]");
    laser->SetGrid();
    I_alpha->GetYaxis()->SetTitle("Corrente misurata [#mu A]");
    I_alpha->Draw();
    malus_fit->Draw("same");
    I_alpha->Write("Verifica legge di Malus");
    malus_fit->Write("Malus Fit");

    //------------------------------- concentrazioni delle soluzioni -----------------------------

    double offset     = 0.0108;
    double media_s_20 = 6.7;
    double media_s_10 = 3.604;
    double media_f_20 = -5.896;
    double s_err_20   = mis_sacc_err_20[0];
    double s_err_10   = mis_sacc_err_10[0];
    double f_err      = mis_frut_err[0];

    /*
    Ho modificato l'errore sull'angolo del fruttosio da 0.7 a 0.5 nel CSV
    poiché Filippo mi ha detto che è stato utilizzato l'errore di sensibilità, 
    che rimane uguale a 0.5 su tutte le misure per coerenza interna dell'esperimento.
    */

    printout << "\n\nmedia_s_20  " << media_s_20 << " +/- " << s_err_20;
    printout << "\nmedia_s_10 " << media_s_10 << " +/- " << s_err_10;
    printout << "\nmedia_f_20" << media_f_20 << " +/- " << f_err;

    double conc_sacc_20 = media_s_20 / (prov_sacc_20 * pot_rot_sacc);
    double conc_sacc_10 = media_s_10 / (prov_sacc_10 * pot_rot_sacc);
    double conc_frut    = media_f_20 / (prov_sacc_20 * pot_rot_frut);


    double conc_sacc_20_err = sqrt (1. / (prov_sacc_20 * pot_rot_sacc) * 1. / (prov_sacc_20 * pot_rot_sacc) * s_err_20 * s_err_20 
                                    + media_s_20 / (prov_sacc_20 * prov_sacc_20 * pot_rot_sacc) * media_s_20 / (prov_sacc_20 * prov_sacc_20 * pot_rot_sacc) * 0.002 * 0.002);

    double conc_sacc_10_err = sqrt (1. / (prov_sacc_10 * pot_rot_sacc) * 1. / (prov_sacc_10 * pot_rot_sacc) * s_err_10 * s_err_10 
                                    + media_s_10 / (prov_sacc_10 * prov_sacc_10 * pot_rot_sacc) * media_s_10 / (prov_sacc_10 * prov_sacc_10* pot_rot_sacc) * 0.002 * 0.002);

    double conc_frut_err    = sqrt (1. / (prov_sacc_20 * pot_rot_frut) * 1. / (prov_sacc_20 * pot_rot_frut) * s_err_20 * s_err_20 
                                    + media_s_20 / (prov_sacc_20 * prov_sacc_20 * pot_rot_frut) * media_s_20 / (prov_sacc_20 * prov_sacc_20 * pot_rot_frut) * 0.002 * 0.002);
    
    printout << "\n\n\n\n Concentrazione saccarosio 20 \t\t\t= "    << conc_sacc_20  << " +/- " << conc_sacc_20_err << "\t g/cm^3";
    printout << "\n\n Concentrazione saccarosio 10 \t\t\t= "        << conc_sacc_10  << " +/- " << conc_sacc_10_err << "\t g/cm^3";
    printout << "\n\n Concentrazione fruttosio  20 \t\t\t= "        << conc_frut     << " +/- " << conc_frut_err    << "\t g/cm^3";
    
    double z = (conc_sacc_20 - conc_sacc_10)/sqrt(conc_sacc_10_err*conc_sacc_10_err+conc_sacc_20_err*conc_sacc_20_err);
    printout << "\n\n Test z di Gauss = \t\t\t"          <<z;
    printout << "\n p-value = \t\t\t\t" << st.pvalZtwotailed(z);



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