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
const char rootfilepath[256]      = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper2.root";
const char outfilepath[256]       = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper2_results.txt";
const char nVetro[256]            = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/INTERF_n_vetro.csv";
const char lLaser[256]            = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/INTERF_lambda_laser.csv";
const char nAria[256]             = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/INTERF_n_aria.csv";

//Percorsi Federica
//const char rootfilepath[256]      = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper2.root";
//const char outfilepath[256]       = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper2_results.txt";
//const char nVetro[256]            = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/INTERF_n_vetro.csv";
//const char lLaser[256]            = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/INTERF_lambda_laser.csv";
//const char nAria[256]             = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/INTERF_n_aria.csv;

//Percorsi Filippo
//const char rootfilepath[256]      = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper2.root";
//const char outfilepath[256]       = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper2_results.txt";
//const char nVetro[256]            = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/INTERF_n_vetro.csv";
//const char lLaser[256]            = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/INTERF_lambda_laser.csv";
//const char nAria[256]             = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 2-3/INTERF_n_aria.csv;

void esper2()
{
    double pi = TMath::Pi();

    cerr << "\nStart esper2.cxx data analysis";
    //Avviamo il timer che misura il tempo di calcolo impiegato dalla CPU per completare l'analisi
    auto start = std::chrono::system_clock::now();

    //Apriamo un nuovo file di output in formato ROOT (.root) dove salvare i plots
    TFile out_file(rootfilepath, "RECREATE");
    
    //Apriamo un nuovo file di testo in cui salvare l'output dell'analisi
    fstream printout;
    printout.open(outfilepath, std::fstream::out);
    printout << setprecision(8);
    

    //Stampiamo l'intestazione 
    printout << endl;
    printout << "  +---------------------------------------------------------------+ \n";
    printout << "  |  Oreglia, Sibilla, Pagliarino - Corso B - C.d.L. in Fisica    | \n";
    printout << "  |           Universita' degli Studi di Torino -  ESP2           | \n";
    printout << "  +---------------------------------------------------------------+ \n";
    printout << "  |             ANALISI DATI  - C++11 + CERN ROOT 6               | \n";
    printout << "  +---------------------------------------------------------------+ \n";
    printout << "  | Experiment:  Ottica - esperienza 2                            | \n";
    printout << "  | Date:        13/12/3030                                       | \n";
    printout << "  | Revision:    2.1.0                                            | \n";
    printout << "  | Description: Polarimetro                                      | \n";
    printout << "  +---------------------------------------------------------------+ \n\n";

    //Apertura files CSV
    csvdata CSVnVetro(nVetro);
    CSVnVetro.setDelimiters(';');
    csvdata CSVlLaser(lLaser);
    CSVlLaser.setDelimiters(';');
    csvdata CSVnAria(nAria);
    CSVnAria.setDelimiters(';');


    //---------------------------- IMPORTAZIONE DATI SPERIMENTALI --------------------------------


    //Misura lunghezza d'onda del laser
    double nFrHeNe          = 80.;      //     Frange di interferenza misurate
    double FrErr            = 0.1;      //     Incertezza sul conteggio frange (misura fatta con software di analisi immagini + webcam)
    double refLambda        = 632.8e-9;   //[m]  Lunghezza d'onda di riferimento laser dal manuale UniPhase 
    double refLambda_err    = 1e-10;     //[m]  Incertezza associata a ^^^
    double opticalPaths[10];            //[m]  Cammini ottici misurati
    double opticalPath_s    = 0.5e-6;   //[um] Sensibilità del micrometro

    for (int i = 1; i < 11; ++i)
    {
        opticalPaths[i-1] = CSVlLaser.getDouble(i, 0);
        //cout << "\n" << CSVlLaser.getDouble(i, 0);
    }
    
    //Misura indice di rifrazione del vetro
    double ppTichness       = 5.38 * 1e-3;   //[m] Spessore del parallel plate ruotabile
    double ppTichness_err   = 0.05 * 1e-3;   //[m] Incertenzza associata a ^^^
    double nFrGlass         = 50;            //    Numero frange di interferenza misurate
    double angles[10];                       //[deg] Angolo di rotazione del parallel plate
    double angles_s         = 0.1;           //[deg] Sensibilità del goniometro
    double zeros[10];                        //[deg] Zero del goniometro, misure ripetute
    for (int i = 1; i < 11; ++i)
    {
        angles[i-1] = CSVnVetro.getDouble(i, 0);
        zeros[i-1]  = CSVnVetro.getDouble(i, 1);
    }

    //Misura indice di rifrazione dell'aria
    double vacLen           = 32.33e-3;      //[m] Lunghezza della camera a vuoto
    double vacLen_err       = 4e-3;       //[m] Incertezza associata a ^^^
    double envPressure      = 987.7;         //[hPa] Pressione ambiente misurata da stazione meteoclimatica di Fisica
    double envPressure_err  = 0.1;           //[hPa] Incertezza associata a ^^^

 
    double nFrAir_s         = 0.3;          //Sensibilità associata al sistema di conteggio delle frange di interferenza
    double deltaP_s         = 100.;         //[mBar] Sensibilità del manometro
    double nFrAir[9];                       //Numero frange di interferenza misurate
    double deltaP[9];                       //[mBar] Variazione di pressione associata allo scorrimento delle frange di interferenza
    double nFrAir_err[9];                   //Errori associati a ^^^
    double deltaP_err[9];                   //[mBar] errori associati a ^^^

    for (int i = 1; i < 10; ++i)
    {
        nFrAir[i-1]         = CSVnAria.getDouble(i, 2); 
        deltaP[i-1]         = CSVnAria.getDouble(i, 0);
    }
    for (int i = 0; i  < 10; ++i)
    {
        nFrAir_err[i]     = nFrAir_s;
        deltaP_err[i]     = deltaP_s;
    }

    //------------------------------- ANALISI DATI SPERIMENTALI ----------------------------------

    auto * lLaser = new TH1D("nVetro", "Variazione di cammino ottico che causa lo spostamento di 80 frange", 6, 24.e-6, 26.5e-6);
    for (int i = 0; i < 10; ++i)
    {
        lLaser->Fill(opticalPaths[i]);
    }
    lLaser->Write("Laser lambda");
    printout << "\n\nLunghezza d'onda del laser:\n";
    printout << "\nNumero di misure = " << lLaser->GetEntries();
    printout << "\nMedia = " << lLaser->GetMean();
    printout << "\nDeviazione standard = " << lLaser->GetStdDev();
    printout << "\nErrore sulla media = " << lLaser->GetMeanError();
    
    //--------------------------MISURA LUNGHEZZA D'ONDA DI LASER He-Ne----------------------------

    double lambdaLaser      = 2 * lLaser->GetMean() / nFrHeNe; 
    double lambdaLaserErr   = sqrt(
        (2 / nFrHeNe) * (2 / nFrHeNe) * lLaser->GetMeanError() * lLaser->GetMeanError() +
        (1 / (nFrHeNe * nFrHeNe)) * 2 * lLaser->GetMean() * (1 / (nFrHeNe * nFrHeNe)) * 2 * lLaser->GetMean() * FrErr * FrErr
    );

    printout << "\n\nLunghezza d'onda misurata per il laser He-Ne = " << lambdaLaser << " +/- " << lambdaLaserErr << " m";
    printout << "\nValore di riferimento = " << refLambda << " +/- " << refLambda_err;

    //------------------------------------- MISURA n VETRO----------------------------------------

    auto * zero_h = new TH1D("zero_h", "Zero del goniometro", 10, -0.8, 0.5);
    for (int i = 0; i < 10; ++i)
    {
        zero_h->Fill(zeros[i]);
    }
    zero_h->Write("Zero_h");
    double zero_mean = zero_h->GetMean();
    double zero_err  = zero_h->GetMeanError();
    
    auto * angles_h = new TH1D("angles_h", "Angoli misurati", 10, 6.5, 7.8);
    for (int i = 0; i < 10; i++)
    {
        angles_h->Fill(angles[i] - zero_mean);
    }
    angles_h->Write("Angles_h");
    double angle_mean = angles_h->GetMean();
    double angleMerr  = angles_h->GetMeanError();
    double angle_err  = sqrt(zero_err * zero_err + angleMerr * angleMerr);
    printout << "\n\n-----------------------------------\n";
    printout << "\nMisura n del vetro:\n";
    printout << "\nZero     = " << zero_mean << " +/- " << zero_err << " deg";
    printout << "\nMedia θ  = " << angle_mean << " +/- " << angle_err << " deg";

    //Calcolo indice di rifrazione del vetro
    double glassN     = 1.58;
    double glassN_err = 0.03;

    printout << "\n\nIndice di rifrazione del vetro n = " << glassN << " +/- " <<  glassN_err;

    printout << "\n\n-----------------------------------\n\n";

    //------------------------------------ MISURA n(p) ARIA---------------------------------------

    printout << "Misura indice di rifrazione dell'aria";

    auto * airN = new TGraphErrors(10, deltaP, nFrAir, deltaP_err, nFrAir_err);
    airN->Sort();
    //Differenza di pressione tra la camera a vuoto e il laboratorio [mbar]
    airN->GetXaxis()->SetTitle("#Delta p [mbar]");
    //Spostamento frange di interferenza rispetto allo zero
    airN->GetYaxis()->SetTitle("#Delta N");
    
    //Fit Lineare
    auto * f1 = new TF1("f1", "pol1", 0., 100.);
    airN->Fit(f1);
    airN->SetTitle("Misurazione indice di rifrazione dell'aria");
    double slope        = f1->GetParameter(1);
    double slope_err    = f1->GetParError(1);
    double bias         = f1->GetParameter(0);
    double bias_err     = f1->GetParError(0);
    
    printout << "\n\nFit lineare per calcolo indice di rifrazione dell'aria  N = a + b*x\n";
    printout << "\nChi2  = " << f1->GetChisquare();
    printout << "\ndof   = " << f1->GetNDF();
    printout << "\np-val = " << f1->GetProb();
    printout << "\nCoff. = " << f1->GetProb();
    printout << "\na     = " << bias << " +/- " << bias_err << " frange";
    printout << "\nb     = " << slope << " +/- " << slope_err << " frange / mbar";
    airN->Write("Fit lineare N aria");

    TCanvas * C1 = new TCanvas("C1","Indice di rigrazione dell'aria", 600,400);
    C1->cd();
    C1->SetGrid();
    airN->Draw();
    C1->Draw();

    double dnAria     = slope * refLambda  / (2 * vacLen);
    double dnAria_err = sqrt(
        (refLambda  / (2 * vacLen)) * (refLambda  / (2 * vacLen)) * slope_err * slope_err +
        (slope  / (2 * vacLen)) * (slope  / (2 * vacLen)) * refLambda_err * refLambda_err +
        (slope * refLambda  / (vacLen * vacLen * 2)) * (slope * refLambda  / (vacLen * vacLen * 2)) * vacLen_err * vacLen_err
    );

    printout << "\n\nΔn/ΔP = " << dnAria << " +/- " << dnAria_err;

    double nAria        = 1 + dnAria * envPressure;
    double nAria_err    = sqrt(dnAria * dnAria * envPressure_err * envPressure_err +
        envPressure * envPressure * dnAria_err * dnAria_err);

    double refPressure  = 1013.25; //Please see @ https://emtoolbox.nist.gov/Wavelength/Documentation.asp
    double refN         = 1.000271800;

    double nAria2        = 1 + dnAria * refPressure;
    double nAria_err2    = refPressure * dnAria_err;

    printout << setprecision(8);
    printout << "\n\nn a 987.7 mbar = " << nAria << " +/- " << nAria_err;
    printout << "\n\nn a 101.325 mbar = " << nAria2 << " +/- " << nAria_err2;

    //--------------------------- Confronto con i valori di riferimento --------------------------

    printout << "\n\n-----------------------------------\n\nConfronto con i valori di riferimento:";
    printout << "\n\nLunghezza d'onda laser";
    double z1    = (refLambda - lambdaLaser) / sqrt(lambdaLaserErr * lambdaLaserErr + refLambda_err * refLambda_err);
    if (z1 > 0) z1 = z1 * -1.;
    double pval1 = (2. * (ROOT::Math::normal_cdf(z1)));
    printout << "\nZ test 2-tailed:  Z = " << z1 << "   pval = " << pval1;

    printout << "\n\nIndice di rifrazione dell'aria alla pressione di riferimento:";
    double z2    = (nAria2 - refN) / nAria_err2;
    if (z2 > 0) z2 = z2 * -1.;
    double pval2 = (2. * (ROOT::Math::normal_cdf(z2)));
    printout << "\nZ test 2-tailed:  Z = " << z2 << "   pval = " << pval2;
    printout << "\nN_ref = " << refN << "\n";
    

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