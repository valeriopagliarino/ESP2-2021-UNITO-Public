#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>  
#include "statlib.cxx"
#include "csvdata.cxx"
#include "espToolkit.cxx"
using namespace std;

statlib st;
double pi = TMath::Pi();

//Percorsi Valerio:
const char rootfilepath[256]  = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper5.root";
const char outfilepath[256]   = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper5_results.txt";
const char diodeTable[256]    = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati elettrotecnica 5 e 6/diodeTable.csv";
const char ledTable[256]      = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati elettrotecnica 5 e 6/ledTable.csv";

//Percorsi Federica:
//const char rootfilepath[256] = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper5.root";
//const char outfilepath[256]  = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper5_results.txt";
//const char diodeTable[256]   = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati elettrotecnica 5 e 6/diodeTable.csv";
//const char ledTable[256]     = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati elettrotecnica 5 e 6/ledTable.csv";


//Percorsi Filippo:
//const char rootfilepath[256] = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper5.root";
//const char outfilepath[256]  = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper5_results.txt";
//const char diodeTable[256]   = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati elettrotecnica 5 e 6/diodeTable.csv";
//const char ledTable[256]     = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati elettrotecnica 5 e 6/ledTable.csv";

void esper5()
{
    cerr << "\nStart esper5.cxx data analysis";
    //Avviamo il timer che misura il tempo di calcolo impiegato dalla CPU per completare l'analisi
    auto start = std::chrono::system_clock::now();

    //Apriamo un nuovo file di output in formato ROOT (.root) dove salvare i plots
    TFile out_file(rootfilepath, "RECREATE");
    
    //Apriamo un nuovo file di testo in cui salvare l'output dell'analisi
    fstream printout;
    printout.open(outfilepath, std::fstream::out);
    

    //Stampiamo l'intestazione 
    printout << endl;
    printout << "  +---------------------------------------------------------------+ \n";
    printout << "  |  Oreglia, Sibilla, Pagliarino - Corso B - C.d.L. in Fisica    | \n";
    printout << "  |           Università degli Studi di Torino -  ESP2            | \n";
    printout << "  +---------------------------------------------------------------+ \n";
    printout << "  |             ANALISI DATI  - C++11 + CERN ROOT 6               | \n";
    printout << "  +---------------------------------------------------------------+ \n";
    printout << "  | Experiment:  Elettrotecnica - esperienza 5                    | \n";
    printout << "  | Date:        15/03/2021                                       | \n";
    printout << "  | Revision:    1.0.0                                            | \n";
    printout << "  | Description: Caratterizzazione di un diodo al silicio         | \n";
    printout << "  +---------------------------------------------------------------+ \n\n";

    //Apriamo il file CSV con dati filtro passa basso
    csvdata diodeIV(diodeTable);
    diodeIV.setDelimiters(';');

    csvdata ledIV(ledTable);
    ledIV.setDelimiters(';');

    //---------------------------- IMPORTAZIONE DATI SPERIMENTALI --------------------------------
    
    // Diodo al Silicio
    double R1_Si            = 1015.2;   //Ohm
    double R1_Si_err        = 2.3;      //Ohm
    double VmaxPSU          = 20;       //V
    double VmaxPSU_err      = 0.1;      //V
    double ImaxPSU          = 0.1;      //A
    double ImaxPSU_err      = 0.001;    //A
    double Rint_Volt        = 9200000;  //Ohm
    double Rint_Volt_err    = 400000;   //Ohm
    double Rint_Amm         = 7.4;      //Ohm
    double Rint_Amm_err     = 0000.05;  //Ohm
    double Temp             = 23.1;     //Kelvin
    double Temp_err         = 1.046;    //Kelvin


    // Diodo LED
    double R1_LED           = 1015.2;   //Ohm
    double R1_LED_err       = 0002.3;   //Ohm
    double Temp_LED         = 23.1;     //Kelvin
    double Temp_LED_err     = 1.046;    //Kelvin

    const int num_Si        = 31;
    const int num_LED       = 30;

    double V_Si[num_Si];
    double V_Si_err[num_Si];
    double I_Si[num_Si];
    double I_Si_err[num_Si];

    double V_LED[num_LED];
    double V_LED_err[num_LED];
    double I_LED[num_LED];
    double I_LED_err[num_LED];

    for (int i = 1; i <= num_Si; ++i)
    {
        V_Si[i-1]       = diodeIV.getDouble (i,  1) / 1000.;
        V_Si_err[i-1]   = diodeIV.getDouble (i,  2) / 1000.;
        I_Si[i-1]       = diodeIV.getDouble (i,  3) / 1000.;
        I_Si_err[i-1]   = diodeIV.getDouble (i,  4) / 1000.;
    }

    for (int i = 1; i <= num_LED; ++i)
    {
        V_LED[i-1]      = ledIV.getDouble (i,  1) / 1000.;
        V_LED_err[i-1]  = ledIV.getDouble (i,  2) / 1000.;
        I_LED[i-1]      = ledIV.getDouble (i,  3) / 1000.;
        I_LED_err[i-1]  = ledIV.getDouble (i,  4) / 1000.;
    }

    //---------------------------------- ANALISI DATI --------------------------------------------
    
    auto * diode   = new TGraphErrors(num_Si,  V_Si,  I_Si,   V_Si_err,  I_Si_err);
    auto * diode2  = new TGraphErrors(num_Si,  V_Si,  I_Si,   V_Si_err,  I_Si_err);
    auto * l_inv   = new TGraphErrors(num_Si,  I_Si,  V_Si,   I_Si_err,  V_Si_err);
    auto * l_inv2   = new TGraphErrors(num_Si,  I_Si,  V_Si,   I_Si_err,  V_Si_err);
    auto * led     = new TGraphErrors(num_LED, V_LED, I_LED, V_LED_err, I_LED_err);

    diode->Sort();
    led->Sort();
    printout << setprecision(4);

    //Parametri di fit:
    /*
    * [0] Is        Corrente di polarizzazione inversa
    * [1] η * Vt    Fattore di idealità * tensione termica
    * [2] R_d       Resistenza della giunzione LED
    */

    //Fit diodo raddrizzatore al silicio
    auto * df = new TF1("df", "[0]*(exp(x/[1])-1)", 0.4637, 10.0);
    diode->SetTitle("I(V) Diodo al silicio");
    df->SetParameter(0, 1e-8);
    df->SetParameter(1, 5e-2);
    diode->Fit("df", "R");
    diode->SetTitle("I(V) Diodo al silicio");
    printout << "\n\nFit diodo raddrizzatore al silicio\n";
    printout << "Chi^2  "   << df->GetChisquare();
    printout << "\np-val  " << df->GetProb();
    printout << "\nndof   " << df->GetNDF();
    printout << "\nrX^2   " << df->GetChisquare() / df->GetNDF();
    printout << "\nIs     " << df->GetParameter(0) << " +/- " << df->GetParError(0) << " A";
    printout << "\nη*Vt   " << df->GetParameter(1) << " +/- " << df->GetParError(1) << " V";
    printout << "\n";
    //----------------------------FATTORE DI IDEALITA' ETA----------------------------------------
    double kb       = 8.617333262e-5;                 //eV/K                                   //adim
    double Vter     = (Temp+273.15)*kb;
    double Vter_err = kb*Temp_err;
    double etas             = df->GetParameter(1)/Vter;      //adim
    double etas_err = sqrt(
        (1/Vter)*(1/Vter)*df->GetParError(1)*df->GetParError(1)+
        (df->GetParameter(1)/(Vter*Vter))*(df->GetParameter(1)/(Vter*Vter))*Vter_err*Vter_err
    );
    printout << "\n\nη diodo silicio   " << etas << " +/- " << etas_err << " adim";
    printout << "\n";

    auto * c1 = new TCanvas("c3","I(V) Diodo al silicio",10,3,600,400);
    c1->cd();
    diode->GetYaxis()->SetTitle("Corrente [A]");
    diode->GetXaxis()->SetTitle("Tensione [V]");
    diode->Draw();
    c1->SetGrid();
    c1->Draw();
    c1->Write("C1");
    

    //Calcolo del chi^2 per i punti eliminati con la restrizione dell'intervallo di fit
    const int fnum = 10;
    
    TF1 * dfr[fnum];
    printout << "\n\nValori di chi^2 aggiungendo i valori esterni alla restrizione di fit\n";
    for (int i = 0; i < fnum; ++i)
    {
        double x1, x2, y1, y2, lowerBound = 0.; 
        diode2->GetPoint(i, x2, y2);
        if (i != 0)
        {
            diode2->GetPoint(i - 1, x1, y1);
            lowerBound = (x1 + x2) / 2.;  
        }
        if (i == 0) lowerBound = (x2 / 2.);
        dfr[i] = new TF1((std::string("dfr") + std::to_string(i)).c_str(), "[0]*(exp(x/[1])-1)", lowerBound, 10.0);
        dfr[i]->FixParameter(0, df->GetParameter(0));
        dfr[i]->FixParameter(1, df->GetParameter(1));
        diode2->Fit((std::string("dfr") + std::to_string(i)).c_str(), "R");
        printout << "\nVmin = " << lowerBound << "\t NDoF = " << dfr[i]->GetNDF() - 2<< "  chi^2 = " << dfr[i]->GetChisquare();
    }
  

    //Fit diodo LED
    auto * lf = new TF1("lf"  , "[0]*(exp(x/[1])-1)", 0, 10.0);
    lf->SetParameter(0, 1e-8);
    lf->SetParameter(1, 5e-2);
    led  ->Fit("lf", "r");
    led->SetTitle("I(V) Diodo LED");
    printout << "\n\nFit diodo LED\n";
    printout << "Chi^2  "   << lf->GetChisquare();
    printout << "\np-val  " << lf->GetProb();
    printout << "\nndof   " << lf->GetNDF();
    printout << "\nrX^2   " << lf->GetChisquare() / lf->GetNDF();
    printout << "\nIs     " << lf->GetParameter(0) << " +/- " << lf->GetParError(0) << " A";
    printout << "\nη*Vt   " << lf->GetParameter(1) << " +/- " << lf->GetParError(1) << " V";
    printout << "\n";

    //Fit diodo LED con funzione invertita e termine di resistenza interna giunzione
    auto * rf = new TF1("rf", "[1]*log(1+x/[0])+ [2]*x", 0.002, 0.02);
    rf->SetParameter(0, 1e-9);
    rf->SetParLimits(0, 1e-10, 1e-6);
    rf->FixParameter(0, 1e-8);
    rf->SetParameter(1, 5e-2);
    rf->SetParameter(2, 50.0);
    l_inv->SetTitle("V(I) Diodo LED");
    l_inv->Fit("rf", "R");
    printout << "\n\nFit diodo LED con f. inv. e resistenza giunzione\n";
    printout << "Chi^2  "   << rf->GetChisquare();
    printout << "\np-val  " << rf->GetProb();
    printout << "\nndof   " << rf->GetNDF();
    printout << "\nrX^2   " << rf->GetChisquare() / rf->GetNDF();
    printout << "\nIs     " << rf->GetParameter(0) << " +/- " << rf->GetParError(0) << " A";
    printout << "\nη*Vt   " << rf->GetParameter(1) << " +/- " << rf->GetParError(1) << " V";
    printout << "\nR_d    " << rf->GetParameter(2) << " +/- " << rf->GetParError(2) << " Ω";
    printout << "\n";
    double etal             = rf->GetParameter(1)/Vter;      //adim
    double etal_err = sqrt(
        (1/Vter)*(1/Vter)*rf->GetParError(1)*rf->GetParError(1)+
        (rf->GetParameter(1)/(Vter*Vter))*(rf->GetParameter(1)/(Vter*Vter))*Vter_err*Vter_err
    );
    printout << "\n\nη diodo led   " << etal << " +/- " << etal_err << " adim";
    printout << "\n";

    l_inv->Write("V(I) Diodo");
    diode->Write("Diodo");
    led->Write("Led");
    auto * c3 = new TCanvas("c3","V(I) Diodo LED fit invertito e completo con resistenza della giunzione",10,3,600,400);
    c3->cd();
    l_inv->GetXaxis()->SetTitle("Corrente [A]");
    l_inv->GetYaxis()->SetTitle("Tensione [V]");
    l_inv->Draw();
    rf->Draw("same");
    c3->SetGrid();
    c3->Draw();
    c3->Write("C3");

    l_inv->Write("V(I) Diodo LED comp.");
    diode->Write("I(V) Diodo");
    led->Write("I(V) Led");

    //Calcolo del chi^2 per i punti eliminati con la restrizione dell'intervallo di fit
    const int lnum = 10;
    
    TF1 * dlr[lnum];
    printout << "\n\nValori di chi^2 aggiungendo i valori esterni alla restrizione di fit\n";
    for (int i = 0; i < lnum; ++i)
    {
        double x1, x2, y1, y2, lowerBound = 0.; 
        l_inv2->GetPoint(i, x2, y2);
        if (i != 0)
        {
            l_inv2->GetPoint(i - 1, x1, y1);
            lowerBound = (x1 + x2) / 2.;  
        }
        if (i == 0) lowerBound = (x2 / 2.);
        dlr[i] = new TF1((std::string("dlr") + std::to_string(i)).c_str(), "[1]*log(1+x/[0]) + [2]*x", lowerBound, 0.02);
        dlr[i]->FixParameter(0, rf->GetParameter(0));
        dlr[i]->FixParameter(1, rf->GetParameter(1));
        dlr[i]->SetParameter(2, rf->GetParameter(2));
        l_inv2->Fit((std::string("dlr") + std::to_string(i)).c_str(), "R");
        printout << "\nImin = " << lowerBound << "\t NDoF = " << dlr[i]->GetNDF() - 2<< "  chi^2 = " << dlr[i]->GetChisquare();
    }




    //------------- SCRITTURA DEL TEMPO DI CALCOLO E CHIUSURA DEI FILES DI OUTPUT ----------------
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
        std::cout << f.rdbuf();

}