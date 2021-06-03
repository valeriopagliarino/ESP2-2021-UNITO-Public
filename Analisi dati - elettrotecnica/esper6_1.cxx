#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>  
#include "statlib.cxx"
#include "csvdata.cxx"
#include "espToolkit.cxx"
#include <string> 
using namespace std;

statlib st;
double pi = TMath::Pi();

//Percorsi Valerio:
//const char rootfilepath[256]       = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper6.root";
//const char outfilepath[256]        = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper6_results.txt";
//const char transistorTable[256]    = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati elettrotecnica 5 e 6/transistorTable.csv";

//Percorsi Federica:
const char rootfilepath[256]      = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper6.root";
const char outfilepath[256]       = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper6_results.txt";
const char transistorTable[256]   = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati elettrotecnica 5 e 6/transistorTable.csv";


//Percorsi Filippo:
//const char rootfilepath[256]      =
//const char outfilepath[256]       =
//const char transistorTable[256]   =


void esper6_1()
{
    cerr << "\nStart esper6.cxx data analysis";
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
    printout << "  | Experiment:  Elettrotecnica - esperienza 6                    | \n";
    printout << "  | Date:        15/03/2021                                       | \n";
    printout << "  | Revision:    1.0.0                                            | \n";
    printout << "  | Description: Caratterizzazione di un transistor BJT NPN CE    | \n";
    printout << "  +---------------------------------------------------------------+ \n\n";

    //Apriamo il file CSV con dati del transistor
    csvdata transistor(transistorTable);
    transistor.setDelimiters(';');

    //---------------------------- IMPORTAZIONE DATI SPERIMENTALI --------------------------------
    
    //Immissione dati
    double ImaxPSU          = 0.1;    //A 
    double R                = 9930.0; //Ohm
    double Rc               = 1015.2; //Ohm
    double ImaxPSU_err      = 0.001;  //A 
    double R_err            = 23.0;   //Ohm
    double Rc_err           = 2.3;    //Ohm
    double Temp             = 23.1;   //C
    double Temp_err         = 1.046;  //C

    //Costante di Boltzmann
    double Kb = 1.380649e-23;    // J/K
    
    //Carica dell'elettrone
    double e  = 1.602176634e-19; // C

    //Conversione temperature in Kelvin
    Temp     = Temp     + 273.15;
    Temp_err = Temp_err + 273.15;

    const int num_Ib = 90;

    double Ib[num_Ib];
    double Ib_err[num_Ib];
    double Vbe[num_Ib];
    double Vbe_err[num_Ib];
    double Vce[num_Ib];
    double Vce_err[num_Ib];
    double Ic[num_Ib];
    double Ic_err[num_Ib];

    for (int i = 1; i <= num_Ib; ++i)
    {
        Ib[i-1]        = transistor.getDouble (i,  0) / 1000000.0;
        Ib_err[i-1]    = transistor.getDouble (i,  1) / 1000000.0;
        Vbe[i-1]       = transistor.getDouble (i,  2) / 1000.0;
        Vbe_err[i-1]   = transistor.getDouble (i,  3) / 1000.0;
        Vce[i-1]       = transistor.getDouble (i,  4) / 1000.0;
        Vce_err[i-1]   = transistor.getDouble (i,  5) / 1000.0;
        Ic[i-1]        = transistor.getDouble (i,  6) / 1000.0;
        Ic_err[i-1]    = transistor.getDouble (i,  7) / 1000.0;
    }


    //---------------------------------- ANALISI DATI --------------------------------------------
    
    //ANALISI CARATTERISTICA DI USCITA DEL TRANSISTOR---------------------------------------------
    const int s = 18;

    printout << "\n\nTABELLA DATI SPERIMENTALI\nVbe\t\tVce\n-----------------------------\n";
    printout << std::setprecision(6);

    TGraphErrors * t[5];
    for (int j = 0; j < 5; ++j)
    {
        t[j] = new TGraphErrors(s);
        t[j]->SetName((std::string("t") + std::to_string(j)).c_str());

        int k = 0;
        for (int i = j*s ; i <= (s*(j+1)-1); i++)
        {
            t[j]->SetPoint(k, Vce[i], Ic[i]);
            t[j]->SetPointError(k, Vce_err[i], Ic_err[i]);
            printout << i << "\t" << Vce[i] << "\t \t" << Ic[i] << std::endl;
            k++;
        }
        printout << std::endl << "-----------------------------" << std::endl;
        t[j]->Sort();
    }


    //INTERPOLAZIONE CARATTERISTICA DI USCITA ----------------------------------------------------

    //Tensione termica
    double Vt     = Kb * Temp / e;
    double Vt_err = (Kb / e) * Temp_err;

    //Parametri di fit:
    /*
    * [0] Is        Corrente di polarizzazione inversa
    * [1] η * Vt    Fattore di idealità * tensione termica
    */

   //I\cdot e^{\frac{V}{t}}\cdot\left(\frac{-x}{A}\right)+I\cdot\left(e^{\frac{-x}{nt}}-1\right)

    TF1 * f[5];
    for (int i = 0; i < 5; ++i) 
        f[i]= new TF1((std::string("f") + std::to_string(i)).c_str(), "-[0]*exp([2]/[1])", 0., 10.);

    TCanvas * cOut = new TCanvas("cOut","Caratteristica in uscita",10,3,600,400);
    cOut->cd();
    cOut->SetTitle("Caratteristica in uscita");
    printout << ("\n\nCARATTERISTICA DI USCITA\n");
    
    for (int j = 0; j < 5; ++j)
    {
        f[j]->SetParameter(1, 5e-2);
        f[j]->SetParameter(0, 1e-8);
        /*t[j]->Fit(f[j], "r");
        printout << "\n\n----- Fit Ib = " << Ib[2+s*j] <<" A -----\n";
        printout << "Chi^2  "   << f[j]->GetChisquare();
        printout << "\np-val  " << f[j]->GetProb();
        printout << "\nndof   " << f[j]->GetNDF();
        printout << "\nrX^2   " << f[j]->GetChisquare() / f[j]->GetNDF();
        printout << "\nIs     " << f[j]->GetParameter(0) << " +/- " << f[j]->GetParError(0);
        printout << "\nη*Vt   " << f[j]->GetParameter(1) << " +/- " << f[j]->GetParError(1);
        */
        //t[j]->GetXaxis()->SetLimits(0.,9.);
        //t[j]->GetYaxis()->SetLimits(0.,3e-2);
        t[j]->Write((std::string("complete_t") + std::to_string(j)).c_str());  
    }

    t[3]->GetXaxis()->SetLimits(-0.5, 9.0);
    t[3]->Draw();
    t[3]->GetXaxis()->SetTitle("Vce [V]");
    t[3]->GetYaxis()->SetTitle("Ic [A]        ");
    t[3]->SetTitle("Caratteristica di uscita");
    for (int j = 0; j < 5; ++j) if (j != 3) t[j]->Draw("same");

    TPaveText * pt;
    pt = new TPaveText(5.99949,0.02420185,6.989971,0.0254392,"br");
    pt->AddText("400 #mu A");
    pt->Draw();
    
    pt = new TPaveText(5.99949,0.01914936,6.989971,0.0203867,"br");
    pt->AddText("350 #mu A");
    pt->Draw();
    
    pt = new TPaveText(5.99949,0.01430308,6.989971,0.01554043,"br");
    pt->AddText("300 #mu A");
    pt->Draw();
    
    pt = new TPaveText(5.977958,0.007703899,6.968438,0.008941246,"br");
    pt->AddText("200 #mu A");
    pt->Draw();
    
    pt = new TPaveText(5.977958,0.002857625,6.968438,0.004094971,"br");
    pt->AddText("100 #mu A");
    pt->Draw();

    cOut->SetGrid();
    cOut->Write("C_Out");
    cOut->Draw();

    //EFFETTO EARLY NELLA CARATTERISTICA DI USCITA -----------------------------------------------

    //Parametri di fit:
    /*
    * [0] Is        Corrente di polarizzazione inversa
    * [1] Vbe/ Vt   Tensione base-emettitore / Tensione termica
    * [2] Va        Tensione di Early
    */


    TCanvas * cEarly = new TCanvas("cEarly","Caratteristica di uscita",10,3,600,400);
    cEarly->cd();
    cEarly->SetTitle("Caratteristica in uscita, effetto Early");
     
    TF1 * a[5];
    TF1 * ad[5];
    double sVam = 0.; //Somma delle tensioni di Early pesate sull'inverso delle varianze
    double rVar = 0.; //Somma degli inversi delle varianze

    for (int j = 0; j < 5; ++j)
    {
        printout << "\n\n----- a[" << j << "] -----";
        const char * name = (std::string("a") + std::to_string(j)).c_str();
        if ((j != 4) && (j != 2)) a[j] = new TF1(name, "-[0]*exp([1]) * (1+x/[2])", 2.0, 7.5);
        if (j == 2)               a[j] = new TF1(name, "-[0]*exp([1]) * (1+x/[2])", 2.0, 6.8);
        if (j == 4)               a[j] = new TF1(name, "-[0]*exp([1]) * (1+x/[2])", 2.4, 6.8);
        ad[j] = new TF1(name, "-[0]*exp([1]) * (1+x/[2])", -30., 10.);
        a[j]->SetParameter(0, 1e-8);
        a[j]->SetParameter(1, (0.10 / 5e-2));
        a[j]->SetParameter(2, -20.);
        t[j]->Fit(a[j], "r");
        ad[j]->SetParameter(0, a[j]->GetParameter(0));
        ad[j]->SetParameter(1, a[j]->GetParameter(1));
        ad[j]->SetParameter(2, a[j]->GetParameter(2));
        ad[j]->SetLineStyle(6);
        printout << "\nIs     " << a[j]->GetParameter(0) << " +/- " << a[j]->GetParError(0);
        printout << "\nVbe/Vt " << a[j]->GetParameter(1) << " +/- " << a[j]->GetParError(1);
        printout << "\nVa     " << a[j]->GetParameter(2) << " +/- " << a[j]->GetParError(2);
        //Media pesata sugli errori delle tensioni di Early Va.
        sVam = sVam + a[j]->GetParameter(2) / (a[j]->GetParError(2) * a[j]->GetParError(2));
        rVar = rVar + 1./(a[j]->GetParError(2) * a[j]->GetParError(2));

        t[j]->Write((std::string("Early_t") + std::to_string(j)).c_str());
    }

    double Va_err = sqrt(1./rVar);       //Tensione di Early ottenuta come media pesata di parametri di fit
    double Va_val = sVam / rVar * (-1.); //Relativa deviazione standard

    printout << "\n\n Tensione di Early = " << Va_val << " +/- " << Va_err << " V\n\n";

    t[3]->GetXaxis()->SetLimits(-34.,9.);
    t[3]->SetTitle("Caratteristica di uscita, effetto Early");
    //t[3]->GetYaxis()->SetLimits(0.0, 0.02);
    t[3]->Draw();
    ad[3]->Draw("same");
    
    auto * l = new TLine(Va_val,-0.001,Va_val,0.005);
    l->SetLineStyle(2);
    l->SetLineWidth(2);
    l->SetLineColor(kBlue);
    l->Draw();
    pt = new TPaveText(1.694484,0.02241189,4.2777,0.02340765,"br");
    pt->AddText("400 #mu A");
    pt->Draw();
   
    pt = new TPaveText(1.717143,0.01758245,4.30036,0.01857821,"br");
    pt->AddText("350 #mu A");
    pt->Draw();
   
    pt = new TPaveText(1.671824,0.01325088,4.25504,0.01424665,"br");
    pt->AddText("300 #mu A");
    pt->Draw();

    pt = new TPaveText(-30.0187,0.00517765,-28.45932,0.006518109,"br");
    pt->AddText("Va");
    pt->Draw();
   
    pt = new TPaveText(1.694484,0.006828223,4.2777,0.007823984,"br");
    pt->AddText("200 #mu A");
    pt->Draw();
   
    pt = new TPaveText(1.671824,0.002546449,4.25504,0.003542211,"br");
    pt->AddText("100 #mu A");
    pt->Draw();

    for (int j = 0; j < 5; ++j) if (j != 3)
    {
        t[j]->Draw("same");
        ad[j]->Draw("same");
    }

    cEarly->SetGrid();
    cEarly->Write("C_Early");
    cEarly->Draw();
    return;

    //CARATTERISTICA DI INGRESSO -----------------------------------------------------------------

    TGraphErrors * fb[s];

    TCanvas * cIn = new TCanvas("cEarly","Caratteristica di ingresso",10,3,600,400);
    cIn->cd();
    cIn->SetTitle("Caratteristica in ingresso");

    for (int j = 0; j < s; ++j)
    {
        fb[j] = new TGraphErrors(5);
        fb[j]->SetName((std::string("fb") + std::to_string(j)).c_str());
        fb[j]->SetTitle((std::string("Vce ~ ") + std::to_string(Vce[j+s*2]) + std::string(" V")).c_str());

        for (int i = 0; i < 5; ++i)
        {
            fb[j]->SetPoint(i, Vbe[j+i*s], Ib[j+i*s]);
            fb[j]->SetPointError(i, Vbe_err[j+i*s], Ib_err[j+i*s]);
        }
        
        //fb[j]->SetMarkerStyle(20 + j);
        fb[j]->Sort();
        fb[j]->Write((std::string("fb") + std::to_string(j)).c_str());
    }

    auto * l1 = new TLegend(0.1,0.7,0.48,0.9);
    l1->AddEntry(fb[0], fb[0]->GetTitle(), "l");
    l1->SetHeader("Legenda");

    //Attenzione che in ogni caso alcune linee non vengono visualizzate e bisogna cambiare colore a manina!!!
    fb[1]->SetLineColor(4);
    fb[2]->SetLineColor(2);
    fb[3]->SetLineColor(3);
    fb[4]->SetLineColor(7);
    fb[5]->SetLineColor(94);
    fb[6]->SetLineColor(65);
    fb[7]->SetLineColor(102);
    fb[8]->SetLineColor(156);
    fb[9]->SetLineColor(51);
    fb[10]->SetLineColor(9);
    fb[11]->SetLineColor(108);
    fb[12]->SetLineColor(159);
    fb[13]->SetLineColor(207);
    fb[14]->SetLineColor(223);
    fb[15]->SetLineColor(214);
    fb[16]->SetLineColor(56);

    //fb[0]->SetLineStyle(2);
    fb[0]->SetLineWidth(1);
    fb[0]->GetXaxis()->SetTitle("Vbe [V]");
    fb[0]->GetYaxis()->SetTitle("Ib [A]");
    fb[0]->SetTitle("Caratteristica in ingresso Ib = f(Vbe)");
    fb[0]->Draw();

    for (int j = 0; j < s; ++j)
    {
        fb[j]->GetXaxis()->SetLimits(0.4 , 1.8);
        fb[j]->GetYaxis()->SetLimits(3e-5, 9e-4);
        if (j != 0) 
        {
            l1->AddEntry(fb[j], fb[j]->GetTitle(), "l");
            //fb[j]->SetLineStyle(2);
            fb[j]->SetLineWidth(1);
            //fb[j]->SetMarkerSize(3);
            fb[j]->Draw("same");
        }
        l1->Draw();
    }


    cIn->SetGrid();
    cIn->Draw();
    cIn->Write("cIn");


    //BETA_F(I_COLLETTORE) -----------------------------------------------------------------------

    double BetaF     [5];
    double BetaF_err [5];
    double Ib_s      [5];
    double Ib_s_err  [5];
    double Ic_s      [5];
    double Ic_s_err  [5];

    //Ib_s e Ic_s contengono i valori al vaiare di Vbe per Vce corrispondente a 6V
    const int u = 15;
    const int z = 5;

    for (int j = 0; j < z; ++j)
    {
        Ib_s[j]      = Ib[s*j + u];
        Ib_s_err[j]  = Ib_err[s*j + u];
        Ic_s[j]      = Ic[s*j + u];
        Ic_s_err[j]  = Ic_err[s*j + u];
        BetaF[j]     = 1. / (Ib_s[j] / Ic_s[j]);
        BetaF_err[j] =  sqrt(
            (1./Ib_s[j]) * (1./Ib_s[j]) * Ic_s_err[j] * Ic_s_err[j] + 
            Ic_s[j] * Ic_s[j] * (1./(Ib_s[j] * Ib_s[j])) * (1./(Ib_s[j] * Ib_s[j])) * Ib_s_err[j] * Ib_s_err[j]
        );
    }

    auto * bf = new TGraphErrors(z, Ic_s, BetaF, Ic_s_err, BetaF_err);
    //auto * bf_fit = new TF1();
    bf->Write("BetaF(Ic)");

    TCanvas * beta = new TCanvas("beta","Bf(Ic)",10,3,600,400);
    beta->cd();
    beta->SetTitle("Bf(Ic)");
    bf->Sort();
    bf->GetXaxis()->SetTitle("Ic [A]");
    bf->GetYaxis()->SetTitle("#beta f");
    bf->SetTitle("Amplificazione al variare di Ic #beta f(Ic)");
    //bf->SetMarkerStyle(21);
    //bf->SetMarkerSize(1);
    bf->Draw();
    beta->SetGrid();
    beta->Draw();

    
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

    std::ifstream ff(outfilepath);
    if (ff.is_open())
        std::cout << ff.rdbuf();

}