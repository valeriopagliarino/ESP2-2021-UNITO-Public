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

//Percorsi Valerio:
//"/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper2.root";
//"/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper2_results.txt";

//Percorsi Federica:
const char rootfilepath[256] = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper2.root";
const char outfilepath[256]  = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper2_results.txt";

//Percorsi Filippo:
//
//

void esper2()
{
    cerr << "\nStart esper2.cxx data analysis";
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
    printout << "  | Experiment:  Elettrotecnica - esperienza 2                    | \n";
    printout << "  | Date:        13/12/2020                                       | \n";
    printout << "  | Revision:    1.1.0                                            | \n";
    printout << "  | Description: Caratterizzazione oscilloscopio e trasformatori  | \n";
    printout << "  +---------------------------------------------------------------+ \n\n";

    //---------------------------- IMPORTAZIONE DATI SPERIMENTALI --------------------------------

    //Trasformatori
    double freq            = 50;        //Hz
    double R1              = 100;       //Ω
    double R1_err          = 5;         //Ω
    double R2              = 485.3e3;   //Ω
    double R2_err          = 2.9e3;     //Ω
    double C               = 5.074e-6;  //F
    double C_err           = 0.013e-6;  //F
    double VinPP           = 325*2;       //V
    double VinPP_err       = 33;        //V
    double VabPP           = 40.4;      //V
    double VabPP_err       = 1.;        //V
    double VxPP            = 36.6;      //V
    double VxPP_err        = 1.;        //V
    double V1PP            = 5.56;      //V
    double V1PP_err        = 0.2;       //V
    double V2PP            = 11.1;      //V
    double V2PP_err        = 0.4;       //V
    double VyPP            = 15.5e-3;   //V
    double VyPP_err        = 0.4e-3;    //V

    //Resistenza interna oscilloscopio
    double Vcg             = V2PP;
    double Vcg_err         = V2PP_err;
    double Vdg             = 7.52;      //V
    double Vdg_err         = 0.20;      //V

    //Impedenza interna dichiarata dal costruttoore
    double R_i_teo         = 1.0e6;     //Ω
    double C_teo           = 20e-12;    //F

    //Prove condensatori integratore
    double C_a             = C;         //F
    double C_b             = 1.024e-6;  //F   
    double C_c             = 143.22e-9; //F
    double C_d             = 46.98e-9;  //F
    double C_a_err         = C_err;     //F
    double C_b_err         = 0.005e-6;  //F
    double C_c_err         = 0.32e-9;   //F
    double C_d_err         = 0.12e-9;   //F

    double Vy_a            = 13.5e-3;   //V
    double Vy_b            = 56.4e-3;   //V
    double Vy_c            = 378e-3;    //V
    double Vy_d            = 1.17;      //V
    double Vy_a_err        = 0.4e-3;    //V
    double Vy_b_err        = 2e-3;      //V
    double Vy_c_err        = 10e-3;     //V
    double Vy_d_err        = 0.04;      //V

    //---------------------------------- ANALISI DATI --------------------------------------------
    
    //Rapporti di trasformazione T1 e T2
    double RaT1 = VabPP / VinPP;
    double RaT1_err = sqrt(((1./VinPP) * (1./VinPP)) * VabPP_err * VabPP_err +  (-1.*VabPP / (VinPP*VinPP)) * (-1.*VabPP / (VinPP*VinPP)) * VinPP_err * VinPP_err);

    double RaT2 = V2PP / V1PP;
    double RaT2_err = sqrt(((1./V1PP) * (1./V1PP)) * V2PP_err * V2PP_err +  (-1.*V2PP / (V1PP*V1PP)) * (-1.*V2PP / (V1PP*V1PP)) * V1PP_err * V1PP_err);

    //Impedenza interna oscilloscopio
    double R_int_osc = (Vdg * R2) / (Vcg - Vdg);
    double R_int_osc_err = sqrt(
        (R2/(Vcg-Vdg)+Vdg*R2/((Vcg-Vdg)*(Vcg-Vdg)))  *  (R2/(Vcg-Vdg)+Vdg*R2/((Vcg-Vdg)*(Vcg-Vdg))) * Vdg_err * Vdg_err +
        (-1*Vdg*R2/((Vcg-Vdg)*(Vcg-Vdg)))  *  (-1*Vdg*R2/((Vcg-Vdg)*(Vcg-Vdg))) * Vcg_err * Vcg_err +
        (Vdg/(Vcg-Vdg))  *  (Vdg/(Vcg-Vdg)) * R2_err * R2_err
    );

    printout << "\n";
    printout << "\nRapporto di trasformazione T1 = " << RaT1 << "  +/- " << RaT1_err;
    printout << "\nRapporto di trasformazione T2 = " << RaT2 << "  +/- " << RaT2_err;
    printout << "\nImpedenza interna oscilloscopio Ri [Ohm] = " << R_int_osc <<  " +/- " << R_int_osc_err;

    //Confronto con test Z resistenza interna oscilloscopio con specifiche del costruttore
    double pi   = 3.1415926535;
    double chiC = (1./(C_teo * 2 * pi * freq));
    double Imp_ref = (chiC*R_i_teo)/sqrt(((R_i_teo)*(R_i_teo)  + (chiC)*(chiC)));

    double z_osc_ires_ref = (Imp_ref - R_int_osc) / R_int_osc_err;
    double pval_osc_ires_ref = st.pvalZtwotailed(z_osc_ires_ref);

    printout << "\n\nImpedenza interna dichiarata dal costruttore a 50Hz [Ohm] = " << Imp_ref;
    printout << "\nTest Z di compatibilità misurata <=> dichiarata, Z = " << z_osc_ires_ref;
    printout << "  pvalue (two-tailed) = " << pval_osc_ires_ref;

    //Calcolo tau e f_c
    double tau = C*R2;
    double tau_err = sqrt(
        R2*R2 * C_err * C_err + 
        C * C * R2_err * R2_err
    );
    printout << "\ntau [s] = " << tau << " +/- " << tau_err;

    double f_c      = 1. / tau;
    double f_c_err  = abs((-1. / (tau*tau)) * tau_err);
    printout << "\nf_c [Hz] = " << f_c << " +/- " << f_c_err;


    //Filtro passivo integratore con diversi condensatori
    TGraphErrors * integrator = new TGraphErrors(4);
    integrator->SetPoint(1, C_a, Vy_a);
    integrator->SetPoint(2, C_b, Vy_b);
    integrator->SetPoint(3, C_c, Vy_c);
    integrator->SetPoint(4, C_d, Vy_d);

    integrator->SetPointError(1, C_a_err, Vy_a_err);
    integrator->SetPointError(2, C_b_err, Vy_b_err);
    integrator->SetPointError(3, C_c_err, Vy_c_err);
    integrator->SetPointError(4, C_d_err, Vy_d_err);
    integrator->SetTitle("Tensione in uscita al variare della capacità dell'integratore RC");
    integrator->GetXaxis()->SetTitle("Capacita' C [F]");
    integrator->GetYaxis()->SetTitle("Tensione in uscita [V]");

    //Attenuazione in dB dell'integratore
    TGraphErrors * integratordB = new TGraphErrors(4);
    integratordB->SetPoint(1, C_a, getDB(V2PP, Vy_a));
    integratordB->SetPoint(2, C_b, getDB(V2PP, Vy_b));
    integratordB->SetPoint(3, C_c, getDB(V2PP, Vy_c));
    integratordB->SetPoint(4, C_d, getDB(V2PP, Vy_d));

    integratordB->SetPointError(1, C_a_err, getDBerr(V2PP, Vy_a, V2PP_err, Vy_a_err));
    integratordB->SetPointError(2, C_b_err, getDBerr(V2PP, Vy_b, V2PP_err, Vy_b_err));
    integratordB->SetPointError(3, C_c_err, getDBerr(V2PP, Vy_c, V2PP_err, Vy_c_err));
    integratordB->SetPointError(4, C_d_err, getDBerr(V2PP, Vy_d, V2PP_err, Vy_d_err));

    printout << "\n\nAttenuazione integratore: condensatore [F], tensione in uscita [V], guadagno in uscita [dB]";
    printout << "\nC_a\t" << C_a << " +/- " << C_a_err << "\t\t" << Vy_a << " +/- " << Vy_a_err << "\t"   << getDB(V2PP, Vy_a) << " +/- " << getDBerr(V2PP, Vy_a, V2PP_err, Vy_a_err);
    printout << "\nC_a\t" << C_b << " +/- " << C_b_err << "\t\t" << Vy_b << " +/- " << Vy_b_err << "\t"   << getDB(V2PP, Vy_b) << " +/- " << getDBerr(V2PP, Vy_b, V2PP_err, Vy_b_err);
    printout << "\nC_a\t" << C_c << " +/- " << C_c_err << "\t\t" << Vy_c << " +/- " << Vy_c_err << "\t\t" << getDB(V2PP, Vy_c) << " +/- " << getDBerr(V2PP, Vy_c, V2PP_err, Vy_c_err);
    printout << "\nC_a\t" << C_d << " +/- " << C_d_err << "\t\t" << Vy_d << " +/- " << Vy_d_err << "\t\t" << getDB(V2PP, Vy_d) << " +/- " << getDBerr(V2PP, Vy_d, V2PP_err, Vy_d_err);


    integratordB->SetTitle("Guadagno al variare della capacita' dell'integratore RC");
    integratordB->GetXaxis()->SetTitle("Capacità C [F]");
    integratordB->GetYaxis()->SetTitle("Guadagno [dB]");
    integrator->SetMarkerSize(2.2);
    integratordB->SetMarkerSize(2.2);
    integrator->SetDrawOption("A*");
    integratordB->SetDrawOption("A*");

    integrator->Write("Filtro integratore Vout(Cap)");
    integratordB->Write("Filtro integratore Gain(Cap)");

    gROOT->Reset();
    auto * c1 = new TCanvas("c1","gerrors2",200,10,700,500);
    TPad * pad = new TPad("pad","",0,0,1,1);
    //pad->SetFillColor(0);
    pad->SetGrid();
    pad->Draw();
    pad->cd();

    // draw a frame to define the range
    TH1F *hr = c1->DrawFrame(2e-8,0.,6e-6,1.4);
    hr->SetTitle("Attenuazione introdotta dall'integratore");
    hr->SetXTitle("Capacita' C [F]");
    hr->SetYTitle("Ampiezza in uscita picco-picco [V]");
    //pad->GetFrame()->SetFillColor(21);
    pad->GetFrame()->SetBorderSize(12);

    integrator->SetMarkerColor(kBlue);
    integrator->SetMarkerStyle(21);
    gPad->SetLogx();
    integrator->SetLineWidth(0);
    integrator->Draw("LP");


    c1->cd();
    TPad *overlay = new TPad("overlay","",0,0,1,1);
    gPad->SetLogx();
    overlay->SetFillColor(0);
    overlay->SetFrameFillStyle(4000);
    overlay->Draw();
    overlay->cd();
    integratordB->SetMarkerColor(kRed);
    integratordB->SetMarkerStyle(20);
    integratordB->SetLineWidth(0);
    Double_t xmin = pad->GetUxmin();
    Double_t ymin = -60;
    Double_t xmax = pad->GetUxmax();
    Double_t ymax = 0;
    TH1F *hframe = overlay->DrawFrame(xmin,ymin,xmax,ymax);
    hframe->GetYaxis()->SetTickLength(0);
    hframe->GetXaxis()->SetLabelOffset(99);
    hframe->GetYaxis()->SetLabelOffset(99);
    gPad->SetLogx();
    integrator->GetXaxis()->SetTitle("Capacita' C [F]");
    integrator->GetYaxis()->SetTitle("Guadagno in uscita [dB]");
    integratordB->Draw("LP");
    gPad->SetLogx();
    
    //Draw an axis on the right side
    TGaxis * axis = new TGaxis(xmax,ymin,xmax, ymax,ymin,ymax,510,"+L");
    axis->SetTitle("Guadagno in uscita [dB]");
    axis->SetLabelSize(0.035);
    axis->SetTextFont(2);
    axis->SetLineColor(kRed);
    axis->SetLabelColor(kRed);
    axis->Draw();

    auto * leg = new TLegend(0.7,0.1,0.9,0.2);
    leg->SetHeader("Legenda");
    leg->AddEntry(integrator,"Ampiezza in uscita");
    leg->AddEntry(integratordB,"Guadagno");
    leg->Draw();

    c1->Write("Canvas integrator");


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


/*  APPUNTI ANALISI DATI 12/01/21

    L'unica parte richiesta nella relazione sono le misure per caratterizzare i rapporti di trasformazione
    con opportune incertezze e le misure di resistenza interna con opportune incertezze.

    Fare attenzione ai fattori 2 e sqrt(2) e a non confrontare valori RMS con valori di picco

    Anche in questo caso la parte di misura del ciclo di isteresi non è richiesta nelle relazioni, ma sarà richiesta nella prova pratica.

*/