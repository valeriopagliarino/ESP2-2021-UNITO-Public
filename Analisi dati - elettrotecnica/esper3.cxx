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
//const char rootfilepath[256]   = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper3.root";
//const char outfilepath[256]    = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper3_results.txt";
//const char filterdpath[256]    = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati elettrotecnica 3 e 4/RCfilter.csv";

//Percorsi Federica:
//const char rootfilepath[256] = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper3.root";
//const char outfilepath[256]  = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper3_results.txt";
//const char filterdpath[256]  = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati elettrotecnica 3 e 4/RCfilter.csv";

//Percorsi Filippo:
const char rootfilepath[256] = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper3.root";
const char outfilepath[256]  = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper3_results.txt";
const char filterdpath[256]  = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati elettrotecnica 3 e 4/RCfilter.csv";

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
    

    //Stampiamo l'intestazione 
    printout << endl;
    printout << "  +---------------------------------------------------------------+ \n";
    printout << "  |  Oreglia, Sibilla, Pagliarino - Corso B - C.d.L. in Fisica    | \n";
    printout << "  |           Università degli Studi di Torino -  ESP2            | \n";
    printout << "  +---------------------------------------------------------------+ \n";
    printout << "  |             ANALISI DATI  - C++11 + CERN ROOT 6               | \n";
    printout << "  +---------------------------------------------------------------+ \n";
    printout << "  | Experiment:  Elettrotecnica - esperienza 3                    | \n";
    printout << "  | Date:        07/01/2021                                       | \n";
    printout << "  | Revision:    1.0.0                                            | \n";
    printout << "  | Description: Studio di filtri passa alto e passa basso        | \n";
    printout << "  +---------------------------------------------------------------+ \n\n";

    //Apriamo il file CSV con dati filtro passa basso
    csvdata filtercsv(filterdpath);
    filtercsv.setDelimiters(';');

    //---------------------------- IMPORTAZIONE DATI SPERIMENTALI --------------------------------

    //Valori dei componenti misurati che verranno poi confrontati con i valori ricavati come parametri dei fit
    double R            =    1058.9;        //Ω
    double R_err        =       2.4;        //Ω
    double C            =  46.78e-9;        //F
    double C_err        =   0.12e-9;        //F

    double Vgen         =        20;        //V
    double Vgen_err     =         1;        //V

    //Misure alla frequenza di taglio
    double Vin_fc_hi_pass      = 20.06;    //V
    double Vin_fc_hi_pass_err  = 1.;       //V
    double Vout_fc_hi_pass     = 14;       //V
    double Vout_fc_hi_pass_err = 1.;       //V
    double Delay_hi_pass       = 39e-6;    //s
    double Delay_hi_pass_err   = 39e-6;    //s

    double Vin_fc_lo_pass      = 19.8;     //V
    double Vin_fc_lo_pass_err  = 1.;       //V
    double Vout_fc_lo_pass     = 14.1;     //V
    double Vout_fc_lo_pass_err = 1.;       //V
    double Delay_lo_pass       = 41.e-6;   //s
    double Delay_lo_pass_err   = 1e-6;     //s

    double Fc                  = 3212.95;  //Hz
    double Fc_err              = 0.64;     //Hz

    //Integratore e derivatore
    double Diff_f             = 20.;       //Hz
    double Diff_f_err         = 0.004;     //Hz
    double Diff_Vout          = 2.98;      //V
    double Diff_Vout_err      = 0.04;      //V
    double Diff_Vin           = 18.9;      //V
    double Diff_Vin_err       = 1.;        //V
    double Diff_Delay         = 12.8e-3;   //s
    double Diff_Delay_err     = 1e-4;      //s

    double Intg_f             = 150.e3;    //Hz
    double Intg_f_err         = 0.03e3;    //Hz
    double Intg_Vout          = 4.36;      //V
    double Intg_Vout_err      = 0.2;       //V
    double Intg_Vin           = 19.7;      //V
    double Intg_Vin_err       = 1.;        //V
    double Intg_Delay         = 1.5e-6;    //s
    double Intg_Delay_err     = 0.1e-6;    //s

    

    //Risposta in frequenza filtro
    int rowNum = 45;
    double Vin_lo_pass_val[rowNum];  //V
    double Vin_lo_pass_err[rowNum];  //V
    double Vin_hi_pass_val[rowNum];  //V
    double Vin_hi_pass_err[rowNum];  //V
    double frequency_val  [rowNum];  //Hz
    double frequency_err  [rowNum];  //Hz
    double Vout_lo_pass_v [rowNum];  //V
    double Vout_lo_pass_e [rowNum];  //V
    double Vout_hi_pass_v [rowNum];  //V
    double Vout_hi_pass_e [rowNum];  //V

    for (int i = 2; i <= rowNum; ++i)
    {        
        Vin_lo_pass_val[i-2] = filtercsv.getDouble (i,  0);
        Vin_lo_pass_err[i-2] = filtercsv.getDouble (i,  3);
        Vin_hi_pass_val[i-2] = filtercsv.getDouble (i,  8);
        Vin_hi_pass_err[i-2] = filtercsv.getDouble (i, 11);
        frequency_val  [i-2] = filtercsv.getDouble (i,  1);
        frequency_err  [i-2] = filtercsv.getDouble (i,  4);
        Vout_lo_pass_v [i-2] = filtercsv.getDouble (i,  2);
        Vout_lo_pass_e [i-2] = filtercsv.getDouble (i,  5);
        Vout_hi_pass_v [i-2] = filtercsv.getDouble (i, 10);
        Vout_hi_pass_e [i-2] = filtercsv.getDouble (i, 13);

    }


    printout << "\n\nParsing completed.";

    //---------------------------------- ANALISI DATI --------------------------------------------
    
    //Calcolo grandezze derivate
    double HiPass_Gain        [rowNum]; //Adimensionale  
    double HiPass_Gain_err    [rowNum]; //Adimensionale  
    double HiPass_Gain_dB     [rowNum]; //dB
    double HiPass_Gain_dB_err [rowNum]; //dB
    double LoPass_Gain        [rowNum]; //Adimensionale  
    double LoPass_Gain_err    [rowNum]; //Adimensionale  
    double LoPass_Gain_dB     [rowNum]; //dB
    double LoPass_Gain_dB_err [rowNum]; //dB

    

    printout << "\n\n\nGUADAGNO FILTRI\n-----------------------\nFrequenza [Hz]\t\t\t\tGain_lo_pass\t\t\t\tGain_lo_pass_dB\t\t\t\tGain_hi_pass\t\t\t\tGain_hi_pass_dB\n";
    for (int i = 0; i < rowNum -1; ++i)
    {        
        HiPass_Gain[i] = Vout_hi_pass_v[i] / Vin_hi_pass_val[i];
        LoPass_Gain[i] = Vout_lo_pass_v[i] / Vin_lo_pass_val[i];

        HiPass_Gain_err[i] = sqrt(
            (1./(Vin_hi_pass_val[i])) * (1./(Vin_hi_pass_val[i])) * Vout_hi_pass_e[i] * Vout_hi_pass_e[i] +
            ((Vout_hi_pass_v[i])/((Vin_hi_pass_val[i])*(Vin_hi_pass_val[i]))) * ((Vout_hi_pass_v[i])/((Vin_hi_pass_val[i])*(Vin_hi_pass_val[i]))) * Vin_hi_pass_err[i]*Vin_hi_pass_err[i]
        );

        LoPass_Gain_err[i] = sqrt(
            (1./(Vin_lo_pass_val[i])) * (1./(Vin_lo_pass_val[i])) * Vout_lo_pass_e[i] * Vout_lo_pass_e[i] +
            ((Vout_lo_pass_v[i])/((Vin_lo_pass_val[i])*(Vin_lo_pass_val[i]))) * ((Vout_lo_pass_v[i])/((Vin_lo_pass_val[i])*(Vin_lo_pass_val[i]))) * Vin_lo_pass_err[i]*Vin_lo_pass_err[i]
        );

        HiPass_Gain_dB[i] = getDB(Vin_hi_pass_val[i], Vout_hi_pass_v[i]);
        LoPass_Gain_dB[i] = getDB(Vin_lo_pass_val[i], Vout_lo_pass_v[i]);

        HiPass_Gain_dB_err[i] = getDBerr(Vin_hi_pass_val[i], Vout_hi_pass_v[i], Vin_hi_pass_err[i], Vout_hi_pass_e[i]);
        LoPass_Gain_dB_err[i] = getDBerr(Vin_lo_pass_val[i], Vout_lo_pass_v[i], Vin_lo_pass_err[i], Vout_lo_pass_e[i]);
        
        printout <<  std::setprecision(2) << "\n" << frequency_val[i] << "\t +/- \t" << frequency_err[i] << "\t\t\t\t" << LoPass_Gain[i] << "\t +/- \t" << LoPass_Gain_err[i] << "\t\t\t\t";
        printout <<  std::setprecision(2) <<  LoPass_Gain_dB[i] << "\t +/- \t" << LoPass_Gain_dB_err[i]   << "\t\t\t\t" << HiPass_Gain[i] << "\t +/- \t" << HiPass_Gain_err[i] << "\t\t\t\t";
        printout <<  std::setprecision(2) <<  HiPass_Gain_dB[i] << "\t +/- \t" << HiPass_Gain_dB_err[i];
    }

    TGraphErrors * lo_pass_gain = new TGraphErrors(44, frequency_val, LoPass_Gain, frequency_err, LoPass_Gain_err);
    TGraphErrors * hi_pass_gain = new TGraphErrors(44, frequency_val, HiPass_Gain, frequency_err, HiPass_Gain_err);
    TGraphErrors * lo_pass_dB   = new TGraphErrors(44, frequency_val, LoPass_Gain_dB, frequency_err, LoPass_Gain_dB_err);
    TGraphErrors * hi_pass_dB   = new TGraphErrors(44, frequency_val, HiPass_Gain_dB, frequency_err, HiPass_Gain_dB_err);
    double vv[2] = {50.,500000};
    double ff[2] = {-3.,-3.};
    TGraph       * refLine3dB   = new TGraph(2, vv, ff);
      
    lo_pass_gain->Write("Passa basso dati");
    hi_pass_gain->Write("Passa alto dati");
    lo_pass_dB->Write("dB - Passa basso dati");
    hi_pass_dB->Write("dB - Passa alto dati");
    refLine3dB->Write("Ref a -3 dB");

    //Fit
    auto * hi_fit = new TF1("lo_fit", "x/(sqrt([0]*[0]+x*x))");
    auto * lo_fit = new TF1("lo_fit", "[0]/(sqrt([0]*[0]+x*x))");
    auto * hi_dB  = new TF1("lo_fit", "20*log(x/(sqrt([0]*[0]+x*x)))/log(10)");
    auto * lo_dB  = new TF1("lo_fit", "20*log([0]/(sqrt([0]*[0]+x*x)))/log(10)");
    hi_fit->SetParameter(0,12000);   
    lo_fit->SetParameter(0,12000);   
    hi_dB ->SetParameter(0,2000);   
    lo_dB ->SetParameter(0,2000); 


    lo_pass_gain->Fit(lo_fit);
    printout << "\n\n\n\nFit passa basso gain";
    printout << "\n---------------------------------------------"; 
    printout << "\nFit chi square = " << lo_fit->GetChisquare() << "  ndof = " << lo_fit->GetNDF() << "   chi/ndf = " << lo_fit->GetChisquare() / lo_fit->GetNDF(); 
    printout << "\np-value = " << lo_fit->GetProb();
    printout << "\nf di taglio = \t" << lo_fit->GetParameter(0) << "\t +/- " << lo_fit->GetParError(0) << "  Hz";
    double ft_lo_fit = lo_fit->GetParameter(0);
    double ft_lo_fit_err = lo_fit->GetParError(0);
    lo_fit->Write("Fit passa basso gain"); //scrivere il fit sul file


    hi_pass_gain->Fit(hi_fit);
    printout << "\n\n\n\nFit passa alto gain";
    printout << "\n---------------------------------------------"; 
    printout << "\nFit chi square = " << hi_fit->GetChisquare() << "  ndof = " << hi_fit->GetNDF() << "   chi/ndf = " << hi_fit->GetChisquare() / hi_fit->GetNDF(); 
    printout << "\np-value = " << hi_fit->GetProb();
    printout << "\nf di taglio = \t" << hi_fit->GetParameter(0) << "\t +/- " << hi_fit->GetParError(0) << "  Hz";
    double ft_hi_fit = hi_fit->GetParameter(0);
    double ft_hi_fit_err = lo_fit->GetParError(0);
    hi_fit->Write("Fit passa alto gain"); //scrivere il fit sul file


    lo_pass_dB->Fit(lo_dB);
    printout << "\n\n\n\nFit passa basso dB";
    printout << "\n---------------------------------------------"; 
    printout << "\nFit chi square = " << lo_dB->GetChisquare() << "  ndof = " << lo_dB->GetNDF() << "   chi/ndf = " << lo_dB->GetChisquare() / lo_dB->GetNDF(); 
    printout << "\np-value = " << lo_dB->GetProb();
    printout << "\nf di taglio = \t" << lo_dB->GetParameter(0) << "\t +/- " << lo_dB->GetParError(0) << "  Hz";
    lo_dB->Write("Fit passa basso dB"); //scrivere il fit sul file


    hi_pass_dB ->Fit(hi_dB);
    printout << "\n\n\n\nFit passa alto dB";
    printout << "\n---------------------------------------------"; 
    printout << "\nFit chi square = " << hi_dB->GetChisquare() << "  ndof = " << hi_dB->GetNDF() << "   chi/ndf = " << hi_dB->GetChisquare() / hi_dB->GetNDF(); 
    printout << "\np-value = " << hi_dB->GetProb();
    printout << "\nf di taglio = \t" << hi_dB->GetParameter(0) << "\t +/- " << hi_dB->GetParError(0) << "  Hz";
    hi_dB->Write("Fit passa alto dB"); //scrivere il fit sul file

    
    
    //_________________________________Grafici__________________________________________________________________________________

    TCanvas * passabasso = new TCanvas("passabasso","Risposta in frequenza filtro passa basso",10,3,600,400);
    passabasso->SetFillColor(0);
    passabasso->cd();
    passabasso->Range(1.255128,-0.1396286,6.238813,1.316306);
    passabasso->SetFillColor(0);
    passabasso->SetBorderMode(0);
    passabasso->SetBorderSize(2);
    passabasso->SetLogx();
    passabasso->SetGridx();
    passabasso->SetGridy();
    passabasso->SetFrameBorderMode(0);
    passabasso->SetFrameBorderMode(0);
    lo_pass_gain->SetMarkerSize(0.6);
    lo_pass_gain->SetMarkerStyle(21);
    lo_pass_gain->SetTitle("Risposta in frequenza filtro passa basso");
    lo_pass_gain->GetXaxis()->SetTitle("Frequenza [Hz]");
    lo_pass_gain->GetYaxis()->SetTitle("Guadagno ");
    lo_pass_gain->Draw("AP");
    lo_fit->Draw("same");
    TLegend * l1 = new TLegend(0.1,0.7,0.48,0.9);
    l1->SetHeader("Legenda","C"); 				// option "C" allows to center the header
    l1->AddEntry(lo_pass_gain, "Dati sperimentali" ,"p");
    l1->AddEntry(lo_fit, "Fit calcolo simbolico","l");
    l1->Draw();
    

    TCanvas * passaalto = new TCanvas("passaalto","Risposta in frequenza filtro passa alto",10,3,600,400);
    passaalto->SetFillColor(0);
    passaalto->cd();
    passaalto->Range(1.255128,-0.1396286,6.238813,1.316306);
    passaalto->SetFillColor(0);
    passaalto->SetBorderMode(0);
    passaalto->SetBorderSize(2);
    passaalto->SetLogx();
    passaalto->SetGridx();
    passaalto->SetGridy();
    passaalto->SetFrameBorderMode(0);
    passaalto->SetFrameBorderMode(0);
    hi_pass_gain->SetMarkerSize(0.6);
    hi_pass_gain->SetMarkerStyle(21);
    hi_pass_gain->SetTitle("Risposta in frequenza filtro passa alto");
    hi_pass_gain->GetXaxis()->SetTitle("Frequenza [Hz]");
    hi_pass_gain->GetYaxis()->SetTitle("Guadagno ");
    hi_pass_gain->Draw("AP");
    hi_fit->Draw("same");
    TLegend * l2 = new TLegend(0.1,0.7,0.48,0.9);
    l2->SetHeader("Legenda","C"); 				// option "C" allows to center the header
    l2->AddEntry(lo_pass_gain, "Dati sperimentali" ,"p");
    l2->AddEntry(lo_fit, "Fit calcolo simbolico","l");
    l2->Draw();

    passaalto->Write("Grafico passa alto");
    passabasso->Write("Grafico passa basso");


    TCanvas * passabasso_dB = new TCanvas("passabasso_dB","Risposta in frequenza filtro passa basso",10,3,600,400);
    passabasso_dB->SetFillColor(0);
    passabasso_dB->cd();
    passabasso_dB->Range(1.255128,-0.1396286,6.238813,1.316306);
    passabasso_dB->SetFillColor(0);
    passabasso_dB->SetBorderMode(0);
    passabasso_dB->SetBorderSize(2);
    passabasso_dB->SetLogx();
    passabasso_dB->SetGridx();
    passabasso_dB->SetGridy();
    passabasso_dB->SetFrameBorderMode(0);
    passabasso_dB->SetFrameBorderMode(0);
    lo_pass_dB->SetMarkerSize(0.6);
    lo_pass_dB->SetMarkerStyle(21);
    lo_pass_dB->SetTitle("Risposta in frequenza filtro passa basso");
    lo_pass_dB->GetXaxis()->SetTitle("Frequenza [Hz]");
    lo_pass_dB->GetYaxis()->SetTitle("Guadagno [dB]");
    lo_pass_dB->Draw("AP");
    lo_fit->Draw("same");
    TLine * line1 = new TLine(0,-3.,550000,-3.);
    line1->SetLineColor(kRed);
    line1->Draw();
    TLine * line1v = new TLine(3200,-48,3200,-3);
    line1v->SetLineColor(kBlue);
    line1v->Draw();
    TLegend * l3 = new TLegend(0.1,0.7,0.48,0.9);
    l3->SetHeader("Legenda","C"); 				// option "C" allows to center the header
    l3->AddEntry(lo_pass_gain, "Dati sperimentali" ,"p");
    l3->AddEntry(lo_fit, "Fit calcolo simbolico","l");
    l3->AddEntry(line1, "-3 dB","l");
    l3->AddEntry(line1v, "Frequenza critica","l");
    l3->Draw();
    passabasso_dB->Draw();
    

    TCanvas * passaalto_dB = new TCanvas("passaalto_dB","Risposta in frequenza filtro passa alto",10,3,600,400);
    passaalto_dB->SetFillColor(0);
    passaalto_dB->cd();
    passaalto_dB->Range(1.255128,-0.1396286,6.238813,1.316306);
    passaalto_dB->SetFillColor(0);
    passaalto_dB->SetBorderMode(0);
    passaalto_dB->SetBorderSize(2);
    passaalto_dB->SetLogx();
    passaalto_dB->SetGridx();
    passaalto_dB->SetGridy();
    passaalto_dB->SetFrameBorderMode(0);
    passaalto_dB->SetFrameBorderMode(0);
    hi_pass_dB->SetMarkerSize(0.6);
    hi_pass_dB->SetMarkerStyle(21);
    hi_pass_dB->SetTitle("Risposta in frequenza filtro passa alto");
    hi_pass_dB->GetXaxis()->SetTitle("Frequenza [Hz]");
    hi_pass_dB->GetYaxis()->SetTitle("Guadagno [dB]");
    hi_pass_dB->Draw("AP");
    hi_fit->Draw("same");
    TLine * line2 = new TLine(0,-3.,550000,-3.);
    line2->SetLineColor(kRed);
    line2->Draw();
    TLine * line2v = new TLine(3100,-35.7,3100,-3);
    line2v->SetLineColor(kBlue);
    line2v->Draw();
    TLegend * l4 = new TLegend(0.1,0.7,0.48,0.9);
    l4->SetHeader("Legenda","C"); 				// option "C" allows to center the header
    l4->AddEntry(lo_pass_gain, "Dati sperimentali" ,"p");
    l4->AddEntry(lo_fit, "Fit calcolo simbolico","l");
    l4->AddEntry(line2, "-3 dB","l");
    l4->AddEntry(line2v, "Frequenza critica","l");
    l4->Draw();
    passaalto_dB->Draw();

    passaalto_dB->Write("dB - Grafico passa alto");
    passabasso_dB->Write("dB - Grafico passa basso");


    //refLine3dB->Draw("same");
    
    //Confronto frequenze di taglio
    double f_RC_id      = 1. / (2 * pi * (R*C)); //Non chiamatela teorica!!
    double f_RC_id_err  = sqrt(
        (1/(C*R*R)) * (1/(C*R*R)) * R_err * R_err +
        (1/(C*C*R)) * (1/(C*C*R)) * C_err * C_err
    ) / (2*pi);
    printout << "\n\nFrequenza di taglio con componenti ideali = " << f_RC_id << " +/- " << f_RC_id_err << "  Hz";

    statlib st;

    //Test Z
    printout << "\n\nTest Z compatibilità f_id e f-taglio passa basso";
    double Z_1 = (f_RC_id - ft_lo_fit) / sqrt(ft_lo_fit_err * ft_lo_fit_err + f_RC_id_err * f_RC_id_err);
    double pval1 = st.pvalZtwotailed(Z_1);
    printout << "\nZ = " << Z_1 << "\t\t p-val 2-tailed = " << pval1;

    printout << "\n\nTest Z compatibilità f_id e f-taglio passa alto";
    double Z_2 = (f_RC_id - ft_hi_fit) / sqrt(ft_hi_fit_err * ft_hi_fit_err + f_RC_id_err * f_RC_id_err);
    double pval2 = st.pvalZtwotailed(Z_2);
    printout << "\nZ = " << Z_2 << "\t\t p-val 2-tailed = " << pval2;

    printout << "\n\nTest Z compatibilità f-taglio passa basso con passa basso";
    double Z_3 = (ft_hi_fit - ft_lo_fit) / sqrt(ft_lo_fit_err * ft_lo_fit_err +  ft_hi_fit_err *ft_hi_fit_err );
    double pval3 = st.pvalZtwotailed(Z_3);
    printout << "\nZ = " << Z_3 << "\t\t p-val 2-tailed = " << pval3;


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

    Confrontare le frequenze di taglio del passa basso e del passa alto, dovrebbero essere compatibili tra loro.
    Non chiamare la frequenza ottenuta dai valori di R e di C "frequenza teorica": si tratta semplicemente della frequenza
    ottenuta approssimando R e C come componenti ideali.

    La parte dell'esperienza su integratori e derivatori è qualitativa e non è richiesto metterla in analisi dati.
    (Verrà richiesto nella prova pratica di realizzare integratori e derivatori)

    Esiste una correlazione tra f_high e f_low (frequenze di taglio dei due filtri), ma la possiamo trascurare e fare comunque
    un confronto mediante test Z.
*/