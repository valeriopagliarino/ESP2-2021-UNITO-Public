#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>  
#include <string>
#include "statlib.cxx"
#include "csvdata.cxx"
#include "espToolkit.cxx"
using namespace std;

statlib st;

//Percorsi Valerio:
const char rootfilepath[256]    = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati/esper4.root";
const char outfilepath[256]     = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati/esper4_results.txt";
const char rclfilepath[256]     = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati esperienze 3 e 4/RCLfilter.csv";
const char vrefilepath[256]     = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati esperienze 3 e 4/voltage_res.csv";
const char bandwdtpath[256]     = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati esperienze 3 e 4/bandwidth_vtable.csv";

//Percorsi Federica:
//const char rootfilepath[256]  = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati/esper4.root";
//const char outfilepath[256]   = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati/esper4_results.txt";
//const char rclfilepath[256]   = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati esperienze 3 e 4/RCLfilter.csv";
//const char vrefilepath[256]   = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati esperienze 3 e 4/voltage_res.csv";
//const char bandwdtpath[256]   = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati esperienze 3 e 4/bandwidth_vtable.csv";


//Percorsi Filippo:
//const char path[256]  = "";


void esper4_1()
{
    cerr << "\nStart esper4.cxx data analysis";
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
    printout << "  | Experiment:  Elettrotecnica - esperienza 4                    | \n";
    printout << "  | Date:        07/01/2021                                       | \n";
    printout << "  | Revision:    1.0.0                                            | \n";
    printout << "  | Description: Studio di RCL e rifasamento                      | \n";
    printout << "  +---------------------------------------------------------------+ \n\n";

    double pi = TMath::Pi();

    //Apriamo il file CSV con dati filtro RCL
    csvdata rclcsv(rclfilepath);
    rclcsv.setDelimiters(';');              //ROWS 2-44
    int rnum_rcl = 43;

    //Apriamo il file CSV con dati tensioni in uscita ai capi dell'induttanza vs R
    csvdata voltagerescsv(vrefilepath);
    voltagerescsv.setDelimiters(';');       //ROWS 2-8
    int rnum_vre = 7;

    //Apriamo il file CSV con dati ampiezza di banda vs R
    csvdata bandwdcsv(bandwdtpath);         
    bandwdcsv.setDelimiters(';');           //ROWS 2-25
    int rnum_bwd = 24;



    //---------------------------- IMPORTAZIONE DATI SPERIMENTALI --------------------------------

    double R        = 262.8     ; //Ohm
    double R_err    = 0.8       ; //Ohm
    double C        = 46.78e-9  ; //F
    double C_err    = 0.12e-9   ; //F  
    double L        = 4.19e-3   ; //H
    double L_err    = 0.04e-3   ; //mH
    double Q_L      = 12.70     ; //Adimensionale
    double Q_L_err  = 0.03      ; //Adimensionale


    //RCL Filter (scan in frequenza)
    double Vinp_RCL[rnum_rcl];
    double freq_RCL[rnum_rcl];
    double Vpas_RCL[rnum_rcl];
    double Vcut_RCL[rnum_rcl];
    double Vinp_RCL_err[rnum_rcl];
    double freq_RCL_err[rnum_rcl];
    double Vpas_RCL_err[rnum_rcl];
    double Vcut_RCL_err[rnum_rcl];

    for (int i = 2; i <= rnum_rcl; ++i)
    {
        Vinp_RCL[i-2]     = rclcsv.getDouble (i,  0);
        freq_RCL[i-2]     = rclcsv.getDouble (i,  1);
        Vpas_RCL[i-2]     = rclcsv.getDouble (i,  2);
        Vcut_RCL[i-2]     = rclcsv.getDouble (i,  7);
        Vinp_RCL_err[i-2] = rclcsv.getDouble (i,  3);
        freq_RCL_err[i-2] = rclcsv.getDouble (i,  4);
        Vpas_RCL_err[i-2] = rclcsv.getDouble (i,  5);
        Vcut_RCL_err[i-2] = rclcsv.getDouble (i,  8);
    }
    
    //Gain ai capi di L al variare di R
    double Rmea_Vre[rnum_vre];
    double Vinp_Vre[rnum_vre];
    double Vout_Vre[rnum_vre];
    double Gain_Vre[rnum_vre];
    double Rmea_Vre_err[rnum_vre];
    double Vinp_Vre_err[rnum_vre];
    double Vout_Vre_err[rnum_vre];
    double Gain_Vre_err[rnum_vre];



    for (int i = 2; i <= rnum_vre; ++i)
    {
        Rmea_Vre[i-2]     = voltagerescsv.getDouble (i,  0);
        Vinp_Vre[i-2]     = voltagerescsv.getDouble (i,  2);
        Vout_Vre[i-2]     = voltagerescsv.getDouble (i,  4);
        Rmea_Vre_err[i-2] = voltagerescsv.getDouble (i,  1);
        Vinp_Vre_err[i-2] = voltagerescsv.getDouble (i,  3);
        Vout_Vre_err[i-2] = voltagerescsv.getDouble (i,  5);
    }


    //Larghezza di banda al variare di R
    double Rmea_bwd[rnum_bwd];
    double Freq_bwd[rnum_bwd];
    double Vinp_bwd[rnum_bwd];
    double Vout_bwd[rnum_bwd];
    double Rmea_bwd_err[rnum_bwd];
    double Freq_bwd_err[rnum_bwd];
    double Vinp_bwd_err[rnum_bwd];
    double Vout_bwd_err[rnum_bwd];

    for (int i = 1; i <= rnum_bwd; ++i)
    {
        Rmea_bwd[i-1]     = bandwdcsv.getDouble (i,  2);
        Freq_bwd[i-1]     = bandwdcsv.getDouble (i,  0);
        Vinp_bwd[i-1]     = bandwdcsv.getDouble (i,  4);
        Vout_bwd[i-1]     = bandwdcsv.getDouble (i,  6);
        Rmea_bwd_err[i-1] = bandwdcsv.getDouble (i,  3);
        Freq_bwd_err[i-1] = bandwdcsv.getDouble (i,  1);
        Vinp_bwd_err[i-1] = bandwdcsv.getDouble (i,  5);
        Vout_bwd_err[i-1] = bandwdcsv.getDouble (i,  7);
    }


    //---------------------------------- ANALISI DATI --------------------------------------------

    const int rowNum = 42;

    double bandPass_Gain        [rowNum]; //Adimensionale  
    double bandPass_Gain_err    [rowNum]; //Adimensionale  
    double bandPass_Gain_dB     [rowNum]; //dB  
    double bandPass_Gain_dB_err [rowNum]; //dB 
    double bandCut_Gain         [rowNum]; //Adimensionale  
    double bandCut_Gain_err     [rowNum]; //Adimensionale  
    double bandCut_Gain_dB      [rowNum]; //dB  
    double bandCut_Gain_dB_err  [rowNum]; //dB

    TGraphErrors * band_cut_gain = new TGraphErrors(rowNum, freq_RCL, bandCut_Gain, freq_RCL_err, bandCut_Gain_err);
    TGraphErrors * band_pass_gain = new TGraphErrors(rowNum, freq_RCL, bandPass_Gain, freq_RCL_err, bandPass_Gain_err);
    TGraphErrors * band_cut_dB   = new TGraphErrors(rowNum, freq_RCL, bandCut_Gain_dB, freq_RCL_err, bandCut_Gain_dB_err);
    TGraphErrors * band_pass_dB   = new TGraphErrors(rowNum, freq_RCL, bandPass_Gain_dB, freq_RCL_err, bandPass_Gain_dB_err);

    printout << "\n\n\n------------------------------\nFrequenza \t\t\tVcut_gain\t\t\tVcut_dB\t\t\tVpass_gain\t\t\tVpass_dB\n";

    for (int i = 0; i < rowNum; ++i)
    {
        bandPass_Gain[i] = Vpas_RCL[i] / Vinp_RCL[i];
        bandCut_Gain[i] = Vcut_RCL[i] / Vinp_RCL[i];

        bandPass_Gain_err[i] = sqrt(
            (1./(Vinp_RCL[i])) * (1./(Vinp_RCL[i])) * Vpas_RCL_err[i] * Vpas_RCL_err[i] +
            ((Vpas_RCL[i])/((Vinp_RCL[i])*(Vinp_RCL[i]))) * ((Vpas_RCL[i])/((Vinp_RCL[i])*(Vinp_RCL[i]))) * Vinp_RCL_err[i]*Vinp_RCL_err[i]
        );

        bandCut_Gain_err[i] = sqrt(
            (1./(Vinp_RCL[i])) * (1./(Vinp_RCL[i])) * Vcut_RCL_err[i] * Vcut_RCL_err[i] +
            ((Vcut_RCL[i])/((Vinp_RCL[i])*(Vinp_RCL[i]))) * ((Vcut_RCL[i])/((Vinp_RCL[i])*(Vinp_RCL[i]))) * Vinp_RCL_err[i]*Vinp_RCL_err[i]
        );

        bandPass_Gain_dB[i] = getDB(Vinp_RCL[i], Vpas_RCL[i]);
        bandCut_Gain_dB[i] = getDB(Vinp_RCL[i], Vcut_RCL[i]);

        bandPass_Gain_dB_err[i] = getDBerr(Vinp_RCL[i], Vpas_RCL[i], Vinp_RCL_err[i], Vpas_RCL_err[i]);
        bandCut_Gain_dB_err[i] = getDBerr(Vinp_RCL[i], Vcut_RCL[i], Vinp_RCL_err[i], Vcut_RCL_err[i]);
        
        printout <<  std::setprecision(2) << "\n" << freq_RCL[i] << "\t +/- \t" << freq_RCL_err[i] << "\t\t\t\t" << bandCut_Gain[i] << "\t +/- \t" << bandCut_Gain_err[i] << "\t\t\t\t";
        printout <<  std::setprecision(2) <<  bandCut_Gain_dB[i] << "\t +/- \t" << bandCut_Gain_dB_err[i]   << "\t\t\t\t" << bandPass_Gain[i] << "\t +/- \t" << bandPass_Gain_err[i] << "\t\t\t\t";
        printout <<  std::setprecision(2) <<  bandPass_Gain_dB[i] << "\t +/- \t" << bandPass_Gain_dB_err[i];

        band_cut_gain ->SetPoint(i, freq_RCL[i], bandCut_Gain[i]);     
        band_pass_gain->SetPoint(i, freq_RCL[i], bandPass_Gain[i]);    
        band_cut_dB->SetPoint(i, freq_RCL[i], bandCut_Gain_dB[i]);
        band_pass_dB->SetPoint(i, freq_RCL[i], bandPass_Gain_dB[i]); 

        band_cut_gain ->SetPointError(i, freq_RCL_err[i], bandCut_Gain_err[i]);     
        band_pass_gain->SetPointError(i, freq_RCL_err[i], bandPass_Gain_err[i]);    
        band_cut_dB->SetPointError(i, freq_RCL_err[i], bandCut_Gain_dB_err[i]);
        band_pass_dB->SetPointError(i, freq_RCL_err[i], bandPass_Gain_dB_err[i]); 
    }

    band_cut_gain->Sort();
    band_pass_gain->Sort();  
    band_cut_dB->Sort();   
    band_pass_dB->Sort();  
    double vv[2] = {100.,300000.};
    double ff[2] = {-3.,-3.};
    //TGraph       * refLine3dB   = new TGraph(2, vv, ff);
      
    band_pass_gain->Write("Passa banda dati");
    band_cut_gain->Write("Taglia banda dati");
    band_pass_dB->Write("dB - Passa banda dati");
    band_cut_dB->Write("dB - Taglia banda dati");
    //refLine3dB->Write("Ref a -3 dB");

    //________________________________________________________FIT BANDA PASSANTE RCL_________________________________________

    printout <<  std::setprecision(5);
   
    //"R/(sqrt(R^2 + (pi*Q*R/(2*f) - Q*R/f)^2))" (Passa Banda)
    TF1 * pass_fit    = new TF1("pass_fit", "[1]/sqrt(1+[2]*[2]*((x*x-[0]*[0])/(x*[0]))*((x*x-[0]*[0])/(x*[0])))",2,50000);
    //TF1 * pass_fit_dB = new TF1("pass_fit_dB", "20*log([3]/(sqrt(([3]+[0])^2+(-[1]*[2]*([3]+[0])/x+[1]*x*([3]+[0])/[2])^2))) / log(10)", 2, 50000);
    //[0] = F0
    //[1] = R/R+Rs
    //[2] = Band with
    pass_fit       ->SetParameter(0, 10000);      //Hz
    pass_fit       ->SetParameter(1, 1);       //Adimensionale
    pass_fit       ->SetParameter(2, 11);   //Ohm/H

    //pass_fit_dB    ->SetParameter(0, 20.);      //Ohm
    //pass_fit_dB    ->SetParameter(1, 5.);       //Adimensionale
    //pass_fit_dB    ->SetParameter(2, 10000.);   //Hz
    //pass_fit_dB    ->FixParameter(3, R);        //Ohm

    //band_pass_dB ->Fit(pass_fit_dB, "R");
    band_pass_gain ->Fit(pass_fit, "R");
    printout << "\n\n";
    printout << "\n\nFit filtro passa banda:\n";
    printout << "-----------------------------";
    printout << "\nchi^2 = " << pass_fit->GetChisquare();    
    printout << "\nndf   = " << pass_fit->GetNDF();           
    printout << "\np-val = " << pass_fit->GetProb();         
    printout << "\nfo    = " << pass_fit->GetParameter(0)   << " +/- " <<  pass_fit->GetParError(0);
    printout << "\nR/Rs     = " << pass_fit->GetParameter(1)   << " +/- " <<  pass_fit->GetParError(1);
    printout << "\nBw   = " << pass_fit->GetParameter(2)   << " +/- " <<  pass_fit->GetParError(2);
   



    //"R/(sqrt(R^2 + (pi*Q*R/(2*f) - Q*R/f)^2))" (Taglia Banda)
    auto * cut_fit    = new TF1("cut_fit", "(sqrt([0]^2+(-[1]*[2]*([3]+[0])/x+[1]*x*([3]+[0])/[2])^2))/(sqrt(([3]+[0])^2+(-[1]*[2]*([3]+[0])/x+[1]*x*([3]+[0])/[2])^2))", 10., 50000.);
    auto * cut_fit_dB = new TF1("cut_fit", "20*log((sqrt([0]^2+(-[1]*[2]*([3]+[0])/x+[1]*x*([3]+[0])/[2])^2))/(sqrt(([3]+[0])^2+(-[1]*[2]*([3]+[0])/x+[1]*x*([3]+[0])/[2])^2))) / log(10)", 10., 50000);
    //[0] = Resistenza Rs parassita in serie all'induttore
    //[1] = Fattore di merito Q
    //[2] = Frequenza di risonanza
    //[3] = Termine noto bloccato, valore della resistenza R
    cut_fit    ->SetParLimits(0, 43., 50.);     //Ohm
    cut_fit    ->SetParLimits(1, 0.8,2.);       //Adimensionale
    cut_fit    ->SetParLimits(2, 8000.,11000); //Hz
    cut_fit    ->FixParameter(3, R);            //Ohm

    cut_fit    ->SetParameter(0, 100.);         //Ohm
    cut_fit    ->SetParameter(1, 1.);           //Adimensionale
    cut_fit    ->SetParameter(2, 10000);        //Hz
    cut_fit    ->FixParameter(3, R);            //Ohm
    
    cut_fit_dB ->SetParameter(0, 100.);         //Ohm
    cut_fit_dB ->SetParameter(1, 1.);           //Adimensionale
    cut_fit_dB ->SetParameter(2, 10000);        //Hz
    cut_fit_dB ->FixParameter(3, R);            //Ohm

    band_cut_dB ->Fit(cut_fit_dB);
    band_cut_gain->Fit(cut_fit); 
    printout << "\n\n";
    printout << "\n\nFit filtro taglia banda:\n";
    printout << "-----------------------------";
    printout << "\nchi^2 = " << cut_fit->GetChisquare();    
    printout << "\nndf   = " << cut_fit->GetNDF();           
    printout << "\np-val = " << cut_fit->GetProb();         
    printout << "\nRs    = " << cut_fit->GetParameter(0)   << " +/- " <<  cut_fit->GetParError(0);
    printout << "\nQ     = " << cut_fit->GetParameter(1)   << " +/- " <<  cut_fit->GetParError(1);
    printout << "\nf_r   = " << cut_fit->GetParameter(2)   << " +/- " <<  cut_fit->GetParError(2);
    printout << "\nR !L! = " << cut_fit->GetParameter(3)   << " +/- " <<  cut_fit->GetParError(3);
    

    //(Canvas plot passa banda)
    TCanvas * band_pass_c = new TCanvas("band_pass_c","Filtro passa banda",10,3,600,400);
    band_pass_c->cd();
    band_pass_c->SetLogx();
    band_pass_gain->SetTitle("Risposta in frequenza filtro passa banda");
    band_pass_gain->Sort();
    band_pass_gain->SetMarkerStyle(23);
    band_pass_gain->GetXaxis()->SetTitle("Freq. [Hz]");
    band_pass_gain->GetYaxis()->SetTitle("Guadagno");
    band_pass_c->SetGrid();
    band_pass_gain->Write("Passa Banda");
    band_pass_gain->Draw();
    band_pass_c->Draw();
    band_pass_c->Write("Canvas - passa Banda");

    //(Canvas plot taglia banda)
    TCanvas * band_cut_c = new TCanvas("band_cut_c","Filtro soppressore di banda",10,3,600,400);
    band_cut_c->cd();
    band_cut_c->SetLogx();
    band_cut_gain->SetTitle("Risposta in frequenza filtro soppressore di banda");
    band_cut_gain->Sort();
    band_cut_gain->SetMarkerStyle(23);
    band_cut_gain->GetXaxis()->SetTitle("Freq. [Hz]");
    band_cut_gain->GetYaxis()->SetTitle("Guadagno");
    band_cut_c->SetGrid();
    band_cut_gain->Write("Taglia Banda");
    band_cut_gain->Draw();
    band_cut_c->Draw();
    band_cut_c->Write("Canvas - taglia Banda");


    

//________________________________________________________Matrice di covarianza_________________________________________

    

    // Get fit results
    TFitResultPtr fitResult = band_pass_gain->Fit(pass_fit,"S");

    // Get total covariance matrix (will be used for integral error)
    TMatrixDSym cov = fitResult->GetCovarianceMatrix();
    //TMatrixD cov = pass_fit->GetCovarianceMatrix();

    double RsR_fit        = pass_fit->GetParameter(1);
    double RsR_fit_var    = cov[1][1];
    double f0_fit       = pass_fit->GetParameter(0);
    double f0_fit_var   = cov[0][0];
    double Bw_fit       = pass_fit->GetParameter(0);
    double Bw_fit_var   = cov[2][2];
 
    double cov_RsR_f0     = cov[1][0];
    double cov_Bw_RsR     = cov[1][2];
    double cov_f0_Bw    = cov[0][2];
 
    //Ricaviamo L e C:
    double L_fit = (RsR_fit/Bw_fit);
    double C_fit = 1/(4*pi*pi*L_fit*f0_fit*f0_fit);
    
    double L_fit_err = sqrt(
        1./(Bw_fit*Bw_fit)*RsR_fit_var+
        RsR_fit*RsR_fit/(Bw_fit*Bw_fit*Bw_fit*Bw_fit)*Bw_fit_var+
        1./(Bw_fit)*RsR_fit/(Bw_fit*Bw_fit)*cov_Bw_RsR
    );

    double C_fit_err = sqrt(
        1./(4*L_fit*L_fit*f0_fit*f0_fit*pi*pi)*1./(4*L_fit*L_fit*f0_fit*f0_fit*pi*pi) * L_fit_err* L_fit_err+
        1./(2*L_fit*f0_fit*f0_fit*f0_fit*pi*pi)*1./(2*L_fit*f0_fit*f0_fit*f0_fit*pi*pi) * f0_fit_var
    );

    //Vedere foglio di calcolo simbolico XCAS


    printout << "\n\n L e C ricavate dai parametri:  " << L_fit*1000 << " +/- " << L_fit_err*1000 << "  mH   " << C_fit*1e9 << " +/- " << C_fit_err*1e9 << "  nF"; 
    
    
    //________________________________________________________Guadagno al variare di R_________________________________________


    for (int i = 0; i < rnum_vre - 1; ++i)
    {
        Gain_Vre[i] = Vout_Vre[i] / Vinp_Vre[i];

        Gain_Vre_err[i] = sqrt(
            (1./(Vinp_Vre[i])) * (1./(Vinp_Vre[i])) * Vout_Vre_err[i] * Vout_Vre_err[i] +
            ((Vout_Vre[i])/((Vinp_Vre[i])*(Vinp_Vre[i]))) * ((Vout_Vre[i])/((Vinp_Vre[i])*(Vinp_Vre[i]))) * Vinp_Vre_err[i]*Vinp_Vre_err[i]
        );
    }

    auto * gainVre = new TGraphErrors(rnum_vre, Rmea_Vre, Gain_Vre, Rmea_Vre_err, Gain_Vre_err);
    gainVre->Sort();
    gainVre->SetMarkerStyle(21);
    gainVre->SetTitle("Guadagno al variare di R");
    gainVre->GetXaxis()->SetTitle("Resistenza R [Ohm]");
    gainVre->GetYaxis()->SetTitle("Guadagno");
    gainVre->GetXaxis()->SetLimits(-900,6000);
    gainVre->RemovePoint(0); 
    auto * gainVreFit = new TF1("gainVreFit" , "[0]*1/(x+[1])",0,5500);
    gainVre->Fit(gainVreFit, "S");
    printout << "\n\nGuadagno vs resistenza R";
    printout << "\n------------------------\n";
    printout << "Chi^2 = " << gainVreFit->GetChisquare() << "\n";
    printout << "N.DOF = " << gainVreFit->GetNDF() << "\n";
    printout << "p-val = " << gainVreFit->GetProb() << "\n";
    gainVre->Write("Gain vs R");
    auto * gainVreC = new TCanvas("gainVreC", "Guadagno al variare di R", 10,3,600,400);
    gainVreC->cd();
    gainVre->Draw();
    TLine *line = new TLine(-900,1.,6000,1.);
    line->SetLineColor(kRed);
    gainVreC->SetGrid();
    line->Draw();
    gainVreC->Write("Canvas - Gain vs R");
    gainVreC->Draw();

    //___________________________________________________ Ampiezza di banda al variare di R _________________________________________

   const int nFrq = 7;

   printout << "\n\n\n\n";
   auto * tR1 = new TGraph(nFrq);
   auto * tR2 = new TGraph(nFrq);
   auto * tR3 = new TGraph(nFrq);

   for (int i = 0; i < nFrq +1; ++i)
   {
       double gain_bwd = Vout_bwd[i] / Vinp_bwd[i];
       tR1->SetPoint(i, Freq_bwd[i], gain_bwd);
   }

   for (int i = nFrq; i < 2*nFrq + 2; ++i)
   {
       double gain_bwd = Vout_bwd[i] / Vinp_bwd[i];
       tR2->SetPoint(i, Freq_bwd[i], gain_bwd);
       
   }

   for (int i = nFrq*2 + 2; i < 3*nFrq + 3; ++i)
   {
       double gain_bwd = Vout_bwd[i] / Vinp_bwd[i];
       tR3->SetPoint(i, Freq_bwd[i], gain_bwd);  
   }

    //tR1->RemovePoint(0);
    //tR2->RemovePoint(0);
    //tR3->RemovePoint(0);

    auto * gainBwdC = new TCanvas("gainBwdC", "Banda passante al variare di R", 10,3,600,400);
    TF1  * tfR1     = new TF1("pass_fit", "[3]/(sqrt(([3]+[0])^2+(-[1]*[2]*([3]+[0])/x+[1]*x*([3]+[0])/[2])^2))",2,50000);
    TF1  * tfR2     = new TF1("pass_fit", "[3]/(sqrt(([3]+[0])^2+(-[1]*[2]*([3]+[0])/x+[1]*x*([3]+[0])/[2])^2))",2,50000);
    TF1  * tfR3     = new TF1("pass_fit", "[3]/(sqrt(([3]+[0])^2+(-[1]*[2]*([3]+[0])/x+[1]*x*([3]+[0])/[2])^2))",2,50000);
    
    tfR1       ->SetParameter(0, 20.);      //Ohm
    tfR1       ->SetParameter(1, 5.);       //Adimensionale
    tfR1       ->SetParameter(2, 10000.);   //Hz
    tfR1       ->SetParameter(3, R);        //Ohm

    tfR2       ->SetParameter(0, 20.);      //Ohm
    tfR2       ->SetParameter(1, 5.);       //Adimensionale
    tfR2       ->SetParameter(2, 10000.);   //Hz
    tfR2       ->SetParameter(3, R);        //Ohm

    tfR3       ->SetParameter(0, 20.);      //Ohm
    tfR3       ->SetParameter(1, 5.);       //Adimensionale
    tfR3       ->SetParameter(2, 10000.);   //Hz
    tfR3       ->SetParameter(3, R);        //Ohm
    
    gainBwdC->cd();
    gainBwdC->SetLogx();
    tR1->GetXaxis()->SetLimits(100,200000);
    tR2->GetXaxis()->SetLimits(100,200000);
    tR3->GetXaxis()->SetLimits(100,200000);
    tfR1->SetLineColor(4);
    tfR2->SetLineColor(5);
    tfR3->SetLineColor(2);
    tR1->GetXaxis()->SetTitle("F [Hz]");
    tR1->GetYaxis()->SetTitle("Guadagno");
    tR2->GetXaxis()->SetTitle("F [Hz]");
    tR2->GetYaxis()->SetTitle("Guadagno");
    tR3->GetXaxis()->SetTitle("F [Hz]");
    tR3->GetYaxis()->SetTitle("Guadagno");
    tR1->SetTitle("Confronto ampiezza di banda al variare di R");
    tR2->SetTitle("Confronto ampiezza di banda al variare di R");
    tR3->SetTitle("Confronto ampiezza di banda al variare di R");
    tR1->Fit(tfR1);
    tR2->Fit(tfR2);
    tR3->Fit(tfR3);
    tR1->Sort();
    tR2->Sort();
    tR3->Sort();
    tR1->SetMarkerStyle(21);
    tR2->SetMarkerStyle(28);
    tR3->SetMarkerStyle(3);
    //tR2->Draw();
    tR3->Draw();
    tR1->Draw("same");
    tR2->Write("bw R2");
    tR3->Write("bw R3");
    tR1->Write("bw R1");
    //Legenda
    TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("Legenda","C");
    legend->AddEntry(tfR1, "Resistenza R1","l");
    legend->AddEntry(tfR3, "Resistenza R2","l");
    gainBwdC->SetTitle("Confronto ampiezza di banda al variare di R");
    legend->Draw();
    gainBwdC->Write("Canvas - bandwidth");
    gainBwdC->Draw();

    //---------------------------TEST Z PER L E C-------------------------------------------------
    double om = 2*3.1415926535*1000;
    double Rs_teo = om*L/Q_L;
    double Rs_teo_err = sqrt(om*om/(Q_L*Q_L)*L_err*L_err+om*om*L*L/(Q_L*Q_L*Q_L*Q_L)*Q_L_err*Q_L_err);
    double fr_teo = 1/(2*3.1415926535*sqrt(L*C));
    double pig = 3.1415926535;
    double fr_teo_err = sqrt(
        C*L/(16*C*C*L*L*L*L*pig*pig) * L_err * L_err +
        C*L/(16*C*C*L*L*C*C*pig*pig) * C_err * C_err 
    );

    /*printout << "\n\nTest Z compatibilità Rs ";
    double Z_1 = (Rs_teo - Rs_fit) / sqrt(Rs_teo_err * Rs_teo_err + Rs_fit_var);
    double pval1 = st.pvalZtwotailed(Z_1);
    printout << "\nZ = " << Z_1 << "\t\t p-val 2-tailed = " << pval1;

    printout << "\n\nTest Z compatibilità fr";
    double Z_2 = (fr_teo - f0_fit) / sqrt(fr_teo_err * fr_teo_err + f0_fit_var);
    double pval2 = st.pvalZtwotailed(Z_2);
    printout << "\nZ = " << Z_2 << "\t\t p-val 2-tailed = " << pval2;

    printout << "\n\nTest Z compatibilità L";
    double Z_3 = (L - L_fit) / sqrt(L_err * L_err + L_fit_err*L_fit_err);
    double pval3 = st.pvalZtwotailed(Z_3);
    printout << "\nZ = " << Z_2 << "\t\t p-val 2-tailed = " << pval3;

    printout << "\n\nTest Z compatibilità C";
    double Z_4 = (C - C_fit / sqrt(C_err*C_err + C_fit_err*C_fit_err));
    double pval4 = st.pvalZtwotailed(Z_4);
    printout << "\nZ = " << Z_4 << "\t\t p-val 2-tailed = " << pval4;

    printout << "\n\nTest Z compatibilità le frequenze di taglio";
    double Z_5 = (cut_fit->GetParameter(2) - pass_fit->GetParameter(2)) / sqrt(cut_fit->GetParError(2) * cut_fit->GetParError(2) + pass_fit->GetParError(2)*pass_fit->GetParError(2));
    double pval5 = st.pvalZtwotailed(Z_5);
    printout << "\nZ = " << Z_5 << "\t\t p-val 2-tailed = " << pval5;*/

    printout << "\nRs_teo= " << Rs_teo << "\t\tfr_teo" << fr_teo << "\nRs_teo_err= " << fr_teo_err;
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

    Quasi sicuramente nel fit per l'induttanza bisognerà prevedere la modellizzazione con la Rs in serie, altrimenti chi2 rimane troppo alto.
    https://tinyurl.com/ModelEsp4
    
    Nella relazione la prima parte sul rifasamento non è necessaria, bisogna saperla fare per la prova pratica.

    Nei fit non si usano mai i valori di Rs,C,L misurati con il ponte RCL, si ricavano invece dai parametri del fit che devono essere facilmente
    inizializzabili e non devono essere necessariamente i valori di L,C, Rs sconosciuti.
    Quando si costruisce la funzione di fit, si prevede un parametro libero da cui ricavare la Rs (resistenza in serie all'induttore dovuta al filo avvolto).
    Si considera invece nota la resistenza R, altrimenti non si potrebbe ricavare Rs.

    Parametri consigliati:
    R/Rs
    Frisonanza
    larghezza di banda
    etc.

    PROBLEMA:       L, C e f0 non formano un set indipendente di parametri (L e C dei componenti dipendono debolmente dalla frequenza).
                    Quindi si tiene conto della matrice di covarianza nel calcolo di L,C, f0 a partire dai parametri di fit.

    A seconda della frequenza utilizzata per la caratterizzazione del filtro e della Rs, il guadagno nella frequenza di risonanza del filtro passa
    banda potrebbe differire da 1.

    SCRIVERE LE FORME FUNZIONALI IN MODO FACILMENTE INTERPRETABILE

    Si utilizzano poi test Z per ricavare L,C, Rs misurati con il ponte e ottenuti dal fit. Confrontando f_ris_ideale e f_ris_fit stiamo valutando
    quanto le approssimazioni ad elementi ideali siano buone, se differiscono vuol dire che ci sono induttanze e capacità parassite.
    Cercare di inferire quali possono essere le cause di una eventuale capacità o induttanza parassita osservata.


*/

/* 
    eMail prof. Bellan del 20/02/21: 
    L'induttore utilizzato manifesta un'induttanza variabile al crescere della frequenza,
    quindi siamo autorizzati a restringere il fit tagliando parte delle alte frequenze.
*/