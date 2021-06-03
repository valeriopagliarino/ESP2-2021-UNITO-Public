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
const char spectrum[256]          = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 1/Spettrofotometro/spec.csv";

//Percorsi Federica
//const char rootfilepath[256]      = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper1.root";
//const char outfilepath[256]       = "/Users/federicasibilla/Documentƒs/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper1_results.txt";
//const char HgLine[256]            = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 1/E1_HgLines.csv";
//const char nGlass[256]            = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 1/E1_nGlass.csv";
//const char spectrum[256]          = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 1/Spettrofotometro/spec.csv";


//Percorsi Filippo:
//const char rootfilepath[256]      = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper1.root";
//const char outfilepath[256]       = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper1_results.txt";
//const char HgLine[256]            = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 1/E1_HgLines.csv";
//const char nGlass[256]            = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 1/E1_nGlass.csv";
//const char spectrum[256]          = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 1/Spettrofotometro/spec.csv";

void esper1()
{
    double pi = TMath::Pi();

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
    csvdata spectrumCsv(spectrum);
    spectrumCsv.setDelimiters(';');
    

    //---------------------------- IMPORTAZIONE DATI SPERIMENTALI --------------------------------

    //Spettrofotometro digitale
    int N_spec = 1345;
    
    double alpha        = 45.;         //deg   Angolo al vertice del prisma
    double grating      = 1. / 300000; //m     Distanza tra le righe del reticolo (passo)
    double lampPower    = 1000;        //W  1 kW @ 1A 3B Scientific Potenza della lampada spettrale (mercurio)

    int N_HeLines       = 8;
    int N_nGlass        = 6;

    double c = pi/180.;

    //Incertezze:
    double  errLetturaAngolare = 0.01;    //deg (mezzo primo)
    double  errLambdaDig       = 1e-9;      //m
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
    double  minDevAngle[N_nGlass];      //deg
    double  minDevAngle_err[N_nGlass];  //deg
    double  lambda_v[N_nGlass];          //m


    TGraph * spec = new TGraph();
    spec->SetTitle("Spettrofotometro digitale");
    spec->GetXaxis()->SetLimits(370e-9, 700e-9);
    double maxVal = 0.;

    for (int i = 1; i < N_spec; ++i)
    {
        double v = spectrumCsv.getDouble(i, 1);
        if (v > maxVal) maxVal = v;
    }
    
    for (int i = 1; i < N_spec; ++i)
    {
        spec->SetPoint(i-1, spectrumCsv.getDouble(i, 0)  * 1e-9, (spectrumCsv.getDouble(i, 1) / maxVal));
    }

    spec->Sort();

    //printout << "\n\nAcquisito lo spettrofotometro digitale, " << spec->GetEntries() << " entries, con media " << spec->GetMean() << "\n\n";

    spec->SetLineColor(kBlack);
    spec->Write("Spettrofotometro");

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
        lambda_v[i-1]        = nGlassCsv.getDouble(i, 2);
        //minDevAngle[i-1]     = nGlassCsv.getDouble(i, 3);
        //minDevAngle_err[i-1] = nGlassCsv.getDouble(i, 4);
    }

    printout << "\n\n--------- errori Δα in radianti ---------\n\n"; 
    for (int i = 1; i < N_nGlass+1; i++)
    {
        minDevAngle[i-1]     = (nGlassCsv.getDouble(i, 5)  + nGlassCsv.getDouble(i, 12) +
        nGlassCsv.getDouble(i, 6)  + nGlassCsv.getDouble(i, 13) +
        nGlassCsv.getDouble(i, 7)  + nGlassCsv.getDouble(i, 14) +
        nGlassCsv.getDouble(i, 8)  + nGlassCsv.getDouble(i, 15) +
        nGlassCsv.getDouble(i, 9)  + nGlassCsv.getDouble(i, 16) +
        nGlassCsv.getDouble(i, 10) + nGlassCsv.getDouble(i, 17) +
        nGlassCsv.getDouble(i, 11) + nGlassCsv.getDouble(i, 18))/14;

        double varianzaDev = 0.;
        for (int j = 5; j < 19; j++)
            varianzaDev += (minDevAngle[i-1] - nGlassCsv.getDouble(i, j))*(minDevAngle[i-1] - nGlassCsv.getDouble(i, j));

        minDevAngle_err[i-1] = sqrt(varianzaDev) * c / sqrt(N_nGlass);  
        printout << "\n" << i << ")" << minDevAngle_err[i-1] << "\n";                     
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

    int colori[12] = {6,9,4,30,3,5,2,46,1,10,11,12};

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

        line1->SetLineColor(colori[i]);
        line2->SetLineColor(colori[i]);
        line3->SetLineColor(colori[i]);
        line4->SetLineColor(colori[i]);
        line5->SetLineColor(colori[i]);
        line6->SetLineColor(colori[i]);
        line7->SetLineColor(colori[i]);
        line8->SetLineColor(colori[i]);

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
    double lambda_r_err[N_HeLines];

    printout << "\nLunghezze d'onda: \n";

    for (int i = 0; i < N_HeLines; ++i)
    {
        double lambda_ord1_deg_semidisp = sqrt(pow(0.5 * (abs(abs(dx1[i]) - abs(sx1[i]))),2) + errLetturaAngolare * errLetturaAngolare); //Somma statistico e sistematico
        double lambda_ord2_deg_semidisp = sqrt(pow(0.5 * (abs(abs(dx2[i]) - abs(sx2[i]))),2) + errLetturaAngolare * errLetturaAngolare);
        double lambda_ord3_deg_semidisp = sqrt(pow(0.5 * (abs(abs(dx3[i]) - abs(sx3[i]))),2) + errLetturaAngolare * errLetturaAngolare);
        if (i == (N_HeLines - 1)) lambda_ord1_deg_semidisp = .1;

        double lambda_ord1_semidisp = cos(abs(abs(dx1[i] * c) + abs(sx1[i] * c)) / 2.) * lambda_ord1_deg_semidisp * c / (1./grating);
        double lambda_ord2_semidisp = cos(abs(abs(dx2[i] * c) + abs(sx2[i] * c)) / 2.) * lambda_ord2_deg_semidisp * c / (2./grating);
        double lambda_ord3_semidisp = cos(abs(abs(dx3[i] * c) + abs(sx3[i] * c)) / 2.) * lambda_ord3_deg_semidisp * c / (3./grating);

        double lambda_ord1 = sin(abs(dx1[i] + sx1[i]) * c / 2.) / (1. * (1./grating));
        double lambda_ord2 = sin(abs(dx2[i] + sx2[i]) * c / 2.) / (2. * (1./grating));
        double lambda_ord3 = sin(abs(dx3[i] + sx3[i]) * c / 2.) / (3. * (1./grating));
        //double lambda_ord4 = sin(abs(dx4[i] + sx4[i]) * c / 2.) / (4. * (1./grating));

        //printout << "\nLunghezze d'onda: \n";
        //st.quickZtwotailed(lambda_ord1, lambda_ord2, lambda_ord1_semidisp, lambda_ord2_semidisp, printout);
        //st.quickZtwotailed(lambda_ord1, lambda_ord3, lambda_ord1_semidisp, lambda_ord3_semidisp, printout);
        //st.quickZtwotailed(lambda_ord2, lambda_ord3, lambda_ord2_semidisp, lambda_ord3_semidisp, printout);

        std::vector<double> vm;
        std::vector<double> vm_err;
        if (lambda_ord1 > 1e-12) vm.push_back(lambda_ord1); 
        if (lambda_ord2 > 1e-12) vm.push_back(lambda_ord2); 
        if (lambda_ord3 > 1e-12) vm.push_back(lambda_ord3); 
        if (lambda_ord1 > 1e-12) vm_err.push_back(lambda_ord1_semidisp); 
        if (lambda_ord2 > 1e-12) vm_err.push_back(lambda_ord2_semidisp); 
        if (lambda_ord3 > 1e-12) vm_err.push_back(lambda_ord3_semidisp); 
        //cout << "\n\n" <<  "err____" << vm_err.size();
        
        lambda_r[i]     =  media_n_variabili(vm);
        lambda_r_err[i] =  media_n_variabili_err(vm, vm_err);
        printout << "\n" << i+1 << "] [metri]  ordine 1 (" << lambda_ord1 << ")  ordine 2 (" << lambda_ord2 << ")  ordine 3 (" << lambda_ord3 << ")";
        printout << "\n         media = " << lambda_r[i] << " +/- " << lambda_r_err[i];
        printout << "\nSemidispersione: ord1 = " << lambda_ord1_semidisp << "  ord2 = " << lambda_ord2_semidisp << "  ord3 = " << lambda_ord3_semidisp << "   v = " << vm.size() << "\n\n\n";
        

    }

    TCanvas * C2 = new TCanvas("C2","Linee spettrali del Mercurio",10,3,600,400);
    C2->cd();
    C2->SetTitle("Linee spettrali del Mercurio");

    spec->GetXaxis()->SetLimits(370e-9, 700e-9); //Fondo scala plot delle righe del mercurio
    spec->GetYaxis()->SetLimits(0, 1.);
    C2->SetGrid();
    spec->GetXaxis()->SetTitle("Lunghezza d'onda #lambda [m]");
    spec->GetYaxis()->SetTitle("Scala di intensita' relativa normalizzata");
    spec->SetTitle("Linee spettrali del Mercurio");
    spec->Draw();
    TLine * line1b;
    TLine * line1b_esup;
    TLine * line1b_einf;
    double lineH = 1.;

    printout << "\n\nLinee spettrali mercurio:\n";

    for (int i = 0; i < N_HeLines; ++i)
    {
        double xx_, yy_;
        spec->GetPoint(i, xx_, yy_);
        double zz, z_err, pvalZ;
        z_err = sqrt(lambda_r_err[i]*lambda_r_err[i] + errLambdaDig * errLambdaDig);
        zz = (lambda_r[i] - ldm[i]) / z_err;
        pvalZ = st.pvalZtwotailed(zz);
        printout << "\nLunghezza d'onda = " << lambda_r[i] << " +/- " << lambda_r_err[i] << " metri  \t|  Dig.: " << ldm[i] << " +/- " << errLambdaDig << " nm" << "   \t|  Z = " << zz << "  p-val =  " << pvalZ;
        line1b = new TLine(lambda_r[i],0.,lambda_r[i], lineH);
        line1b_esup = new TLine(lambda_r[i] + lambda_r_err[i], 0., lambda_r[i] + lambda_r_err[i], lineH);
        line1b_einf = new TLine(lambda_r[i] - lambda_r_err[i], 0., lambda_r[i] - lambda_r_err[i], lineH);
        line1b_esup->SetLineWidth(2);
        line1b_einf->SetLineWidth(2);
        switch(i)
        {
            case 0:
                line1b->SetLineColor(6);
                line1b_esup->SetLineColor(6);
                line1b_einf->SetLineColor(6);
                break;
            case 1:
                line1b->SetLineColor(9);
                line1b_esup->SetLineColor(9);
                line1b_einf->SetLineColor(9);
                break;

            case 2:
                line1b->SetLineColor(4);
                line1b_esup->SetLineColor(4);
                line1b_einf->SetLineColor(4);
                break;

            case 3:
                line1b->SetLineColor(41);
                line1b_esup->SetLineColor(41);
                line1b_einf->SetLineColor(41);
                break;

            case 4:
                line1b->SetLineColor(3);
                line1b_esup->SetLineColor(3);
                line1b_einf->SetLineColor(3);
                break;

            case 5:
                line1b->SetLineColor(92);
                line1b_esup->SetLineColor(92);
                line1b_einf->SetLineColor(92);
                break;

            case 6:
                line1b->SetLineColor(97);
                line1b_esup->SetLineColor(97);
                line1b_einf->SetLineColor(97);
                break;

            case 7:
                line1b->SetLineColor(100);
                line1b_esup->SetLineColor(100);
                line1b_einf->SetLineColor(100);
                break;
        }
        line1b->SetLineWidth(4);
        line1b->Draw();
        line1b_einf->Draw();
        line1b_esup->Draw();
    }
    //spec->Draw("same");
    C2->Draw();
    C2->Write("Lambda mercurio");
    printout << "\n\n";


    //-------------------------MISURA DI n(lambda) NEL PRISMA IN VETRO----------------------------

    //http://www.arass-brera.org/it/indice-analitico/60-inventario/ottica/rifrazione/257-prisma
    
    double glass_indices[N_nGlass];
    double indices_err  [N_nGlass];
    double lambda_n_err [N_nGlass];
    for (int i = 1; i < N_nGlass+1; i++)
    {
        lambda_n_err[i-1]=lambda_r_err[i];
        //printout << "\n" << i  << "] err = " << lambda_r_err[i]; //Assumiamo l'errore dello spettrofotometro 1 nm
        printout << "\n" << i  << "] err = " << errLambdaDig << "  nm";
    }

    printout << "\n\nIndice di rifrazione del vetro:";
    
    for (int i = 0; i < N_nGlass; ++i)
    {
        glass_indices[i] = sin(c * ((180 - (minDevAngle[i] - 180.) + alpha)) / 2.) / (sin(c * (alpha / 2.)));
        indices_err[i]   = 0.5 / (sin(c * (alpha / 2.))) * cos(c * ((180 - (minDevAngle[i] - 180.) + alpha)) / 2.) * minDevAngle_err[i];
        printout << "\n" << i + 1 << "] n = " << glass_indices[i] << "\t\t at  " << lambda_v[i]*1e9 << "nm";
    }
    
    auto * t1 = new TGraphErrors(N_nGlass);
    for (int i = 0; i < N_nGlass; ++i)
    {
        t1->SetPoint(i, lambda_v[i], glass_indices[i]);
        t1->SetPointError(i, errLambdaDig, indices_err[i]);
        //t1->SetPointError(i, lambda_n_err[i], indices_err[i]); Utilizziamo lo spettrofotometro digitale, quindi assumiamo 1nm, pari alla sensibilità del cursore che abbiamo utilizzato.
    }

    //Fit con legge di Cauchi n = A + B/λ
    auto * cauchy = new TF1("cauchy", "[0]+[1]/(x*x)", 0., 100.);
    cauchy->SetParameter(0, 1.5);
    cauchy->SetParameter(1, 1e-14);
    t1->Fit("cauchy", "R");

    cauchy->SetParameter(0,1.8);
    cauchy->SetParameter(1,1.3e-14);
    
    printout << "\n\nFit indice di rifrazione del prisma n = A + B/(λ^2)\n";
    printout << "Chi^2  "   << cauchy->GetChisquare();
    printout << "\np-val  " << cauchy->GetProb();
    printout << "\nndof   " << cauchy->GetNDF();
    printout << "\nrX^2   " << cauchy->GetChisquare() / cauchy->GetNDF();
    printout << "\nA      " << cauchy->GetParameter(0) << " +/- " << cauchy->GetParError(0) << " ";
    printout << "\nB      " << cauchy->GetParameter(1) << " +/- " << cauchy->GetParError(1) << " m";
    printout << "\n";

    TCanvas * C3 = new TCanvas("C3","Indice di rifrazione del vetro",10,3,600,400);
    C3->cd();
    C3->SetTitle("Indice di rifrazione del vetro al variare della lunghezza d'onda");
    t1->SetTitle("Indice di rifrazione del vetro al variare della lunghezza d'onda");
    t1->GetXaxis()->SetTitle("Lunghezza d'onda #lambda [m]");
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


}