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
const char rootfilepath[256]       = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper4.root";
const char outfilepath[256]        = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper4_results.txt";
const char biconvessa[256]         = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 4-5/E4_biconvessa.csv";
const char biconcava[256]          = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 4-5/E4_biconcava.csv";



//Percorsi Federica
//const char rootfilepath[256]      = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper4.root";
//const char outfilepath[256]       = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper4_results.txt";
//const char biconvessa[256]        = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 4-5/E4_biconvessa.csv";
//const char biconcava[256]         = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 4-5/E4_biconcava.csv";

//Percorsi Filippo
//const char rootfilepath[256]      = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper4.root";
//const char outfilepath[256]       = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - ottica/esper4_results.txt";
//const char biconvessa[256]        = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 4-5/E4_biconvessa.csv";
//const char biconcava[256]         = "/home/filippo/uni/esp2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati ottica 4-5/E4_biconcava.csv";

void esper4_1()
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
    csvdata biconvessaCsv(biconvessa); // Convergente
    biconvessaCsv.setDelimiters(';');
    csvdata biconcavaCsv(biconcava);   // Divergente
    biconcavaCsv.setDelimiters(';');

    //---------------------------- IMPORTAZIONE DATI SPERIMENTALI --------------------------------

    double zeroCondensatore = 281.e-3;   //m
    double spessoreFinestra = 3.94e-3;   //m

    double zeroCondensatore_err = 0.003; //m 
    double spessoreFinestra_err = 0.002; //m

    //Offset addizionali
    double offset_p_div = 0.277; //[m] Offset sulla posizione "p" con il sistema lente convergente + divergente
    double lens_thk_div = 0.010; //[m] Metà spessore medio del sistema lente convergente + divergente
    double offset_pconv = 0.277; //[m] Offset sulla posizione "p" con la lente convergente
    double lens_thkconv = 0.005; //[m] Metà spessore medio della lente convergente
    double offsetDi_err = 0.006; //[m] Errore sull'offset sistema di lenti convergente + divergente
    double offsetCo_err = 0.003; //[m] Errore sull'offset lente divergente

    //Calcolo dell'errore sistematico totale sulla posizione
    double p_err     = 1e-3;     //m
    double errSistDi = sqrt(zeroCondensatore_err * zeroCondensatore_err + spessoreFinestra_err * spessoreFinestra_err + offsetDi_err * offsetDi_err);
    double errSistCo = sqrt(zeroCondensatore_err * zeroCondensatore_err + spessoreFinestra_err * spessoreFinestra_err + offsetCo_err * offsetCo_err);

    double ref0 = zeroCondensatore - spessoreFinestra;


    //____________________________________________________________________________________
    printout << "\n\n\n--------------- MISURE RIPETUTE LENTE CONVERGENTE DIRITTA A 0.490";
    int N1 = 30;

    auto * convDirH490 = new TH1D("convDirH490", "Misure ripetute lente convergente diritta 490 mm", 17, 0.955, 0.972);
    for (int i = 1; i < N1; ++i)
    {
        convDirH490->Fill(biconvessaCsv.getDouble(i, 2));
        //cout << "\n" << biconvessaCsv.getDouble(i, 2);
    }
    convDirH490->Write("convDirH490");

    printout << "\n\n Media \t\t\t= "        << convDirH490->GetMean() << " m";
    printout << "\n Dev. Std. \t\t= "        << convDirH490->GetStdDev() << " m";
    printout << "\n Errore sulla media \t= " << convDirH490->GetMeanError() << " m";
    printout << "\n Numero di bins \t= "     << convDirH490->GetNbinsX() << "  (0.955, 0.972) m";
    printout << "\n Numero entries \t= "     << convDirH490->GetEntries();


    //____________________________________________________________________________________
    printout << "\n\n\n--------------- MISURE RIPETUTE LENTE CONVERGENTE RUOTATA A 0.490";
    int N2 = 30;

    auto * convRotH490 = new TH1D("convRotH490", "Misure ripetute lente convergente ruotata 490 mm", 30, 0.950, 0.980);
    for (int i = 1; i < N1; ++i)
    {
        convRotH490->Fill(biconvessaCsv.getDouble(i, 3));
        //cout << "\n" << biconvessaCsv.getDouble(i, 3);
    }
    convRotH490->Write("convRotH490");

    printout << "\n\n Media \t\t\t= "        << convRotH490->GetMean() << " m";
    printout << "\n Dev. Std. \t\t= "        << convRotH490->GetStdDev() << " m";
    printout << "\n Errore sulla media \t= " << convRotH490->GetMeanError() << " m";
    printout << "\n Numero di bins \t= "     << convRotH490->GetNbinsX() << "  (0.950, 0.980) m";
    printout << "\n Numero entries \t= "     << convRotH490->GetEntries();


    //____________________________________________________________________________________
    printout << "\n\n\n--------------- MISURE RIPETUTE LENTE CONVERGENTE DIRITTA A 0.530";
    int N3 = 30;

    auto * convDirH530 = new TH1D("convDirH530", "Misure ripetute lente convergente diritta 530 mm", 24, 0.882, 0.906);
    for (int i = 1 + N1; i < N3 + 1 + N1; ++i)
    {
        convDirH530->Fill(biconvessaCsv.getDouble(i, 2));
        //cout << "\n" << biconvessaCsv.getDouble(i, 3);
    }
    convDirH530->Write("convDirH530");

    printout << "\n\n Media \t\t\t= "        << convDirH530->GetMean() << " m";
    printout << "\n Dev. Std. \t\t= "        << convDirH530->GetStdDev() << " m";
    printout << "\n Errore sulla media \t= " << convDirH530->GetMeanError() << " m";
    printout << "\n Numero di bins \t= "     << convDirH530->GetNbinsX() << "  (0.882, 0.906) m";
    printout << "\n Numero entries \t= "     << convDirH530->GetEntries();


    //____________________________________________________________________________________
    printout << "\n\n\n--------------- MISURE RIPETUTE LENTE CONVERGENTE RUOTATA A 0.530";
    int N4 = 100;

    auto * convRotH530 = new TH1D("convRotH530", "Misure ripetute lente convergente ruotata 530 mm", 17, 0.886, 0.902);
    for (int i = 1 + N1; i < N3 + 1 + N1; ++i)
    {
        convRotH530->Fill(biconvessaCsv.getDouble(i, 3));
        //cout << "\n" << biconvessaCsv.getDouble(i, 3);
    }
    convRotH530->Write("convRotH530");

    printout << "\n\n Media \t\t\t= "        << convRotH530->GetMean() << " m";
    printout << "\n Dev. Std. \t\t= "        << convRotH530->GetStdDev() << " m";
    printout << "\n Errore sulla media \t= " << convRotH530->GetMeanError() << " m";
    printout << "\n Numero di bins \t= "     << convRotH530->GetNbinsX() << "  (0.886, 0.902) m";
    printout << "\n Numero entries \t= "     << convRotH530->GetEntries();

    //____________________________________________________________________________________
    printout << "\n\n\n--------------- MISURE RIPETUTE LENTE DIVERGENTE DIRITTA A 0.620";
    int N5   = 30;
    int N5_0 = 29;

    auto * divDirH620 = new TH1D("divDirH620", "Misure ripetute lente divergente diritta 620 mm", 15, 1.150, 1.170);
    for (int i = N5_0; i < N5_0 + N5; ++i)
    {
        divDirH620->Fill(biconcavaCsv.getDouble(i, 2));
    }
    divDirH620->Write("divDirH620");

    printout << "\n\n Media \t\t\t= "        << divDirH620->GetMean() << " m";
    printout << "\n Dev. Std. \t\t= "        << divDirH620->GetStdDev() << " m";
    printout << "\n Errore sulla media \t= " << divDirH620->GetMeanError() << " m";
    printout << "\n Numero di bins \t= "     << divDirH620->GetNbinsX() << "  (1.150, 1.170) m";
    printout << "\n Numero entries \t= "     << divDirH620->GetEntries();

    //____________________________________________________________________________________
    //printout << "\n\n\n--------------- REGRESSIONE LINEARE LENTE CONVERGENTE";
    int start_l = 61;
    int end_l   = 59;
    int n_mea   = 59;

    //----------------------------- STAMPA DEGLI ISTOGRAMMI -------------------------------
    auto * histogramC = new TCanvas("histogramC", "Istogrammi", 800,800);
    histogramC->Divide(2,3);
    histogramC->cd(1);
    convDirH490->GetXaxis()->SetTitle("Distanza q [m]");
    convDirH490->Draw();
    histogramC->cd(2);
    convRotH490->GetXaxis()->SetTitle("Distanza q [m]");
    convRotH490->Draw();
    histogramC->cd(3);
    convDirH530->GetXaxis()->SetTitle("Distanza q [m]");
    convDirH530->Draw();
    histogramC->cd(4);
    convRotH530->GetXaxis()->SetTitle("Distanza q [m]");
    convRotH530->Draw();
    histogramC->cd(5);
    divDirH620->GetXaxis()->SetTitle("Distanza q [m]");
    divDirH620->Draw();
    histogramC->SetGrid();
    histogramC->Update();
    histogramC->Write("istogrammi");

    
    auto * cvCanvas = new TCanvas("cvCanvas", "Lente convergente", 800, 400);
    cvCanvas->cd();
    auto * cvGraph  = new TGraphErrors();

    std::vector<double> conv_p;
    std::vector<double> conv_q;
    conv_p.push_back(0.490);
    conv_q.push_back((convDirH490->GetMean() + convRotH490->GetMean())/2.);
    conv_p.push_back(0.530);
    conv_q.push_back((convDirH530->GetMean() + convRotH530->GetMean())/2.);

    conv_p.push_back(biconvessaCsv.getDouble(61, 0));
    conv_q.push_back((biconvessaCsv.getDouble(61, 2) + biconvessaCsv.getDouble(62, 2) + biconvessaCsv.getDouble(63, 2)) / 3.);
    conv_p.push_back(biconvessaCsv.getDouble(64, 0));
    conv_q.push_back((biconvessaCsv.getDouble(64, 2) + biconvessaCsv.getDouble(65, 2) + biconvessaCsv.getDouble(66, 2) + biconvessaCsv.getDouble(67, 2)) / 4.);
    conv_p.push_back(biconvessaCsv.getDouble(68, 0));
    conv_q.push_back((biconvessaCsv.getDouble(68, 2) + biconvessaCsv.getDouble(69, 2) + biconvessaCsv.getDouble(70, 2) + biconvessaCsv.getDouble(71, 2)) / 4.);
    conv_p.push_back(biconvessaCsv.getDouble(72, 0));
    conv_q.push_back((biconvessaCsv.getDouble(72, 2) + biconvessaCsv.getDouble(73, 2) + biconvessaCsv.getDouble(74, 2) + biconvessaCsv.getDouble(75, 2)) / 4.);
    conv_p.push_back(biconvessaCsv.getDouble(76, 0));
    conv_q.push_back((biconvessaCsv.getDouble(76, 2) + biconvessaCsv.getDouble(77, 2) + biconvessaCsv.getDouble(78, 2) + biconvessaCsv.getDouble(79, 2)) / 4.);

    for (int k = 0; k < conv_p.size(); ++k)
    {
        cvGraph->SetPoint(k, 1. / (conv_p[k] - offset_pconv + lens_thkconv), 1. / (conv_q[k] - conv_p[k] - lens_thkconv));
        cvGraph->SetPointError(k, p_err, sqrt(convDirH490->GetStdDev() * convDirH490->GetStdDev() + errSistCo * errSistCo));
    }
    cvGraph->Sort(); 
    cvGraph->RemovePoint(3);

    auto * fcv  = new TF1("fcv" , "[0]*x+ (1./[1])", -2000., 2000.);
    fcv->SetParameter(1, 1.);
    cvGraph->Fit(fcv);
    
    cvGraph->GetXaxis()->SetTitle("1/p [1/m]");
    cvGraph->GetYaxis()->SetTitle("1/q [1/m]");
    cvGraph->SetTitle("Lente convergente");
    cvGraph->Write("Graph Convergente");
    cvGraph->Draw();
    cvCanvas->SetGrid();
    cvCanvas->Write("Conv. canvas");
    cvCanvas->Draw();

    //____________________________________________________________________________________
    //printout << "\n\n\n--------------- REGRESSIONE LINEARE LENTE DIVERGENTE";

    auto * divCanvas = new TCanvas("divCanvas", "Sistema di lenti convergente + divergente", 800, 400);
    divCanvas->cd();
    auto * divGraph  = new TGraphErrors();

    std::vector<double> div_p;
    std::vector<double> div_q;
    div_p.push_back(0.620);
    div_q.push_back(divDirH620->GetMean());

    div_p.push_back(biconcavaCsv.getDouble(1, 0));
    div_q.push_back((biconcavaCsv.getDouble(1, 2) + biconcavaCsv.getDouble(2, 2) + biconcavaCsv.getDouble(3, 2) + biconcavaCsv.getDouble(4, 2)) / 4.);
    div_p.push_back(biconcavaCsv.getDouble(5, 0));
    div_q.push_back((biconcavaCsv.getDouble(5, 2) + biconcavaCsv.getDouble(6, 2) + biconcavaCsv.getDouble(7, 2) + biconcavaCsv.getDouble(8, 2)) / 4.);
    div_p.push_back(biconcavaCsv.getDouble(9, 0));
    div_q.push_back((biconcavaCsv.getDouble(9, 2) + biconcavaCsv.getDouble(10, 2) + biconcavaCsv.getDouble(11, 2) + biconcavaCsv.getDouble(12, 2)) / 4.);
    div_p.push_back(biconcavaCsv.getDouble(13, 0));
    div_q.push_back((biconcavaCsv.getDouble(13, 2) + biconcavaCsv.getDouble(14, 2) + biconcavaCsv.getDouble(15, 2) + biconcavaCsv.getDouble(16, 2)) / 4.);
    div_p.push_back(biconcavaCsv.getDouble(17, 0));
    div_q.push_back((biconcavaCsv.getDouble(17, 2) + biconcavaCsv.getDouble(18, 2) + biconcavaCsv.getDouble(19, 2) + biconcavaCsv.getDouble(20, 2)) / 4.);
    div_p.push_back(biconcavaCsv.getDouble(25, 0));
    div_q.push_back((biconcavaCsv.getDouble(25, 2) + biconcavaCsv.getDouble(26, 2) + biconcavaCsv.getDouble(27, 2) + biconcavaCsv.getDouble(28, 2)) / 4.);


    for (int k = 0; k < div_p.size(); ++k)
    {
        if ((div_p[k]  > 1e-5) && (div_q[k] > 1e-5))
        {
            divGraph->SetPoint(k, 1. / (div_p[k] - offset_p_div + lens_thk_div), 1. / (div_q[k] - div_p[k] - lens_thk_div));
            divGraph->SetPointError(k, p_err, sqrt(divDirH620->GetStdDev() * divDirH620->GetStdDev() / 3. + errSistDi * errSistDi));
        }
    }
    divGraph->Sort(); 
    divGraph->Write("Graph Divergente");

    //____________________________________________________________________________________
    printout << "\n\n\n--------------- FIT TGraph";

    auto * fdiv = new TF1("fdiv", "[0]*x+ (1./[1])", -2000., 2000.);
    fdiv->SetParameter(1, 1.);
    divGraph->Fit(fdiv);

    divGraph->SetTitle("Sistema lente convergente + divergente");
    divGraph->GetXaxis()->SetTitle("1/p [1/m]");
    divGraph->GetYaxis()->SetTitle("1/q [1/m]");
    divGraph->Draw();
    divCanvas->SetGrid();
    divCanvas->Write("Div. canvas");
    divCanvas->Draw();



    printout << "\n\nFit Lente Divergente\n";
    printout << "Chi^2  "   << fdiv->GetChisquare();
    printout << "\np-val  " << fdiv->GetProb();
    printout << "\nndof   " << fdiv->GetNDF();
    printout << "\nrX^2   " << fdiv->GetChisquare() / fdiv->GetNDF();
    printout << "\nslope  " << fdiv->GetParameter(0) << " +/- " << fdiv->GetParError(0) << " ";
    printout << "\nf      " << fdiv->GetParameter(1) << " +/- " << fdiv->GetParError(1) << " m";
    double pp1 = st.pvalZtwotailed((fdiv->GetParameter(0) -1) / fdiv->GetParError(0));
    printout << "\np_p_0  " << pp1 << "\n";

    printout << "\n\nFit Lente Convergente\n";
    printout << "Chi^2  "   << fcv->GetChisquare();
    printout << "\np-val  " << fcv->GetProb();
    printout << "\nndof   " << fcv->GetNDF();
    printout << "\nrX^2   " << fcv->GetChisquare() / fcv->GetNDF();
    printout << "\nslope  " << fcv->GetParameter(0) << " +/- " << fcv->GetParError(0) << " ";
    printout << "\nf      " << fcv->GetParameter(1) << " +/- " << fcv->GetParError(1) << " m";
    double pp2 = st.pvalZtwotailed((fcv->GetParameter(0) -1) / fcv->GetParError(0));
    printout << "\np_p_0  " << pp2 << "\n";
    
    printout << "\n\nLa lunghezza focale della lente convergente è  = " << fcv->GetParameter(1)  << " +/- " << fcv->GetParError(1) << " m";
    printout << "\n\nLa lunghezza focale del sistema di lenti c+d è = " << fdiv->GetParameter(1) << " +/- " << fdiv->GetParError(1) << " m";

    //Calcolo della lunghezza focale della lente divergente:
    double lf_div       = 1./((1./fdiv->GetParameter(1)) - (1./(fcv->GetParameter(1))));
    double tt = fdiv->GetParameter(1);
    double ff = fcv->GetParameter(1);
    double tt_err = fdiv->GetParError(1);
    double ff_err = fcv->GetParError(1);

    double lf_div_err   = sqrt(
        (1/(((1/tt-1)/ff) * ((1/tt-1)/ff) * tt * tt)) * ((1/(((1/tt-1)/ff) * ((1/tt-1)/ff) * tt * tt))) * tt_err * tt_err +
        (1/(((1/tt-1)/ff) * ((1/tt-1)/ff) * ff * ff))  *  (1/(((1/tt-1)/ff) * ((1/tt-1)/ff) * ff * ff))  * ff_err * ff_err);

    printout << "\n\nLa lunghezza focale della lente divergente è = " << lf_div << " +/- " << lf_div_err;  

    printout << "\n\n\nTest di compatibilità per le misure diritto e rovescio:";

    double z490 = (convDirH490->GetMean() - convRotH490->GetMean()) / sqrt(convRotH490->GetMeanError() * convRotH490->GetMeanError() + convDirH490->GetMeanError() * convDirH490->GetMeanError());
    double z530 = (convDirH530->GetMean() - convRotH530->GetMean()) / sqrt(convRotH530->GetMeanError() * convRotH530->GetMeanError() + convDirH530->GetMeanError() * convDirH530->GetMeanError());
    double p490 = st.pvalTtwotailed(z490, 30);
    double p530 = st.pvalTtwotailed(z530, 30);

    printout << "\n490\t\t z = " << z490 << "  p = " << p490; 
    printout << "\n530\t\t z = " << z530 << "  p = " << p530 << "\n\n"; 
    
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
