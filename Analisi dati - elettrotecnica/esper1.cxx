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

//Percorsi Valerio:
//const char rootfilepath[256]  = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper1.root";
//const char outfilepath[256]   = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnicaesper1_results.txt";
//const char lampcsvpath[256]   = "/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati elettrotenica 1 e 2/lamp_IV.csv";
//Percorsi Federica:
const char rootfilepath[256] = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper1.root";
const char outfilepath[256]  = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Analisi dati - elettrotecnica/esper1_results.txt";
const char lampcsvpath[256]  = "/Users/federicasibilla/Documents/UNI-2/Università-ESP2/ESP2-GR1-FISICA-UNITO-20-21/Fisica - Esperimentazioni 2/Misure e dati sperimentali/Presa dati elettrotecnica 1 e 2/lamp_IV.csv";

//Percorsi Filippo:
//const char rootfilepath[256] =
//const char outfilepath[256]  =
//const char lampcsvpath[256]  =


void esper1()
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
    printout << "  | Experiment:  Elettrotecnica - esperienza 1                    | \n";
    printout << "  | Date:        13/12/3030                                       | \n";
    printout << "  | Revision:    2.1.0                                            | \n";
    printout << "  | Description: Caratterizzazione strumenti + IV Lamp            | \n";
    printout << "  +---------------------------------------------------------------+ \n\n";


    //Apriamo il file CSV esistente in cui è stata salvata la tabella con i dati IV della lampadina
    csvdata csv(lampcsvpath);
    csv.setDelimiters(';');

    //---------------------------- IMPORTAZIONE DATI SPERIMENTALI --------------------------------

    //Dichiariamo dei vettori dinamici e poi leggiamo i dati IV della lampadina dal file CSV

    std::vector<double> Vlamp;      //V
    std::vector<double> Ilamp;      //mA
    std::vector<double> Vlamp_err;  //V
    std::vector<double> Ilamp_err;  //mA

    for (int i = 2; i < 26; ++i)
    {
        Vlamp.push_back(csv.getDouble(i, 1));
        Ilamp.push_back(csv.getDouble(i, 3));  // (riga, colonna)  -!- Gli indici sono base-1

        Vlamp_err.push_back(csv.getDouble(i, 2));
        Ilamp_err.push_back(csv.getDouble(i, 4));
    }
    cerr << "\nInput data parsing completed.";

    printout << "\nDati sperimentali lampadina\n\nV [V]\t\tVerr [V]\t\tI [mA]\t\tIerr [mA]";

    for (int i = 0; i < Vlamp.size(); ++i)
    {
        printout << "\n" << Vlamp[i] << "\t\t" << Vlamp_err[i] << "\t\t" << Ilamp[i] << "\t\t" << Ilamp_err[i];
    }

    printout << "\n\n";

    //Resistenza interna voltmetro analogico
    double Vag_a      = 7.9;      //V
    double Vag_err_a  = 0.1;      //V
    double Vbg_a      = 2.2;      //V
    double Vbg_err_a  = 0.1;      //V
    double R470K      = 485.3e3;  //Ω
    double R470K_err  = 2.9e3;    //Ω

    //Resistenza interna volmetro digitale
    double Vag_d      = 7.994;     //V
    double Vag_err_d  = 0.013;     //V
    double Vbg_d      = 7.593;     //V
    double Vbg_err_d  = 0.013;     //V

    //Resistenza interna amperometro analogico
    double Vtot_a     = 7.992;     //V
    double Vtot_err_a = 0.013;     //V
    double Vamp_a     = 556.3e-3;  //V
    double Vamp_err_a = 1.1e-3;    //V
    double I_a        = 49.0e-3;   //A
    double I_err_a    = 0.5e-3;    //A
    double R150       = 155.76;    //Ω
    double R150_err   = 0.34;      //Ω

    //Resistenza interna amperometro digitale
    double Vtot_d     = 7.991;     //V
    double Vtot_err_d = 0.013;     //V
    double Vamp_d     = 340.7e-3;  //V
    double Vamp_err_d = 0.8e-3;    //V
    double I_d        = 47.70e-3;  //A
    double I_err_d    = 0.29e-3;   //A

    //Dati costruttore strumenti di misura da confrontare
    double  Rint_volt_a_prod        = 0.;       //Ω
    double  Rint_volt_a_prod_err    = 0.;       //Ω
    double  Rint_volt_d_prod        = 0.;       //Ω
    double  Rint_volt_d_prod_err    = 0.;       //Ω
    double  Rint_amp_a_prod         = 0.;       //Ω
    double  Rint_amp_a_prod_err     = 0.;       //Ω
    double  Rint_amp_d_prod         = 0.;       //Ω
    double  Rint_amp_d_prod_err     = 0.;       //Ω

    //Lamp
    double Rfreddo_lamp      = 50.57;           //Ω 
    double Rfreddo_lamp_err  = 0.13;            //Ω
    

    //---------------------------------- ANALISI DATI --------------------------------------------

    //------- Caratterizzazione strumenti ------------
    
    // Resistenza interna voltmetro analogico
    double Rint_volt_a     = R470K * Vbg_a / (Vag_a - Vbg_a);
    double Rint_volt_err_a = sqrt(
        ((Vbg_a / (Vag_a - Vbg_a))*(Vbg_a / (Vag_a - Vbg_a))) * R470K_err * R470K_err +
        (R470K * Vbg_a / ((Vag_a - Vbg_a) * (Vag_a - Vbg_a))) * (R470K * Vbg_a / ((Vag_a - Vbg_a) * (Vag_a - Vbg_a))) * Vag_err_a * Vag_err_a + 
        ((R470K * Vag_a) / ((Vag_a - Vbg_a) * (Vag_a - Vbg_a))) * ((R470K * Vag_a) / ((Vag_a - Vbg_a) * (Vag_a - Vbg_a))) * Vbg_err_a * Vbg_err_a
    );
    printout << "\nResistenza interna voltmetro analogico = \t\t(" << Rint_volt_a << "\t\t +/- " << Rint_volt_err_a << ")\tΩ";

    // Resistenza interna voltmetro digitale
    double Rint_volt_d     = R470K * Vbg_d / (Vag_d - Vbg_d);
    double Rint_volt_err_d = sqrt(
        ((Vbg_d / (Vag_d - Vbg_d))*(Vbg_d / (Vag_d - Vbg_d))) * R470K_err * R470K_err +
        (R470K * Vbg_d / ((Vag_d - Vbg_d) * (Vag_d - Vbg_d))) * (R470K * Vbg_d / ((Vag_d - Vbg_d) * (Vag_d - Vbg_d))) * Vag_err_d * Vag_err_d + 
        ((R470K * Vag_d) / ((Vag_d - Vbg_d) * (Vag_d - Vbg_d))) * ((R470K * Vag_d) / ((Vag_d - Vbg_d) * (Vag_d - Vbg_d))) * Vbg_err_d * Vbg_err_d
    );
    printout << "\nResistenza interna voltmetro digitale = \t\t(" << Rint_volt_d << "\t +/- " << Rint_volt_err_d << ")\tΩ";


    // Resistenza interna amperometro analogico (metodo partitore) - Si tiene conto della resistenza interna del voltmetro digitale impiegato
    double Rint_amp_a_part = Vamp_a * Rint_volt_d * R150 / (Vtot_a * Rint_volt_d - Vamp_a * R150 - Vamp_a * Rint_volt_d);
    
    double diffVamp = Rint_volt_d*R150/(Vtot_a*Rint_volt_d-Vamp_a*R150-Vamp_a*Rint_volt_d)+Vamp_a*Rint_volt_d*R150*(R150+Rint_volt_d)/((Vtot_a*Rint_volt_d-Vamp_a*R150-Vamp_a*Rint_volt_d)*(Vtot_a*Rint_volt_d-Vamp_a*R150-Vamp_a*Rint_volt_d));
    double diffRint = Vamp_a*R150/(Vtot_a*Rint_volt_d-Vamp_a*R150-Vamp_a*Rint_volt_d)+Vamp_a*Rint_volt_d*R150*(-Vtot_a+Vamp_a)/((Vtot_a*Rint_volt_d-Vamp_a*R150-Vamp_a*Rint_volt_d)*(Vtot_a*Rint_volt_d-Vamp_a*R150-Vamp_a*Rint_volt_d));
    double diffR    = Vamp_a*Rint_volt_d/(Vtot_a*Rint_volt_d-Vamp_a*R150-Vamp_a*Rint_volt_d)+Vamp_a*Rint_volt_d*R150*Vamp_a/((Vtot_a*Rint_volt_d-Vamp_a*R150-Vamp_a*Rint_volt_d)*(Vtot_a*Rint_volt_d-Vamp_a*R150-Vamp_a*Rint_volt_d));
    double diffVtot = -Vamp_a*Rint_volt_d*R150*Rint_volt_d/((Vtot_a*Rint_volt_d-Vamp_a*R150-Vamp_a*Rint_volt_d)*(Vtot_a*Rint_volt_d-Vamp_a*R150-Vamp_a*Rint_volt_d));
    
    double Rint_amp_a_err_part = sqrt( 
        diffVamp * diffVamp * Vamp_err_a * Vamp_err_a +
        diffRint * diffRint * Rint_volt_err_d * Rint_volt_err_d +
        diffR * diffR * R150_err * R150_err +
        diffVtot * diffVtot * Vtot_err_a * Vtot_err_a
    );
    printout << "\nResistenza amperometro analogico con metodo partitore = (" <<  Rint_amp_a_part << "\t +/- " << Rint_amp_a_err_part << ")\tΩ";

    // Resistenza interna amperometro analogico (metodo diretto)
    double Rint_amp_a_dir = 1. / (I_a / Vamp_a - 1. / Rint_volt_d);

    double diffI_a     = -1./((I_a/Vamp_a-1./Rint_volt_d)*(I_a/Vamp_a-1./Rint_volt_d)*Vamp_a);
    double diffVamp_a  = I_a/((I_a/Vamp_a-1./Rint_volt_d)*(I_a/Vamp_a-1./Rint_volt_d)*Vamp_a*Vamp_a);
    double diffRintvd  = -1./((I_a/Vamp_a-1/Rint_volt_d)*(I_a/Vamp_a-1./Rint_volt_d)*Rint_volt_d*Rint_volt_d);

    double Rint_amp_a_err_dir = sqrt(
        diffI_a * diffI_a * I_err_a * I_err_a +
        diffVamp_a * diffVamp_a * Vamp_err_a * Vamp_err_a +
        diffRintvd * diffRintvd * Rint_volt_err_d * Rint_volt_err_d
    );

    printout << "\nResistenza amperometro analogico con metodo diretto = \t(" <<  Rint_amp_a_dir << "\t +/- " << Rint_amp_a_err_dir << ")\tΩ";


  // Resistenza interna amperometro digitale (metodo partitore) - Si tiene conto della resistenza interna del voltmetro digitale impiegato
    double Rint_amp_d_part = Vamp_d * Rint_volt_d * R150 / (Vtot_d * Rint_volt_d - Vamp_d * R150 - Vamp_d * Rint_volt_d);
    
    double diffVamp_d = Rint_volt_d*R150/(Vtot_d*Rint_volt_d-Vamp_d*R150-Vamp_d*Rint_volt_d)+Vamp_d*Rint_volt_d*R150*(R150+Rint_volt_d)/((Vtot_d*Rint_volt_d-Vamp_d*R150-Vamp_d*Rint_volt_d)*(Vtot_d*Rint_volt_d-Vamp_d*R150-Vamp_d*Rint_volt_d));
    double diffRint_d = Vamp_d*R150/(Vtot_d*Rint_volt_d-Vamp_d*R150-Vamp_d*Rint_volt_d)+Vamp_d*Rint_volt_d*R150*(-Vtot_d+Vamp_d)/((Vtot_d*Rint_volt_d-Vamp_d*R150-Vamp_d*Rint_volt_d)*(Vtot_d*Rint_volt_d-Vamp_d*R150-Vamp_d*Rint_volt_d));
    double diffR_d    = Vamp_d*Rint_volt_d/(Vtot_d*Rint_volt_d-Vamp_d*R150-Vamp_d*Rint_volt_d)+Vamp_d*Rint_volt_d*R150*Vamp_d/((Vtot_d*Rint_volt_d-Vamp_d*R150-Vamp_d*Rint_volt_d)*(Vtot_d*Rint_volt_d-Vamp_d*R150-Vamp_d*Rint_volt_d));
    double diffVtot_d = -Vamp_d*Rint_volt_d*R150*Rint_volt_d/((Vtot_d*Rint_volt_d-Vamp_d*R150-Vamp_d*Rint_volt_d)*(Vtot_d*Rint_volt_d-Vamp_d*R150-Vamp_d*Rint_volt_d));
    
    double Rint_amp_d_err_part = sqrt( 
        diffVamp_d * diffVamp_d * Vamp_err_d * Vamp_err_d +
        diffRint_d * diffRint_d * Rint_volt_err_d * Rint_volt_err_d +
        diffR_d * diffR_d * R150_err * R150_err +
        diffVtot_d * diffVtot_d * Vtot_err_d * Vtot_err_d
    );
    printout << "\nResistenza amperometro digitale con metodo partitore = \t(" <<  Rint_amp_d_part << "\t +/- " << Rint_amp_d_err_part << ")\tΩ";

    // Resistenza interna amperometro digitale (metodo diretto)
    double Rint_amp_d_dir = 1. / (I_d / Vamp_d - 1. / Rint_volt_d);

    double diffI_d     = -1./((I_d/Vamp_d-1./Rint_volt_d)*(I_d/Vamp_d-1./Rint_volt_d)*Vamp_d);
    double diffVamp_d_dir = I_d/((I_d/Vamp_d-1./Rint_volt_d)*(I_d/Vamp_d-1./Rint_volt_d)*Vamp_d*Vamp_d);
    double diffRintvd_dir  = -1./((I_d/Vamp_d-1/Rint_volt_d)*(I_d/Vamp_d-1./Rint_volt_d)*Rint_volt_d*Rint_volt_d);

    double Rint_amp_d_err_dir = sqrt(
        diffI_d * diffI_d * I_err_d * I_err_d +
        diffVamp_d_dir * diffVamp_d_dir * Vamp_err_d * Vamp_err_d +
        diffRintvd_dir * diffRintvd_dir * Rint_volt_err_d * Rint_volt_err_d
    );

    printout << "\nResistenza amperometro digitale con metodo diretto = \t(" <<  Rint_amp_d_dir << "\t +/- " << Rint_amp_d_err_dir << ")\tΩ";

    //Potenza dissipata dai voltmetri
    //-----------------------------------------
    //P = P(Va,Vb,R470K) 
    double P_a_Va = 1000* Vag_a * Vag_a / (R470K * Vbg_a / (Vag_a - Vbg_a)); //mW
    double P_a_Vb = 1000* Vbg_a * Vbg_a / (R470K * Vbg_a / (Vag_a - Vbg_a)); //mW
    double P_d_Va = 1000* Vag_d * Vag_d / (R470K * Vbg_d / (Vag_d - Vbg_d)); //mW
    double P_d_Vb = 1000* Vbg_d * Vbg_d / (R470K * Vbg_d / (Vag_d - Vbg_d)); //mW

    double P_a_Va_err = sqrt(
        (3*Vag_a*Vag_a-2*Vag_a*Vbg_a)/(R470K*Vbg_a) * (3*Vag_a*Vag_a-2*Vag_a*Vbg_a)/(R470K*Vbg_a) * Vag_err_a * Vag_err_a +
        (Vag_a*Vag_a*Vag_a/(R470K*Vbg_a*Vbg_a)) * (Vag_a*Vag_a*Vag_a/(R470K*Vbg_a*Vbg_a)) * Vbg_err_a * Vbg_err_a +
        ((Vag_a*Vag_a*Vag_a+Vag_a*Vag_a*Vbg_a)/(R470K*R470K*Vbg_a)) * ((Vag_a*Vag_a*Vag_a+Vag_a*Vag_a*Vbg_a)/(R470K*R470K*Vbg_a)) * R470K_err * R470K_err
    ) * 1000; //mW

    printout << "\n\n\nPotenza dissipata dai voltmetri durante la misurazione:";
    printout << "\nVolmetro analogico, Va, P = (" << P_a_Va << " +/- " << P_a_Va_err << ") mW";

    double P_a_Vb_err = sqrt(
        (Vbg_a/R470K)*(Vbg_a/R470K) * Vag_err_a * Vag_err_a +
        ((Vag_a-2*Vbg_a)/R470K) * ((Vag_a-2*Vbg_a)/R470K) * Vbg_err_a * Vbg_err_a +
        ((-Vag_a*Vbg_a+Vbg_a*Vbg_a)/(R470K*R470K))*((-Vag_a*Vbg_a+Vbg_a*Vbg_a)/(R470K*R470K))* R470K_err * R470K_err
    ) * 1000; //mW

    printout << "\nVolmetro analogico, Vb, P = (" << P_a_Vb << " +/- " << P_a_Vb_err << ") mW";

    //-----------------------------

    double P_d_Va_err = sqrt(
        (3*Vag_d*Vag_d-2*Vag_d*Vbg_d)/(R470K*Vbg_d) * (3*Vag_d*Vag_d-2*Vag_d*Vbg_d)/(R470K*Vbg_d) * Vag_err_d * Vag_err_d +
        (Vag_d*Vag_d*Vag_d/(R470K*Vbg_d*Vbg_d)) * (Vag_d*Vag_d*Vag_d/(R470K*Vbg_d*Vbg_d)) * Vbg_err_d * Vbg_err_d +
        ((Vag_d*Vag_d*Vag_d+Vag_d*Vag_d*Vbg_d)/(R470K*R470K*Vbg_d)) * ((Vag_d*Vag_d*Vag_d+Vag_d*Vag_d*Vbg_d)/(R470K*R470K*Vbg_d)) * R470K_err * R470K_err
    ) * 1000; //mW

    printout << "\nVolmetro digitale, Va, P = (" << P_d_Va << " +/- " << P_d_Va_err << ")  mW";

    double P_d_Vb_err = sqrt(
        (Vbg_d/R470K)*(Vbg_d/R470K) * Vag_err_d * Vag_err_d +
        ((Vag_d-2*Vbg_d)/R470K) * ((Vag_d-2*Vbg_d)/R470K) * Vbg_err_d * Vbg_err_d +
        ((-Vag_d*Vbg_d+Vbg_d*Vbg_d)/(R470K*R470K))*((-Vag_d*Vbg_d+Vbg_d*Vbg_d)/(R470K*R470K))* R470K_err * R470K_err
    ) * 1000; //mW

    printout << "\nVolmetro digitale, Vb, P = (" << P_d_Vb << " +/- " << P_d_Vb_err << ")  mW";

    //------- Relazioni lampadina ------------

    //Array per potenza e resistenza
    std::vector<double> P_R;        //mW
    std::vector<double> R_V;        //KΩ
    std::vector<double> P_R_err;    //mW
    std::vector<double> R_V_err;    //KΩ

    //Calcolo potenza e resistenza
    printout << "\n\nCalcolo potenza (K) e resistenza (m)\n Potenza (m) \t\t Errore \t\t Resistenza (K) \t\t Errore";
    for (int i = 0; i < Vlamp_err.size() ; i++){
        P_R.push_back(Vlamp[i]*Ilamp[i]);
        R_V.push_back(Vlamp[i]/Ilamp[i]);
        P_R_err.push_back(sqrt(Vlamp[i]*Vlamp[i]*Ilamp_err[i]*Ilamp_err[i]+Ilamp[i]*Ilamp[i]*Vlamp_err[i]*Vlamp_err[i]));
        R_V_err.push_back(sqrt(Vlamp_err[i]*Vlamp_err[i]/(Ilamp[i]*Ilamp[i]) + (Vlamp[i]*Vlamp[i]/(Ilamp[i]*Ilamp[i]*Ilamp[i]*Ilamp[i])) * Ilamp_err[i]*Ilamp_err[i]));
    }
    for (int i = 0; i < Vlamp_err.size() ; i++){
        printout << "\n" << P_R[i] << "\t\t" << P_R_err[i] << "\t\t" << R_V[i] << "\t\t" << R_V_err[i];
    }

    TGraphErrors * currVsVolt = new TGraphErrors(Ilamp.size());
    TGraphErrors * resVsVolt  = new TGraphErrors(Ilamp.size());
    TGraphErrors * potVsRes   = new TGraphErrors(Ilamp.size());

    for (int i = 1; i < Vlamp.size(); ++i)
    {
        currVsVolt->SetPoint(i, Vlamp[i], Ilamp[i]);
        currVsVolt->SetPointError(i, Vlamp_err[i], Ilamp_err[i]);
        
        resVsVolt->SetPoint(i, Vlamp[i], R_V[i]);
        resVsVolt->SetPointError(i, Vlamp_err[i], R_V_err[i]);

        potVsRes->SetPoint(i, R_V[i], P_R[i]);
        potVsRes->SetPointError(i, R_V_err[i], P_R_err[i]);
    }

    //Titoli degli assi, griglie e impostazioni grafiche
    currVsVolt->SetTitle("I(V) filamento di tungsteno");
    resVsVolt ->SetTitle("R(V) filamento di tungsteno");
    potVsRes  ->SetTitle("P(R) filamento di tungsteno");
    
    currVsVolt->GetXaxis()->SetTitle("Tensione [V]");
    currVsVolt->GetYaxis()->SetTitle("Corrente [mA]");
    currVsVolt->SetMarkerSize(2);
    currVsVolt->SetMarkerStyle('o');

    resVsVolt->GetXaxis()->SetTitle("Tensione [V]");
    resVsVolt->GetYaxis()->SetTitle("Resistenza [KΩ]");
    resVsVolt->SetMarkerSize(2);
    resVsVolt->SetMarkerStyle('o');

    potVsRes->GetXaxis()->SetTitle("Resistenza [KΩ]");
    potVsRes->GetYaxis()->SetTitle("Potenza [mW]");
    resVsVolt->SetMarkerSize(2);
    resVsVolt->SetMarkerStyle('o');

    printout << setprecision(4);
    //__________________________________________________________________________________________________________
    //Regressione con la (4) della relazione P(R) - Modello mr^q
    //__________________________________________________________________________________________________________
    TF1 * f4 = new TF1("f4", "[0]*x^[1]", 0., 0.5);
    f4->SetParameter(0, 10000.);
    f4->SetParameter(1, 4.);
    TVirtualFitter::SetMaxIterations(10000);
    potVsRes->Fit(f4);
    printout << "\n\n\n\nRegressione con la (4) della relazione P(R)";
    printout << "\n---------------------------------------------"; 
    printout << "\nFit chi square = " << f4->GetChisquare() << "  ndof = " << f4->GetNDF() << "   chi/ndf = " << f4->GetChisquare() / f4->GetNDF(); 
    printout << "\np-value = " << f4->GetProb();
    printout << "\nm = \t" << f4->GetParameter(0) << "\t +/- " << f4->GetParError(0);
    printout << "\nq = \t" << f4->GetParameter(1) << "\t +/- " << f4->GetParError(1);
    potVsRes->Write("P(R) con f4 - Modello mr^q"); //scrivere il fit sul file

    //__________________________________________________________________________________________________________
    //1A.a //Regressione con la (4) con gradi di libertà aggiuntivi - modello completo
    //__________________________________________________________________________________________________________

    TF1 * f11 = new TF1("f11", "[0]*x^[1]+[2]*x^([1]/4)+[3]", 0.00, 0.5);
    f11->SetParameter(0, 10000.);
    f11->SetParameter(1, 4.);
    TVirtualFitter::SetMaxIterations(10000);
    potVsRes->Fit(f11, "R");
    printout << "\n\n\n\nRegressione con la (4) della relazione P(R) con dof aggiuntivi m*x^q+b*x^(q/4)+c";
    printout << "\n---------------------------------------------"; 
    printout << "\nFit chi square = " << f11->GetChisquare() << "  ndof = " << f11->GetNDF() << "   chi/ndf = " << f11->GetChisquare() / f11->GetNDF(); 
    printout << "\np-value = " << f11->GetProb();
    printout << "\nm = \t" << f11->GetParameter(0) << "\t +/- " << f11->GetParError(0);
    printout << "\nq = \t" << f11->GetParameter(1) << "\t +/- " << f11->GetParError(1);
    printout << "\nb = \t" << f11->GetParameter(2) << "\t +/- " << f11->GetParError(2);
    printout << "\nc = \t" << f11->GetParameter(3) << "\t +/- " << f11->GetParError(3);
    potVsRes->Write("P(R) con f11 - Modello completo"); //scrivere il fit sul file

    //__________________________________________________________________________________________________________
        //Discussione del range di validità del modello semplificato mr^q, poiché
        // non essendo trascendente è utile perché è facilmente invertibile.
    //__________________________________________________________________________________________________________
    
    /*printout << "\n\nStudiamo l'intervallo di validità del modello semplificato mr^q";
    for (int i = 0; i < 8; ++i)
    {
        TF1 * f4r = new TF1("f4r", "[0]*x^[1]", (0.005 + R_V[i]), 0.5);
        f4r->SetParameter(0, 10000.);
        f4r->SetParameter(1, 4.);
        TVirtualFitter::SetMaxIterations(10000);
        potVsRes->Fit(f4r, "R");
        printout << "\n\nRestrizione togliendo " << i << " punti. (" << (0.005 + R_V[i]) << ")" ;
        printout << "\nFit chi square = " << f4r->GetChisquare() << "  ndof = " << f4r->GetNDF() << "   chi/ndf = " << f4r->GetChisquare() / f4r->GetNDF(); 
    }*/

    const int firstV = 1;
    printout << "\n\nIntervallo ristretto: [" << firstV << ";end]";

    //__________________________________________________________________________________________________________
    //Regressione con la (4) di P(R) - Modello semplificato mr^q con intervallo ristretto
    //__________________________________________________________________________________________________________
    TF1 * f4s = new TF1("f4s", "[0]*x^[1]", 0., 0.5);
    f4s->SetParameter(0, 10000.);
    f4s->SetParameter(1, 4.);
    TVirtualFitter::SetMaxIterations(10000);
    potVsRes->Fit(f4s);
    printout << "\n\n\n\nRegressione con la (4) di P(R) - Modello semplificato mr^q con intervallo ristretto";
    printout << "\n---------------------------------------------"; 
    printout << "\nFit chi square = " << f4s->GetChisquare() << "  ndof = " << f4s->GetNDF() << "   chi/ndf = " << f4s->GetChisquare() / f4s->GetNDF(); 
    printout << "\np-value = " << f4s->GetProb();
    printout << "\nm = \t" << f4s->GetParameter(0) << "\t +/- " << f4s->GetParError(0);
    printout << "\nq = \t" << f4s->GetParameter(1) << "\t +/- " << f4s->GetParError(1);
    potVsRes->Write("P(R) con f4 RISTRETTO - Modello mr^q"); //scrivere il fit sul file
    
    //__________________________________________________________________________________________________________
    //Regressione P(R) modello mr^q + termine convettivo 
    //__________________________________________________________________________________________________________
    TF1 * f12 = new TF1("f12", "[0]*x^[1]+[2]*x^([1]/4)", 0, 0.5);
    f12->SetParameter(0, 10000.);
    f12->SetParameter(1, 4.);
    TVirtualFitter::SetMaxIterations(10000);
    potVsRes->Fit(f12);
    printout << "\n\n\n\nRegressione P(R) modello mr^q + termine convettivo ";
    printout << "\n---------------------------------------------"; 
    printout << "\nFit chi square = " << f12->GetChisquare() << "  ndof = " << f12->GetNDF() << "   chi/ndf = " << f12->GetChisquare() / f12->GetNDF(); 
    printout << "\np-value = " << f12->GetProb();
    printout << "\nm = \t" << f12->GetParameter(0) << "\t +/- " << f12->GetParError(0);
    printout << "\nq = \t" << f12->GetParameter(1) << "\t +/- " << f12->GetParError(1);
    printout << "\nb = \t" << f12->GetParameter(2) << "\t +/- " << f12->GetParError(2);
    potVsRes->Write("Regressione P(R) modello mr^q + termine convettivo");

    //__________________________________________________________________________________________________________
    //Regressione P(R) modello mr^q + termine a freddo, costante
    //__________________________________________________________________________________________________________
    TF1 * f13 = new TF1("f13", "[0]*x^[1]+[2]", 0, 0.5);
    f13->SetParameter(0, 10000.);
    f13->SetParameter(1, 4.);
    TVirtualFitter::SetMaxIterations(10000);
    potVsRes->Fit(f13);
    printout << "\n\n\n\nRegressione P(R) modello mr^q + termine a freddo, costante";
    printout << "\n---------------------------------------------"; 
    printout << "\nFit chi square = " << f13->GetChisquare() << "  ndof = " << f13->GetNDF() << "   chi/ndf = " << f13->GetChisquare() / f13->GetNDF(); 
    printout << "\np-value = " << f13->GetProb();
    printout << "\nm = \t" << f13->GetParameter(0) << "\t +/- " << f13->GetParError(0);
    printout << "\nq = \t" << f13->GetParameter(1) << "\t +/- " << f13->GetParError(1);
    printout << "\nc = \t" << f13->GetParameter(2) << "\t +/- " << f13->GetParError(2);
    potVsRes->Write("Regressione P(R) modello mr^q + termine a freddo, costante");

    //_________________________________________MODELLI VINCOLATI CON q=4________________________________________

    //__________________________________________________________________________________________________________
    //Modello mr^q con q=4 (3A)
    //__________________________________________________________________________________________________________
    TF1 * f14 = new TF1("f14", "[0]*x^4", 0, 0.5);
    f14->SetParameter(0, 10000.);
    f14->SetParameter(1, 4.);
    TVirtualFitter::SetMaxIterations(10000);
    potVsRes->Fit(f14);
    printout << "\n\n\n\nModello mr^q con q=4 (3A)";
    printout << "\n---------------------------------------------"; 
    printout << "\nFit chi square = " << f14->GetChisquare() << "  ndof = " << f14->GetNDF() << "   chi/ndf = " << f14->GetChisquare() / f14->GetNDF(); 
    printout << "\np-value = " << f14->GetProb();
    printout << "\nm = \t" << f14->GetParameter(0) << "\t +/- " << f14->GetParError(0);
    potVsRes->Write("Modello mr^q con q=4 (3A)");


    //__________________________________________________________________________________________________________
    //Modello mr^q con q=4 + termine a freddo, costante (3B)
    //__________________________________________________________________________________________________________
    TF1 * f15 = new TF1("f15", "[0]*x^4+[1]", 0, 0.5);
    f15->SetParameter(0, 10000.);
    f15->SetParameter(1, 4.);
    TVirtualFitter::SetMaxIterations(10000);
    potVsRes->Fit(f15);
    printout << "\n\n\n\nModello mr^q con q=4 + termine a freddo, costante (3B)";
    printout << "\n---------------------------------------------"; 
    printout << "\nFit chi square = " << f15->GetChisquare() << "  ndof = " << f15->GetNDF() << "   chi/ndf = " << f15->GetChisquare() / f15->GetNDF(); 
    printout << "\np-value = " << f15->GetProb();
    printout << "\nm = \t" << f15->GetParameter(0) << "\t +/- " << f15->GetParError(0);
    printout << "\nc = \t" << f15->GetParameter(1) << "\t +/- " << f15->GetParError(1);
    potVsRes->Write("Modello mr^q con q=4 + termine a freddo, costante (3B)");


    //__________________________________________________________________________________________________________
    //Modello mr^q con q^4 + termine costante + termine convettivo (4)
    //__________________________________________________________________________________________________________
    TF1 * f16 = new TF1("f16", "[0]*x^4+[1]*x+[2]", 0, 0.5);
    f16->SetParameter(0, 10000.);
    f16->SetParameter(1, 4.);
    TVirtualFitter::SetMaxIterations(10000);
    potVsRes->Fit(f16);
    printout << "\n\n\n\nModello mr^q con q^4 + termine costante + termine convettivo (4)";
    printout << "\n---------------------------------------------";
    printout << "\nFit chi square = " << f16->GetChisquare() << "  ndof = " << f16->GetNDF() << "   chi/ndf = " << f16->GetChisquare() / f16->GetNDF(); 
    printout << "\np-value = " << f16->GetProb();
    printout << "\nm = \t" << f16->GetParameter(0) << "\t +/- " << f16->GetParError(0);
    printout << "\nb = \t" << f16->GetParameter(1) << "\t +/- " << f16->GetParError(1);
    printout << "\nc = \t" << f16->GetParameter(2) << "\t +/- " << f16->GetParError(2);
    potVsRes->Write("Modello mr^q con q^4 + termine costante + termine convettivo (4)");


    //________________________________ALTRI_PLOT_I(V)_R(V)_Pol4_________________________________________________

    //__________________________________________________________________________________________________________
    //1B //Regressione con la (5) della relazione I(V)
    //__________________________________________________________________________________________________________
    TF1 * f5 = new TF1("f5", "([0]^(1./([1]+1)))*x^(([1]-1)/([1]+1))", (Vlamp[firstV]+0.1), 12.5);
    f5->SetParameter(0, 10000.);
    f5->SetParameter(1, 3.5);
    TVirtualFitter::SetMaxIterations(10000);
    currVsVolt->Fit(f5, "R");
    printout << "\n\n\n\nRegressione con la (5) della relazione I(V)";
    printout << "\n---------------------------------------------"; 
    printout << "\nFit chi square = " << f5->GetChisquare() << "  ndof = " << f5->GetNDF() << "   chi/ndf = " << f5->GetChisquare() / f5->GetNDF(); 
    printout << "\np-value = " << f5->GetProb();
    printout << "\nm = \t" << f5->GetParameter(0) << "\t +/- " << f5->GetParError(0);
    printout << "\nq = \t" << f5->GetParameter(1) << "\t +/- " << f5->GetParError(1);
    currVsVolt->Write("I(V) con f5"); //scrivere il fit sul file

    //__________________________________________________________________________________________________________
    //1C //Regressione relazione R(V)
    //__________________________________________________________________________________________________________
    TF1 * f6 = new TF1("f6", "((([0]^(2.0/([1]+1))) * x^(2*([1]-1)/([1]+1)))/[0])^(1.0/([1]-1))", (Vlamp[firstV]+0.1), 12.);
    f6->SetParameter(0, 14000);
    f6->SetParameter(1, 3.5);
    resVsVolt->Fit(f6, "R");
    printout << "\n\n\n\n1C Regressione relazione R(V) - Ricavato dal modello semplificato mr^q";
    printout << "\n---------------------------------------------"; 
    printout << "\nFit chi square = " << f6->GetChisquare() << "  ndof = " << f6->GetNDF() << "   chi/ndf = " << f6->GetChisquare() / f6->GetNDF(); 
    printout << "\np-value = " << f6->GetProb();
    printout << "\nm = \t" << f6->GetParameter(0) << "\t +/- " << f6->GetParError(0);
    printout << "\nq = \t" << f6->GetParameter(1) << "\t +/- " << f6->GetParError(1);
    resVsVolt->Write("R(V) con f6 - Ricavato dal modello semplificato mr^q"); //scrivere il fit sul file

    //__________________________________________________________________________________________________________
    //2A //Regressione polinomiale di grado 4 della P(R)
    //__________________________________________________________________________________________________________
    TF1 * f7 = new TF1("f7", "pol4", 0, 0.5);
    TVirtualFitter::SetMaxIterations(10000);
    potVsRes->Fit(f7);
    //potVsRes->Draw();
    printout << "\n\n\n\n2A Regressione polinomiale di grado 4 della P(R) a+bx+cx^2+dx^3+ex^4";
    printout << "\n---------------------------------------------"; 
    printout << "\nFit chi square = " << f7->GetChisquare() << "  ndof = " << f7->GetNDF() << "   chi/ndf = " << f7->GetChisquare() / f7->GetNDF(); 
    printout << "\np-value = " << f7->GetProb();
    printout << "\na = \t" << f7->GetParameter(0) << "\t +/- " << f7->GetParError(0);
    printout << "\nb = \t" << f7->GetParameter(1) << "\t +/- " << f7->GetParError(1);
    printout << "\nc = \t" << f7->GetParameter(2) << "\t +/- " << f7->GetParError(2);
    printout << "\nd = \t" << f7->GetParameter(3) << "\t +/- " << f7->GetParError(3);
    printout << "\ne = \t" << f7->GetParameter(4) << "\t +/- " << f7->GetParError(4);
    potVsRes->Write("P(R) con f7"); //scrivere il fit sul file


    //________________________________Canvas con fit di P(R)____________________________________________________

    const int nmisure = Vlamp.size();
    TCanvas * cPofR = new TCanvas("cPofR","P(R)",10,3,600,400);
    cPofR->SetFillColor(0);
    cPofR->cd();
    potVsRes->SetMarkerSize(0.6);
    potVsRes->SetMarkerStyle(21);
    potVsRes->SetTitle("P(R)");
    potVsRes->GetXaxis()->SetTitle("R [k#Omega]");
    potVsRes->GetYaxis()->SetTitle("P [mW]");
    potVsRes->Draw("AP");
        f4->SetLineColor(kRed+1);                                          //fit Rosso con mr^q
            f4->Draw("same");
        f4s->SetLineColor(kOrange+8);                                      //fit Arancio con mr^q ristretta 
            f4s->Draw("same");
        f11->SetLineColor(kGreen+4);                                       //fit Verde Scuro con modello completo 
            f11->Draw("same");
        f12->SetLineColor(kGreen+2);                                       //fit Verde medio con mr^q più convettivo 
            f12->Draw("same");
        f13->SetLineColor(kGreen-9);                                       //fit Verde chiaro mr^q più costante 
            f13->Draw("same");
        f16->SetLineColor(kBlue+2);                                        //fit Blu medio con mr^4 completo 
            f16->Draw("same");
        f15->SetLineColor(kAzure-4);                                       //fit Azzurro medio con mr^4 più costante 
            f15->Draw("same");
        f14->SetLineColor(kAzure+8);                                       //fit Azzurro chiaro con mr^4 
            f14->Draw("same");
        f7->SetLineColor(kBlack);                                          //fit nero chiaro con polinomiale
            f7->Draw("same");
    //creo una legenda
    TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("Legenda","C"); 				// option "C" allows to center the header
    legend->AddEntry(f4, "(1) mr^q","l");
    legend->AddEntry(f4s,"(2) mr^q ristretto","l");
    legend->AddEntry(f11,"(3) modello completo","l");
    legend->AddEntry(f12,"(4) mr^q con termine convettivo","l");
    legend->AddEntry(f13,"(5) mr^q con termine a freddo","l");
    legend->AddEntry(f16,"(6) mr^4 completo","l");
    legend->AddEntry(f15,"(7) mr^4 con termine a freddo","l");
    legend->AddEntry(f14,"(8) mr^4","l");
    legend->AddEntry(f7, "(9) polinomio quarto grado","l");
    legend->Draw();

    cPofR->Write("Grafico P(R) con fit");

    //__________________________________________________________________________________________________________
    //-------------- TEST STATISTICI ---------------------------------------------------------------------------
    //__________________________________________________________________________________________________________

    //Sommario:
    // - Fisher tra modello completo e mr^q + convettivo
    // - Fisher tra modello completo e mr^q + costante
    // - Fisher tra modello completo e modello completo ma con q=4
    // - Gauss di q tra P(R) e I(V) ristretto 
    // - Gauss di q tra R(V) e I(V) ristretti
    // - Gauss di q tra P(R) e R(V) ristretti
    // - Gauss di compatibilità tra q di P(R), R(V), I(R) e 4 
    // - Gauss di compatibilità tra q di P(R) completo e 4 


    // - Test di Fisher tra modello completo e mr^q + convettivo
    //-------------------------------------------------------------------------------
    double fvar_a = (f11->GetChisquare()/f11->GetNDF())/(f12->GetChisquare()/f12->GetNDF());
    double fpvalue_a = ROOT::Math::fdistribution_cdf(fvar_a,  f11->GetNDF(), f12->GetNDF());
    printout << "\n\n\n+-------------------------------------------------------------+";
    printout << "\nTest di Fisher tra modello completo e mr^q + convettivo";
    printout << "\n" << "Fisher var (chi2/ndf  /  chi2/ndf)   f = " << fvar_a;
    printout << "\n" << "p-value = " << fpvalue_a;

    // - Fisher tra modello completo e mr^q + costante
    //-------------------------------------------------------------------------------
    double fvar_b = (f11->GetChisquare()/f11->GetNDF())/(f13->GetChisquare()/f13->GetNDF());
    double fpvalue_b = ROOT::Math::fdistribution_cdf(fvar_b,  f11->GetNDF(), f13->GetNDF());
    printout << "\n\n\n+-------------------------------------------------------------+";
    printout << "\nFisher tra modello completo e mr^q + costante";
    printout << "\n" << "Fisher var (chi2/ndf  /  chi2/ndf)   f = " << fvar_b;
    printout << "\n" << "p-value = " << fpvalue_b;

    // - Fisher tra modello completo e modello completo ma con q=4
    //-------------------------------------------------------------------------------
    double fvar_c = (f11->GetChisquare()/f11->GetNDF())/(f16->GetChisquare()/f16->GetNDF());
    double fpvalue_c = ROOT::Math::fdistribution_cdf(fvar_c,  f11->GetNDF(), f16->GetNDF());
    printout << "\n\n\n+-------------------------------------------------------------+";
    printout << "\nFisher tra modello completo e modello completo ma con q=4";
    printout << "\n" << "Fisher var (chi2/ndf  /  chi2/ndf)   f = " << fvar_c;
    printout << "\n" << "p-value = " << fpvalue_c;

    // - TEST DI GAUSS LAMP - ATTENZIONE!! Sono problematici perché le grandezze confrontate non sono indipendenti! 
    //                        Quindi non si mettono nella relazione.
    /*
    printout << "\n\n+-------------------------------------------------------------+";

    printout << "\n\nTest Z di q tra P(R) e I(V) ristretto";
    quickZtwotailed(f4s->GetParameter(1), f5->GetParameter(1), f4s->GetParError(1), f5->GetParError(1), printout);

    printout << "\n\nGauss di q tra R(V) e I(V) ristretti";
    quickZtwotailed(f6->GetParameter(1),  f5->GetParameter(1), f6->GetParError(1), f5->GetParError(1), printout);

    printout << "\n\nGauss di q tra R(V) e P(R) ristretti";
    quickZtwotailed(f6->GetParameter(1), f4s->GetParameter(1), f6->GetParError(1), f4s->GetParError(1), printout);

    printout << "\n\nGauss di q tra R(V) 4";
    quickZtwotailed(f6->GetParameter(1), 4.,  f6->GetParError(1), 0., printout);

    printout << "\n\nGauss di q tra I(V) 4";
    quickZtwotailed(f5->GetParameter(1), 4.,  f5->GetParError(1), 0., printout);

    printout << "\n\nGauss di q tra P(R) 4";
    quickZtwotailed(f4s->GetParameter(1), 4., f4s->GetParError(1), 0., printout);

    printout << "\n\nGauss di q tra P(R) completo e 4";
    quickZtwotailed(f11->GetParameter(1), 4., f11->GetParError(1), 0., printout);
    */

    // - TEST DI GAUSS CARATTERIZZAZIONE STRUMENTI
    printout << "\n\n+-------------------------------------------------------------+";

    printout << "\n\nTest Z Voltmetro analogico - Misure <=> Costruttore";
    quickZtwotailed(Rint_volt_a, Rint_volt_a_prod, Rint_volt_err_a, Rint_volt_a_prod_err, printout);

    printout << "\n\nTest Z Voltmetro digitale - Misure <=> Costruttore";
    quickZtwotailed(Rint_volt_d, Rint_volt_d_prod, Rint_volt_err_d, Rint_volt_d_prod_err, printout);

    printout << "\n\nTest Z Amperometro analogico - Metodo diretto <=> Costruttore";
    quickZtwotailed(Rint_amp_a_dir, Rint_amp_a_prod, Rint_amp_a_err_dir, Rint_amp_a_prod_err, printout);

    printout << "\n\nTest Z Amperometro analogico - Metodo partitore <=> Costruttore";
    quickZtwotailed(Rint_amp_a_part, Rint_amp_a_prod, Rint_amp_a_err_part, Rint_amp_a_prod_err, printout);

    printout << "\n\nTest Z Amperometro digitale - Metodo diretto <=> Costruttore";
    quickZtwotailed(Rint_amp_d_dir, Rint_amp_d_prod, Rint_amp_d_err_dir, Rint_amp_d_prod_err, printout);

    printout << "\n\nTest Z Amperometro digitale - Metodo partitore <=> Costruttore";
    quickZtwotailed(Rint_amp_d_part, Rint_amp_d_prod, Rint_amp_d_err_part, Rint_amp_d_prod_err, printout);    

    printout << "\n\n+-------------------------------------------------------------+";

    //Scrittura oggetti rimanenti TGraphErrors su file ROOT
    currVsVolt->Write("Lamp I(V)");
    resVsVolt->Write("Lamp R(V)");
    potVsRes->Write("Lamp P(R)");

    //------------- SCRITTURA DEL TEMPO MACCHINA E CHIUSURA DEI FILES DI OUTPUT ----------------

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

    Fare attenzione a non dimenticare alcun dettaglio richiesto nelle schede di laboratorio.
    Riportare con incertezze le resistenze interne e i valori di Vag e Vbg.
    Calcolare la POTENZA DISSIPATA DAI MULTIMETRI.

    La tabella dei valori è superflua, al massimo si può mettere in appendice, non è richiesta nel corpo, dove invece
    fa fede il grafico.

    Sovrapporre sullo stesso grafico tutti i diversi tentativi di regressione effettuati con diverse curve (diversi colori) e LEGENDA.

    Per quasi tutti quelli che hanno i due multimetri digitali, facilmente servirà il modello completo di P(R) per riuscire a fittare correttamente,
    perché gli errori sono piccoli. Si è utilizzata l'approssimazione per poter avere una funzione invertibile da mettere in P(V) e R(V).
    In questo caso facilmente il modello semplificato invertito non funzionerà su tutti i dati, ma bisognerà (A PARTIRE DA P(R)) rimuovere 
    i primi punti fino a quando la temperatura diventa abbastanza elevata da consentire l'appossimazione di corpo nero.

    Riportare le funzioni, in particolare la R(V) e le altre ricavate in autonomia. Riportare i grafici I(V) e R(V), oltre a P(R).

    I valori di q si ricavano da numerosi fit, ATTENZIONE: i valori di q ottenuti non sono indipendenti tra loro, perché sono fortemente correlati!!
    Non si possono fare test Z! Di conseguenza non si possono confrontare tra di loro, ma 
*/