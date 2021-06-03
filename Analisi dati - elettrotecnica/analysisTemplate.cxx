#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>    
//#include "statlib.cxx"
#include "csvdata.cxx"
using namespace std;

//EXPERIMENTAL DATA
std::vector<long int> data1;
std::vector<double> data2;
std::vector<double> data3;

//ANALYSIS OBJECTS

//TGraphError * h = new TGraphError();



void parseData();
void analysis1();
void analysis2();
void analysis3();
void exportValues();
void plotResults();

void analysisTemplate()
{
    auto start = std::chrono::system_clock::now();
    cerr << endl;
    cerr << "  +-------------------------------------------------+ \n";
    cerr << "  |   V. Pagliarino - Corso B - C.d.L. in Fisica    | \n";
    cerr << "  |    Università degli Studi di Torino -  ESP2     | \n";
    cerr << "  +-------------------------------------------------+ \n";
    cerr << "  | Experiment:  ********************************** | \n";
    cerr << "  | Date:        **/**/****                         | \n";
    cerr << "  | Revision:    ***.***.***.***                    | \n";
    cerr << "  | Description: *********************************  | \n";
    cerr << "  +-------------------------------------------------+ \n\n";
    parseData();
    analysis1();
    analysis2();
    analysis3();
    exportValues();
    plotResults();
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    std::cout << "\nFinished computation at " << std::ctime(&end_time)
          << "elapsed time: " << elapsed_seconds.count() << "s\n";
}

void parseData()
{
    csvdata csv("/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/Fisica - Esperimentazioni 2/dataex.csv");
    csv.setDelimiters(';');

    for (int i = 1; i < 10; ++i)
    {
        data2.push_back(csv.getDouble(i, 0));
        data3.push_back(csv.getDouble(i, 1));
        cout << "\n" <<  data2[data2.size() - 1];
    }

    cerr << "\nParsing completed.";

    
    
    csv.flush();
}

void analysis1()
{

}

void analysis2()
{

}

void analysis3()
{

}

void exportValues()
{
    cerr << "\nWriting analysis results to CSV+ROOT files...";
    system("touch \"/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/Fisica - Esperimentazioni 2/dataexOutput.txt\"");
    fstream printout;
    printout.open("/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/Fisica - Esperimentazioni 2/dataexOutput.csv");
    TFile out_file("/Users/ValerioPagliarino/Desktop/Esperimentazioni II/Workspace/Fisica - Esperimentazioni 2/dataexOutGraphics.root", "RECREATE");
//.....................................
    

    //h->Write();

    printout << "name;value";






//.....................................
    // Close the file
    out_file.Close(); 
    printout.close();
    cerr << "\nFiles closed successfully.";
}

void plotResults()
{
    TBrowser * tb = new TBrowser();
}