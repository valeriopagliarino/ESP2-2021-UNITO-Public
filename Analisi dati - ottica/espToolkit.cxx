#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
//Include CERN ROOT LIBRARIES etc.

//Prototipi static functions
double getDB            (double Vin,    double Vout);
double getDBerr         (double Vin,    double Vout,    double Vin_err, double Vout_err);
double getDBPower       (double Pin,    double Pout);
double getDBerrPower    (double Pin,    double Pout,    double Pin_err, double Pout_err);
double media_n_variabili(std::vector<double> v);


//Calcolo del guadagno in uscita da un quadripolo date le tensioni in ingresso e in uscita
double getDB(double Vin, double Vout)
{
    return 20. * log10(Vout / Vin);
}

double getDBerr(double Vin, double Vout, double Vin_err, double Vout_err)
{
    return sqrt(
        (20. / (log(10) * Vout)) * (20. / (log(10) * Vout)) * Vout_err * Vout_err + 
        (-20. / (log(10) * Vin)) * (-20. / (log(10) * Vin))  * Vin_err  * Vin_err );
}

//Calcolo del guadagno in uscita da un quadripolo date le potenze in ingresso e in uscita
double getDBPower(double Pin, double Pout)
{
    return 10. * log10(Pout / Pin);
}

double getDBerrPower(double Pin, double Pout, double Pin_err, double Pout_err)
{
    return sqrt(
        (10. / (log(10) * Pout)) * (10. / (log(10) * Pout)) * Pout_err * Pout_err + 
        (-10. / (log(10) * Pin)) * (-10. / (log(10) * Pin))  * Pin_err  * Pin_err );
}

double media_n_variabili(std::vector<double> v)
{
    double sum = 0.;
    for (int i = 0; i < v.size(); ++i)
    {
        sum = sum + v[i];
    }
    sum = sum / v.size();
    //cout << "\n\nMedia di " << v.size() << " elementi\n\n";
    return sum;
}

double media_n_variabili_err(std::vector<double> v, std::vector<double> v_err)
{
    // Please see https://drive.google.com/file/d/1JkzlBUEWPx4YleAo2tmC6fC5VdFpLw35/view

    double s;
    
    for (int i = 0; i < v.size(); ++i)
    {
        s = s + 1 / (v_err[i] * v_err[i]);
    }

    s = 1. / s;

    return sqrt(s);   
}