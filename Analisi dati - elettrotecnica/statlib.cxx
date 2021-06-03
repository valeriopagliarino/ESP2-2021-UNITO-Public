#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
//Include CERN ROOT LIBRARIES ROOStat etc.


class statlib
{

public:
//only public functions, the object itself is not to be used.

double pvalChi2left(double chi2, int ndof);
double pvalChi2right(double chi2, int ndof);
double pvalChi2twotailed(double chi2, int ndof);
double pvalZtwotailed(double z);
double pvalTtwotailed(double t, int ndof);
double pvalFleft(double f, int ndof_num, int ndof_den);
double pvalFright(double f, int ndof_num, int ndof_den);
double quickZtwotailed(double val1, double val2, double err1, double err2, std::fstream &streambuffer);

};


//Implementation


double statlib::pvalChi2left(double chi2, int ndof)
{
    return ROOT::Math::chisquared_cdf(chi2, (double)ndof);
}

double statlib::pvalChi2right(double chi2, int ndof)
{
    return ROOT::Math::chisquared_cdf(chi2, (double)ndof);
}

double statlib::pvalChi2twotailed(double chi2, int ndof)
{
    double pchi2l = ROOT::Math::chisquared_cdf(chi2, (double)ndof);
    double pchi2r = ROOT::Math::chisquared_cdf(chi2, (double)ndof);
    pchi2l = pchi2l * 2;
    pchi2r = pchi2r * 2;
    double pchi2 = pchi2l;
    if (pchi2 >= 1) pchi2 = pchi2r;
    return pchi2;
}

double statlib::pvalZtwotailed(double z)
{
   if (z > 0) z = z * -1.;
   return (2. * (ROOT::Math::normal_cdf(z)));
}

double statlib::pvalTtwotailed(double t, int ndof)
{
    if (t > 0) t = t * -1.;
    return ROOT::Math::tdistribution_cdf(t, (double)ndof);
}

double statlib::pvalFleft(double f, int ndof_num, int ndof_den)
{
    return ROOT::Math::fdistribution_cdf(f, (double)ndof_num, (double)ndof_den);
}

double statlib::pvalFright(double f, int ndof_num, int ndof_den)
{
    return 1 - (ROOT::Math::fdistribution_cdf(f, (double)ndof_num, (double)ndof_den));
}

double quickZtwotailed(double val1, double val2, double err1, double err2, std::fstream &streambuffer)
{
    double z    = (val1 - val2) / sqrt(err1 * err1 + err2 * err2);
    if (z > 0) z = z * -1.;
    double pval = (2. * (ROOT::Math::normal_cdf(z)));

    streambuffer << "\nZ test 2-tailed:  Z = " << z << "   pval = " << pval;
    return pval;
}