#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

bool debug = false;
bool debug2 = false;

class csvdata
{

public:
    csvdata();
    csvdata(const char * filename);
    long int    getlInt     (int row, int column);
    long double getlDouble  (int row, int column);
    int         getInt      (int row, int column);
    double      getDouble   (int row, int column);
    std::string getString   (int row, int column);

    void replace     (int row, int coumn, const char * content);
    void appendline  (const char * strline);
    void setDelimiters(char d);
    void flush();

private:
    char sep = ',';
    std::string line;
    std::string filepath;
    std::fstream rfile;
    std::vector<std::string> filecontent;
};

csvdata::csvdata(const char * filename)
{
    //if (debug) cout << endl << endl << filename << endl << endl;
    //if (debug2) system(strcat((strcat("more", " ")), filename));
    std::cerr << "\n csvData lib is opening " << filename;
    std::string str(filename);
    filepath = str;
    rfile.open(filename, std::fstream::in | std::fstream::out);
    if (rfile.is_open())
    {
        std::string line;
        while (std::getline(rfile, line))
        {
            filecontent.push_back(line);
        }
        std::cerr << "\n Done.\n";
        //__________DEBUG
            if (debug == true)
            {
                cout <<  "\nline_no: " << filecontent.size();
            }
        //__________END_DEBUG
    }
    else
        std::cerr << "\n Err: unable to open the file.\n";
}

void csvdata::setDelimiters(char d)
{
    sep = d;
}

long int csvdata::getlInt(int row, int column)
{
    return std::stol(getString(row, column));
}

long double csvdata::getlDouble(int row, int column)
{
    return std::stold(getString(row, column));
}

int csvdata::getInt(int row, int column)
{
    return std::stoi(getString(row, column));
}

double csvdata::getDouble(int row, int column)
{
    return std::stod(getString(row, column));
}

std::string csvdata::getString(int row, int column)
{
    if (row >= filecontent.size())
    {
        std::cerr << "\n csvData lib: Err: out of bounds. File has " << filecontent.size() << " rows, requested " << row;
        std::string snull("");
        return snull;
    }

    std::string currRow = filecontent[row];
    std::vector<std::string> sarray;

    istringstream f(currRow);
    string s;
    while (getline(f, s, sep))
        sarray.push_back(s);
    return sarray[column];
}

void csvdata::replace(int row, int column, const char * content)
{
    if (row >= filecontent.size())
    {
        std::cerr << "\n csvData lib: Err: out of bounds. File has " << filecontent.size() << " rows, requested " << row;
        return;
    }

    std::string currRow = filecontent[row];
    std::vector<std::string> sarray;

    istringstream f(currRow);
    string s;
    while (getline(f, s, sep))
        sarray.push_back(s);

    std::string rep(content);
    sarray[column] = rep;
    //WRITE FILE
}

void csvdata::appendline(const char * strline)
{
    //cout << std::string(strline).c_str();
    rfile << std::string(strline);
}

void csvdata::flush()
{
    rfile.close();
}
