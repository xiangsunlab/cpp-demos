#ifndef AUX_PARAM_h
#define AUX_PARAM_h

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <random>
#include <complex>
#include <iomanip>
#include <stdexcept>
#include <memory>
#include <sys/time.h>
using namespace std;

//--------------- TYPES -------------------
typedef double Real;
typedef std::complex<Real>  Complex;
typedef std::vector<vector<Complex> > Complex_Matrix;
typedef std::vector<vector<vector<Complex> > > Complex_3D_Array;
typedef std::vector<vector<Real> > Real_Matrix;
typedef std::vector<vector<vector<Real> > > Real_3D_Array;


//---------------- CONSTANTS ----------------
const Real hbar = 1;
const Real pi = std::acos(-1.0);//3.14159265358979;
const Complex I(0.0,1.0);  //the imaginary I



//---------------- PARAMENTRY CLASS ----------------
// ParamEntry for a single entry of parameter
class ParamEntry {
public:
    string name;
    string str_value;
    int    int_value;
    Real   real_value;
    string type;
    
    ParamEntry(){}
    ParamEntry(string n, string v);
    ParamEntry(string n, int v);
    ParamEntry(string n, Real v);
    ~ParamEntry(){}
};


//---------------- PARAMETERS CLASS ----------------
// Parameters is a singleton class, only 1 instance, stores global parameters
class Parameters {
public:
    map<string, ParamEntry> pm;

    static Parameters& get() {
        static Parameters instance;
        return instance;
    }

    void addParam(string key, ParamEntry entry);
    
private:
    Parameters(){}
    ~Parameters(){}

friend void loadParam(string inputfile, Parameters& para);
};


//--------------- Param loading -----------------------
void init(int argc, char *argv[], Parameters& par);

void readInputControlFile(Parameters& par, string inputfile);

void readInputControlCommandline(Parameters& par, string inputline);

void load1Ddata(string inputfile, vector<double>& data, int header=0);

string load2Ddata(string inputfile, Real_Matrix& data, const int ncol, char delimiter=',', int header=0); //returning header (column names) if exists

void save1Ddata(string outputfile, vector<double>& data);

//--------------- AUXILIARY FUNCTIONS -----------------
Real delta(int a, int b);

Real heaviside(Real x);

Real get_wall_time();

Real get_cpu_time();

int Job_finished(int &jobdone, int count, int total);

int Job_finished_percent(int &jobdone, int count, int total);

void split_string(const string& str, vector<string>& list);

vector<string> parse_csv_line(const string& str, char dlim);

double** Create_matrix_pointer(int row, int col);

void FFT(int dir, int m, Real *x, Real *y);

void FFT(int dir, int m, vector<Real> &x, vector<Real> &y);

void CorrFunc(const int trun, const int tcor, double A[], double B[], double Corr[]);

void CorrFunc(const int trun, const int tcor, vector<double>& A, vector<double>& B, vector<double>& Corr);

double Sum(double *data, int n);

double Sum(vector<double> data);

double Mean(vector<double> data);

//Note template must be defined in header file!
template <typename Container>
void split_template(const std::string& str, Container& cont) {
    std::istringstream iss(str);
    std::copy(std::istream_iterator<std::string>(iss),
         std::istream_iterator<std::string>(),
         std::back_inserter(cont)  );
}


//LAPACK routine (use -llapack -lblas -lgfortran flags when compile)
extern "C" {
    void dsyev_(const char &JOBZ, const char &UPLO, const int &N, double *A,
                const int &LDA, double *W, double *WORK, const int &LWORK,
                int &INFO);
}







//-----------------------------------------------------------------

#endif
