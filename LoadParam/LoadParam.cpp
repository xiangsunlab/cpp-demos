/* This code illustrates how to load paramters from file and commandline
   To compile: g++ -Wall -g -O3 -o LoadParam LoadParam.cpp aux_param.cpp
   To run:     ./LoadParam <control.inp> [--Parame1=value] ...
   (c) Xiang Sun 2020
 */

#include "aux_param.h"

int main (int argc, char *argv[]) {
    
    double start_time = get_wall_time();

    Parameters& par = Parameters::get();

    init(argc, argv, par);
    
    
    //--------------- do something -----------------
    cout.precision(14);
    int header(0);
    int j, n;
    
    //Load 1D data from file (no header)
    vector<double> data1;
    string filename = par.get().pm["load_file1"].str_value;
    
    if (filename != "") {
        load1Ddata(filename, data1);
    }
    cout << data1[0] << ",  "<< data1[1] << ",  ..." << endl << endl;
    
    
    //Load 1D data from file (with header)
    filename = "1Ddata-header.dat";
    header = 1;
    load1Ddata(filename, data1, header);
    cout << data1[0] << ",  "<< data1[1] << ",  ..." << endl << endl;
    
    
    //Load 2D data from file (no header)
    Real_Matrix data2;
    filename = par.get().pm["load_file2"].str_value;
    int numColumn = 4;
    
    if (filename != "") {
        load2Ddata(filename, data2, numColumn);
    }
    for (j = 0; j < numColumn; j++) cout << data2[0][j]  << "\t";
    cout << endl;
    for (j = 0; j < numColumn; j++) cout << data2[1][j]  << "\t";
    cout << endl;
    n = data2.size();
    cout << "..." << endl;
    for (j = 0; j < numColumn; j++) cout << data2[n-1][j]  << "\t";
    cout << endl;
    
    
    //Load 2D data from file (with header)
    filename = "2Ddata-header.dat";
    header = 1; //there is a header line
    string colnames;
    
    if (filename != "") {
        colnames = load2Ddata(filename, data2, numColumn, ',', header);
    }
    cout << "column names = " << colnames << endl;
    for (j = 0; j < numColumn; j++) cout << data2[0][j]  << "\t";
    cout << endl;
    for (j = 0; j < numColumn; j++) cout << data2[1][j]  << "\t";
    cout << endl;
    n = data2.size();
    cout << "..." << endl;
    for (j = 0; j < numColumn; j++) cout << data2[n-1][j]  << "\t";
    cout << endl;
    
    //-----------------------------------------
    
    cout << "\n--> Finished successfully!\n    Time used "
         << get_wall_time() - start_time << " sec." << endl;

    return 0;
}
