#include "aux_param.h"


// --------------------------- ParamEntry ------------------------
ParamEntry::ParamEntry(string n, string v) {
    name = n;
    str_value = v;
    int_value = 0;
    real_value = 0.0;
    type = "str";
}

ParamEntry::ParamEntry(string n, int v) {
    stringstream ss;
    ss << v;
    name = n;
    str_value = ss.str();
    int_value = v;
    real_value = v;
    type = "int";
}

ParamEntry::ParamEntry(string n, Real v) {
    stringstream ss;
    ss << v;
    name = n;
    str_value = ss.str();
    int_value = 0;
    real_value = v;
    type = "real";
}



// --------------------------- Parameters ------------------------
void Parameters::addParam(string key, ParamEntry entry) {
    pm[key] = entry;
}

void loadParam(string inputfile, Parameters& para) {
    //reading data type <str>, <int>, <real> labels from input control file
    ParamEntry mypar;
    stringstream ss;
    string line;
    string key;
    int i;
    Real r;
    string s;
    vector<string> list;
    ifstream input;
    input.open(inputfile);
    if (!input.is_open()) {
        throw std::runtime_error("Error: unable to open input file: "+inputfile);
    }

    //mypar = ParamEntry(key, value);
    //para.pm[inputfile] = mypar;

    while (getline(input, line)) {
        ss.str("");
        ss.clear();
        ss << line;
        ss >> key;
        if (key[0] == '#') continue;
        list.clear();
        split_string(line, list);
        
        if (list[2] == "<str>") {
            s = list[1];
            mypar = ParamEntry(key, s);
        }
        else if (list[2] == "<int>") {
            ss >> i;
            mypar = ParamEntry(key, i);
        }
        else if (list[2] == "<real>") {
            ss >> r;
            mypar = ParamEntry(key, r);
        }
        else {
            throw std::runtime_error("Error: unkown parameter type: " + list[2] + " for " + key + " in " + inputfile);
        }
        //add one parameter
        para.pm[key] = mypar;
    }
}


//---------------------------- LOADING PARAMETERS ------------------------------

void init(int argc, char *argv[], Parameters& par) {
    ifstream infile;
    stringstream ss;
    //ss.str("");
    //ss.clear();
    string varstr("");
    string inputfile("control.inp");
    cout.precision(12);

    //cout << ">>> # of command-line argument: " << argc-1 << endl;
    if (argc == 2) {
        ss.str("");
        ss << argv[1];
        inputfile = ss.str();
        ss.clear();
        readInputControlFile(par, inputfile);
    }
    else if (argc >= 3) {
        ss.str("");
        ss << argv[1];
        inputfile = ss.str();
        ss.clear();
        readInputControlFile(par, inputfile);

        //load overwritten parameters from command line like "--jobid=123"
        cout << "--> Overwrote the following parameters from command line: \n" << endl;
        for (int n = 2; n < argc; n++) {
            ss.str("");
            ss << argv[n];
            varstr = ss.str();
            readInputControlCommandline(par, varstr);
            ss.clear();
        }
        cout << endl;
    }
    else {
        throw std::runtime_error("Usage: ./EXECUTABLE  <control.inp>  [--Param=Value ...]");
    }
    return;
}




void readInputControlFile(Parameters& par, string inputfile) {
    int fw = 14; //field width for output
    //load control parameters from input file
    loadParam(inputfile, par);

    //add command-line parameters
    //ParamEntry ent=ParamEntry("jobid", jobid);
    //par.get().addParam("jobid", ent);

    map<string, ParamEntry>::iterator it; //parameter map iterator
    
    //output the loaded control parameters
    cout << "--> Loaded " << par.get().pm.size() << " parameters from input file <"
    << inputfile << ">:" << endl;
    cout << endl;
    for (it = par.get().pm.begin(); it != par.get().pm.end(); it++) {
        if (it->second.type == "str") {
            cout << left << "    " << setw(fw) << it ->first << left << " : " << setw(fw) << it->second.str_value  << "   <str>" << endl;
        }
    }
    for (it = par.get().pm.begin(); it != par.get().pm.end(); it++) {
        if (it->second.type == "int") {
            cout << left << "    " << setw(fw) << it ->first << left << " : " << setw(fw) << it->second.int_value << "   <int>"<< endl;
        }
    }
    for (it = par.get().pm.begin(); it != par.get().pm.end(); it++) {
        if (it->second.type == "real") {
            cout << left << "    " << setw(fw) << it ->first << left << " : " << setw(fw) << it->second.real_value <<  "   <real>"<< endl;
        }
    }
    cout << endl;
    return;
}



void readInputControlCommandline(Parameters& par, string inputline) {
    //overwrite parameters from commandline of "--id=123" style
    int fw = 14; //field width for output
    if (inputline.substr(0,2) != "--")
        throw std::runtime_error("Error: commandline overwrite parameter should be like --Param=Value, starting with --");

    string shortline;
    shortline = inputline.substr(2);
    //parse "--id=123" style into key and value strings
    vector<string> list;
    list = parse_csv_line(shortline, '=');
    if (list.size() != 2)
        throw std::runtime_error("Error: commandline overwrite parameter should be like --Param=Value, indicating both key and value!");

    string keystr = list[0];
    string valstr = list[1];
    if ( par.get().pm.find(keystr) == par.get().pm.end() ) {
        // not found key
        cout << "Error: parameter " << keystr << " was not preset in control file." << endl;
    } else {
        // found key
        if (par.get().pm[keystr].type == "str") {
            ParamEntry ent=ParamEntry(keystr, valstr);
            par.get().addParam(keystr, ent);
        }
        else if (par.get().pm[keystr].type == "int") {
            int i = stoi(valstr);
            ParamEntry ent=ParamEntry(keystr, i);
            par.get().addParam(keystr, ent);
        }
        else if (par.get().pm[keystr].type == "real") {
            double d = stod(valstr);
            ParamEntry ent=ParamEntry(keystr, d);
            par.get().addParam(keystr, ent);
        }
        cout << left << "    " << setw(fw) << keystr << left << " : " << setw(fw) << valstr
        << "   <" << setw(5) << par.get().pm[keystr].type + ">" << "  [Overwrite]" << endl;
        
    }
    return;
}




void load1Ddata(string inputfile, vector<double>& data, int header) {
    //load 1D colum data, default setting "header = 0" means no header
    ifstream infile;
    string line;
    vector<string> list;
    double x;
    int count(0);
    int orig_size;
    orig_size = data.size();
    
    infile.open(inputfile);
    if (!infile.is_open()) {
        throw std::runtime_error("Error: cannot open input file: "+inputfile);
    }
    //first line (if there is header(!=0), then skip the first line)
    if (header != 0) {
        std::getline(infile, line);
    }
    //the data starts
    count = 0;
    while(std::getline(infile, line)) {
        list = parse_csv_line(line, ',');
        x = stod(list[0]);  //col index: j
        if (count < orig_size) {
            data[count] = x;
        }
        else {
            data.push_back(x);
        }
        count++;
    }
    data.erase(data.begin()+count, data.end());
    cout << ">>  Loaded 1D data (" << data.size() << " rows) from file: <" << inputfile << ">" << endl;
    return;
}




string load2Ddata(string inputfile, Real_Matrix& data, const int ncol, char delimiter,
                int header) {
    //load 2D colum data, default setting: delimiter=',' and "header = 0" means no header
    ifstream infile;
    string line;
    vector<string> list;
    string colnames;
    double x;
    int j;
    int count(0);
    int orig_size;
    orig_size = data.size();
    vector<double> emptyrow(ncol, 0);
    
    infile.open(inputfile);
    if (!infile.is_open()) {
        throw std::runtime_error("Error: cannot open input file: "+inputfile);
    }
    //first line (if there is header(!=0), then save to colnames)
    if (header != 0) {
        std::getline(infile, line);
        colnames = line;
    }
    //the data starts
    count = 0;
    while(std::getline(infile, line)) {
        list = parse_csv_line(line, delimiter);
        if (count < orig_size) {
            data[count].resize(ncol,0);
        }
        else {
            data.push_back(emptyrow);
        }
        for (j = 0; j < ncol; j++) {
            x = stod(list[j]);  //col index: j
            data[count][j] = x;
        }
        count++;
    }
    data.erase(data.begin()+count, data.end());
    cout << ">>  Loaded 2D data (" << data.size() << " rows, " << ncol << " cols) from file: <" << inputfile << ">" << endl;


    return colnames;
}



void save1Ddata(string outputfile, vector<double>& data) {
    int i;
    int n = data.size();
    ofstream outfile;
    outfile.open(outputfile);
    for (i = 0; i < n; i++) {
        outfile << data[i] << endl;
    }
    outfile.close();
    outfile.clear();
    cout << "<<  Saved  1D data (" << data.size() << " rows) to file: <" << outputfile << ">" << endl;
    return;
}



//---------------------------- AUXILIARY FUNCTIONS -----------------------------
void split_string(const string& str, vector<string>& list) {
    std::istringstream iss(str);
    std::copy(std::istream_iterator<std::string>(iss),
         std::istream_iterator<std::string>(),
         std::back_inserter(list)  );
}

vector<string> parse_csv_line(const string& line, char dlim){
    string field;
    vector<string> list;
    istringstream s(line);
    while (getline(s, field, dlim)) {
        list.push_back(field);
    }
    return list;
}

Real delta(int a, int b) {
    if (a == b) return 1;
    else return 0;
}

Real heaviside(Real x) {
    if (x >= 0) return 1;
    else return 0;
}

Real get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (Real)time.tv_sec + (Real)time.tv_usec * .000001;
}

Real get_cpu_time(){
    return (Real)clock() / CLOCKS_PER_SEC;
}

int Job_finished(int &jobdone, int count, int total) {
    int tenpercent;
    tenpercent = static_cast<int> (10 * static_cast<Real> (count)/ static_cast<Real> (total) );
    if ( tenpercent > jobdone ) {
        jobdone = tenpercent;
        time_t rawtime;
        time (&rawtime);
        cout << "    Job finished "<< jobdone <<"0%.  Local time: " << ctime(&rawtime) << endl;
    }
    return tenpercent;
}

int Job_finished_percent(int &jobdone, int count, int total) {
    int percent;
    percent = static_cast<int> (100 * static_cast<Real> (count)/ static_cast<Real> (total) );
    if ( percent > jobdone ) {
        jobdone = percent;
        time_t rawtime;
        time (&rawtime);
        cout << "    Job finished "<< jobdone <<"%.  Local time: " << ctime(&rawtime) << endl;
    }
    return percent;
}


double** Create_matrix_pointer(int row, int col) {
    double **matrix = new double* [col];
    matrix[0] = new double [col*row];
    for (int i = 1; i < col; ++i)
        matrix[i] = matrix[i-1] + row;
    return matrix;
}


void FFT(int dir, int m, Real *x, Real *y) {
    /* This code computes an in-place complex-to-complex FFT Written by Paul Bourke
      x and y are the real and imaginary arrays of N=2^m points.
      dir =  1 gives forward transform
      dir = -1 gives reverse transform
      Formula: forward
                  N-1
                  ---
              1   \           - i 2 pi k n / N
      X(n) = ----  >   x(k) e                       = forward transform
              1   /                                    n=0..N-1
                  ---
                  k=0
      
      Formula: reverse
                  N-1
                  ---
               1  \           i 2 pi k n / N
      X(n) =  ---  >   x(k) e                  = reverse transform
               N  /                               n=0..N-1
                  ---
                  k=0
      */
        
    int n,i,i1,j,k,i2,l,l1,l2;
    double c1,c2,tx,ty,t1,t2,u1,u2,z;
    
    // Calculate the number of points
    n = 1;
    for (i=0;i<m;i++) n *= 2;
    
    // Do the bit reversal
    i2 = n >> 1; //i2 = (010 ...0)_2,second highest bit of n=(100 ...0)_2
    j = 0; //reversely bit accumulater from the second highest bit, i2.
    for (i=0;i<n-1;i++) {
        if (i < j) {
            tx = x[i]; //swap(i,j)
            ty = y[i];
            x[i] = x[j];
            y[i] = y[j];
            x[j] = tx;
            y[j] = ty;
        }
        //to find the highest non-one bit, k, from the second highest bit
        k = i2;
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k; //add 1 reversly
    }
    
    // Compute the Radix-2 FFT: Cooley-Tukey Algorithm
    c1 = -1.0; // c1+i*c2 = -1 = c^(i 2Pi/2) = W_2, def W_N^j = e^(i 2j*Pi/N)
    c2 = 0.0;
    l2 = 1;
    for (l=0;l<m;l++) {
        l1 = l2;
        l2 <<= 1;
        u1 = 1.0;
        u2 = 0.0;
        for (j=0;j<l1;j++) {
            for (i=j;i<n;i+=l2) {
                //Butterfly calculation of x,y[i] and x,y[i1]:
                //t1+i*t2 =(u1+i*u2)(x[i1]+i*y[i2]) where u1+i*u2=W_N^j=e^(i 2j*Pi/N)
                i1 = i + l1;
                t1 = u1 * x[i1] - u2 * y[i1];
                t2 = u1 * y[i1] + u2 * x[i1];
                x[i1] = x[i] - t1;
                y[i1] = y[i] - t2;
                x[i] += t1;
                y[i] += t2;
            }
            // i1+i*u2 *= c1+i*c2, or W_N
            z =  u1 * c1 - u2 * c2;
            u2 = u1 * c2 + u2 * c1;
            u1 = z;
        }
        //c1+i*c2 = sqrt(c1+i*c2) eg. W_2 --> W_4 ...
        c2 = sqrt((1.0 - c1) / 2.0);
        if (dir == 1)
            c2 = -c2;
        c1 = sqrt((1.0 + c1) / 2.0);
    }
    
    // times STEPS*DeltaT forward FFT (time --> freq)
    /*if (dir == 1) {
        for (i=0; i<n; i++) {
        x[i] *= 1;//DeltaT;
        y[i] *= 1;//DeltaT;
        }
        }*/
    
    // Scaling for inverse transform
    
    if (dir == -1) {
        for (i=0;i<n;i++) {
            x[i] /= n;
            y[i] /= n;
        }
    }
    
    /*
    //for symmetrical FT,
    double sqn;
    sqn = sqrt(n);
    for (i=0;i<n;i++) {
        x[i] /= sqn;
        y[i] /= sqn;
    }
    */
    return;
}


void FFT(int dir, int m, vector<Real> &x, vector<Real> &y) {
    /* This code computes an in-place complex-to-complex FFT Written by Paul Bourke
      x and y are the real and imaginary arrays of N=2^m points.
      dir =  1 gives forward transform
      dir = -1 gives reverse transform
      Formula: forward
                  N-1
                  ---
              1   \           - i 2 pi k n / N
      X(n) = ----  >   x(k) e                       = forward transform
              1   /                                    n=0..N-1
                  ---
                  k=0
      
      Formula: reverse
                  N-1
                  ---
               1  \           i 2 pi k n / N
      X(n) =  ---  >   x(k) e                  = reverse transform
               N  /                               n=0..N-1
                  ---
                  k=0
      */
        
    int n,i,i1,j,k,i2,l,l1,l2;
    double c1,c2,tx,ty,t1,t2,u1,u2,z;
    
    // Calculate the number of points
    n = 1;
    for (i=0;i<m;i++) n *= 2;
    
    // Do the bit reversal
    i2 = n >> 1; //i2 = (010 ...0)_2,second highest bit of n=(100 ...0)_2
    j = 0; //reversely bit accumulater from the second highest bit, i2.
    for (i=0;i<n-1;i++) {
        if (i < j) {
            tx = x[i]; //swap(i,j)
            ty = y[i];
            x[i] = x[j];
            y[i] = y[j];
            x[j] = tx;
            y[j] = ty;
        }
        //to find the highest non-one bit, k, from the second highest bit
        k = i2;
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k; //add 1 reversly
    }
    
    // Compute the Radix-2 FFT: Cooley-Tukey Algorithm
    c1 = -1.0; // c1+i*c2 = -1 = c^(i 2Pi/2) = W_2, def W_N^j = e^(i 2j*Pi/N)
    c2 = 0.0;
    l2 = 1;
    for (l=0;l<m;l++) {
        l1 = l2;
        l2 <<= 1;
        u1 = 1.0;
        u2 = 0.0;
        for (j=0;j<l1;j++) {
            for (i=j;i<n;i+=l2) {
                //Butterfly calculation of x,y[i] and x,y[i1]:
                //t1+i*t2 =(u1+i*u2)(x[i1]+i*y[i2]) where u1+i*u2=W_N^j=e^(i 2j*Pi/N)
                i1 = i + l1;
                t1 = u1 * x[i1] - u2 * y[i1];
                t2 = u1 * y[i1] + u2 * x[i1];
                x[i1] = x[i] - t1;
                y[i1] = y[i] - t2;
                x[i] += t1;
                y[i] += t2;
            }
            // i1+i*u2 *= c1+i*c2, or W_N
            z =  u1 * c1 - u2 * c2;
            u2 = u1 * c2 + u2 * c1;
            u1 = z;
        }
        //c1+i*c2 = sqrt(c1+i*c2) eg. W_2 --> W_4 ...
        c2 = sqrt((1.0 - c1) / 2.0);
        if (dir == 1)
            c2 = -c2;
        c1 = sqrt((1.0 + c1) / 2.0);
    }
    
    // times STEPS*DeltaT forward FFT (time --> freq)
    /*if (dir == 1) {
        for (i=0; i<n; i++) {
        x[i] *= 1;//DeltaT;
        y[i] *= 1;//DeltaT;
        }
        }*/
    
    // Scaling for inverse transform
    
    if (dir == -1) {
        for (i=0;i<n;i++) {
            x[i] /= n;
            y[i] /= n;
        }
    }
    
    /*
    //for symmetrical FT,
    double sqn;
    sqn = sqrt(n);
    for (i=0;i<n;i++) {
        x[i] /= sqn;
        y[i] /= sqn;
    }
    */
    return;
}


void CorrFunc(const int trun, const int tcor, double A[], double B[], double Corr[]) {
    //calculate C_AB(t)=<A(0) B(t)> equilibrium correlation function
    //Input: trun numbers, A[trun], B[trun]
    //Output:tcor numbers, Corr[tcor]
    double norm[tcor];
    int t, t0, tt0, tt0max;
    double A0(0);
    
    for (t = 0; t < tcor; t++){
        Corr[t]=0;
        norm[t]=0;
    }
    
    for (t0 = 0; t0 < trun; t0++) {
        A0=A[t0];
        tt0max = min(trun, t0+tcor);
        
        for (tt0 = t0; tt0 < tt0max; tt0++) {
            t=tt0-t0;
            Corr[t]+=A0*B[tt0];
            norm[t]++;
        }
    }
    for (t=0; t< tcor; t++){
        Corr[t]/=norm[t];
    }
    return;
}


void CorrFunc(const int trun, const int tcor, vector<double>& A, vector<double>& B, vector<double>& Corr) {
    //calculate C_AB(t)=<A(0) B(t)> equilibrium correlation function
    //Input: trun numbers, A[trun], B[trun]
    //Output:tcor numbers, Corr[tcor]
    vector<double> norm;
    norm.resize(tcor,0);
    int t, t0, tt0, tt0max;
    double A0(0);
    
    for (t = 0; t < tcor; t++){
        Corr[t]=0;
        norm[t]=0;
    }
    
    for (t0 = 0; t0 < trun; t0++) {
        A0=A[t0];
        tt0max = min(trun, t0+tcor);
        
        for (tt0 = t0; tt0 < tt0max; tt0++) {
            t=tt0-t0;
            Corr[t]+=A0*B[tt0];
            norm[t]++;
        }
    }
    for (t=0; t< tcor; t++){
        Corr[t]/=norm[t];
    }
    return;
}


double Sum(double *data, int n){
    double S = 0;
    for (int i = 0; i< n; i++) {
        S += data[i];
    }
    return S;
}



double Sum(vector<double> data) {
    int n = data.size();
    double S = 0;
    for(int i = 0; i < n; i++) {
        S += data[i];
    }
    return S;
}


double Mean(vector<double> data) {
    int n = data.size();
    double S = 0;
    for(int i = 0; i < n; i++) {
        S += data[i];
    }
    return S/n;
}
