/*
 * Implementations of the commonly-used functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "CommFunc.h"

double CommFunc::Abs(const double &x)
{
	complex<double> cld(x);
	double ldAbs = abs(cld);
	return(ldAbs);
}

double CommFunc::sum(const vector<double> &x)
{
    int size = x.size();
    int i=0;
    double d_buf=0.0;
    for(i=0; i<size; i++) d_buf+=x[i];
    return (double)d_buf;
}

double CommFunc::mean(const vector<double> &x)
{
    int size = x.size();
    int i=0;
    double d_buf=0.0;
    for(i=0; i<size; i++) d_buf+=x[i];
    d_buf/=(double)size;
    return (double)d_buf;
}

double CommFunc::median(const vector<double> &x)
{
    vector<double> b(x);
    int size = b.size();
    if(size==1) return b[0];
    stable_sort(b.begin(), b.end());
    if(size%2==1) return b[(size-1)/2];
    else return (b[size/2]+b[size/2-1])/2;
}

double CommFunc::var(const vector<double> &x)
{
    int size = x.size();
    if(size<=1) return(0.0);
    int i=0;
    double mu=0.0, s2=0.0;
    for(i=0; i<size; i++) mu+=x[i];
    mu/=(double)size;
    for(i=0; i<size; i++) s2+=(x[i]-mu)*(x[i]-mu);
    s2/=(double)(size-1);
    return (double)s2;
}

double CommFunc::cov(const vector<double> &x, const vector<double> &y)
{
    int size = x.size();
    int i=0;
    double mu1=0.0, mu2=0.0, c=0.0;
    for(i=0; i<size; i++){
        mu1+=x[i];
        mu2+=y[i];
    }
    mu1/=(double)size;
    mu2/=(double)size;

    for(i=0; i<size; i++) c+=(x[i]-mu1)*(y[i]-mu2);
    c/=(double)(size-1);
    return c;
}

bool CommFunc::FloatEqual(double lhs, double rhs)
{
	if (Abs(lhs - rhs) < FloatErr) return true;
	return false;
}

bool CommFunc::FloatNotEqual(double lhs, double rhs)
{
	if (Abs(lhs - rhs) >= FloatErr) return true;
	return false;
}

const double CommFunc::Sqr(const double &a)
{
	return a*a;
}

const double CommFunc::Max(const double &a, const double &b)
{
	return b > a ? (b) : (a);
}

const double CommFunc::Min(const double &a, const double &b)
{
	return b < a ? (b) : (a);
}

const double CommFunc::Sign(const double &a, const double &b)
{
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

int CommFunc::rand_seed()
{
  	stringstream str_strm;
	str_strm<<time(NULL);
	string seed_str=str_strm.str();
    reverse(seed_str.begin(), seed_str.end());
    seed_str.erase(seed_str.begin()+7, seed_str.end());
	return(abs(atoi(seed_str.c_str())));
}

void CommFunc::FileExist(string filename)
{
    ifstream ifile(filename.c_str());
    if(!ifile) throw("Error: can not open the file ["+filename+"] to read.");
}
int CommFunc::max_abs_id(vector<double> &zsxz)
{
    int id=0;
    double tmpVal, cmpVal=abs(zsxz[0]);
    for( int i=1;i<zsxz.size();i++)
    {
        tmpVal=abs(zsxz[i]);
        if( cmpVal-tmpVal < 1e-6)
        {
            cmpVal=tmpVal;
            id=i;
        }
    }
    return(id);
}

int CommFunc::max_abs_id(VectorXd &zsxz)
{
    int id=0;
    double tmpVal, cmpVal=fabs(zsxz[0]);
    for( int i=1;i<zsxz.size();i++)
    {        
        tmpVal=fabs(zsxz[i]);
        if( cmpVal-tmpVal < 1e-6)
        {
            cmpVal=tmpVal;
            id=i;
         }
    }
    return(id);
}


void CommFunc::getRank(vector<double> &a, vector<int> &b)
{
    b.resize(a.size());
    for (int i = (int)a.size()-1; i >= 0; i--)
    {
        int count = 0;
        for (int j = 0; j < a.size(); j++) if (a[j] < a[i]) count++;
        b[i] = count;
    }
}
void CommFunc::getRank(vector<int> &a, vector<int> &b)
{
    b.resize(a.size());
    for (int i = (int)a.size()-1; i >= 0; i--)
    {
        int count = 0;
        for (int j = 0; j < a.size(); j++) if (a[j] < a[i]) count++;
        b[i] = count;
    }
}
void CommFunc::getRank_norep(vector<int> &a, vector<int> &b)
{
    b.resize(a.size());
    map<int,int> rep_chck;
    long mapsize=0;
    for (int i = (int)a.size()-1; i >= 0; i--)
    {
        int count = 0;
        for (int j = 0; j < a.size(); j++) if (a[j] < a[i]) count++;
        rep_chck.insert(pair<int,int>(count,i));
        while(rep_chck.size()==mapsize)
        {
            count++;
            rep_chck.insert(pair<int,int>(count,i));
        }
        mapsize=rep_chck.size();
        b[i] = count;
    }
}
void CommFunc::getUnique(vector<uint32_t> &a)
{
    sort(a.begin(),a.end());
    vector<uint32_t> ::iterator it=unique(a.begin(),a.end());
    a.erase(it,a.end());
}
void CommFunc::getUnique(vector<int> &a)
{
    sort(a.begin(),a.end());
    vector<int> ::iterator it=unique(a.begin(),a.end());
    a.erase(it,a.end());
}
void CommFunc::match(const vector<uint32_t> &VecA, const vector<uint32_t> &VecB, vector<int> &VecC)
{
    int i=0;
    map<uint32_t, int> id_map;
    map<uint32_t, int>::iterator iter;
    VecC.clear();
    for(i=0; i<VecB.size(); i++) id_map.insert(pair<uint32_t,int>(VecB[i], i));
    for(i=0; i<VecA.size(); i++){
        iter=id_map.find(VecA[i]);
        if(iter==id_map.end()) VecC.push_back(-9);
        else VecC.push_back(iter->second);
    }    
}

void CommFunc::match_only(const vector<uint32_t> &VecA, const vector<uint32_t> &VecB, vector<int> &VecC)
{
    int i=0;
    map<uint32_t, int> id_map;
    map<uint32_t, int>::iterator iter;
    VecC.clear();
    for(i=0; i<VecB.size(); i++) id_map.insert(pair<uint32_t,int>(VecB[i], i));
    for(i=0; i<VecA.size(); i++){
        iter=id_map.find(VecA[i]);
        if(iter!=id_map.end()) VecC.push_back(iter->second);    
    }    
}

void CommFunc::match_only(const vector<int> &VecA, const vector<int> &VecB, vector<int> &VecC)
{
    int i=0;
    map<int, int> id_map;
    map<int, int>::iterator iter;
    VecC.clear();
    for(i=0; i<VecB.size(); i++) id_map.insert(pair<int,int>(VecB[i], i));
    for(i=0; i<VecA.size(); i++){
        iter=id_map.find(VecA[i]);
        if(iter!=id_map.end()) VecC.push_back(iter->second);    
    }    
}

void CommFunc::strcpy2(char** to, string from)
{
    char* tmp=new char[from.size() + 1];
    copy(from.begin(), from.end(), tmp);
    tmp[from.size()]='\0';
    *to=tmp;
}
float CommFunc::readfloat(FILE *f) {
    float v;
    fread((void*)(&v), sizeof(v), 1, f);
    return v;
}
uint64_t CommFunc::readuint64(FILE *f) {
    uint64_t v;
    fread((void*)(&v), sizeof(v), 1, f);
    return v;
}
uint32_t CommFunc::readuint32(FILE *f) {
    uint32_t v;
    fread((void*)(&v), sizeof(v), 1, f);
    return v;
}
int CommFunc::readint(FILE *f) {
    int v;
    fread((void*)(&v), sizeof(v), 1, f);
    return v;
}
double  CommFunc::cor(VectorXd &Y, VectorXd &X, bool centered)
{
    if(Y.size()!= X.size())
    {
        printf("The lenght of vectors not match.\n");
        exit(EXIT_FAILURE);
    }
    double ld=-9;
    long n=Y.size();
    if(centered)
    {
        double xx=X.dot(X), yy=Y.dot(Y), xy=X.dot(Y);
        ld=xy/(sqrt(xx*yy));
    }
    else
    {
        double ysum=Y.sum();
        double xsum=X.sum();
        double xx=X.dot(X), yy=Y.dot(Y), xy=X.dot(Y);
        ld=(n*xy-xsum*ysum)/(sqrt(n*xx-xsum*xsum)*sqrt(n*yy-ysum*ysum));
    }
    return ld;
}

double CommFunc::cor(vector<double> &y, vector<double> &x)
{
    long N = x.size();
    if (N != y.size() || N < 1) {
        printf("Error: The lengths of x and y do not match.\n");
        exit(EXIT_FAILURE);
    }
    
    int i = 0;
    double d_buf = 0.0, y_mu = 0.0, x_mu = 0.0, x_var = 0.0, y_var = 0.0, cov = 0.0;
    for (i = 0; i < N; i++) {
        x_mu += x[i];
        y_mu += y[i];
    }
    x_mu /= (double) N;
    y_mu /= (double) N;
    for (i = 0; i < N; i++) {
        d_buf = (x[i] - x_mu);
        x_var += d_buf*d_buf;
        d_buf = (y[i] - y_mu);
        y_var += d_buf*d_buf;
    }
    x_var /= (double) (N - 1.0);
    y_var /= (double) (N - 1.0);
    for (i = 0; i < N; i++) cov += (x[i] - x_mu)*(y[i] - y_mu);
    cov /= (double) (N - 1);
    double a = 0.0, b = 0.0, sse = 0.0, a_se = 0.0, b_se = 0.0, r = 0.0;
    if (x_var > 0.0) b = cov / x_var;
    a = y_mu - b*x_mu;
    for (i = 0; i < N; i++) {
        d_buf = y[i] - a - b * x[i];
        sse += d_buf*d_buf;
    }
    if (x_var > 0.0) {
        a_se = sqrt((sse / (N - 2.0))*(1.0 / N + x_mu * x_mu / (x_var * (N - 1.0))));
        b_se = sqrt(sse / x_var / (N - 1.0) / (N - 2.0));
    }
    if (x_var > 0.0 && y_var > 0.0) {
        r = cov / sqrt(y_var * x_var);
    }
    
    return (r);
}

