/*
 * Interface to the commonly-used functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#ifndef _COMMFUNC_H
#define _COMMFUNC_H

#define MAX_BUF_SIZE 0x40000000
#if defined _WIN64 || defined _WIN32
#define MAX_LINE_SIZE 0x10000
#else
#define MAX_LINE_SIZE 0x80000
#endif

//#define __DBL_MIN__ 2.225e-308
#define MAX_LINE_BUF 0x1000000

#define MAXSNPNUMPERPROBEINSPARSE 0x300000
#define MAX_PROBE_NUM 0xF0000
#define MAX_SNP_NAME 64
#define DENSE_FILE_TYPE_1 0  // uint32 + floats
#define SPARSE_FILE_TYPE_3F 0x40400000 // uint32 + uint64_t + uint64_ts + uint32_ts + floats
#define SPARSE_FILE_TYPE_3 3 // 10*uint32s + uint64_t + uint64_ts + uint32_ts + floats (indicator+samplesize+snpnumber+probenumber+ 6*-9s +valnumber+cols+rowids+betases) [default]
#define DENSE_FILE_TYPE_3 5  // 10*uint32s + floats (indicator+samplesize+snpnumber+probenumber+ 6*-9s + values) [default]
#define RESERVEDUNITS 16
#define FNAMESIZE 4096

//#define BEST_NUM_HEIDI 41
#define MAX_NUM_LD 500

#define MIN_PVAL_ADJUSTED 1e-150
//typedef unsigned long long         uint64_t;
//typedef unsigned int         uint32_t;

#include <limits>
#include <complex>
#include <vector>
#include <algorithm>
#include <ctime>
#include <fstream>
#include "StrFunc.h"
#include <string.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/Eigenvalues>
#include <queue>
//#include <cmath> //already in Eigen

using namespace Eigen;
using namespace std;

namespace CommFunc
{
	const double FloatErr=numeric_limits<double>::epsilon();
    template <typename T>
    inline T ABS(T const& a)
    {
       return (T{} < a) ? a : -a;
    }
    template <typename T>
    extern void free2(T** to)
    {
        if(*to)
        {
            delete(*to);
            *to=NULL;
        }
    }
    template <typename T>
    extern inline string atosm (T const& a)
    {
        if(a==-9) return("NA");
        stringstream ss;
        ss << a;
        return(ss.str());
    }
	double Abs(const double &x);
	double sum(const vector<double> &x);
    double mean(const vector<double> &x);
    double median(const vector<double> &x);
    double var(const vector<double> &x);
    double cov(const vector<double> &x, const vector<double> &y);
	bool FloatEqual(double lhs, double rhs);
	bool FloatNotEqual(double lhs, double rhs);
	const double Sqr(const double &a);
	const double Max(const double &a, const double &b);
	const double Min(const double &a, const double &b);
	const double Sign(const double &a, const double &b);
	int rand_seed(); //positive value, return random seed using the system time
    void FileExist(string filename);  
    int max_abs_id(VectorXd &zsxz);
    int max_abs_id(vector<double> &zsxz);
    void getRank(vector<double> &a, vector<int> &b);
    void getRank(vector<int> &a, vector<int> &b);
    void getRank_norep(vector<int> &a, vector<int> &b);
    void getUnique(vector<uint32_t> &a);
    void getUnique(vector<int> &a);
    void match(const vector<uint32_t> &VecA, const vector<uint32_t> &VecB, vector<int> &VecC);
    void match_only(const vector<uint32_t> &VecA, const vector<uint32_t> &VecB, vector<int> &VecC);
    void match_only(const vector<int> &VecA, const vector<int> &VecB, vector<int> &VecC);
    static inline unsigned int fputs_checked(const char* ss, FILE* outfile) {
        fputs(ss, outfile);
        return ferror(outfile);
    }
    void strcpy2(char** to, string from);
    float readfloat(FILE *f);
    uint64_t readuint64(FILE *f);
    uint32_t readuint32(FILE *f);
    int readint(FILE *f);
    double cor(vector<double> &y, vector<double> &x);
    double cor(VectorXd &Y, VectorXd &X, bool centered=false);
}

#endif
