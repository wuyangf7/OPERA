/*
 * Interface to the string-operation functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#ifndef _STRFUNC_H
#define _STRFUNC_H

#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
typedef unsigned int         uint32_t;
using namespace std;

namespace StrFunc
{
	bool str_within_quto(const string &str, string &str_buf);
	int split_string(const string &str, vector<string> &vec_str, string separator=" ,\t;\n");
	string first_string(const string &str, const char separator);
	string last_string(const string &str, const char separator);
    void to_upper(char* str, int len);
	void to_upper(string &str);
	void to_lower(string &str);
	string get_sub_str(const string & rst, int pos);
	bool StrEqual(const string &StrA, const string &StrB, bool NoCaseSens=true);
	bool StrVecEqual(const vector<string> &VsBufA, const vector<string> &VsBufB, int Pos);

	// find a string in a string vector ignoring upper or lower case
	vector<string>::iterator find(vector<string> &target_vs, const string &target_str);

	// find a char in a string ignoring upper or lower case
	string::iterator find(string &target_str, const char target_ch);

	// go to the postion of a give string in a stream ignoring upper or lower case
	bool goto_str(std::istream &in_file, const string &str);

	// rewind a stream
	void rewind_if(std::istream &in_file);

	// match two vectors
	void match(const vector<string> &VecA, const vector<string> &VecB, vector<int> &VecC);
	void match_only(const vector<string> &VecA, const vector<string> &VecB, vector<int> &VecC);
   
    
    void set_complement(const vector<string> &VecA, const vector<string> &VecB,const vector<int> &tmp, vector<int> &VecC);
    void set_complement(const vector<string> &VecA, const vector<string> &VecB,const vector<int> &tmp, vector<uint32_t> &VecC);
    int split_string_skip(const string &str, vector<string> &vec_str, string separator, int num2skip);
    void match_only(const vector<string> &VecA, const vector<string> &VecB, vector<uint32_t> &VecC);
    bool has_suffix(const std::string &str, const std::string &suffix);
    void set_intersect(const vector<string> &VecA, const vector<string> &VecB, vector<string> &VecC);
    void set_intersect(const vector<int> &VecA, const vector<int> &VecB, vector<int> &VecC);
    void set_complement(const vector<int> &toRm, const vector<int> &source, vector<int> &VecC);
    bool stringNumCheck(string a, int num);

     /* Permute vector elements to all the combinations
	 * @param elements: vector of each possible values
	 * @param combines: out: combinations of these vectors;
	 * @param temp and col:  only for internal use, don't change these values;
	 * @example
	 *   #include <vector>
	     vector<vector<int>> elements = {
	           {0, 1, 2},
	           {10, 11, 12},
	           {20, 21, 22}
	     };
	     vector<vector<int>> combines;
	     permute_vector(elements, combines);
	 * @depends vector c++11
	 * @Bug report: zhili<zhilizheng@outlook.com>
	 */
	template <typename T>
	void permute_vector(const std::vector<std::vector<T>> &elements, std::vector<std::vector<T>> &combines, std::vector<T> temp = {}, size_t col = 0){
	    size_t e_size = elements.size();
	    temp.resize(e_size);
	    if(col < e_size) {
	        for(size_t i = 0; i < elements[col].size(); i++) {
	            temp[col] = elements[col][i];
	            permute_vector(elements, combines, temp, col + 1);
	        }
	    } else {
	        combines.emplace_back(temp);
	   }
	}
}

#endif
