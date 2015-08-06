//
//  Global.h
//  GetBaseCounts
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 7/31/15.
//  Copyright (c) 2015 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//

#ifndef __GetBaseCounts__Global__
#define __GetBaseCounts__Global__

#include <stdio.h>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include "omp.h"
using namespace std;

void split(const string& line, char delim, vector<std::string>& parsed_item, bool ignore_empty_item);

bool is_number(const string& s);

void output_reference_sequence(string output_fastafile, map<string, string>& ref_seq, vector<string>& orignal_header);

void load_reference_sequence_speedup(string fasta_filename, map<string, string>& ref_seq, int num_thread);


#endif /* defined(__GetBaseCounts__Global__) */
