//
//  Global.cpp
//  GetBaseCounts
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 7/31/15.
//  Copyright (c) 2015 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//

#include "Global.h"

using namespace std;

void split(const string& line, char delim, vector<std::string>& parsed_item, bool ignore_empty_item = false)
{
    stringstream my_ss(line);
    string item;
    while (getline(my_ss, item, delim))
    {
#ifdef _PARSING_DEBUG
        cout << "[DEBUG] Parsed item: " << item << endl;
#endif
        if (ignore_empty_item && item.empty())
            continue;    // whether to skip empty item
        parsed_item.push_back(item);
    }
}


void output_reference_sequence(string output_fastafile, map<string, string>& ref_seq, vector<string>& orignal_header)
{
    size_t base_per_line = 80;
    ofstream out_fs(output_fastafile.c_str());
    if(!out_fs)
    {
        cerr << "[ERROR] Fail to open output fasta file: " << output_fastafile << endl;
        exit(1);
    }
    for(size_t i = 0; i < orignal_header.size(); i++)
    {
        map<string, string>::iterator it = ref_seq.find(orignal_header[i]);
        size_t chrom_len = it->second.length();
#ifdef _FASTA_DEBUG
        cout << "[DEBUG] Output reference sequence: " << it->first << ": " << chrom_len << endl;
#endif
        size_t current_len = 0;
        out_fs << ">" << it->first << endl;
        while(current_len < chrom_len)
        {
            size_t output_len = min(base_per_line, chrom_len - current_len);
            out_fs << it->second.substr(current_len, output_len) << endl;
            current_len += output_len;
        }
    }
    out_fs.close();
}


void load_reference_sequence_speedup(string fasta_filename, map<string, string>& ref_seq, int num_thread)
{
    cout << "[INFO] Loading reference sequence: " << fasta_filename << endl;
    string fasta_index_filename = fasta_filename + ".fai";
    ifstream index_fs(fasta_index_filename.c_str());
    if(!index_fs)
    {
        cerr << "[ERROR] Fail to open reference fasta index file: " << fasta_index_filename << endl;
        exit(1);
    }
#pragma omp parallel num_threads(num_thread)
    {
        string line;
        ifstream ref_fs(fasta_filename.c_str());
        if(!ref_fs)
        {
            cerr << "[ERROR] Fail to open reference fasta file: " << fasta_filename << endl;
            exit(1);
        }
        while(!index_fs.eof())
        {
#pragma omp critical(read_index)
            {
                getline(index_fs, line);
            }
            if(!index_fs.eof())
            {
                vector<string> index_items;
                split(line, '\t', index_items);
                string chrom_name = index_items[0];
                int chrom_len = atoi(index_items[1].c_str());
                long long int chrom_offset = atoll(index_items[2].c_str());
                int base_per_line = atoi(index_items[3].c_str());
                int byte_per_line = atoi(index_items[4].c_str());
                int byte_len = chrom_len + (chrom_len / base_per_line) * (byte_per_line - base_per_line);
                string new_seq;
#pragma omp critical(access_ref_seq_map)
                {
                    ref_seq.insert(make_pair(chrom_name, new_seq));
                    ref_seq[chrom_name].resize(chrom_len);
                }
                ref_fs.seekg(chrom_offset);
                char* seq_buff = new char[byte_len];
                ref_fs.read(seq_buff, byte_len);
                string::iterator it_target;
#pragma omp critical(access_ref_seq_map)
                {
                    it_target = ref_seq[chrom_name].begin();
                }
                char* it_source = seq_buff;
                for(int i = 0; i < byte_len; i++)
                {
                    if(!isspace(*it_source))
                    {
                        *it_target = toupper(*it_source);
                        it_target ++;
                    }
                    it_source ++;
                }
                delete[] seq_buff;
            }
        }
        ref_fs.close();
    }
    index_fs.close();
    cout << "[INFO] Finished loading reference sequence" << endl;
}



bool is_number(const string& s)
{
    string::const_iterator it = s.begin();
    while (it != s.end() && isdigit(*it))
    {
        it++;
    }
    return !s.empty() && it == s.end();
}


